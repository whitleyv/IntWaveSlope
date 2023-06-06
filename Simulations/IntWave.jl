ENV["GKSwstype"] = "nul"

using Statistics
using Printf
using CUDA: has_cuda_gpu
using ArgParse

using Oceananigans

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

# internal wave parameters
const gausW_width = 500/3
const gausW_center = 4500                            # y position 4500 to center forced wave

# sponge regions
const Sp_Region_right = 500                               # size of sponge region on RHS
const Sp_Region_left = 500
const Sp_Center = 20
const Sp_Region_z = 50

# dye parameters
const dye_height = 20                               # height of initial dye
const dye_smoothing = 10                            # tanh initial dye smoothing
const dye_top = 300                                 # tanh initial dye going negative at
const dye_centz = -250                              # center z value for gaussian IC

# drag specifications (Monin-Obukhov drag coefficient)
const z₀ = 1e-2 # Charnock roughness
const κᵥₖ = 0.4 # Von Karman constant

# assuming isotropic
const dz = 2.0 # if we use dy=4 Cd drops from 7.5 × 10⁻³ to 5.6 × 10⁻³
const Cd = (κᵥₖ / log(dz / 2z₀))^2

function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
      "paramset"
          help = "sets which parameters to use"
          default = "U100N100Lz100g100"
      "--fruntime", "-T"
          help = "how many inertial periods sim runs"
          arg_type = Float64
          default = 5.0
      "--resScale"
     	  help = "scale of resolution from default dh = 4 and dz=2"
	        arg_type = Float64
	        default = 1.0  
      "path"
          help = "pathname to save data under"
          default = ""
      "--dflag", "-d"
          help = "determines if you run full diagnostics"
          action = :store_true
  end
  return parse_args(s)
end

args=parse_commandline()

@info("command line args:")
for (arg,val) in args
  @info("   $arg => $val")
end

path_name = args["path"]
setname = args["paramset"]

###########-------- Parameters ----------------#############

@info "Loading in parameters..."

resS = args["resScale"]

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(setname))

dzr = pm.dz * resS
dhr = pm.dh * resS

pm = merge(pm, (;dzr=dzr, dhr=dhr, Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/dzr),
		            nx = round(Int, pm.Lx/dhr),
                m = -π/pm.Lzˢ,
                l = sqrt(((π/pm.Lzˢ)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ

Sp_extra = ifelse(z_start>0, 250.0, 0.0)
Ly = pm.Lyˢ+Sp_extra
ny = round(Int,Ly/pm.dhr)
slope_end = pm.Lzˢ/pm.Tanα

pm = merge(pm, (;Ly=Ly,ny=ny, slope_end=slope_end, Sp_extra=Sp_extra))

# if slope is in different spot than usual, need to move the curved part too!
const zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame

const dye_centy = dye_centz/-pm.Tanα

stop_time = pm.Tfˢ*args["fruntime"]
 
###########-------- GRID SET UP ----------------#############

@info "Setting up the grid..."
# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()

underlying_grid = RectilinearGrid(arch; size = (pm.nx, pm.ny, pm.nz), 
                                            x = (0, pm.Lx), y = (0, pm.Ly),
                                            z = (-pm.Lzˢ, z_start),
                                            halo = (4, 4, 4),
                                            topology = (Periodic, Bounded, Bounded))

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

@inline is_immersed(x, y, z) = z < curvedslope(y)

immersed_grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(is_immersed))

###########-------- FORCING WAVE AND SPONGE FUNCTIONS ----------------#############

@info "Including the forcing for internal waves and sponge layers..."

# internal wave forcing
@inline gaus(y) = exp( -(y - gausW_center)^2 / (2 * gausW_width^2))
@inline Fv_wave(x, y, z, t, p) =   p.U₀ * p.σ *               sin(p.l * y + p.m * z - p.σ * t) * gaus(y)

#sponge regions
@inline mask2nd(X) = heaviside(X)* X^2
@inline right_mask(y, p) = mask2nd((y-p.Lyˢ+Sp_Region_right+p.Sp_extra)/(Sp_Region_right+p.Sp_extra))

@inline left_mask(y) = mask2nd((Sp_Region_left-y)/Sp_Region_left)
@inline top_mask(z)= 0.5*(tanh((z+Sp_Center)/Sp_Region_z) + 1)
@inline corner_mask(y, z) = top_mask(z) * left_mask(y)

@inline u_sponge(x, y, z, t, u, p) = - 0.001 * right_mask(y, p) * u
@inline v_sponge(x, y, z, t, v, p) = - 0.001 * right_mask(y, p) * v
@inline w_sponge(x, y, z, t, w, p) = - 0.001 * right_mask(y, p) * w
@inline b_sponge(x, y, z, t, b, p) =   0.001 * right_mask(y, p) * (p.Ñ^2 * z - b)

@inline vel_sponge_corner(y, z, u) = - 0.002 * corner_mask(y, z) * u
@inline b_sponge_corner(y, z, b, p) =   0.002 * corner_mask(y, z) * (p.Ñ^2 * z - b)

@inline force_u(x, y, z, t, u, p) = u_sponge(x, y, z, t, u, p) + vel_sponge_corner(y, z, u)
@inline force_v(x, y, z, t, v, p) = Fv_wave(x, y, z, t, p) + v_sponge(x, y, z, t, v, p) + vel_sponge_corner(y, z, v)
@inline force_w(x, y, z, t, w, p) = w_sponge(x, y, z, t, w, p) + vel_sponge_corner(y, z, w)
@inline force_b(x, y, z, t, b, p) = b_sponge(x, y, z, t, b, p) + b_sponge_corner(y, z, b, p)

force_params = (; pm.U₀, pm.f, pm.σ, pm.Ñ,
                  pm.l, pm.m, pm.Lyˢ, pm.Sp_extra)

u_forcing = Forcing(force_u, field_dependencies = :u, parameters = force_params)
v_forcing = Forcing(force_v, field_dependencies = :v, parameters = force_params)
w_forcing = Forcing(force_w, field_dependencies = :w, parameters = force_params)
b_forcing = Forcing(force_b, field_dependencies = :b, parameters = force_params)

###########-------- BOUNDARY CONDITIONS ---------------#############

@info "Boundary conditions..."

# BC functions flat in the bounded directions
@inline bottom_drag_u(x, y, t, u, v, w) = - Cd * u * sqrt(u^2 + v^2 + w^2)
@inline bottom_drag_v(x, y, t, u, v, w) = - Cd * v * sqrt(u^2 + v^2 + w^2)

# BC functions not flat in the bounded directions for IB case
@inline bottom_drag_u(x, y, z, t, u, v, w) = - Cd * u * sqrt(u^2 + v^2 + w^2)
@inline bottom_drag_v(x, y, z, t, u, v, w) = - Cd * v * sqrt(u^2 + v^2 + w^2)

u_drag_bc = FluxBoundaryCondition(bottom_drag_u, field_dependencies=(:u, :v, :w))
v_drag_bc = FluxBoundaryCondition(bottom_drag_v, field_dependencies=(:u, :v, :w))

u_bcs = FieldBoundaryConditions(bottom=u_drag_bc, immersed=u_drag_bc)
v_bcs = FieldBoundaryConditions(bottom=v_drag_bc, immersed=v_drag_bc)

boundary_conditions = (; u = u_bcs, v = v_bcs)

###########-------- STARTING UP MODEL/ ICs ---------------#############

start_time = time_ns() # time to start the clock for how long it's been running

@info "Setting up the model..."

kwargs_model = (; grid = immersed_grid,
                closure = SmagorinskyLilly(),
              advection = WENO(),
               coriolis = FPlane(f=pm.f),
                forcing = (v=v_forcing, w=w_forcing, u=u_forcing, b=b_forcing),
    boundary_conditions = boundary_conditions,
               buoyancy = BuoyancyTracer(),
            timestepper = :RungeKutta3)

model = NonhydrostaticModel(; tracers = (:b, :Cg, :Cs,), kwargs_model...)

run(`nvidia-smi`)

@info "Setting initial conditions..."
# linear stratification
b₀(x, y, z) = pm.Ñ^2 * z  

# along the slope
@inline dye(x, y, z, dH, sm, kt) = 0.5*( tanh( (dH-z) / sm ) + tanh(kt - z) )
# only dye in the fluid
@inline above_slope(y, z) = ifelse(z >= curvedslope(y), 1.0, 0.0)
# only as far as the slope goes
@inline slope_side(y) = ifelse(y <= pm.slope_end, 1.0, 0.0)

# gaussian in center
@inline exparg(y, z) = ((y - dye_centy)^2 + (z - dye_centz)^2) / (2*dye_height^2)
@inline dye_gauss(y, z) = exp(-1.0 * exparg(y,z))

# gaussian IC
@inline c₀1(x, y, z) = above_slope(y, z) * dye_gauss(y, z)
# tanh IC
@inline c₀2(x, y, z) = slope_side(y) * above_slope(y, z) * dye(x, y, z, curvedslope(y)+dye_height, dye_smoothing, dye_top)

set!(model, u=0, v=0, w=0, b=b₀, Cg=c₀1, Cs=c₀2)
 
###########-------- SIMULATION SET UP ---------------#############
@info "Simultion set up..."

name_prefix = "IntWave_" * setname
@show data_path = path_name * name_prefix

Δt = 0.8 * underlying_grid.Δyᵃᶜᵃ
simulation = Simulation(model, Δt = Δt, stop_time = stop_time)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=10.0, min_Δt=0.0001) # dec cfl 0.9 -> .5
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) # dec. int 500 -> 5 (max)

progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9))

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(1200.0))

###########-------- DIAGNOSTICS --------------#############

@info "Adding Diagnostics..."
include("diagnostics.jl")
outputs=get_outputs_tuple(model, Fulld=args["dflag"])

# don't want to overwrite existing if checkpointed
simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs;
                                                          schedule = TimeInterval(600.0),
                                                          filename = data_path,
                                                          overwrite_existing = true)


###########-------- RUN! --------------#############

run(`nvidia-smi`)

@info "Run...."

run!(simulation)

