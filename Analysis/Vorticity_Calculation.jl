using Statistics
using Printf
using Measures
using JLD2
using Oceananigans

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"
apath = path_name * "Analysis/"
sn = "U250N100Lz100g100"

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ

Sp_extra = ifelse(z_start>0, 250.0, 0.0)
Ly = pm.Lyˢ+Sp_extra
ny = round(Int,Ly/4)
slope_end = pm.Lzˢ/pm.Tanα

pm = merge(pm, (;Ly=Ly,ny=ny, slope_end=slope_end, Sp_extra=Sp_extra))

# if slope is in different spot than usual, need to move the curved part too!
const zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

name_prefix = "IntWave_mp_noS_" * sn
filepath = path_name * name_prefix * ".jld2"

u_timeseries = FieldTimeSeries(filepath,"u");
v_timeseries = FieldTimeSeries(filepath,"v");
w_timeseries = FieldTimeSeries(filepath,"w");
b_timeseries = FieldTimeSeries(filepath, "b"); # only using for centered nodes

ylength = 880
ui = interior(u_timeseries)[:,1:ylength,:,:];
vi = interior(v_timeseries)[:,1:ylength+1,:,:];
wi = interior(w_timeseries)[:,1:ylength,:,:];

xb, yb, zb = nodes(b_timeseries) #CFC

# centering data to middle
u_ccc = 0.5 .* (ui[1:end-1,:,:,:] .+ ui[2:end,:,:,:]);
v_ccc = 0.5 .* (vi[:,1:end-1,:,:] .+ vi[:,2:end,:,:]);
w_ccc = 0.5 .* (wi[:,:,1:end-1,:] .+ wi[:,:,2:end,:]);

Δx = zb[2]-xb[1]
Δy = yb[2]-yb[1]
Δz = zb[2]-zb[1]
# if x deriv [1,2,3,4], then [1,2,3,4]
# if y deriv [2,1,3,4], then [2,1,3,4]
# if z deriv [3,2,1,4], then [3,2,1,4]
function deriv(A, h, dim_order)
    Ap = permutedims(A, dim_order)
    step = 1/h
    step2 = 1/(h*2)
    da = zeros(size(Ap));

    da[1, :, :, :] = (Ap[2,:,:,:] .- Ap[1,:,:,:]) .* step
    da[end, :, :, :] = (Ap[end, :,:,:] .- Ap[end-1,:,:,:]) .* step
    da[2:end-1,:, :, :] = (Ap[3:end,:,:,:] .- Ap[1:end-2,:,:,:]) .* step2
        
    dap = permutedims(da, dim_order)
    return dap
end

# xlength - 1
du_dy = deriv(u_ccc, Δy, [2,1,3,4]);
du_dz = deriv(u_ccc, Δz, [3,2,1,4]);

# 
dv_dx = deriv(v_ccc, Δx, [1,2,3,4]);
dv_dz = deriv(v_ccc, Δz, [3,2,1,4]);

# zlength - 1
dw_dx = deriv(w_ccc, Δx, [1,2,3,4]);
dw_dy = deriv(w_ccc, Δy, [2,1,3,4]);

# boolean to avoid near the wall
zlength = length(zb)
xlength = length(xb)

Ygrid = reshape(repeat(yb[1:ylength], (xlength-1)*zlength), ylength, zlength, xlength-1)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb';  # all the values greater than slope

boolZYX = permutedims(boolZY, [3, 1, 2])

ω_x = (dw_dy[2:end,:,:,:] .- dv_dz[2:end,:,:,:]) .* boolZYX;
ω_y = (du_dz[:,:,:,:] .- dw_dx[2:end,:,:,:]) .* boolZYX;
ω_z = (dv_dx[2:end,:,:,:] .- du_dy[:,:,:,:]) .* boolZYX;

ω_mag = sqrt.(ω_x.^2 .+ ω_y.^2 .+ ω_z.^2);
savename = apath * "vorticity_terms_" * sn * ".jld2"

jldsave(savename;
ω_x, ω_y, ω_z,
xc=xb[2:end], 
yc=yb[1:ylength], 
zc=zb,
times=b_timeseries.times)
