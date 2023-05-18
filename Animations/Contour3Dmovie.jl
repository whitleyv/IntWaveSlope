using CairoMakie
using GLMakie
using Measures
using Oceananigans
using Printf

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
setname = "U250N100Lz100g100"

@info "getting data from: " * setname

name_prefix = "IntWave_" * setname
filepath = path_name * name_prefix * ".jld2"

c1_timeseries = FieldTimeSeries(filepath, "Cg");
c2_timeseries = FieldTimeSeries(filepath, "Cs");

xc, yc, zc = nodes(c_timeseries) #CCC

include("parameters.jl")

pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSame = 1336.6                           # point where planar and curved corner math up the best

ΔySlopeSame = 0

@inline heaviside(X) = ifelse(X < 0, 0., 1.) # heaviside returns 1 if x >= 0
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lz + pm.Lz * exp( - (y - ystar)^2 / (2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα * y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = @. linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

land = [curvedslope(y) for x in xc, y in yc];

f = Figure(resolution = (1000, 700),fontsize=26) 
ga = f[1, 1] = GridLayout()

ax1 = Axis3(ga[1, 1] , azimuth = π/8, # rotation of plot
            elevation = 0.15, 
            aspect =(1,3,1), xtickwidth = 0,
            xticks = ([0,152], ["152", "0"]), xticklabelpad = -25, 
            zticks = ([-500, -250, 0]),
            xlabeloffset = 10, ylabeloffset = 20, zlabeloffset = 100,
            xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",)
limits!((0,152), (0,5500), (-500,0))
surface!(ax1, -5:4:155, yc, land, shading = true, color=fill((:black, .7),100,100),)

# topography yz cut lower and uppbounds
lower = Point3f.(152, yc, -500);
upper = Point3f.(152, yc, curvedslope.(yc));
band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

n = Observable(1)

c = @lift log10.(clamp.(interior(c_timeseries[$n], :, :, :), 7*1e-5,1))


cg = contour!(ax1, xc, yc, zc, c, levels = (-4:0.05:0), colormap = :thermal)

title = @lift "Internal Wave Breaking, t = " * string(round(c_timeseries.times[$n], digits=2)) * ", Tσ = "*string(round(c_timeseries.times[$n]/pm.Tσ, digits=2))

cb = Colorbar(ga[1, 2], cg, ticks = (-4:1:0, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰"] ),
    label = "Tracer Concentration",)

savename = "gauss3dwaveseries_" * setname
apath  = path_name * "Analysis/"
frames = 1:length(c_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=24) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end
