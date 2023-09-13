using GLMakie
using Measures
using Printf
using JLD2

#setname = "U350N100Lz100g100"
setname = "U350N100Lz130g100Lx300"

pmsetname = "U350N100Lz100g100"
include("../parameters.jl")

pm = getproperty(SimParams(), Symbol(pmsetname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSame = 1336.6                           # point where planar and curved corner math up the best

ΔySlopeSame = 0

@inline heaviside(X) = ifelse(X < 0, 0., 1.) # heaviside returns 1 if x >= 0
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp( - (y - ystar)^2 / (2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα * y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = @. linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

@info "getting data from: " * setname

path_name = "Data/"
filepath1 = path_name * "DyesOnly_" * setname * ".jld2"
filepath2 = path_name * "DyegOnly_" * setname * ".jld2"

f1 = jldopen(filepath1)
#f2 = jldopen(filepath2)

c_timeseries1 = f1["c2i"];
#c_timeseries2 = f2["c1i"];

(cLx, cLy, cLz, tlength) = size(c_timeseries1)

xc = 2:4:(cLx*4)
yc = 2:4:(cLy*4)
zc = -499:2:-1

LxCex = 460 #152 # 460

land = [curvedslope(y) for x in -5:4:LxCex, y in 0:4:3500];
#8.67, 17.33, 17
ctime_idx1 = round(Int, 17.33*3600/600)
#ctime_idx2 = round(Int, 17.67*3600/600)
ctime1 = ctime_idx1*600.0
#ctime2 = ctime_idx1*600.0

title = "Passive Tracer During Internal Wave Breaking, t = " * string(round(ctime1/3600, digits=2)) * " hrs, Tσ = "*string(round(ctime1/pm.Tσ, digits=2))

@info "Plotting..."

f = Figure(resolution = (1650, 600), fontsize=26) 
ga = f[1, 1] = GridLayout()

ax1 = Axis3(ga[1, 1], azimuth = π/8, # rotation of plot
            elevation = 0.15, 
            aspect =(1,3,1), xtickwidth = 0,
            xticks = ([2,LxCex], [string.(cLx*4), "0"]), xticklabelpad = -25, 
            zticks = ([-500, -250, 0]),
            yticks = ([1000, 2000, 3000]),
            xlabeloffset = 10, ylabeloffset = 30, zlabeloffset = 100,
            xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",
            title = title)
limits!((0,LxCex), (0,3500), (-500,0))
surface!(ax1, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100),)

# topography yz cut lower and uppbounds
lower = Point3f.(LxCex, yc, -500);
upper = Point3f.(LxCex, yc, curvedslope.(yc));
band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

c = log10.(clamp.(c_timeseries1[:,:,:,ctime_idx1], 7e-5,1));
c_flat = log10.(clamp.(c_timeseries1[cLx,:,:,ctime_idx1], 7e-5,1));

cg = contour!(ax1, xc, yc, zc, c, levels = (-4:0.05:0), colormap = (:thermal), alpha = 0.10)
contourf!(ax1, yc, zc, c_flat; levels = (-4:0.05:0), colormap = :thermal,
    transformation=(:yz, LxCex))

cb = Colorbar(ga[1, 2], limits = (-4,-1), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
    label = "Tracer Concentration",lowclip = :white, size = 25)

savename = "slope3d" * string.(ctime_idx1) * "_" * setname
apath  = "Analysis/Plots/"

save(apath * savename * ".jpg", f)

cg = false

if cg == true
    title = "Passive Tracer During Internal Wave Breaking, t = " * string(round(ctime2/3600, digits=2)) * " hrs, Tσ = "*string(round(ctime2/pm.Tσ, digits=2))

    f = Figure(resolution = (1700, 600), fontsize=26) 
    ga = f[1, 1] = GridLayout()

    ax1 = Axis3(ga[1, 1], azimuth = π/8, # rotation of plot
                elevation = 0.15, 
                aspect =(1,3,1), xtickwidth = 0,
                xticks = ([2,LxCex], [string.(cLx*4), "0"]), xticklabelpad = -25, 
                zticks = ([-500, -250, 0]),
                yticks = ([1000, 2000, 3000]),
                xlabeloffset = 10, ylabeloffset = 30, zlabeloffset = 100,
                xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",
                title = title)
    limits!((0,LxCex), (0,3500), (-500,0))
    surface!(ax1, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100),)

    # topography yz cut lower and uppbounds
    lower = Point3f.(LxCex, yc, -500);
    upper = Point3f.(LxCex, yc, curvedslope.(yc));
    band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

    c = log10.(clamp.(c_timeseries2[:,:,:,ctime_idx2], 7e-5,1));
    c_flat = log10.(clamp.(c_timeseries2[cLx,:,:,ctime_idx2], 7e-5,1));

    cg = contour!(ax1, xc, yc, zc, c, levels = (-4:0.05:-1), colormap = (:thermal), alpha = 0.15)
    contourf!(ax1, yc, zc, c_flat; levels = (-4:0.05:-1), colormap = :thermal,
        transformation=(:yz, LxCex))

    cb = Colorbar(ga[1, 2], limits = (-4,-1), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
        label = "Tracer Concentration",lowclip = :white, size = 25)

    savename = "gaus3d" * string.(ctime_idx2) * "_" * setname
    apath  = "Analysis/Plots/"

    save(apath * savename * ".png", f)
end