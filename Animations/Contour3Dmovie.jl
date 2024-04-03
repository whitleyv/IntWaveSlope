using GLMakie
using Measures
using Printf
using JLD2

#path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

#setname = "U250N100Lz100g100"
setname = "U250N100Lz100g100" #Lx300"
pmsetname = "U250N100Lz100g100"
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
#filepath1 = path_name * "DyesOnly_" * setname * ".jld2" 
#filepath2 = path_name * "DyegOnly_" * setname * ".jld2"
filepath1 = path_name * "Cs_mp_noS_" * setname * ".jld2"
#name_prefix = "IntWave_" * setname
#filepath = path_name * name_prefix * ".jld2"

#c1_timeseries = FieldTimeSeries(filepath, "Cg");
#c2_timeseries = FieldTimeSeries(filepath, "Cs");

#c1i = interior(c1_timeseries)[:, 1:900, 1:250, 1:2:end]
#c2i = interior(c2_timeseries)[:, 1:900, 1:250, 1:2:end]

f1 = jldopen(filepath1)
#f2 = jldopen(filepath2)

ycut = 875

c_timeseries1 = f1["ci"][:,1:ycut,:,:];
(cLx, cLy, cLz, tlength) = size(c_timeseries1)

xc = f1["xc"];
yc = f1["yc"][1:ycut];
zc = f1["zc"];
ctimes = f1["times"];
#c_timeseries2 = f2["c1i"];

#xc = 2:4:(cLx*4)
#yc = 2:4:(cLy*4)
#zc = -499:2:-1

#LxCex = 460
LxCex = 152

land = [curvedslope(y) for x in -5:4:LxCex, y in yc];

#xc, yc, zc = nodes(c_timeseries) #CCC

#land = [curvedslope(y) for x in xc, y in yc];

n = Observable(1)
#ctimes = 0:600:tlength*600

#title = @lift "Internal Wave Breaking, t = " * string(round(c_timeseries.times[$n], digits=2)) * ", Tσ = "*string(round(c_timeseries.times[$n]/pm.Tσ, digits=2))
title = @lift "Tracer Concentration, t = " * string(round(ctimes[$n]/3600, digits=1)) * " hrs, Tσ = "*string(round(ctimes[$n]/pm.Tσ, digits=2))

#f = Figure(resolution = (1500, 550),fontsize=26) 
f = Figure(resolution = (1600, 600),fontsize=26) 
ga = f[1, 1] = GridLayout()

ax1 = Axis3(ga[1, 1], azimuth = π/8, # rotation of plot
            elevation = 0.15, 
            aspect =(1,3,1), xtickwidth = 0,
            xticks = ([2,LxCex], [string.(cLx*4), "0"]), xticklabelpad = -25, 
            zticks = ([-500, -250, 0]),
            yticks = ([1000, 2000, 3000]),
            xlabeloffset = 10, ylabeloffset = 25, zlabeloffset = 100,
            xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",
            title = title)
limits!((0,LxCex), (0,3500), (-500,0))
surface!(ax1, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100),)

# topography yz cut lower and uppbounds
lower = Point3f.(LxCex, yc, -500);
upper = Point3f.(LxCex, yc, curvedslope.(yc));
band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

c = @lift log10.(clamp.(c_timeseries1[:,:,:,$n], 7e-5,1));
c_flat = @lift log10.(clamp.(c_timeseries1[cLx,:,:,$n], 7e-5,1));

#c = @lift log10.(clamp.(interior(c_timeseries[$n], :, :, :), 7*1e-5,1))

cg = contour!(ax1, xc, yc, zc, c, levels = (-4:0.05:0), colormap = (:thermal), alpha = 0.15)
contourf!(ax1, yc, zc, c_flat; levels = (-4:0.05:0), colormap = :thermal,
    transformation=(:yz, LxCex))

cb = Colorbar(ga[1, 2], limits = (-4,0), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
    label = "Tracer Concentration",lowclip = :white, size = 25)
#cb.alignmode = Mixed(right = 0)
colgap!(ga, -50)

savename = "paper_3D_tracer_" * setname
#savename = "NOPP_slope3dwaveseries_" * setname
#apath  = path_name * "Analysis/"
apath  =  "Analysis/PaperFigures/FinalPaperFigures/FinalFiguresUsed/SupplementaryInfo/"

frames = 1:tlength
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end

#=
fb = Figure(resolution = (1500, 550),fontsize=26) 
ga = fb[1, 1] = GridLayout()

ax1 = Axis3(ga[1, 1], azimuth = π/8, # rotation of plot
            elevation = 0.15, 
            aspect =(1,3,1), xtickwidth = 0,
            xticks = ([2,LxCex], [string.(cLx*4), "0"]), xticklabelpad = -25, 
            zticks = ([-500, -250, 0]),
            yticks = ([1000, 2000, 3000]),
            xlabeloffset = 10, ylabeloffset = 25, zlabeloffset = 100,
            xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",
            title = title)
limits!((0,LxCex), (0,3500), (-500,0))
surface!(ax1, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100),)

# topography yz cut lower and uppbounds
lower = Point3f.(LxCex, yc, -500);
upper = Point3f.(LxCex, yc, curvedslope.(yc));
band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

c = @lift log10.(clamp.(c_timeseries2[:,:,:,$n], 7e-5,1));
#c = @lift log10.(clamp.(interior(c_timeseries[$n], :, :, :), 7*1e-5,1))

global cg = contour!(ax1, xc, yc, zc, c, levels = (-4:0.05:0), colormap = (:thermal), alpha = 0.15)

cb = Colorbar(ga[1, 2], limits = (-4,0), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
    label = "Tracer Concentration",lowclip = :white, size = 25)
#cb.alignmode = Mixed(right = 0)

savename = "gauss3dwaveseries_" * setname
#apath  = path_name * "Analysis/"
#apath  =  "../Analysis/Plots/"

frames = 1:tlength
record(fb, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end
=#