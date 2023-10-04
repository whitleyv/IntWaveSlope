using Measures
using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie
using ArgParse

ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
      "paramset"
          help = "sets which parameters to use"
          default = "U100N100Lz100g100"
      "--resScale"
     	  help = "scale of resolution from default dh = 4 and dz=2"
	        arg_type = Float64
	        default = 1.0  
     "path"
          help = "pathname to save data under"
          default = ""
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

resS = args["resScale"]
@info "Loading in parameters..."

include("../parameters.jl")
pm = getproperty(SimParams(), Symbol(setname))

dzr = pm.dz * resS
dhr = pm.dh * resS

pm = merge(pm, (; dzr=dzr, dhr=dhr, Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
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

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

@info "Pulling Data..."
path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

name1_prefix = "IntWave_" * setname
name2_prefix = "IntWave_smdt_" * setname

filepath1 = path_name * name1_prefix * ".jld2"
filepath2 = path_name * name2_prefix * ".jld2"

# late releases of dye
filepath3 = path_name * "cgl_" * setname * ".jld2"
filepath4 = path_name * "cgl2_" * setname * ".jld2"

f3 = jldopen(filepath3)
CGli = f3["Cgli"];

f4 = jldopen(filepath4)
CGl2i = f4["CGl2i"];

# gaussian at slope but with increased time resolution
Cg_timeseries = FieldTimeSeries(filepath2,"Cg");
Cs_timeseries = FieldTimeSeries(filepath2,"Cs");

spinlength = size(CGli)[4]
spinlength2 = size(CGl2i)[4]
spinupL = 161 - spinlength
spinupL2 = 161 - spinlength2

v_timeseries = FieldTimeSeries(filepath2,"v");
b_timeseries = FieldTimeSeries(filepath2,"b");
ε_timeseries = FieldTimeSeries(filepath2, "ϵ");

v_timeseries = FieldTimeSeries(filepath1,"v");
b_timeseries = FieldTimeSeries(filepath1,"b");
ε_timeseries = FieldTimeSeries(filepath1, "ϵ");
Cgr_timeseries = FieldTimeSeries(filepath1,"Cgr");

xc, yc, zc = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC

land = curvedslope.(yc) 

delta = pm.U₀/pm.Ñ

##############
#   CAIRO PLOT MOVIE
##############

n = Observable(1)
xlocat = 19
telngth = 161
zlength = 250
ylength = 1000

v = @lift interior(v_timeseries[$n],xlocat, :, :,);
b = @lift interior(b_timeseries[$n], xlocat, :, :);
cg = @lift log10.(clamp.(interior(Cg_timeseries[2 * $n - 1],xlocat, :, :), 1e-8,1));
cs = @lift log10.(clamp.(interior(Cs_timeseries[2 * $n - 1], xlocat, :, :), 1e-8,1));
cgr = @lift log10.(clamp.(interior(Cgr_timeseries[$n], xlocat, :, :),  1e-8,1));

function cglf(n)
    if n > spinupL
        cgl = log10.(clamp.(CGli[xlocat, :, :, n - spinupL],  1e-8,1));
    else
        cgl = log10.(clamp.(zeros(1439, 325),  1e-8,1));
    end
    return cgl
end

cgl = @lift(cglf($n))

function cgl2f(n)
    if n > spinupL2
        cgl2 = log10.(clamp.(CGl2i[xlocat, :, :, n - spinupL2],  1e-8,1));
    else
        cgl2 = log10.(clamp.(zeros(1439, 325),  1e-8,1));
    end
    return cgl2
end
cgl2 = @lift(cgl2f($n))

ε = @lift log10.(clamp.(interior(ε_timeseries[$n], xlocat, :, :), 1e-10,1));

title = @lift "Internal Wave Breaking, δ=85.7 m, t = " * string(round(b_timeseries.times[$n]/3600, digits=2)) * " hrs, Tσ = "*string(round(b_timeseries.times[$n]/pm.Tσ, digits=2))

f = Figure(resolution = (1450, 1000),fontsize=26) 

ga = f[1, 1] = GridLayout()
axv = Axis(ga[1, 1], ylabel = "z [m]")
axb = Axis(ga[1, 3])

axcg = Axis(ga[2, 1], ylabel = "z [m]")
axcs = Axis(ga[2, 3])

axcgl = Axis(ga[3, 1], ylabel = "z [m]")
axcgl2 = Axis(ga[3, 3])

axe = Axis(ga[4, 1], ylabel = "z [m]", xlabel = "y [m]")
axn = Axis(ga[4, 3], xlabel = "y [m]")

axe.xticks = 500:1000:1500
axn.xticks = 500:1000:1500

axv.yticks = [-250, 0]
axcg.yticks = [-250, 0]
axcgl.yticks = [-250, 0]
axe.yticks = [-250, 0]

limits!(axv, 0, 2000, -500, 0)
limits!(axcg, 0, 2000, -500, 0)
limits!(axcs, 0, 2000, -500, 0)
limits!(axcgl, 0, 2000, -500, 0)
limits!(axcgl2, 0, 2000, -500, 0)
limits!(axb, 0, 2000, -500, 0)
limits!(axe, 0, 2000, -500, 0)
limits!(axn, 0, 2000, -500, 0)

hidedecorations!(axb)
hidedecorations!(axcs)
hidexdecorations!(axv)
hidexdecorations!(axcg)
hidexdecorations!(axcgl)
hidedecorations!(axcgl2)
hideydecorations!(axn)

# inlcude lines where split happens in bflux:
# z = -334 and -168
z1 = -334
z2 = -168
y1 = -z1/pm.Tanα
y2 = -z2/pm.Tanα

land_pdel = (curvedslope.(yc) .+ delta)[1:382]

global hmv = heatmap!(axv, yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
lines!(axv, yc, land, color=:black, linewidth = 3)
lines!(axv, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axv, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axv, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

global hmb = heatmap!(axb, yc, zc, b, colormap= :thermal, colorrange = (-0.007, 0),)
contour!(axb, yc, zc, b, color = :black, linewidth = 2 , levels = -0.006:0.0005:0, alpha = 0.5)
lines!(axb, yc, land, color=:black, linewidth = 3)
lines!(axb, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axb, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axb, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

global hmcg = heatmap!(axcg, yc, zc, cg, colormap = :thermal, colorrange = (-4, 0))
lines!(axcg, yc, land, color=:black, linewidth = 3)
lines!(axcg, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axcg, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axcg, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

global hmcs = heatmap!(axcs, yc, zc, cs, colormap = :thermal, colorrange = (-4, 0))
lines!(axcs, yc, land, color=:black, linewidth = 3)
lines!(axcs, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axcs, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axcs, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

global hmcgl = heatmap!(axcgl, yc[1:1000], zc, cgl, colormap = :thermal, colorrange = (-4, 0))
lines!(axcgl, yc, land, color=:black, linewidth = 3)
lines!(axcgl, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axcgl, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axcgl, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

global hmcgl2 = heatmap!(axcgl2, yc[1:1000], zc, cgl2, colormap = :thermal, colorrange = (-4, 0))
lines!(axcgl2, yc, land, color=:black, linewidth = 3)
lines!(axcgl2, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axcgl2, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axcgl2, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

global hme = heatmap!(axe, yc, zc, ε, colormap = :thermal, colorrange = (-9, -5))
lines!(axe, yc, land, color=:black, linewidth = 3)
lines!(axe, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axe, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axe, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

global hmn= heatmap!(axn, yc, zc, cgr, colormap = :thermal, colorrange = (-4, 0))
lines!(axn, yc, land, color=:black, linewidth = 3)
lines!(axn, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axn, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axn, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

Label(ga[1, 1:4, Top()], title,
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)

cb1 = Colorbar(ga[1,2], hmv, ticks = (-0.2:.1:0.2), size =35)
cb2 = Colorbar(ga[1,4], hmb, ticks = (-5e-3:2e-3:-1e-3, ["-0.005", "-0.003", "-0.001"] ), size = 35)
cb3 = Colorbar(ga[2:3,2], hmcg, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size = 35)
cb4 = Colorbar(ga[2:4,4], hmcs, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size = 35)
cb5 = Colorbar(ga[4,2], hme, ticks =  (-9:1:-6, ["10⁻⁹", "10⁻⁸", "10⁻⁷", "10⁻⁶"] ), size = 35)

colsize!(ga, 2, Relative(0.04))
colsize!(ga, 4, Relative(0.04))

colgap!(ga, 15)
rowgap!(ga, 5)

savename = "contours_allg_" * setname
apath  = path_name * "Analysis/"

frames = 1:length(b_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end

##############
#   MOVIE OF BETTER TIME RESOLUTION
##############

n = Observable(1)
xlocat = 19
zlength = 250
ylength = 1000

v = @lift interior(v_timeseries[$n],xlocat, :, :,);
b = @lift interior(b_timeseries[$n], xlocat, :, :);
cg = @lift log10.(clamp.(interior(Cg_timeseries[$n],xlocat, :, :), 1e-8,1));
cs = @lift log10.(clamp.(interior(Cs_timeseries[$n],xlocat, :, :), 1e-8,1));

function titlef(n)
    if in(n, negleaps .+ 1)
        title  = "Internal Wave Breaking, δ=85.7 m, Tσ = " * string(round(b_timeseries.times[n]/pm.Tσ, digits=2) ) * " (b̄↓)"
    else
        title  = "Internal Wave Breaking, δ=85.7 m, Tσ = "*string(round(b_timeseries.times[n]/pm.Tσ, digits=2))
    end
    return title
end

title = @lift(titlef($n))

f = Figure(resolution = (1450, 800),fontsize=26) 

ga = f[1, 1] = GridLayout()
axv = Axis(ga[1, 1], ylabel = "z [m]")
axb = Axis(ga[1, 3])

axcg = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]")
axcs = Axis(ga[2, 3], xlabel = "y [m]")

axcg.xticks = 500:1000:1500
axcs.xticks = 500:1000:1500

axcs.yticks = [-250, 0]
axv.yticks = [-250, 0]

limits!(axv, 0, 2000, -500, 0)
limits!(axcg, 0, 2000, -500, 0)
limits!(axcs, 0, 2000, -500, 0)
limits!(axb, 0, 2000, -500, 0)

hidedecorations!(axb)
hidexdecorations!(axv)
hideydecorations!(axcs)

# inlcude lines where split happens in bflux:
# z = -334 and -168
z1 = -334
z2 = -168
y1 = -z1/pm.Tanα
y2 = -z2/pm.Tanα

land_pdel = (curvedslope.(yc) .+ delta)[1:382]

global hmv = heatmap!(axv, yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
lines!(axv, yc, land, color=:black, linewidth = 3)
lines!(axv, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axv, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axv, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

global hmb = heatmap!(axb, yc, zc, b, colormap= :thermal, colorrange = (-0.007, 0),)
contour!(axb, yc, zc, b, color = :black, linewidth = 2 , levels = -0.006:0.0005:0, alpha = 0.5)
lines!(axb, yc, land, color=:black, linewidth = 3)
lines!(axb, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axb, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axb, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

global hmcg = heatmap!(axcg, yc, zc, cg, colormap = :thermal, colorrange = (-4, 0))
lines!(axcg, yc, land, color=:black, linewidth = 3)
lines!(axcg, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axcg, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axcg, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

global hmcs = heatmap!(axcs, yc, zc, cs, colormap = :thermal, colorrange = (-4, 0))
lines!(axcs, yc, land, color=:black, linewidth = 3)
lines!(axcs, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axcs, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axcs, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

Label(ga[1, 1:2, Top()], title,
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)

cb1 = Colorbar(ga[1,2], hmv, ticks = (-0.2:.1:0.2), size =35)
cb2 = Colorbar(ga[1,4], hmb, ticks = (-5e-3:2e-3:-1e-3, ["-0.005", "-0.003", "-0.001"] ), size = 35)
cb3 = Colorbar(ga[2,2], hmcg, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size = 35)
cb4 = Colorbar(ga[2,4], hmcs, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size = 35)

colsize!(ga, 2, Relative(0.04))
colsize!(ga, 4, Relative(0.04))

colgap!(ga, 15)
rowgap!(ga, 5)

savename = "contours_smdt_" * setname
apath  = path_name * "Analysis/"

frames = 1:length(b_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end
