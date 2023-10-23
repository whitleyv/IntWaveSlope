using Measures
using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

setname = "U300N100Lz100g100"

resS = 1.0

@info "Loading in parameters..."

include("parameters.jl")
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

name_mp_prefix = "IntWave_mp_nExD_" * setname
name_reg_prefix = "IntWave_" * setname 

name_mp_prefix = "IntWave_mp_" * setname
name_reg_prefix = "IntWave_smdt_" * setname 

filepath_mp = path_name * name_mp_prefix * ".jld2"
filepath_reg = path_name * name_reg_prefix * ".jld2"

b_mp_timeseries = FieldTimeSeries(filepath_mp,"b");
b_reg_timeseries = FieldTimeSeries(filepath_reg,"b");

cg_reg_timeseries = FieldTimeSeries(filepath_reg, "Cg");
cg_mp_timeseries = FieldTimeSeries(filepath_mp, "Cg");

cs_reg_timeseries = FieldTimeSeries(filepath_reg, "Cs");
cs_mp_timeseries = FieldTimeSeries(filepath_mp, "Cs");

v_mp_timeseries = FieldTimeSeries(filepath_mp,"v");
v_reg_timeseries = FieldTimeSeries(filepath_reg,"v");

xc, yc, zc = nodes(b_mp_timeseries) #CCC
xv, yv, zv = nodes(v_mp_timeseries) #CCC

xc2, yc2, zc2 = nodes(b_reg_timeseries) #CCC

land = curvedslope.(yc) 

delta = pm.U₀/pm.Ñ

##############
#   MOVIE OF BETTER TIME RESOLUTION
##############

n = Observable(1)
xlocat = 19

v_mp = @lift interior(v_mp_timeseries[$n],xlocat, :, :,);
b_mp = @lift interior(b_mp_timeseries[$n], xlocat, :, :);
cg_mp = @lift log10.(clamp.(interior(cg_mp_timeseries[$n],xlocat, :, :), 1e-8,1));
cs_mp = @lift log10.(clamp.(interior(cs_mp_timeseries[$n],xlocat, :, :), 1e-8,1));

#v_reg = @lift interior(v_reg_timeseries[2 * $n - 1],xlocat, :, :);
#b_reg = @lift interior(b_reg_timeseries[2 * $n - 1], xlocat, :, :);
#cg_reg = @lift log10.(clamp.(interior(cg_reg_timeseries[2 * $n - 1],xlocat, :, :), 1e-8,1));
#cs_reg = @lift log10.(clamp.(interior(cs_reg_timeseries[2 * $n - 1],xlocat, :, :), 1e-8,1));

v_reg = @lift interior(v_reg_timeseries[$n],xlocat, :, :);
b_reg = @lift interior(b_reg_timeseries[$n], xlocat, :, :);
cg_reg = @lift log10.(clamp.(interior(cg_reg_timeseries[$n],xlocat, :, :), 1e-8,1));
cs_reg = @lift log10.(clamp.(interior(cs_reg_timeseries[$n],xlocat, :, :), 1e-8,1));


title  = @lift "Internal Wave Breaking, δ=85.7 m, Tσ = "*string(round(b_mp_timeseries.times[$n]/pm.Tσ, digits=2))

f = Figure(resolution = (1800, 800),fontsize=26) 

ga = f[1, 1] = GridLayout()
axv = Axis(ga[1, 1], ylabel = "m < 0, z [m]")
axb = Axis(ga[1, 3])
axc = Axis(ga[1, 5])
axcb = Axis(ga[1,7])

axv2 = Axis(ga[2, 1], ylabel = "m > 0, z [m]", xlabel = "y [m]")
axb2 = Axis(ga[2, 3], xlabel = "y [m]")
axc2 = Axis(ga[2, 5], xlabel = "y [m]")
axc2b = Axis(ga[2, 7], xlabel = "y [m]")

axc2b.xticks = 500:1000:1500
axc2.xticks = 500:1000:1500
axv2.xticks = 500:1000:1500
axb2.xticks =  500:1000:1500

axv2.yticks = [-250, 0]
axv.yticks = [-250, 0]

limits!(axv, 0, 2000, -500, 0)
limits!(axv2, 0, 2000, -500, 0)

limits!(axc, 0, 2000, -500, 0)
limits!(axc2, 0, 2000, -500, 0)

limits!(axb, 0, 2000, -500, 0)
limits!(axb2, 0, 2000, -500, 0)

limits!(axcb, 0, 2000, -500, 0)
limits!(axc2b, 0, 2000, -500, 0)

hidedecorations!(axb)
hidedecorations!(axcb)
hidedecorations!(axc)
hidexdecorations!(axv)
hideydecorations!(axc2)
hideydecorations!(axb2)
hideydecorations!(axc2b)

# inlcude lines where split happens in bflux:
# z = -334 and -168
z1 = -334
z2 = -168
y1 = -z1/pm.Tanα
y2 = -z2/pm.Tanα

land_pdel = (curvedslope.(yc) .+ delta)[1:382]

global hmv = heatmap!(axv, yv, zv, v_reg, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
lines!(axv, yc, land, color=:black, linewidth = 3)
lines!(axv, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axv, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axv, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

global hmb = heatmap!(axb, yc2, zc2, b_reg, colormap= :thermal, colorrange = (-0.007, 0),)
contour!(axb, yc2, zc2, b_reg, color = :black, linewidth = 2 , levels = -0.006:0.0005:0, alpha = 0.5)
lines!(axb, yc, land, color=:black, linewidth = 3)
lines!(axb, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axb, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axb, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

global hmcg = heatmap!(axc, yc, zc, cg_reg, colormap = :thermal, colorrange = (-4, 0))
lines!(axc, yc, land, color=:black, linewidth = 3)
lines!(axc, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axc, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axc, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

global hmcs = heatmap!(axcb, yc, zc, cs_reg, colormap = :thermal, colorrange = (-4, 0))
lines!(axcb, yc, land, color=:black, linewidth = 3)
lines!(axcb, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axcb, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axcb, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

heatmap!(axv2, yv, zv, v_mp, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
lines!(axv2, yc, land, color=:black, linewidth = 3)
lines!(axv2, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axv2, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axv2, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

heatmap!(axb2, yc, zc, b_mp, colormap= :thermal, colorrange = (-0.007, 0),)
contour!(axb2, yc, zc, b_mp, color = :black, linewidth = 2 , levels = -0.006:0.0005:0, alpha = 0.5)
lines!(axb2, yc, land, color=:black, linewidth = 3)
lines!(axb2, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)
lines!(axb2, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray30, linewidth = 2)
lines!(axb2, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray30, linewidth = 2)

heatmap!(axc2, yc, zc, cg_mp, colormap = :thermal, colorrange = (-4, 0))
lines!(axc2, yc, land, color=:black, linewidth = 3)
lines!(axc2, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axc2, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axc2, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

heatmap!(axc2b, yc, zc, cs_mp, colormap = :thermal, colorrange = (-4, 0))
lines!(axc2b, yc, land, color=:black, linewidth = 3)
lines!(axc2b, yc[1:382], land_pdel, color=:white, linewidth = 2, linestyle = :dash)
lines!(axc2b, 1:1:y1, ones(round(Int,y1)).* z1, color = :gray70, linewidth = 2)
lines!(axc2b, 1:1:y2, ones(round(Int,y2)).* z2, color = :gray70, linewidth = 2)

wtimes = b_mp_timeseries.times./ pm.Tσ

Label(ga[1, 1:2, Top()], title,
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)

cb1 = Colorbar(ga[1:2,2], hmv, ticks = (-0.2:.1:0.2), size =35)
cb2 = Colorbar(ga[1:2,4], hmb, ticks = (-5e-3:2e-3:-1e-3, ["-0.005", "-0.003", "-0.001"] ), size = 35)
cb3 = Colorbar(ga[1:2,6], hmcg, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size = 35)

colsize!(ga, 2, Relative(0.04))
colsize!(ga, 4, Relative(0.04))
colsize!(ga, 6, Relative(0.04))

colgap!(ga, 15)
rowgap!(ga, 5)

savename = "contours_mtest_" * setname
apath  = path_name * "Analysis/"


frames = 1:length(b_mp_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    @info msg
    n[] = j
end
