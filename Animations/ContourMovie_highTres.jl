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

name2_prefix = "IntWave_smdt_" * setname

filepath2 = path_name * name2_prefix * ".jld2"

# gaussian at slope but with increased time resolution
Cg_timeseries = FieldTimeSeries(filepath2,"Cg");
Cs_timeseries = FieldTimeSeries(filepath2,"Cs");
v_timeseries = FieldTimeSeries(filepath2,"v");
b_timeseries = FieldTimeSeries(filepath2,"b");

xc, yc, zc = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC

land = curvedslope.(yc) 

delta = pm.U₀/pm.Ñ

##############
#   TRACER WEIGHTED CALCULATION
##############
ylength = 1000
CGi = interior(Cg_timeseries)[:, 1:ylength,:,:];
Bi = interior(b_timeseries)[:, 1:ylength,:,:];

csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

@info "Calculating first moment..."
# ∫ cb dV
cb = CGi .* Bi;
cbsum = sum(cb, dims=(1,2,3))[1,1,1,:];

# b̄ =  ∫ cb dV / ∫ c dV
b̄_cg = cbsum ./ csum;
Δb̄_cg = b̄_cg[2:end] .- b̄_cg[1:end-1]
neg_indices = findall(Δb̄_cg .< 0)

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
cwtd = @lift b̄_cg[$n]
tim = @lift b_timeseries.times[$n] ./ pm.Tσ;
tcwtd = @lift Point(b_timeseries.times[$n] ./ pm.Tσ, b̄_cg[$n]);

function titlef(n)
    if in(n, neg_indices .+ 1)
        title  = "Internal Wave Breaking, δ=85.7 m, Tσ = " * string(round(b_timeseries.times[n]/pm.Tσ, digits=2) ) * " (b̄↓)"
    else
        title  = "Internal Wave Breaking, δ=85.7 m, Tσ = "*string(round(b_timeseries.times[n]/pm.Tσ, digits=2))
    end
    return title
end

title = @lift(titlef($n))

f = Figure(resolution = (1650, 800),fontsize=26) 

ga = f[1, 1] = GridLayout()
axv = Axis(ga[1, 1], ylabel = "z [m]")
axb = Axis(ga[1, 3])

axcg = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]")
axcs = Axis(ga[2, 3], xlabel = "y [m]")

axwtd = Axis(ga[3, 1:3], xlabel = "WavePeriods [Tσ]", ylabel = "b̄ [ms⁻²]")

axcg.xticks = 500:1000:1500
axcs.xticks = 500:1000:1500
axwtd.xticks = 1:1:10

axcs.yticks = [-250, 0]
axv.yticks = [-250, 0]

limits!(axv, 0, 2000, -500, 0)
limits!(axcg, 0, 2000, -500, 0)
limits!(axcs, 0, 2000, -500, 0)
limits!(axb, 0, 2000, -500, 0)
limits!(axwtd, 0, 11, -3.1e-3, -2.85e-3)

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

wtimes = b_timeseries.times./ pm.Tσ
global wtdp = lines!(axwtd, wtimes, b̄_cg , color = :firebrick2, linewidth = 5)
scatter!(axwtd, tcwtd, markersize = 25, marker = :circle, 
color =:dodgerblue2, strokewidth = 1, strokecolor = :black)

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
    @info msg
    n[] = j
end
