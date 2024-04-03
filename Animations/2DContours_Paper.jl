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

path_name = "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"

apath = path_name * "Analysis/"

setname = "U250N100Lz100g100"

@info "Loading in parameters..."

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
m = -π/pm.Lz,
l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
Tf = 2*π/pm.f, 
Tσ = 2*π/pm.σ))
zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
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
name1_prefix = "IntWave_mp_noS_" * setname 

filepath1 = path_name * name1_prefix * ".jld2"

v_timeseries = FieldTimeSeries(filepath1,"v");
b_timeseries = FieldTimeSeries(filepath1,"b");
Cs_timeseries = FieldTimeSeries(filepath1,"Cs");

xc, yc, zc = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC

land = curvedslope.(yc) 

delta = pm.U₀/pm.Ñ

##############
#   CAIRO PLOT MOVIE
##############
 
n = Observable(1)
xlocat = 19
tlength = 161
zlength = 250
ylength = 1000

n = Observable(1)
xlocat = 19

v = @lift interior(v_timeseries[$n],xlocat, :, :,);
b = @lift interior(b_timeseries[$n], xlocat, :, :);

title = @lift @sprintf( "Horizontal Velocity (v), t = %.2f hrs, Tσ = %.2f", b_timeseries.times[$n]/3600, b_timeseries.times[$n]/pm.Tσ)

f = Figure(resolution = (2000, 1000),fontsize=35) 

    ga = f[1, 1] = GridLayout()
    axv = Axis(ga[1, 1], ylabel = "z [m]")

    axv.xticks = 500:500:3000
    axv.yticks = [-500, -250, 0]

    limits!(axv, 0, 3000, -500, 0)

    land_pdel = (linslope.(yc) .+ delta)[1:382]
    blevels  = range( -0.006,0; length = 75)
    global hmv = heatmap!(axv, yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
    contour!(axv, yc, zc, b, color = :gray30, linewidth = 0.5, levels=blevels, alpha = 0.5)
    band!(axv, yc, land, -500, color=:black)
    lines!(axv, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)

    Label(ga[1, 1:2, Top()], title,
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :center)

    cb1 = Colorbar(ga[1,2], hmv, ticks = (-0.2:.1:0.2), size =35, label = "v [ms⁻¹]")

    colsize!(ga, 2, Relative(0.05))

savename = "presentation_vcontours_mp_" * setname

frames = 1:length(b_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=6) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end

n = Observable(1)

title = @lift @sprintf( "Tracer Concentration, t = %.2f hrs, Tσ = %.2f", b_timeseries.times[$n]/3600, b_timeseries.times[$n]/pm.Tσ)

b = @lift interior(b_timeseries[$n], xlocat, :, :);
cs = @lift log10.(clamp.(interior(Cs_timeseries[$n], xlocat, :, :), 1e-8,1));

f = Figure(resolution = (2000, 1000),fontsize=35) 

    ga = f[1, 1] = GridLayout()
    axv = Axis(ga[1, 1], ylabel = "z [m]")

    axv.xticks = 500:500:3000
    axv.yticks = [-500, -250, 0]

    limits!(axv, 0, 3000, -500, 0)

    land_pdel = (linslope.(yc) .+ delta)[1:382]
    blevels  = range( -0.006,0; length = 75)
    global hmv = heatmap!(axv, yc, zc, cs, colormap = :thermal, colorrange = (-4, 0))
    contour!(axv, yc, zc, b, color = :gray80, linewidth = 0.5, levels=blevels, alpha = 0.5)
    band!(axv, yc, land, -500, color=:black)
    lines!(axv, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)

    Label(ga[1, 1:2, Top()], title,
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :center)

    cb1 = Colorbar(ga[1,2], hmv, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size =35, label = "Tracer, c")

    colsize!(ga, 2, Relative(0.05))

savename = "presentation_Cscontours_mp_" * setname

frames = 1:length(b_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=6) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end

n = Observable(1)

titlev = @lift @sprintf( "Horizontal Velocity (v), t = %.2f hrs, Tσ = %.2f", b_timeseries.times[$n]/3600, b_timeseries.times[$n]/pm.Tσ)
titlec = @sprintf( "Tracer Concentration")

b = @lift interior(b_timeseries[$n], xlocat, :, :);
cs = @lift log10.(clamp.(interior(Cs_timeseries[$n], xlocat, :, :), 1e-8,1));
v = @lift interior(v_timeseries[$n],xlocat, :, :,);

f = Figure(resolution = (2000, 1500),fontsize=35) 

    ga = f[1, 1] = GridLayout()
    axv = Axis(ga[1, 1], ylabel = "z [m]")
    axc = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]")

    axv.xticks = 500:500:3000
    axv.yticks = [-500, -250, 0]
    axc.xticks = 500:500:3000
    axc.yticks = [-500, -250, 0]

    limits!(axv, 0, 3000, -500, 0)
    limits!(axc, 0, 3000, -500, 0)

    hidexdecorations!(axv, grid = false)

    land_pdel = (linslope.(yc) .+ delta)[1:382]
    blevels  = range( -0.006,0; length = 75)

    global hmc = heatmap!(axc, yc, zc, cs, colormap = :thermal, colorrange = (-4, 0))
    contour!(axc, yc, zc, b, color = :gray80, linewidth = 0.5, levels=blevels, alpha = 0.5)
    band!(axc, yc, land, -500, color=:black)
    lines!(axc, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)

    global hmv = heatmap!(axv, yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
    contour!(axv, yc, zc, b, color = :gray30, linewidth = 0.5, levels=blevels, alpha = 0.5)
    band!(axv, yc, land, -500, color=:black)
    lines!(axv, yc[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)

    Label(ga[1, 1:2, Top()], titlev,
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :center)

    Label(ga[2, 1:2, Top()], titlec,
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :center)

    cb1 = Colorbar(ga[1,2], hmv, ticks = (-0.2:.1:0.2), size =35, label = "v [ms⁻¹]")

    cb2 = Colorbar(ga[2,2], hmc, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size =35, label = "Tracer, c")

    colsize!(ga, 2, Relative(0.05))

savename = "presentation_vCscontours_mp_" * setname

frames = 1:length(b_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=6) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end

