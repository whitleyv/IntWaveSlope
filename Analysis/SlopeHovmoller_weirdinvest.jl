using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

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

name_mp_prefix = "IntWave_mp_" * setname
name_reg_prefix = "IntWave_" * setname
name_reg_prefix = "IntWave_smdt_" * setname 
name_reg_prefix = "IntWave_wraised_" * setname

filepath_mp = path_name * name_mp_prefix * ".jld2"
filepath_reg = path_name * name_reg_prefix * ".jld2"

b_mp_timeseries = FieldTimeSeries(filepath_mp,"b");
b_reg_timeseries = FieldTimeSeries(filepath_reg,"b");

c_reg_timeseries = FieldTimeSeries(filepath_reg, "Cg");
c_mp_timeseries = FieldTimeSeries(filepath_mp, "Cg");

v_mp_timeseries = FieldTimeSeries(filepath_mp,"v");
v_reg_timeseries = FieldTimeSeries(filepath_reg,"v");

xc, yc, zc = nodes(b_mp_timeseries) #CCC
xv, yv, zv = nodes(v_mp_timeseries) #CCC

xc2, yc2, zc2 = nodes(b_reg_timeseries) #CCC


xlocat = 19
b_mpi = interior(b_mp_timeseries)[xlocat,1:334,:,:];
v_mpi = interior(v_mp_timeseries)[xlocat,1:334,:,:];

c_mpi = interior(c_mp_timeseries)[xlocat,1:334,:,:];
c_regi = interior(c_reg_timeseries)[xlocat,1:334,:,1:2:end];

b_regi = interior(b_reg_timeseries)[xlocat,1:334,:,1:2:end];
v_regi = interior(v_reg_timeseries)[xlocat,1:334,:,1:2:end];

lin_land = linslope.(yc) # the actual z value the land hits for each value of y.
delta = pm.U₀/pm.Ñ

##########
#    DATA for HOVMOLLER at z* m above slope line 
#  \
#\  \
# \  \
#  \  \
#   \|z*\
#########
tlength = 161
ylength = 334

function b_slope_line(zval, lin_land, zc, bi, tlength, ylength)

    land_pz = lin_land .+ zval # where we want the hovmoller data to come from

    # for each y value, find the first z value above that number
    gtZ(x) = x > 0
    # num of indices that are too high =  total length of pushed of slope line - number of indices below zero 
    start_index  = length(land_pz) - sum(land_pz .< zc[end]) + 1
    indexLength = length(start_index:ylength)
    depth_indices = zeros(Int, indexLength);

    for (zdx, i) in enumerate(start_index:ylength)
        depth_indices[zdx] = sum(zc .< land_pz[i]) +1
    end

    b_slopelines = zeros(indexLength, tlength)

    # at each y value
    for (zdx, j) in enumerate(start_index:ylength)
        b_slopelines[zdx, :] = bi[j, depth_indices[zdx], :]
    end

    return b_slopelines

end

b_slopelines_10_mp = b_slope_line(10, lin_land, zc, b_mpi, tlength, ylength)
b_slopelines_30_mp = b_slope_line(30, lin_land, zc, b_mpi, tlength, ylength)
b_slopelines_60_mp = b_slope_line(60, lin_land, zc, b_mpi, tlength, ylength)

b_slopelines_10_reg = b_slope_line(10, lin_land, zc2, b_regi, tlength, ylength)
b_slopelines_30_reg = b_slope_line(30, lin_land, zc2, b_regi, tlength, ylength)
b_slopelines_60_reg = b_slope_line(60, lin_land, zc2, b_regi, tlength, ylength)

v_slopelines_10_mp = b_slope_line(10, lin_land, zc, v_mpi, tlength, ylength)
v_slopelines_30_mp = b_slope_line(30, lin_land, zc, v_mpi, tlength, ylength)
v_slopelines_60_mp = b_slope_line(60, lin_land, zc, v_mpi, tlength, ylength)

v_slopelines_10_reg = b_slope_line(10, lin_land, zc, v_regi, tlength, ylength)
v_slopelines_30_reg = b_slope_line(30, lin_land, zc, v_regi, tlength, ylength)
v_slopelines_60_reg = b_slope_line(60, lin_land, zc, v_regi, tlength, ylength)

# number of data points cut off bc how far out from the slope you are:
ycut10 = ylength - size(b_slopelines_10_mp)[1]
ycut30 = ylength - size(b_slopelines_30_mp)[1]
ycut60 = ylength - size(b_slopelines_60_mp)[1]

b_slopelines_10_mp = vcat(zeros(ycut10, tlength), b_slopelines_10_mp)
b_slopelines_30_mp = vcat(zeros(ycut30, tlength), b_slopelines_30_mp)
b_slopelines_60_mp = vcat(zeros(ycut60, tlength), b_slopelines_60_mp)
v_slopelines_10_mp = vcat(zeros(ycut10, tlength), v_slopelines_10_mp)
v_slopelines_30_mp = vcat(zeros(ycut30, tlength), v_slopelines_30_mp)
v_slopelines_60_mp = vcat(zeros(ycut60, tlength), v_slopelines_60_mp)

ycut102 = ylength - size(b_slopelines_10_reg)[1]
ycut302 = ylength - size(b_slopelines_30_reg)[1]
ycut602 = ylength - size(b_slopelines_60_reg)[1]

b_slopelines_10_reg = vcat(zeros(ycut102, tlength), b_slopelines_10_reg)
b_slopelines_30_reg = vcat(zeros(ycut302, tlength), b_slopelines_30_reg)
b_slopelines_60_reg = vcat(zeros(ycut602, tlength), b_slopelines_60_reg)
v_slopelines_10_reg = vcat(zeros(ycut10, tlength), v_slopelines_10_reg)
v_slopelines_30_reg = vcat(zeros(ycut30, tlength), v_slopelines_30_reg)
v_slopelines_60_reg = vcat(zeros(ycut60, tlength), v_slopelines_60_reg)

c_slopelines_10_mp = b_slope_line(10, lin_land, zc, c_mpi, tlength, ylength)
c_slopelines_30_mp = b_slope_line(30, lin_land, zc, c_mpi, tlength, ylength)
c_slopelines_60_mp = b_slope_line(60, lin_land, zc, c_mpi, tlength, ylength)

c_slopelines_10_reg = b_slope_line(10, lin_land, zc, c_regi, tlength, ylength)
c_slopelines_30_reg = b_slope_line(30, lin_land, zc, c_regi, tlength, ylength)
c_slopelines_60_reg = b_slope_line(60, lin_land, zc, c_regi, tlength, ylength)

############
#       PLOTTING Hovmollers
#   ______________ ___
#   | b10 | b30  | | |
#   |_____|______| |_|
#   | e10 |  e30 | | |
#   |_____|______| |_|
#   | v10 |  v30 | | |
#   |_____|______| |_|
#       t
###########

wave_times = b_mp_timeseries.times./ pm.Tσ
ycut = yc[1:ylength]

f = Figure(resolution = (2000, 2000), fontsize=26)
    ga = f[1, 1] = GridLayout()
    gcb = f[1, 2] = GridLayout()

    ax1 = Axis(ga[1, 1], ylabel = "y' [m]")
    ax2 = Axis(ga[1, 2])
    ax2b = Axis(ga[1, 3])
    ax3 = Axis(ga[2, 1], ylabel = "y' [m]")
    ax4 = Axis(ga[2, 2])
    ax4b = Axis(ga[2, 3])
    ax5 = Axis(ga[3, 1], ylabel = "y' [m]")
    ax6 = Axis(ga[3, 2],)
    ax6b = Axis(ga[3, 3])
    ax7 = Axis(ga[4, 1], ylabel = "y' [m]", xlabel = "Tσ")
    ax8 = Axis(ga[4, 2], xlabel = "Tσ")
    ax8b = Axis(ga[4, 3], xlabel = "Tσ")

    limits!(ax1, 0, 8, 0, 1300)
    limits!(ax2, 0, 8, 0, 1300)
    limits!(ax3, 0, 8, 0, 1300)
    limits!(ax4, 0, 8, 0, 1300)
    limits!(ax5, 0, 8, 0, 1300)
    limits!(ax6, 0, 8, 0, 1300)
    limits!(ax6b, 0, 8, 0, 1300)
    limits!(ax4b, 0, 8, 0, 1300)
    limits!(ax2b, 0, 8, 0, 1300)
    limits!(ax7, 0, 8, 0, 1300)
    limits!(ax8, 0, 8, 0, 1300)
    limits!(ax8b, 0, 8, 0, 1300)

    ax1.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
    ax3.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
    ax5.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
    ax7.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )

    ax7.xticks = 1:1:11
    ax8.xticks = 1:1:11
    ax8b.xticks = 1:1:11

    hidedecorations!(ax2)
    hidedecorations!(ax2b)
    hidedecorations!(ax4)
    hidedecorations!(ax4b)
    hidedecorations!(ax6)
    hidedecorations!(ax6b)
    hidexdecorations!(ax1)
    hidexdecorations!(ax3)
    hidexdecorations!(ax5)
    hideydecorations!(ax8)
    hideydecorations!(ax8b)

    reverse!(b_slopelines_10_mp, dims = 1)
    reverse!(b_slopelines_30_mp, dims = 1)
    reverse!(b_slopelines_60_mp, dims = 1)
    reverse!(v_slopelines_10_mp, dims = 1)
    reverse!(v_slopelines_30_mp, dims = 1)
    reverse!(v_slopelines_60_mp, dims = 1)

    reverse!(b_slopelines_10_reg, dims = 1)
    reverse!(b_slopelines_30_reg, dims = 1)
    reverse!(b_slopelines_60_reg, dims = 1)
    reverse!(v_slopelines_10_reg, dims = 1)
    reverse!(v_slopelines_30_reg, dims = 1)
    reverse!(v_slopelines_60_reg, dims = 1)

    reverse!(c_slopelines_10_mp, dims = 1)
    reverse!(c_slopelines_30_mp, dims = 1)
    reverse!(c_slopelines_60_mp, dims = 1)

    reverse!(c_slopelines_10_reg, dims = 1)
    reverse!(c_slopelines_30_reg, dims = 1)
    reverse!(c_slopelines_60_reg, dims = 1)


    # axes 1 and 2 for Buoynacy 
    hm = heatmap!(ax1, wave_times, ycut, b_slopelines_10_reg', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax1, wave_times, ycut, b_slopelines_10_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax1, Point.(0.25, 1200), text = "m < 0, 10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    contour!(ax1, wave_times, ycut, c_slopelines_10_reg', color = (:white, 0.5), linewidth = 3, levels = 1e-2:1:1.5e-2)
    contour!(ax1, wave_times, ycut, c_slopelines_10_reg', color = (:deepskyblue1, 0.5), linewidth = 3, levels = 1e-3:1:1.5e-3)

    heatmap!(ax2, wave_times, ycut, b_slopelines_30_reg', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax2, wave_times, ycut, b_slopelines_30_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax2, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    contour!(ax2, wave_times, ycut, c_slopelines_30_reg', color = (:white, 0.5), linewidth = 3, levels = 1e-2:1:1.5e-2)
    contour!(ax2, wave_times, ycut, c_slopelines_30_reg', color = (:deepskyblue1, 0.5), linewidth = 3, levels = 1e-3:1:1.5e-3)

    heatmap!(ax2b, wave_times, ycut, b_slopelines_60_reg', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax2b, wave_times, ycut, b_slopelines_60_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax2b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    contour!(ax2b, wave_times, ycut, c_slopelines_60_reg', color = (:white, 0.5), linewidth = 3, levels = 1e-2:1:1.5e-2)
    contour!(ax2b, wave_times, ycut, c_slopelines_60_reg', color = (:deepskyblue1, 0.5), linewidth = 3, levels = 1e-3:1:1.5e-3)

    # axes 3 and 4 for Dissipation 
    hm2 = heatmap!(ax3, wave_times, ycut, b_slopelines_10_mp', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax3, wave_times, ycut, b_slopelines_10_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax3, Point.(0.25, 1200), text = "m > 0, 10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    contour!(ax3, wave_times, ycut, c_slopelines_10_mp', color = (:white, 0.5), linewidth = 3, levels = 1e-2:1:1.5e-2)
    contour!(ax3, wave_times, ycut, c_slopelines_10_mp', color = (:deepskyblue1, 0.5), linewidth = 3, levels = 1e-3:1:1.5e-3)

    heatmap!(ax4, wave_times, ycut, b_slopelines_30_mp', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax4, wave_times, ycut, b_slopelines_30_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax4, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    contour!(ax4, wave_times, ycut, c_slopelines_30_mp', color = (:white, 0.5), linewidth = 3, levels = 1e-2:1:1.5e-2)
    contour!(ax4, wave_times, ycut, c_slopelines_30_mp', color = (:deepskyblue1, 0.5), linewidth = 3, levels = 1e-3:1:1.5e-3)

    heatmap!(ax4b, wave_times, ycut, b_slopelines_60_mp', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax4b, wave_times, ycut, b_slopelines_60_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax4b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    contour!(ax4b, wave_times, ycut, c_slopelines_60_mp', color = (:white, 0.5), linewidth = 3, levels = 1e-2:1:1.5e-2)
    contour!(ax4b, wave_times, ycut, c_slopelines_60_mp', color = (:deepskyblue1, 0.5), linewidth = 3, levels = 1e-3:1:1.5e-3)

    # axes 5 and 6 for Velocity 
    hm3 = heatmap!(ax5, wave_times, ycut, v_slopelines_10_reg', colormap = :balance, colorrange = (-0.3, 0.3))
    contour!(ax5, wave_times, ycut, b_slopelines_10_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax5, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax6, wave_times, ycut, v_slopelines_30_reg', colormap = :balance, colorrange = (-0.3, 0.3))
    contour!(ax6, wave_times, ycut, b_slopelines_30_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax6, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax6b, wave_times, ycut, v_slopelines_60_reg', colormap = :balance, colorrange = (-0.3, 0.3))
    contour!(ax6b, wave_times, ycut, b_slopelines_60_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax6b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    # axes 5 and 6 for Velocity 
    heatmap!(ax7, wave_times, ycut, v_slopelines_10_mp', colormap = :balance, colorrange = (-0.3, 0.3))
    contour!(ax7, wave_times, ycut, b_slopelines_10_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax7, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax8, wave_times, ycut, v_slopelines_30_mp', colormap = :balance, colorrange = (-0.3, 0.3))
    contour!(ax8, wave_times, ycut, b_slopelines_30_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax8, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax8b, wave_times, ycut, v_slopelines_60_mp', colormap = :balance, colorrange = (-0.3, 0.3))
    contour!(ax8b, wave_times, ycut, b_slopelines_60_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax8b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    colgap!(ga, 20)
    rowgap!(ga, 20)

    cb = Colorbar(gcb[1:2,1], hm, ticks = (-0.005:0.001:-0.001, ["-5×10⁻³", "-4×10⁻³", "-3×10⁻³", "-2×10⁻³", "-1×10⁻³"] ), size = 35, label = "b")
    cb = Colorbar(gcb[3:4,1], hm3, ticks = (-0.3:.15:0.3, ["-0.3", "-0.15", "0", "0.15", "0.3"] ), size = 35, label = "v")

    colsize!(f.layout, 2, Relative(0.04))

    savename = "bvc_mp_Hovmol_" * setname
    apath  = path_name * "Analysis/"

save(apath * savename * ".png", f)

############
#       PLOTTING Hovmollers
#   ______________ ___
#   | b10 | b30  | | |
#   |_____|______| |_|
#   | b10 |  b30 | | |
#   |_____|______| |_|

#       t
###########

wave_times = b_mp_timeseries.times./ pm.Tσ
ycut = yc[1:ylength]

f = Figure(resolution = (2000, 1500), fontsize=26)
    ga = f[1, 1] = GridLayout()
    gcb = f[1, 2] = GridLayout()

    ax1 = Axis(ga[1, 1], ylabel = "y' [m]")
    ax2 = Axis(ga[1, 2])
    ax2b = Axis(ga[1, 3])
    ax3 = Axis(ga[2, 1], ylabel = "y' [m]", xlabel = "Tσ")
    ax4 = Axis(ga[2, 2], xlabel = "Tσ")
    ax4b = Axis(ga[2, 3], xlabel = "Tσ")

    limits!(ax1, 0, 8, 0, 1300)
    limits!(ax2, 0, 8, 0, 1300)
    limits!(ax3, 0, 8, 0, 1300)
    limits!(ax4, 0, 8, 0, 1300)
    limits!(ax4b, 0, 8, 0, 1300)
    limits!(ax2b, 0, 8, 0, 1300)

    ax1.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
    ax3.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )

    ax3.xticks = 1:1:11
    ax4.xticks = 1:1:11
    ax4b.xticks = 1:1:11

    hidedecorations!(ax2)
    hidedecorations!(ax2b)
    hidexdecorations!(ax1)
    hideydecorations!(ax4)
    hideydecorations!(ax4b)

    reverse!(b_slopelines_10_mp, dims = 1)
    reverse!(b_slopelines_30_mp, dims = 1)
    reverse!(b_slopelines_60_mp, dims = 1)

    reverse!(b_slopelines_10_reg, dims = 1)
    reverse!(b_slopelines_30_reg, dims = 1)
    reverse!(b_slopelines_60_reg, dims = 1)

    # axes 1 and 2 for Buoynacy 
    hm = heatmap!(ax1, wave_times, ycut, b_slopelines_10_reg', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax1, wave_times, ycut, b_slopelines_10_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax1, Point.(0.25, 1200), text = "m > 0, no Ex z, 10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax2, wave_times, ycut, b_slopelines_30_reg', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax2, wave_times, ycut, b_slopelines_30_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax2, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax2b, wave_times, ycut, b_slopelines_60_reg', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax2b, wave_times, ycut, b_slopelines_60_reg', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax2b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    # axes 3 and 4 for Dissipation 
    hm2 = heatmap!(ax3, wave_times, ycut, b_slopelines_10_mp', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax3, wave_times, ycut, b_slopelines_10_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax3, Point.(0.25, 1200), text = "m > 0, 10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax4, wave_times, ycut, b_slopelines_30_mp', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax4, wave_times, ycut, b_slopelines_30_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax4, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax4b, wave_times, ycut, b_slopelines_60_mp', colormap = :thermal, colorrange = (-6e-3, 0))
    contour!(ax4b, wave_times, ycut, b_slopelines_60_mp', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
    text!(ax4b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    colgap!(ga, 20)
    rowgap!(ga, 20)

    cb = Colorbar(gcb[1:2,1], hm, ticks = (-0.005:0.001:-0.001, ["-5×10⁻³", "-4×10⁻³", "-3×10⁻³", "-2×10⁻³", "-1×10⁻³"] ), size = 35, label = "b")

    colsize!(f.layout, 2, Relative(0.04))

    savename = "b_mpnoExZ_Hovmol_" * setname
    apath  = path_name * "Analysis/"

save(apath * savename * ".png", f)

#############
# Gaussian trcaer release
#############

ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)

c_mp_Fi = interior(c_mp_timeseries)[:,1:ylength,:,:];
c_reg_Fi = interior(c_reg_timeseries)[:,1:ylength,:,:];

b_mp_Fi = interior(b_mp_timeseries)[:,1:ylength,:,:];
b_reg_Fi = interior(b_reg_timeseries)[:,1:ylength,:,:];

@info "Tracer Weighted Calculations.."


function calculate_tracer_weight(CGi, Bi)
    
    @info "Calculating Denominator of tracer wtd average..."
    # ∫ c dV
    csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

    @info "Calculating first moment..."
    # ∫ cb dV
    cb = CGi .* Bi;
    cbsum = sum(cb, dims=(1,2,3))[1,1,1,:];

    # b̄ =  ∫ cb dV / ∫ c dV
    c_weighted_bavg = cbsum ./ csum;

    return c_weighted_bavg
end

b̄_Cg_mp = calculate_tracer_weight(c_mp_Fi, b_mp_Fi)
b̄_Cg_reg = calculate_tracer_weight(c_reg_Fi, b_reg_Fi)


######################
#                  b CENTROID PLOT FOR ALL GAUSSIAN
######################
f = Figure(resolution = (1000, 1400), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "Δb̄ [ms⁻²]")
limits!(ax1, 0, 11, -3.4e-4, 3.4e-4)

ax2 = Axis(ga[2, 1], ylabel = "b̄ [ms⁻²]", xlabel = "Wave Periods [Tσ]")
ax2.xticks =  2:2:10
limits!(ax2, 0, 11, -3.1e-3, -2.5e-3)

hidexdecorations!(ax1)

dbmp = lines!(ax1, wave_times, b̄_Cg_mp .- b̄_Cg_mp[1] , color = :dodgerblue2, linewidth = 5)
dbmn = lines!(ax1, wave_times, b̄_Cg_reg .- b̄_Cg_reg[1] , color = :firebrick2, linewidth = 5)

lines!(ax2, wave_times, b̄_Cg_mp, color = :dodgerblue2, linewidth = 5)
lines!(ax2, wave_times, b̄_Cg_reg, color = :firebrick2, linewidth = 5)

axislegend(ax1, [dbmp, dbmn], ["m > 0", "m < 0"], position = :lt)

savename = "cwtd_bbar_mtest_" * setname

save(apath * savename * ".png", f)

################################ INTRUSIONS

# averaging over all times and space vs delta
Cheight_havg_tavg = zeros(Lvals)

################################  THORPE SCALE
# mean and max over all values for each delta
thorpe_hmax_tavg = zeros(Lvals)

# how many times were saved?
zlength = pm.nz
ylength = 425 # roughly out to y = 1700 (is this enough still??)
xlength = pm.nx

# cutting off top and bottom to avoid boundaries:

lastH = 450 #start at z = -50 (used to be -250)
firstH = 50 #end at z = -450
# indices to start and stop 
z_st = round(Int, lastH/2) # start at z = -50
z_en = round(Int, firstH/2) # end at z = -450 (smaller indices = deeper)
y_st = round(Int, ((pm.Lz-lastH)/pm.Tanα)/4) # find corresponding y value on slope when z = 250
y_en = round(Int, 2500/4) # just choosing this y value to include most of dye excursions

zlength_sm = length(z_en:z_st)
ylength_sm = length(y_st:y_en)

cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns
isemp(A) = length(A) < 1
isNemp(A) = length(A) > 1
thresh(z) = z>10^(-6)

# rolling averages width
rolWidy = 20
rolWidz = 3

# number of individual results to log at each (y,t) location
numRes = 25

# total number of data points
Wtlength = wave_info.Wl * wave_info.nTσ
# 4 waves long
cutWtlength = 4*wave_info.Wl
# starts at 7Tσ
W7length = wave_info.Wl * 7  + 1 
# thorpe starts at 3Tσ
W3length = wave_info.Wl * 3  + 1 
# ends at 10.9Tσ
W11length = wave_info.Wl * 11
LtcutWtlength = 8*wave_info.Wl
#for rolling wave avg, can't use the last half of wave if no more waves after that...
hWl = floor(Int64, wave_info.Wl/2)

# - 1/2 wave if doesn't fit since rolling 
if Wtlength > (W11length + hWl)
    # if room to go to end of "last" wave, then go!
    rWtlength = cutWtlength
else 
    rWtlength = cutWtlength - hWl
end


    ### Calculation ###########################################################################

    All_Cheights = zeros(cutWtlength,ylength_sm, numRes)

    @info "Calculating Dye Heights at every time..."
    for (j, i) in enumerate(wave_info.WavePeriods[W7length:W11length])

        # pulling c values within range
        c = c_timeseries[i]
        ci = mean(interior(c)[:,:, z_en:z_st], dims =1)[1,:,:]

        # averaging rolling window
        cmi = zeros(ylength_sm, zlength_sm)

        # yi = y index in cut range
        # y_md y index from total range
        for (yi, y_md) in enumerate(y_st:y_en)
            cmi[yi,:] = mean(ci[y_md-rolWidy:rolWidy+y_md, :], dims = 1)
        end

        Δc = cmi[:,1:end-1]-cmi[:,2:end]
        dcdz = Δc /-2

        # finding all "roots"
        for (yi, y_md) in enumerate(y_st:y_en)
            roots = []
            for k = 2:size(dcdz)[2]-2
                if dcdz[yi,k] <= 0 && dcdz[yi,k+1] >0
                    roots = [roots;k]
                end
            end

            if isNemp(roots)

                top_cond(z) = z>=curvedslope(yc[y_md])
                # first index above the slope, fy : end in the fluid
                fy = findfirst(top_cond, zc)
                # cutting off any roots that are within 4 indices of bottom (6m)
                cutbot = sum( roots .< fy + 3)
                roots = roots[cutbot+1:end]

                rt_widths = []

                w=1

                for p = 2:length(roots)
                    maxc = maximum(cmi[yi,roots[p-1]:roots[p]])
                    minr = maxc*.5 
                    if maxc > 10^(-4) && cmi[yi, roots[p]] <= minr && cmi[yi, roots[p-1]] <= minr
                        rng = roots[p-1]:roots[p]
                        fr = findfirst(thresh, cmi[yi,rng])
                        lr = findfirst(thresh, reverse(cmi[yi,rng]))
                        newr1 = roots[p-1] + fr -1
                        newr2 = roots[p] - lr
                        z1 = zc[newr1 + z_en]
                        z2 = zc[newr2 + z_en]
                        new_wid = z2-z1
                        rt_widths=[rt_widths; new_wid]
                        All_Cheights[j, yi, w] = new_wid
                        w +=1
                    end
                end
            end    
            
        end

    end

    @info "Calculating Thorpe Displacements..."
    # calculating RMS value over only the over turn at every time and y position
    All_Toverturns = zeros(LtcutWtlength, ylength-1, numRes)

    for (ni, i) in enumerate(wave_info.WavePeriods[W3length:W11length])

        Displace_Lt = zeros(ylength-1, zlength)
        #@info "time $i/$tlength..."
        b = b_timeseries[i]
        bi = interior(b)[:, 1:ylength, :]
    
        # average in x you get in terms of y and z
        b_xavg = mean(bi, dims = 1)[1,:,:]
    
        for (j,y) in enumerate(yb[2:ylength])

            top_cond(z) = z>=curvedslope(y)
            # first index above the slope, fy : end in the fluid
            fy = findfirst(top_cond, zb)
            # concatenate vertical profile and z values
            bz = cat(b_xavg[j+1,:], zb, dims=2)
            # cut arrays to only include relevant z values after topo removes
            bz_cut = bz[fy:zlength,:]
    
            # sort both rows by buoyancy values
            bz_std = bz_cut[sortperm(bz_cut[:,1]),:]
    
            # z in bz_cut is actually their new loc, while bz_std stored old val. at idx
            # taking the new poisition of the parcel - old
            # + displacement means a downward move to stabilize
            #              d_o[k] =  z new[k] .- z old[k]  for k = slope:sfc
            Displace_Lt[j,fy:zlength] =  bz_cut[:,2] .- bz_std[:,2]
    
            # now take the sum of these displacement from the sfc to each depth
            Sum_Disp = zeros(zlength)
            for k = fy:zlength
                Sum_Disp[k] = sum(Displace_Lt[j,k:end])
            end

            # find all the values above threshold, defining overturns
            O_idxs = findall(cond_gt4, Sum_Disp)
            
            if isNemp(O_idxs)
                # find all the jumps in values 
                sO_idxs = findall(cond_gt1, O_idxs[2:end] - O_idxs[1:end-1])
                # find idx for 1st value that changed, (beginning of overturn)
                b_over_idxs = isNemp(sO_idxs) ? [O_idxs[1];O_idxs[sO_idxs .+ 1]] : O_idxs[1]
                e_over_idxs =  isNemp(sO_idxs) ? [O_idxs[sO_idxs]; O_idxs[end]] :  O_idxs[end]
    
                # if nothing then just set it to 1 so that it takes an index at all
                Overturns_rmsL = ones(length(b_over_idxs))
                for o in 1:length(b_over_idxs)
                    bo = b_over_idxs[o]
                    eo = e_over_idxs[o]
                    overturn = Displace_Lt[j,bo:eo]
                    overturnL = length(overturn)
                    #Overturns_rmsL[o] = sqrt(sum((overturn).^2)/overturnL)
                    All_Toverturns[ni, j, o] = sqrt(sum((overturn).^2)/overturnL)

                end
            end

        end
        
    end
