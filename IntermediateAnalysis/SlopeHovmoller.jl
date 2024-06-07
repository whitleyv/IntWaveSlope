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

name1_prefix = "IntWave_" * setname
name1_prefix = "IntWave_smdt_U300N100Lz100g100"
filepath1 = path_name * name1_prefix * ".jld2"
b_timeseries = FieldTimeSeries(filepath1,"b");
c_timeseries = FieldTimeSeries(filepath1,"Cg");
e_timeseries = FieldTimeSeries(filepath1,"ϵ");
v_timeseries = FieldTimeSeries(filepath1,"v");

filepath2 = path_name * "cgl_" * setname * ".jld2"
f3 = jldopen(filepath2);
CGli = f3["Cgli"];
spinlength = size(CGli)[4]
spinuplength = tlength - spinlength;

xc, yc, zc = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CCC

xlocat = 19
bi = interior(b_timeseries)[xlocat,1:334,:,:];
ci = interior(c_timeseries)[xlocat,1:334,:,:];
ei = interior(e_timeseries)[xlocat,1:334,:,:];
vi = interior(v_timeseries)[xlocat,1:334,:,:];
cLi =  CGli[xlocat, 1:334, :,:];

lin_land = linslope.(yc) # the actual z value the land hits for each value of y.
delta = pm.U₀/pm.Ñ

##########
#    DATA for (z',t) HOVMOLLER at y* = 703 m  
#  
#\ |
# \|z'  
#  \  
#   \
# y* \
#########
tlength = 161
ystart = 250/pm.Tanα
yidx_start = sum(yc .< ystart) + 1
yidx_end = 350
zidx_start = 125
yidx_length = yidx_end - yidx_start + 1

zidx_end = zidx_start + 125
zidx_length = zidx_end - zidx_start + 1

Idx_z = zidx_start:zidx_end

Idx_y = zeros(Int, zidx_length)

# z = -250 + (1/ta)(y - ystart)
# ta*(z + 250) + ystart = y
perp_slopelinez(z) = ystart + pm.Tanα * (z + 250)

# for the y value greater than y* value
for (j, z) in enumerate(zc[Idx_z])
    @inline abovept(t) = t > perp_slopelinez(z)

    Idx_y[j] = findfirst(abovept, yc)
end

b_perpline = zeros(zidx_length, tlength)
for k = 1:zidx_length
    b_perpline[k,:] = bi[Idx_y[k],Idx_z[k], :]
end

############
#       PLOTTING Hovmollers
#   ______________ ___
#   |  10 |  30  | | |
#   |_____|______| | |
# y |  60 |   85 | | |
#   |_____|______| |_|
#       t
###########
wave_times = b_timeseries.times./ pm.Tσ
zcut = zc[Idx_z]

contline_b_start = (-210) .* pm.Ñ^2
contline_b_end = -250 .* pm.Ñ^2
Δb = contline_b_start - contline_b_end

f = Figure(resolution = (1100, 600), fontsize=26)
ga = f[1, 1] = GridLayout()
gcb = f[1, 2] = GridLayout()

ax3 = Axis(ga[1, 1], ylabel = "z' [m]", xlabel = "Tσ")

limits!(ax3, 0, 8, -250, -150)

ax3.yticks = -250:25:0
ax3.xticks = 0:1:11

hm = heatmap!(ax3, wave_times, zcut, b_perpline', colormap = :thermal, colorrange = (-4e-3, -1e-3))
contour!(ax3, wave_times, zcut, b_perpline', color = :black, linewidth = 3, levels = -0.005:0.001:0, alpha = 0.5)
contour!(ax3, wave_times, zcut, b_perpline', color = :black, linewidth = 0.5, levels = -0.004:0.00025:-0.001, alpha = 0.5)
contour!(ax3, wave_times, zcut, b_perpline', color = (:green, 0.25), linewidth = 1, levels = contline_b_end:Δb/100:contline_b_start)

cb = Colorbar(gcb[1,1], hm, ticks = (-0.005:0.001:0, ["-5×10⁻³", "-4×10⁻³", "-3×10⁻³", "-2×10⁻³", "-1×10⁻³", "0"] ), size = 35)

colsize!(f.layout, 2, Relative(0.04))

savename = "bzHovmol_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f)


##########
#    DATA for HOVMOLLER at z* m above slope line 
#  \
#\  \
# \  \
#  \  \
#   \|z*\
#########
tlength = 161
tlength = 321
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

function b_slope_line(zval, lin_land, zc, bi, tlength, ylength, spin_start)

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
        b_slopelines[zdx, spin_start:end] = bi[j, depth_indices[zdx], :]
    end

    return b_slopelines

end

b_slopelines_10 = b_slope_line(10, lin_land, zc, bi, tlength, ylength)
b_slopelines_30 = b_slope_line(30, lin_land, zc, bi, tlength, ylength)
b_slopelines_60 = b_slope_line(60, lin_land, zc, bi, tlength, ylength)
b_slopelines_85 = b_slope_line(85, lin_land, zc, bi, tlength, ylength)

v_slopelines_10 = b_slope_line(10, lin_land, zc, vi, tlength, ylength)
v_slopelines_30 = b_slope_line(30, lin_land, zc, vi, tlength, ylength)
v_slopelines_60 = b_slope_line(60, lin_land, zc, vi, tlength, ylength)
v_slopelines_85 = b_slope_line(85, lin_land, zc, vi, tlength, ylength)

c_slopelines_10 = b_slope_line(10, lin_land, zc, ci, tlength, ylength)
c_slopelines_30 = b_slope_line(30, lin_land, zc, ci, tlength, ylength)
c_slopelines_60 = b_slope_line(60, lin_land, zc, ci, tlength, ylength)
c_slopelines_85 = b_slope_line(85, lin_land, zc, ci, tlength, ylength)

cL_slopelines_10 = b_slope_line(10, lin_land, zc, cLi, tlength, ylength, spinuplength+1)
cL_slopelines_30 = b_slope_line(30, lin_land, zc, cLi, tlength, ylength, spinuplength+1)
cL_slopelines_60 = b_slope_line(60, lin_land, zc, cLi, tlength, ylength, spinuplength+1)
cL_slopelines_85 = b_slope_line(85, lin_land, zc, cLi, tlength, ylength, spinuplength+1)

e_slopelines_10 = b_slope_line(10, lin_land, zc, ei, tlength, ylength)
e_slopelines_30 = b_slope_line(30, lin_land, zc, ei, tlength, ylength)
e_slopelines_60 = b_slope_line(60, lin_land, zc, ei, tlength, ylength)
e_slopelines_85 = b_slope_line(85, lin_land, zc, ei, tlength, ylength)

b_slopelines_10_full = reverse(b_slopelines_10, dims = 1)
b_slopelines_30_full = reverse(b_slopelines_30, dims = 1)
b_slopelines_60_full = reverse(b_slopelines_60, dims = 1)
b_slopelines_85_full = reverse(b_slopelines_85, dims = 1)

c_slopelines_10_fullNL = reverse(c_slopelines_10, dims = 1)
c_slopelines_30_fullNL = reverse(c_slopelines_30, dims = 1)
c_slopelines_60_fullNL = reverse(c_slopelines_60, dims = 1)
c_slopelines_85_fullNL = reverse(c_slopelines_85, dims = 1)

cL_slopelines_10_fullNL = reverse(cL_slopelines_10, dims = 1)
cL_slopelines_30_fullNL = reverse(cL_slopelines_30, dims = 1)
cL_slopelines_60_fullNL = reverse(cL_slopelines_60, dims = 1)
cL_slopelines_85_fullNL = reverse(cL_slopelines_85, dims = 1)

v_slopelines_10_full = reverse(v_slopelines_10, dims = 1)
v_slopelines_30_full = reverse(v_slopelines_30, dims = 1)
v_slopelines_60_full = reverse(v_slopelines_60, dims = 1)
v_slopelines_85_full = reverse(v_slopelines_85, dims = 1)

e_slopelines_10_fullNL = reverse(e_slopelines_10, dims = 1)
e_slopelines_30_fullNL = reverse(e_slopelines_30, dims = 1)
e_slopelines_60_fullNL = reverse(e_slopelines_60, dims = 1)
e_slopelines_85_fullNL = reverse(e_slopelines_85, dims = 1)

e_slopelines_10_full = log10.(clamp.(e_slopelines_10_fullNL, 1e-16, 1e-2))
e_slopelines_30_full = log10.(clamp.(e_slopelines_30_fullNL, 1e-16, 1e-2))
e_slopelines_60_full = log10.(clamp.(e_slopelines_60_fullNL, 1e-16, 1e-2))
e_slopelines_85_full = log10.(clamp.(e_slopelines_85_fullNL, 1e-16, 1e-2))

c_slopelines_10_full = log10.(clamp.(c_slopelines_10_fullNL, 1e-8, 1))
c_slopelines_30_full = log10.(clamp.(c_slopelines_30_fullNL, 1e-8, 1))
c_slopelines_60_full = log10.(clamp.(c_slopelines_60_fullNL, 1e-8, 1))
c_slopelines_85_full = log10.(clamp.(c_slopelines_85_fullNL, 1e-8, 1))

cL_slopelines_10_full = log10.(clamp.(cL_slopelines_10_fullNL, 1e-8, 1))
cL_slopelines_30_full = log10.(clamp.(cL_slopelines_30_fullNL, 1e-8, 1))
cL_slopelines_60_full = log10.(clamp.(cL_slopelines_60_fullNL, 1e-8, 1))
cL_slopelines_85_full = log10.(clamp.(cL_slopelines_85_fullNL, 1e-8, 1))

# average rolling in time?
roll_window = 1
c_slopelines_10_full_rA = zeros(ylength, tlength - (2*roll_window))
c_slopelines_30_full_rA = zeros(ylength, tlength - (2*roll_window))
c_slopelines_60_full_rA = zeros(ylength, tlength - (2*roll_window))
c_slopelines_85_full_rA = zeros(ylength, tlength - (2*roll_window))

for i = (roll_window+1):(tlength - roll_window)
    c_slopelines_10_full_rA[:,i-roll_window] = mean(c_slopelines_10_full[:,i-roll_window:i+roll_window], dims = 2)
    c_slopelines_30_full_rA[:,i-roll_window] = mean(c_slopelines_30_full[:,i-roll_window:i+roll_window], dims = 2)
    c_slopelines_60_full_rA[:,i-roll_window] = mean(c_slopelines_60_full[:,i-roll_window:i+roll_window], dims = 2)
    c_slopelines_85_full_rA[:,i-roll_window] = mean(c_slopelines_85_full[:,i-roll_window:i+roll_window], dims = 2)

end

roll_window = 1
cL_slopelines_10_full_rA = zeros(ylength, tlength - (2*roll_window))
cL_slopelines_30_full_rA = zeros(ylength, tlength - (2*roll_window))
cL_slopelines_60_full_rA = zeros(ylength, tlength - (2*roll_window))
cL_slopelines_85_full_rA = zeros(ylength, tlength - (2*roll_window))

for i = (roll_window+1):(tlength - roll_window)
    cL_slopelines_10_full_rA[:,i-roll_window] = mean(cL_slopelines_10_full[:,i-roll_window:i+roll_window], dims = 2)
    cL_slopelines_30_full_rA[:,i-roll_window] = mean(cL_slopelines_30_full[:,i-roll_window:i+roll_window], dims = 2)
    cL_slopelines_60_full_rA[:,i-roll_window] = mean(cL_slopelines_60_full[:,i-roll_window:i+roll_window], dims = 2)
    cL_slopelines_85_full_rA[:,i-roll_window] = mean(cL_slopelines_85_full[:,i-roll_window:i+roll_window], dims = 2)

end

############
#       PLOTTING Hovmollers
#   ______________ ___
#   |  10 |  30  | | |
#   |_____|______| | |
# y |  60 |   85 | | |
#   |_____|______| |_|
#       t
###########
wave_times1 = b_timeseries.times./ pm.Tσ
ycut = yc[1:ylength]

f = Figure(resolution = (1600, 1000), fontsize=26)
ga = f[1, 1] = GridLayout()
gcb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "y' [m]")
ax2 = Axis(ga[1, 2])
ax3 = Axis(ga[2, 1], ylabel = "y' [m]", xlabel = "Tσ")
ax4 = Axis(ga[2, 2], xlabel = "Tσ")

limits!(ax1, 0, 8, 0, 1300)
limits!(ax2, 0, 8, 0, 1300)
limits!(ax3, 0, 8, 0, 1300)
limits!(ax4, 0, 8, 0, 1300)

ax1.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
ax3.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )

ax3.xticks = 0:1:11
ax4.xticks = 0:1:11

hidedecorations!(ax2)
hidexdecorations!(ax1)
hideydecorations!(ax4)

tlength1 = 161

hm = heatmap!(ax1, wave_times1, ycut, b_slopelines_10_full', colormap = :thermal, colorrange = (-6e-3, 1e-3))
contour!(ax1, wave_times1, ycut, b_slopelines_10_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax1, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
contour!(ax1, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_10_full_rA', color = (:white, 0.5), linewidth = 3, levels = -2:1:-1.1)
contour!(ax1, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_10_full_rA', color = (:deepskyblue1, 0.5), linewidth = 3, levels = -3:1:-2.1)
#contour!(ax1, wave_times, ycut, b_slopelines_10_full', color = (:green, 0.25), linewidth = 1, levels = contline_b_end:Δb/100:contline_b_start)

heatmap!(ax2, wave_times1, ycut, b_slopelines_30_full', colormap = :thermal, colorrange = (-6e-3, 1e-3))
contour!(ax2, wave_times1, ycut, b_slopelines_30_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax2, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
contour!(ax2, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_30_full_rA', color = (:white, 0.5), linewidth = 3, levels = -2:1:-1.1)
contour!(ax2, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_30_full_rA', color = (:deepskyblue1, 0.5), linewidth = 3, levels = -3:1:-2.1)
#contour!(ax2, wave_times, ycut, b_slopelines_30_full', color = (:green, 0.25), linewidth = 1, levels = contline_b_end:Δb/100:contline_b_start)

heatmap!(ax3, wave_times1, ycut, b_slopelines_60_full', colormap = :thermal, colorrange = (-6e-3, 1e-3))
contour!(ax3, wave_times1, ycut, b_slopelines_60_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax3, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
contour!(ax3, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_60_full_rA', color = (:white, 0.5), linewidth = 3, levels = -2:1:-1.1)
contour!(ax3, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_60_full_rA', color = (:deepskyblue1, 0.5), linewidth = 3, levels = -3:1:-2.1)
#contour!(ax3, wave_times, ycut, b_slopelines_60_full', color = (:green, 0.25), linewidth = 1, levels = contline_b_end:Δb/100:contline_b_start)

heatmap!(ax4, wave_times1, ycut, b_slopelines_85_full', colormap = :thermal, colorrange = (-6e-3, 1e-3))
contour!(ax4, wave_times1, ycut, b_slopelines_85_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax4, Point.(0.25, 1200), text = "85 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
contour!(ax4, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_85_full_rA', color = (:white, 0.5), linewidth = 3, levels = -2:1:-1.1)
contour!(ax4, wave_times[roll_window+1:tlength1-roll_window], ycut, cL_slopelines_85_full_rA', color = (:deepskyblue1, 0.5), linewidth = 3, levels = -3:1:-2.1)
#contour!(ax4, wave_times, ycut, b_slopelines_85_full', color = (:green, 0.25), linewidth = 1, levels = contline_b_end:Δb/100:contline_b_start)

colgap!(ga, 20)
rowgap!(ga, 20)

cb = Colorbar(gcb[1:2,1], hm, ticks = (-0.005:0.001:0, ["-5×10⁻³", "-4×10⁻³", "-3×10⁻³", "-2×10⁻³", "-1×10⁻³", "0"] ), size = 35)

colsize!(f.layout, 2, Relative(0.04))

savename = "bcdtdyeHovmol_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f)


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

wave_times = b_timeseries.times./ pm.Tσ
ycut = yc[1:ylength]

f = Figure(resolution = (2000, 1600), fontsize=26)
ga = f[1, 1] = GridLayout()
gcb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "y' [m]")
ax2 = Axis(ga[1, 2])
ax2b = Axis(ga[1, 3])
ax3 = Axis(ga[2, 1], ylabel = "y' [m]")
ax4 = Axis(ga[2, 2])
ax4b = Axis(ga[2, 3])
ax5 = Axis(ga[3, 1], ylabel = "y' [m]", xlabel = "Tσ")
ax6 = Axis(ga[3, 2], xlabel = "Tσ")
ax6b = Axis(ga[3, 3], xlabel = "Tσ")

limits!(ax1, 0, 8, 0, 1300)
limits!(ax2, 0, 8, 0, 1300)
limits!(ax3, 0, 8, 0, 1300)
limits!(ax4, 0, 8, 0, 1300)
limits!(ax5, 0, 8, 0, 1300)
limits!(ax6, 0, 8, 0, 1300)
limits!(ax6b, 0, 8, 0, 1300)
limits!(ax4b, 0, 8, 0, 1300)
limits!(ax2b, 0, 8, 0, 1300)

ax1.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
ax3.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
ax5.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )

ax5.xticks = 1:1:11
ax6.xticks = 1:1:11
ax6b.xticks = 1:1:11

hidedecorations!(ax2)
hidedecorations!(ax2b)
hidedecorations!(ax4)
hidedecorations!(ax4b)
hidexdecorations!(ax1)
hidexdecorations!(ax3)
hideydecorations!(ax6)
hideydecorations!(ax6b)

# axes 1 and 2 for Buoynacy 
hm = heatmap!(ax1, wave_times, ycut, b_slopelines_10_full', colormap = :thermal, colorrange = (-6e-3, 0))
contour!(ax1, wave_times, ycut, b_slopelines_10_full', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax1, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
#contour!(ax1, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_10_full_rA', color = (:white, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax1, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_10_full_rA', color = (:deepskyblue1, 0.5), linewidth = 5, levels = -3:1:-2.1)

heatmap!(ax2, wave_times, ycut, b_slopelines_30_full', colormap = :thermal, colorrange = (-6e-3, 0))
contour!(ax2, wave_times, ycut, b_slopelines_30_full', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax2, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
#contour!(ax2, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_30_full_rA', color = (:white, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax2, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_30_full_rA', color = (:deepskyblue1, 0.5), linewidth = 5, levels = -3:1:-2.1)

heatmap!(ax2b, wave_times, ycut, b_slopelines_60_full', colormap = :thermal, colorrange = (-6e-3, 0))
contour!(ax2b, wave_times, ycut, b_slopelines_60_full', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax2b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
#contour!(ax2b, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_60_full_rA', color = (:white, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax2b, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_60_full_rA', color = (:deepskyblue1, 0.5), linewidth = 5, levels = -3:1:-2.1)

# axes 3 and 4 for Dissipation 
hm2 = heatmap!(ax3, wave_times, ycut, e_slopelines_10_fullNL', colormap = :thermal, colorrange = (0, 1e-6))
contour!(ax3, wave_times, ycut, b_slopelines_10_full', color = :gray, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax3, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :white, font = :bold, fontsize = 26)
#contour!(ax3, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_10_full_rA', color = (:white, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax3, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_10_full_rA', color = (:deepskyblue1, 0.5), linewidth = 5, levels = -3:1:-2.1)

heatmap!(ax4, wave_times, ycut, e_slopelines_30_fullNL', colormap = :thermal, colorrange = (0, 1e-6))
contour!(ax4, wave_times, ycut, b_slopelines_30_full', color = :gray, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax4, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :white, font = :bold, fontsize = 26)
#contour!(ax4, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_30_full_rA', color = (:white, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax4, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_30_full_rA', color = (:deepskyblue1, 0.5), linewidth = 5, levels = -3:1:-2.1)

heatmap!(ax4b, wave_times, ycut, e_slopelines_60_fullNL', colormap = :thermal, colorrange = (0, 1e-6))
contour!(ax4b, wave_times, ycut, b_slopelines_60_full', color = :gray, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax4b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :white, font = :bold, fontsize = 26)
#contour!(ax4b, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_60_full_rA', color = (:white, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax4b, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_60_full_rA', color = (:deepskyblue1, 0.5), linewidth = 5, levels = -3:1:-2.1)

# axes 5 and 6 for Velocity 
hm3 = heatmap!(ax5, wave_times, ycut, v_slopelines_10_full', colormap = :balance, colorrange = (-0.3, 0.3))
contour!(ax5, wave_times, ycut, b_slopelines_10_full', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax5, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
#contour!(ax5, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_10_full_rA', color = (:gray30, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax5, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_10_full_rA', color = (:navy, 0.5), linewidth = 5, levels = -3:1:-2.1)

heatmap!(ax6, wave_times, ycut, v_slopelines_30_full', colormap = :balance, colorrange = (-0.3, 0.3))
contour!(ax6, wave_times, ycut, b_slopelines_30_full', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax6, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
#contour!(ax6, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_30_full_rA', color = (:gray30, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax6, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_30_full_rA', color = (:navy, 0.5), linewidth = 5, levels = -3:1:-2.1)

heatmap!(ax6b, wave_times, ycut, v_slopelines_60_full', colormap = :balance, colorrange = (-0.3, 0.3))
contour!(ax6b, wave_times, ycut, b_slopelines_60_full', color = :black, linewidth = 2, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax6b, Point.(0.25, 1200), text = "60 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
#contour!(ax6b, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_60_full_rA', color = (:gray30, 0.5), linewidth = 5, levels = -2:1:-1.1)
#contour!(ax6b, wave_times[roll_window+1:tlength-roll_window], ycut, c_slopelines_60_full_rA', color = (:navy, 0.5), linewidth = 5, levels = -3:1:-2.1)

colgap!(ga, 20)
rowgap!(ga, 20)

cb = Colorbar(gcb[1,1], hm, ticks = (-0.005:0.001:-0.001, ["-5×10⁻³", "-4×10⁻³", "-3×10⁻³", "-2×10⁻³", "-1×10⁻³"] ), size = 35, label = "b")
cb = Colorbar(gcb[2,1], hm2, ticks = (1e-7:2e-7:1e-6, ["1×10⁻⁷", "3×10⁻⁷", "5×10⁻⁷", "7×10⁻⁷", "9×10⁻⁹"] ), size = 35, label = "ε")
cb = Colorbar(gcb[3,1], hm3, ticks = (-0.3:.15:0.3, ["-0.3", "-0.15", "0", "0.15", "0.3"] ), size = 35, label = "v")

colsize!(f.layout, 2, Relative(0.04))

savename = "bevHovmol_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f)

############
#       PLOTTING Hovmollers
#   ______________ ___
#   | b10 | b30  | | |
#   |_____|______| | |
# y | e10 |  e30 | | |
#   |_____|______| |_|
#       t
###########

wave_times = b_timeseries.times./ pm.Tσ
ycut = yc[1:ylength]

f = Figure(resolution = (1600, 1000), fontsize=26)
ga = f[1, 1] = GridLayout()
gcb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "y' [m]")
ax2 = Axis(ga[1, 2])
ax3 = Axis(ga[2, 1], ylabel = "y' [m]", xlabel = "Tσ")
ax4 = Axis(ga[2, 2], xlabel = "Tσ")

limits!(ax1, 0, 8, 0, 1300)
limits!(ax2, 0, 8, 0, 1300)
limits!(ax3, 0, 8, 0, 1300)
limits!(ax4, 0, 8, 0, 1300)

ax1.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
ax3.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )

ax3.xticks = 1:1:11
ax4.xticks = 1:1:11

hidedecorations!(ax2)
hidexdecorations!(ax1)
hideydecorations!(ax4)

hm = heatmap!(ax1, wave_times, ycut, b_slopelines_10_full', colormap = :thermal, colorrange = (-6e-3, 0))
contour!(ax1, wave_times, ycut, b_slopelines_10', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax1, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

heatmap!(ax2, wave_times, ycut, b_slopelines_30_full', colormap = :thermal, colorrange = (-6e-3, 0))
contour!(ax2, wave_times, ycut, b_slopelines_30_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax2, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

hm2 = heatmap!(ax3, wave_times, ycut, e_slopelines_10_fullNL', colormap = :thermal, colorrange = (0, 1e-6))
contour!(ax3, wave_times, ycut, b_slopelines_10_full', color = :gray, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax3, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :white, font = :bold, fontsize = 26)

heatmap!(ax4, wave_times, ycut, e_slopelines_30_fullNL', colormap = :thermal, colorrange = (0, 1e-6))
contour!(ax4, wave_times, ycut, b_slopelines_30_full', color = :gray, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax4, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :white, font = :bold, fontsize = 26)

colgap!(ga, 20)
rowgap!(ga, 20)

cb = Colorbar(gcb[1,1], hm, ticks = (-0.005:0.001:-0.001, ["-5×10⁻³", "-4×10⁻³", "-3×10⁻³", "-2×10⁻³", "-1×10⁻³"] ), size = 35, label = "b")
cb = Colorbar(gcb[2,1], hm2, ticks = (1e-7:2e-7:1e-6, ["1×10⁻⁷", "3×10⁻⁷", "5×10⁻⁷", "7×10⁻⁷", "9×10⁻⁹"] ), size = 35, label = "ε")

colsize!(f.layout, 2, Relative(0.04))

savename = "beNLHovmol_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f)


############
#       PLOTTING Hovmollers
#   ______________ ___
#   | b10 | b30  | | |
#   |_____|______| | |
# y | v10 |  v30 | | |
#   |_____|______| |_|
#       t
###########

wave_times = b_timeseries.times./ pm.Tσ
ycut = yc[1:ylength]

f = Figure(resolution = (1600, 1000), fontsize=26)
ga = f[1, 1] = GridLayout()
gcb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "y' [m]")
ax2 = Axis(ga[1, 2])
ax3 = Axis(ga[2, 1], ylabel = "y' [m]", xlabel = "Tσ")
ax4 = Axis(ga[2, 2], xlabel = "Tσ")

limits!(ax1, 0, 8, 0, 1300)
limits!(ax2, 0, 8, 0, 1300)
limits!(ax3, 0, 8, 0, 1300)
limits!(ax4, 0, 8, 0, 1300)

ax1.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )
ax3.yticks = (0:260:1300, ["1300", "1040", "780", "520", "260", "0"] )

ax3.xticks = 1:1:11
ax4.xticks = 1:1:11

hidedecorations!(ax2)
hidexdecorations!(ax1)
hideydecorations!(ax4)

hm = heatmap!(ax1, wave_times, ycut, b_slopelines_10', colormap = :thermal, colorrange = (-6e-3, 0))
contour!(ax1, wave_times, ycut, b_slopelines_10', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax1, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

heatmap!(ax2, wave_times, ycut, b_slopelines_30_full', colormap = :thermal, colorrange = (-6e-3, 0))
contour!(ax2, wave_times, ycut, b_slopelines_30_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax2, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

hm2 = heatmap!(ax3, wave_times, ycut, v_slopelines_10_full', colormap = :balance, colorrange = (-0.3, 0.3))
contour!(ax3, wave_times, ycut, b_slopelines_10_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax3, Point.(0.25, 1200), text = "10 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

heatmap!(ax4, wave_times, ycut, v_slopelines_30_full', colormap = :balance, colorrange = (-0.3, 0.3))
contour!(ax4, wave_times, ycut, b_slopelines_30_full', color = :black, linewidth = 1, levels = -0.005:0.001:0, alpha = 0.5)
text!(ax4, Point.(0.25, 1200), text = "30 m", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

colgap!(ga, 20)
rowgap!(ga, 20)

cb = Colorbar(gcb[1,1], hm, ticks = (-0.005:0.001:-0.001, ["-5×10⁻³", "-4×10⁻³", "-3×10⁻³", "-2×10⁻³", "-1×10⁻³"] ), size = 35, label = "b")
cb = Colorbar(gcb[2,1], hm2, ticks = (-0.3:.15:0.3, ["-0.3", "-0.15", "0", "0.15", "0.3"] ), size = 35, label = "v")

colsize!(f.layout, 2, Relative(0.04))

savename = "bvHovmol_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f)
