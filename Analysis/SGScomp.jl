using Statistics
using Printf
using Oceananigans
using Measures
using JLD2

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
sn = "U300N100Lz100g100"

ENV["GKSwstype"] = "nul" # if on remote HPC

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
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

name_prefix = "IntWave_" * sn
filepath = path_name * name_prefix * ".jld2"

@info "getting data from: " * sn

SGS_timeseries = FieldTimeSeries(filepath, "SGS∇κ∇b");

xb, yb, zb = nodes(SGS_timeseries) #CCC

ylength = 880 
tlength = 161

SGSi = interior(SGS_timeseries)[:,1:ylength,:,:];

include("WaveValues.jl")
wave_info=get_wave_indices(SGS_timeseries, pm, tlength)
Wl = wave_info.Wl

end5waves = wave_info.WavePeriods[:,7:end] # t = 7T : 11 T
beg3waves = wave_info.WavePeriods[:,2:4] # t = 1T : t = 4T
trans3waves = wave_info.WavePeriods[:,5:7] # t = 4T : t = 7T

SGS_beg = SGSi[:,:,:,beg3waves]
SGS_mid = SGSi[:,:,:,trans3waves]
SGS_end = SGSi[:,:,:,end5waves]

SGS_beg_Wavg = mean(SGS_beg, dims = (1,4,5))[1,:,:,1,1];
SGS_mid_Wavg = mean(SGS_mid, dims = (1,4,5))[1,:,:,1,1];
SGS_end_Wavg = mean(SGS_end, dims = (1,4,5))[1,:,:,1,1];


f1 = Figure(resolution = (1000, 1200), fontsize=26)
ga = f1[1, 1] = GridLayout()
gb = f1[2, 1] = GridLayout()
gc = f1[3, 1] = GridLayout()
gcb = f1[1:3, 2] = GridLayout()

axb = Axis(ga[1, 1], ylabel = "z [m]") #v
axm = Axis(gb[1, 1], ylabel = "z [m]") #vh
axe = Axis(gc[1, 1], ylabel = "z [m]", xlabel = "y [m]") #vh

axe.xticks = 500:500:1500

axb.yticks = [-250, 0]
axm.yticks = [-250, 0]
axe.yticks = [-250, 0]


limits!(axb, 0, 2000, -500, 0)
limits!(axm, 0, 2000, -500, 0)
limits!(axe, 0, 2000, -500, 0)

hidexdecorations!(axb)
hidexdecorations!(axm)

colsize!(f1.layout, 2, Relative(0.05))

smax = 2e-8

hv = heatmap!(axb, yb[1:ylength], zb, SGS_beg_Wavg, colormap = :balance, colorrange = (-smax, smax))
    text!(axb, Point.(50, -450), text = "1-4 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)

hvhat = heatmap!(axm, yb[1:ylength], zb, SGS_mid_Wavg, colormap = :balance, colorrange = (-smax, smax))
    text!(axm, Point.(50, -450), text = "4-7 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)

hnab = heatmap!(axe,yb[1:ylength], zb, SGS_end_Wavg, colormap = :balance, colorrange = (-smax, smax))
    text!(axe, Point.(50, -450), text = "7-11 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
cb1 = Colorbar(gcb[1,1], hnab, ticks = (-1e-8:1e-8:1e-8), size =35, label = "⟨∇⋅κ∇b⟩")

colgap!(ga, 15)
colgap!(gb, 15)
colgap!(gc, 15)

save(path_name * "Analysis/SGS_WaveAverage.png", f1)

#####################
# Taking a Volume Integral
#####################

zlength = length(zb)
yslopelength = 334

Ygrid = reshape(repeat(yb[1:ylength], zlength), ylength, zlength);
Zgrid = permutedims(reshape(repeat(zb[1:zlength], ylength), zlength, ylength), [2,1]);
SlopeGridY = curvedslope.(Ygrid);
LinSlopeGridY = linslope.(Ygrid);
################
# \   \
#  \   \
#   \   \
#    \   \
#     \   \
#      \   \
#       \   \
#________\ __\_
#################
endslopedel = (500 +  pm.U₀/pm.Ñ)/pm.Tanα
endslope2del = (500 +  2*pm.U₀/pm.Ñ)/pm.Tanα

# all the values greater than slope
boolZY_aboveslope = (SlopeGridY.<= zb[1:zlength]') .& (Ygrid .<= endslopedel) .& (Zgrid .<= 0) ;
boolZY_belowdelta = boolZY_aboveslope .& ((LinSlopeGridY .+ pm.U₀/pm.Ñ) .>= zb[1:end]');
boolZY_abovedelta = (((LinSlopeGridY .+ pm.U₀/pm.Ñ) .<= zb[1:end]') .& ((LinSlopeGridY .+ 2*pm.U₀/pm.Ñ) .>= zb[1:end]')) .& (Ygrid .<= endslope2del) .& (Zgrid .<= 0);

# on the slope
int_SGSWavg_belowdelta_beg = sum(SGS_beg_Wavg[boolZY_belowdelta])
int_SGSWavg_belowdelta_mid = sum(SGS_mid_Wavg[boolZY_belowdelta])
int_SGSWavg_belowdelta_end = sum(SGS_end_Wavg[boolZY_belowdelta])
# in the interior on delta out
int_SGSWavg_abovedelta_beg = sum(SGS_beg_Wavg[boolZY_abovedelta])
int_SGSWavg_abovedelta_mid = sum(SGS_mid_Wavg[boolZY_abovedelta])
int_SGSWavg_abovedelta_end = sum(SGS_end_Wavg[boolZY_abovedelta])

################
# \  \
#  \  \
#   \__\
#    \  \
#-250 \  \** Gaussian Released
#      \__\
#       \  \
#________\ _\_
#################

# on the slope

boolZY_bd_bot = boolZY_belowdelta .& (Zgrid .< -334) 
boolZY_bd_mid = boolZY_belowdelta .& ((Zgrid .>= -334) .& (Zgrid .< -168))
boolZY_bd_top = boolZY_belowdelta .& ((Zgrid .>= -168) .& (Zgrid .< -0))

# off teh slope
boolZY_ad_bot = boolZY_abovedelta .& (Zgrid .< -334) 
boolZY_ad_mid = boolZY_abovedelta .& ((Zgrid .>= -334) .& (Zgrid .< -168))
boolZY_ad_top = boolZY_abovedelta .& ((Zgrid .>= -168) .& (Zgrid .< -0))

function volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_xpertxavgWavg)
    int_∇_xpertxavgWavg_bot = sum(∇_xpertxavgWavg[boolZY_bd_bot])
    int_∇_xpertxavgWavg_mid = sum(∇_xpertxavgWavg[boolZY_bd_mid])
    int_∇_xpertxavgWavg_top = sum(∇_xpertxavgWavg[boolZY_bd_top])

    boolZY_all = (boolZY_bd_top .* int_∇_xpertxavgWavg_top) .+ (boolZY_bd_mid .* int_∇_xpertxavgWavg_mid) .+ (boolZY_bd_bot .* int_∇_xpertxavgWavg_bot)

    vert_split_int_∇_xpertxavgWavg = (; int_∇_xpertxavgWavg_bot, int_∇_xpertxavgWavg_mid, int_∇_xpertxavgWavg_top)
    return vert_split_int_∇_xpertxavgWavg, boolZY_all
end


vert_split_int_∇_xpertxavgWavg_begT_bd, boolZY_all_begT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, SGS_beg_Wavg);
vert_split_int_∇_xpertxavgWavg_midT_bd, boolZY_all_midT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, SGS_mid_Wavg);
vert_split_int_∇_xpertxavgWavg_endT_bd, boolZY_all_endT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, SGS_end_Wavg);

vert_split_int_∇_xpertxavgWavg_begT_ad, boolZY_all_begT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, SGS_beg_Wavg);
vert_split_int_∇_xpertxavgWavg_midT_ad, boolZY_all_midT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, SGS_mid_Wavg);
vert_split_int_∇_xpertxavgWavg_endT_ad, boolZY_all_endT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, SGS_end_Wavg);

boolZY_all_begT = boolZY_all_begT_bd .+ boolZY_all_begT_ad
boolZY_all_midT = boolZY_all_midT_bd .+ boolZY_all_midT_ad
boolZY_all_endT = boolZY_all_endT_bd .+ boolZY_all_endT_ad

divmax = maximum((maximum(boolZY_all_begT), maximum(boolZY_all_midT), maximum(boolZY_all_endT)))
divmin = minimum((minimum(boolZY_all_begT), minimum(boolZY_all_midT), minimum(boolZY_all_endT)))

divext = maximum((divmax, abs(divmin)))


f1 = Figure(resolution = (1000, 1200), fontsize=26)
ga = f1[1, 1] = GridLayout()
gb = f1[2, 1] = GridLayout()
gc = f1[3, 1] = GridLayout()
gcb = f1[1:3, 2] = GridLayout()

axb = Axis(ga[1, 1], ylabel = "z [m]") #v
axm = Axis(gb[1, 1], ylabel = "z [m]") #vh
axe = Axis(gc[1, 1], ylabel = "z [m]", xlabel = "y [m]") #vh

axe.xticks = 500:500:1500

axb.yticks = [-250, 0]
axm.yticks = [-250, 0]
axe.yticks = [-250, 0]


limits!(axb, 0, 2000, -500, 0)
limits!(axm, 0, 2000, -500, 0)
limits!(axe, 0, 2000, -500, 0)

hidexdecorations!(axb)
hidexdecorations!(axm)

colsize!(f1.layout, 2, Relative(0.05))

hv = heatmap!(axb, yb[1:ylength], zb[1:zlength], boolZY_all_begT, colormap = :balance, colorrange = (-divext, divext))
    text!(axb, Point.(50, -450), text = "1-4 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)

hvhat = heatmap!(axm, yb[1:ylength], zb[1:zlength], boolZY_all_midT, colormap = :balance, colorrange = (-divext, divext))
    text!(axm, Point.(50, -450), text = "4-7 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)

hnab = heatmap!(axe, yb[1:ylength], zb[1:zlength], boolZY_all_endT, colormap = :balance, colorrange = (-divext, divext))
    text!(axe, Point.(50, -450), text = "7-11 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
cb1 = Colorbar(gcb[1,1], hnab, ticks = (-2e-6:2e-6:2e-6), size =35, label = "∫⟨∇⋅κ∇b⟩dV")

colgap!(ga, 15)
colgap!(gb, 15)
colgap!(gc, 15)

save(path_name * "Analysis/IntegratesSGS.png", f1)