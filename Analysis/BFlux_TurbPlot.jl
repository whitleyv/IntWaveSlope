using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using ArgParse
using CairoMakie

dpath = "Data/"

sn = "U300N100Lz100g100"

filescalename1 = dpath * "BFluxTurb_beg_" * sn * ".jld2"
filescalename2 = dpath * "BFluxTurb_mid_" * sn * ".jld2"
filescalename3 = dpath * "BFluxTurb_end_" * sn * ".jld2"

include("../parameters.jl") 
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       nz = round(Int,pm.Lz/2),
                       m = -π/pm.Lz,
                       l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                       Tf = 2*π/pm.f, 
                       Tσ = 2*π/pm.σ))

scale_file1 = jldopen(filescalename1, "r+")
scale_file2 = jldopen(filescalename2, "r+")
scale_file3 = jldopen(filescalename3, "r+")

skeys = keys(scale_file1)

FluxDiv_Wavg_beg  = scale_file1[skeys[1]];
FluxDiv_Phavg_beg  = scale_file1[skeys[2]];
PhaseAveragedVals_beg  = scale_file1[skeys[3]];
WaveAveragedVals_beg  = scale_file1[skeys[4]];

yb = scale_file1[skeys[5]];
zb = scale_file1[skeys[6]];
phase_times = scale_file1[skeys[7]];

FluxDiv_Wavg_mid  = scale_file2[skeys[1]];
FluxDiv_Phavg_mid = scale_file2[skeys[2]];
PhaseAveragedVals_mid  = scale_file2[skeys[3]];
WaveAveragedVals_mid  = scale_file2[skeys[4]];

FluxDiv_Wavg_end  = scale_file3[skeys[1]];
FluxDiv_Phavg_end  = scale_file3[skeys[2]];
PhaseAveragedVals_end  = scale_file3[skeys[3]];
WaveAveragedVals_end  = scale_file3[skeys[4]];

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

land = curvedslope.(yb)
land_pdel = (curvedslope.(yb) .+ pm.U₀/pm.Ñ)[1:382]

zlength = length(zb)
ylength = 880

vmax =  round(maximum(WaveAveragedVals_beg.v_xavgWavg)*100)*1e-2
vbmax = round(maximum(WaveAveragedVals_beg.vb_xpertxavgWavg)*1e6)*1e-6*0.5
wbmax = round(maximum(abs.(WaveAveragedVals_beg.wb_xpertxavgWavg))*1e6)*1e-6*0.5

nabmax = round(maximum(abs.(FluxDiv_Wavg_beg.vb_xpertxavgWavg_dy))*1e7)*1e-7*0.25
dymax = round(maximum(abs.(FluxDiv_Wavg_beg.vb_xpertxavgWavg_dy))*1e7)*1e-7*0.25
dzmax = round(maximum(FluxDiv_Wavg_beg.wb_xpertxavgWavg_dz)*1e7)*1e-7*0.25

# wave average dye:
c_xavgWavg_beg = mean(PhaseAveragedVals_beg.c_xavgPhavg, dims = 3)[:,:,1]
c_xavgWavg_mid = mean(PhaseAveragedVals_mid.c_xavgPhavg, dims = 3)[:,:,1]
c_xavgWavg_end = mean(PhaseAveragedVals_end.c_xavgPhavg, dims = 3)[:,:,1]

function WaveAverageFluxPlot(savename, apath , ylength, v_Wavg_xavg, vb_WpertWavg_xavg, wb_WpertWavg_xavg, v̂_Wavg_xavg, 
    v̂b_WpertWavg_xavg, ŵb_WpertWavg_xavg, vb_WpertWavg_xavg_dy, wb_WpertWavg_xavg_dz, ∇WpertWavg_xavg,c_xavgWavg,
    vmax, vbmax, wbmax, nabmax, dymax, dzmax)

    vblims = (-vbmax,vbmax)
    wblims = (-wbmax,wbmax)
    gradblims = (-nabmax, nabmax)
    dvblims = (-dymax, dymax)
    dwblims = (-dzmax, dzmax)

    yb_cut = yb[1:ylength]

    clog = log10.(clamp.(c_xavgWavg, 1e-8, 1))

    f1 = Figure(resolution = (1300, 1100), fontsize=26)
    ga = f1[2:4, 1] = GridLayout()
    gb = f1[2:4, 2] = GridLayout()
    gc = f1[2:4, 3] = GridLayout()

    gcb1 = f1[1, 1] = GridLayout()
    gcb2 = f1[1, 2] = GridLayout()       
    gcb3 = f1[1, 3] = GridLayout()       

    axv = Axis(ga[1, 1], ylabel = "z [m]") #v
    axvh = Axis(ga[2, 1], ylabel = "z [m]") #v
    axgrad = Axis(ga[3, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
    axgrad.xticks = 500:500:1000
    axv.yticks = [-250, 0]
    axvh.yticks = [-250, 0]
    axgrad.yticks = [-250, 0]

    axvb = Axis(gb[1, 1]) #v'b'
    axvbh = Axis(gb[2, 1]) #vh'b'
    axdyvb = Axis(gb[3, 1], xlabel = "y [m]") #dy v'b'

    axdyvb.xticks = 500:500:1000

    axwb = Axis(gc[1, 1]) #w'b'
    axwbh = Axis(gc[2, 1]) #wh'b'
    axdzwb = Axis(gc[3, 1], xlabel = "y [m]") #dz w'b'
    axdzwb.xticks = 500:500:1000

    limits!(axv, 0, 1500, -500, 0)
    limits!(axvh, 0, 1500, -500, 0)
    limits!(axvb, 0, 1500, -500, 0)
    limits!(axvbh, 0, 1500, -500, 0)
    limits!(axwb, 0, 1500, -500, 0)
    limits!(axwbh, 0, 1500, -500, 0)
    limits!(axgrad, 0, 1500, -500, 0)
    limits!(axdyvb, 0, 1500, -500, 0)
    limits!(axdzwb, 0, 1500, -500, 0)

    hidedecorations!(axvb)
    hidedecorations!(axwb)
    hidedecorations!(axvbh)
    hidedecorations!(axwbh)
    hidexdecorations!(axv)
    hidexdecorations!(axvh)
    hideydecorations!(axdyvb)
    hideydecorations!(axdzwb)

    rowsize!(f1.layout, 1, Relative(0.05))

    hv = heatmap!(axv, yb_cut, zb, v_Wavg_xavg, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axv, yb, land, color=:black, lw = 4)
        lines!(axv, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axv, Point.(50, -450), text = "v", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hvhat = heatmap!(axvh, yb_cut, zb, v̂_Wavg_xavg, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axvh, yb, land, color=:black, lw = 4)
        lines!(axvh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvh, Point.(50, -450), text = "v̂", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hnab = heatmap!(axgrad, yb_cut[2:end], zb[2:end], ∇WpertWavg_xavg, colormap = :balance, colorrange = gradblims)
        lines!(axgrad, yb, land, color=:black, lw = 4)
        lines!(axgrad, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axgrad, Point.(50, -450), text = "∇⋅⟨u⃗'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hvb = heatmap!(axvb, yb_cut, zb, vb_WpertWavg_xavg, colormap = :balance, colorrange = vblims)
        lines!(axvb, yb, land, color=:black, lw = 4)
        lines!(axvb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvb, Point.(50, -450), text = "⟨v'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        contour!(axvb, yb_cut, zb, clog, color = :gray30, lw = 6, levels = -5:1:0, labels=true)

    heatmap!(axvbh, yb_cut, zb, v̂b_WpertWavg_xavg, colormap = :balance, colorrange = vblims)
        lines!(axvbh, yb, land, color=:black, lw = 4)
        lines!(axvbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvbh, Point.(50, -450), text = "⟨v̂'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        

    hdyvb = heatmap!(axdyvb, yb_cut, zb, vb_WpertWavg_xavg_dy, colormap = :balance, colorrange = dvblims)
        lines!(axdyvb, yb, land, color=:black, lw = 4)
        lines!(axdyvb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdyvb, Point.(50, -450), text = "∂y⟨v'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hwb = heatmap!(axwb, yb_cut, zb, wb_WpertWavg_xavg, colormap = :balance, colorrange = wblims)
        lines!(axwb, yb, land, color=:black, lw = 4)
        lines!(axwb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axwb, Point.(50, -450), text = "⟨w'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        contour!(axwb, yb_cut, zb, clog, color = :gray30, lw = 6, levels = -5:1:0)


    heatmap!(axwbh, yb_cut, zb, ŵb_WpertWavg_xavg, colormap = :balance, colorrange = vblims)
        lines!(axwbh, yb, land, color=:black, lw = 4)
        lines!(axwbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axwbh, Point.(50, -450), text = "⟨ŵ'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hdzwb = heatmap!(axdzwb, yb_cut, zb, wb_WpertWavg_xavg_dz, colormap = :balance, colorrange = dwblims)
        lines!(axdzwb, yb, land, color=:black, lw = 4)
        lines!(axdzwb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdzwb, Point.(50, -450), text = "∂z⟨w'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax:vmax:vmax), size =35,vertical = false)
    cb2 = Colorbar(gcb2[1,1], hvb, ticks = (-vbmax/2:vbmax/2:vbmax/2), size =35,vertical = false, )
    cb3 = Colorbar(gcb3[1,1], hwb, ticks = (-wbmax/2:wbmax/2:wbmax/2), size =35,vertical = false, )

    cb3 = Colorbar(gcb1[1,1], hnab, ticks = (-nabmax/2:nabmax/2:nabmax/2),
            size =35,vertical = false, flipaxis=false)
    cb4 = Colorbar(gcb2[1,1], hdyvb, ticks = (-dymax/2:dymax/2:dymax/2),
            size =35,vertical = false, flipaxis=false)
    cb5 = Colorbar(gcb3[1,1], hdzwb, ticks = (-dzmax/2:dzmax/2:dzmax/2),
            size =35,vertical = false, flipaxis=false)

    colgap!(ga, 15)
    rowgap!(ga, 5)
    colgap!(gb, 15)
    rowgap!(gb, 5)
    colgap!(gc, 15)
    rowgap!(gc, 5)

    save(apath * savename * ".png", f1)
end

WaveAverageFluxPlot("blfuxturb_Wavg_beg" * sn,  "Analysis/Plots/", ylength, 
WaveAveragedVals_beg.v_xavgWavg, WaveAveragedVals_beg.vb_xpertxavgWavg,
WaveAveragedVals_beg.wb_xpertxavgWavg, WaveAveragedVals_beg.v̂_xavgWavg, 
WaveAveragedVals_beg.v̂b_xpertxavgWavg, WaveAveragedVals_beg.ŵb_xpertxavgWavg, 
FluxDiv_Wavg_beg.vb_xpertxavgWavg_dy, FluxDiv_Wavg_beg.wb_xpertxavgWavg_dz, FluxDiv_Wavg_beg.∇_xpertxavgWavg,
c_xavgWavg_beg, vmax, vbmax, wbmax, nabmax, dymax, dzmax)

WaveAverageFluxPlot("blfuxturb_Wavg_end" * sn,  "Analysis/Plots/", ylength, 
WaveAveragedVals_end.v_xavgWavg, WaveAveragedVals_end.vb_xpertxavgWavg,
WaveAveragedVals_end.wb_xpertxavgWavg, WaveAveragedVals_end.v̂_xavgWavg, 
WaveAveragedVals_end.v̂b_xpertxavgWavg, WaveAveragedVals_end.ŵb_xpertxavgWavg, 
FluxDiv_Wavg_end.vb_xpertxavgWavg_dy, FluxDiv_Wavg_end.wb_xpertxavgWavg_dz, FluxDiv_Wavg_end.∇_xpertxavgWavg,
c_xavgWavg_end, vmax, vbmax, wbmax, nabmax, dymax, dzmax)

WaveAverageFluxPlot("blfuxturb_Wavg_mid" * sn,  "Analysis/Plots/", ylength, 
WaveAveragedVals_mid.v_xavgWavg, WaveAveragedVals_mid.vb_xpertxavgWavg,
WaveAveragedVals_mid.wb_xpertxavgWavg, WaveAveragedVals_mid.v̂_xavgWavg, 
WaveAveragedVals_mid.v̂b_xpertxavgWavg, WaveAveragedVals_mid.ŵb_xpertxavgWavg, 
FluxDiv_Wavg_mid.vb_xpertxavgWavg_dy, FluxDiv_Wavg_mid.wb_xpertxavgWavg_dz, FluxDiv_Wavg_mid.∇_xpertxavgWavg,
c_xavgWavg_mid, vmax, vbmax, wbmax, nabmax, dymax, dzmax)


#####################
# Taking a Volume Integral
#####################

zlength = length(zb)
yslopelength = 334

Ygrid = reshape(repeat(yb[2:ylength], zlength-1), ylength-1, zlength-1);
Zgrid = permutedims(reshape(repeat(zb[2:zlength], ylength-1), zlength-1, ylength-1), [2,1]);
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
boolZY_aboveslope = ((SlopeGridY .+ 4) .<= zb[2:end]') .& (Ygrid .<= endslopedel) .& (Zgrid .<= 0) ;
boolZY_belowdelta = boolZY_aboveslope .& ((LinSlopeGridY .+ pm.U₀/pm.Ñ) .>= zb[2:end]');
boolZY_abovedelta = (((LinSlopeGridY .+ pm.U₀/pm.Ñ) .<= zb[2:end]') .& ((LinSlopeGridY .+ 2*pm.U₀/pm.Ñ) .>= zb[2:end]')) .& (Ygrid .<= endslope2del) .& (Zgrid .<= 0);

# on the slope
int_∇_xpertxavgWavg_belowdelta_beg = sum(FluxDiv_Wavg_beg.∇_xpertxavgWavg[boolZY_belowdelta])
int_∇_xpertxavgWavg_belowdelta_mid = sum(FluxDiv_Wavg_mid.∇_xpertxavgWavg[boolZY_belowdelta])
int_∇_xpertxavgWavg_belowdelta_end = sum(FluxDiv_Wavg_end.∇_xpertxavgWavg[boolZY_belowdelta])
# in the interior on delta out
int_∇_xpertxavgWavg_abovedelta_beg = sum(FluxDiv_Wavg_beg.∇_xpertxavgWavg[boolZY_abovedelta])
int_∇_xpertxavgWavg_abovedelta_mid = sum(FluxDiv_Wavg_mid.∇_xpertxavgWavg[boolZY_abovedelta])
int_∇_xpertxavgWavg_abovedelta_end = sum(FluxDiv_Wavg_end.∇_xpertxavgWavg[boolZY_abovedelta])

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


vert_split_int_∇_xpertxavgWavg_begT_bd, boolZY_all_begT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, FluxDiv_Wavg_beg.∇_xpertxavgWavg);
vert_split_int_∇_xpertxavgWavg_midT_bd, boolZY_all_midT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, FluxDiv_Wavg_mid.∇_xpertxavgWavg);
vert_split_int_∇_xpertxavgWavg_endT_bd, boolZY_all_endT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, FluxDiv_Wavg_end.∇_xpertxavgWavg);

vert_split_int_∇_xpertxavgWavg_begT_ad, boolZY_all_begT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, FluxDiv_Wavg_beg.∇_xpertxavgWavg);
vert_split_int_∇_xpertxavgWavg_midT_ad, boolZY_all_midT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, FluxDiv_Wavg_mid.∇_xpertxavgWavg);
vert_split_int_∇_xpertxavgWavg_endT_ad, boolZY_all_endT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, FluxDiv_Wavg_end.∇_xpertxavgWavg);

boolZY_all_begT = boolZY_all_begT_bd .+ boolZY_all_begT_ad
boolZY_all_midT = boolZY_all_midT_bd .+ boolZY_all_midT_ad
boolZY_all_endT = boolZY_all_endT_bd .+ boolZY_all_endT_ad

divmax = maximum((maximum(boolZY_all_begT), maximum(boolZY_all_midT), maximum(boolZY_all_endT)))
divmin = minimum((minimum(boolZY_all_begT), minimum(boolZY_all_midT), minimum(boolZY_all_endT)))

divext = maximum((divmax, abs(divmin)))

clog_beg = log10.(clamp.(c_xavgWavg_beg, 1e-8, 1))
clog_mid = log10.(clamp.(c_xavgWavg_mid, 1e-8, 1))
clog_end = log10.(clamp.(c_xavgWavg_end, 1e-8, 1))

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

hv = heatmap!(axb, yb[2:ylength], zb[2:zlength], boolZY_all_begT, colormap = :balance, colorrange = (-divext, divext))
    text!(axb, Point.(50, -450), text = "1-4 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)
    contour!(axb, yb[1:ylength], zb, clog_beg, color = :black, linewidth = 3, levels = -5:1:0)

hvhat = heatmap!(axm, yb[2:ylength], zb[2:zlength], boolZY_all_midT, colormap = :balance, colorrange = (-divext, divext))
    text!(axm, Point.(50, -450), text = "4-7 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)
    contour!(axm, yb[1:ylength], zb, clog_mid, color = :black, linewidth = 3, levels = -5:1:0)

hnab = heatmap!(axe, yb[2:ylength], zb[2:zlength], boolZY_all_endT, colormap = :balance, colorrange = (-divext, divext))
    text!(axe, Point.(50, -450), text = "7-11 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)
    contour!(axe, yb[1:ylength], zb, clog_end, color = :black, linewidth = 3, levels = -5:1:0)

    # create colorbars the size of the whole data set
cb1 = Colorbar(gcb[1,1], hnab, ticks = (-8e-6:4e-6:8e-6), size =35, label = "∫∇⋅⟨u'b'⟩dV")

colgap!(ga, 15)
colgap!(gb, 15)
colgap!(gc, 15)

save("Analysis/Plots/IntegratesFluxDiv.png", f1)