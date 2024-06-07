using Statistics
using Printf
using Measures
using JLD2
using CairoMakie
using Oceananigans

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"
sn = "U350N100Lz100g100"

ENV["GKSwstype"] = "nul" # if on remote HPC

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = π/pm.Lz,
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

name_prefix = "IntWave_mp_noS_" * sn
filepath = path_name * name_prefix * ".jld2"

b_timeseries = FieldTimeSeries(filepath, "b");
xb, yb, zb = nodes(b_timeseries) #CCC
 
ylength = 880 
tlength = 161
land = curvedslope.(yb[1:ylength])

bi = interior(b_timeseries)[:,1:ylength,:,:];

δ = pm.U₀/pm.Ñ
α = atan(pm.Tanα)
δp = δ/cos(α)
habrat = 5#1.1
hab = habrat* δ
δy = hab/pm.Tanα

HABidxL = round(Int, δy/4)
yidx_start = 1
yidx_end = 334
yidx_length = length(yidx_start:yidx_end)

land = curvedslope.(yb[yidx_start:yidx_end])
zidx = sum( zb .< land', dims = 1)[1,:]
Idx_Array_z = reshape(repeat(zidx, HABidxL), yidx_length, HABidxL)
Idx_Array_y = reshape(repeat(yidx_start:yidx_end, HABidxL), yidx_length, HABidxL) .+ (0:HABidxL-1)'

xlength = 38
b_RotatedArray_BelowHAB = zeros(xlength, yidx_length, HABidxL, tlength);

for k = 1:yidx_length
    zdx = zidx[k] # each z index 
    for j = 1:HABidxL
        ydx = Idx_Array_y[k,j]
        b_RotatedArray_BelowHAB[:,k,j, :] = bi[:,ydx, zdx, :]
    end
end

ini_Nisos = 25 #60#50
Δb = -500*pm.Ñ^2/ini_Nisos
V_inb_waves = zeros(ini_Nisos, tlength);
# for each buoyancy class:
for i = 1:tlength
    bt = b_RotatedArray_BelowHAB[:,:, :, i];

    for n = 1:ini_Nisos
        # (1) CCC locations in the density class
        boolB = (bt.< Δb*(n-1)) .& (bt.>= Δb *n)
        # volume integrated dye concentration:
        V_inb_waves[n,i] = sum(boolB) * 32
    end
end

ΔV_inb_waves = V_inb_waves .- V_inb_waves[:,1]

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl

# small wave averages to see full effect:
#1-3
waves25 = wave_info.WavePeriods[:,3:5] # t = 7T : 11 T
waves58 = wave_info.WavePeriods[:,6:8] # t = 7T : 11 T
waves811 = wave_info.WavePeriods[:,9:11] # t = 7T : 11 T

# putting these in wave form
b_waves25 = bi[:,:,:,waves25];
b_waves58 = bi[:,:,:,waves58];
b_waves811 = bi[:,:,:,waves811];

# averaging over waves

b_Wavg25_xavg = mean(b_waves25, dims = (1,4,5))[1,:,:,1,1];
b_Wavg58_xavg = mean(b_waves58, dims = (1,4,5))[1,:,:,1,1];
b_Wavg811_xavg = mean(b_waves811, dims = (1,4,5))[1,:,:,1,1];

b₀ = bi[19,:,:, 1];

b_Wavg25 = mean(b_waves25, dims = (4,5))[:,:,:,1,1];
b_Wavg58 = mean(b_waves58, dims = (4,5))[:,:,:,1,1];
b_Wavg811 = mean(b_waves811, dims = (4,5))[:,:,:,1,1];

b_Rotated_Wavg25 = mean(b_RotatedArray_BelowHAB[:,:,:,waves25], dims = (4,5))[:,:,:,1,1];
b_Rotated_Wavg58 = mean(b_RotatedArray_BelowHAB[:,:,:,waves58], dims = (4,5))[:,:,:,1,1];
b_Rotated_Wavg811 = mean(b_RotatedArray_BelowHAB[:,:,:,waves811], dims = (4,5))[:,:,:,1,1];

ini_Nisos = 25 #60#50
Δb = -500*pm.Ñ^2/ini_Nisos
V_inb_waves25 = zeros(ini_Nisos);
V_inb_waves58 = zeros(ini_Nisos);
V_inb_waves811 = zeros(ini_Nisos);
V_inb_waves0 = zeros(ini_Nisos);

b0 = b_RotatedArray_BelowHAB[:,:,:,1]#bi[:,:,:, 1]
# for each buoyancy class:
for n = 1:ini_Nisos
    # (1) CCC locations in the density class
    # starting with b < 0 should be almost the whole domain
    boolB25 = (b_Rotated_Wavg25.< Δb*(n-1)) .& (b_Rotated_Wavg25.>= Δb *n)
    boolB58 = (b_Rotated_Wavg58.< Δb*(n-1)) .& (b_Rotated_Wavg58.>= Δb *n)
    boolB811 = (b_Rotated_Wavg811.< Δb*(n-1)) .& (b_Rotated_Wavg811.>= Δb *n)
    #boolB25 = (b_Wavg25.< Δb*(n-1)) .& (b_Wavg25.>= Δb *n)
    #boolB58 = (b_Wavg58.< Δb*(n-1)) .& (b_Wavg58.>= Δb *n)
    #boolB811 = (b_Wavg811.< Δb*(n-1)) .& (b_Wavg811.>= Δb *n)

    boolB0 = (b0.< Δb*(n-1)) .& (b0.>= Δb *n)

    V_inb_waves25[n] = sum(boolB25) * 32
    V_inb_waves58[n] = sum(boolB58) * 32
    V_inb_waves811[n] = sum(boolB811) * 32
    V_inb_waves0[n] = sum(boolB0) * 32
end

# range to exlude values with samll initial bins:
inirange = 2:22
b_bins = (Δb .* (1:ini_Nisos))[inirange]

ΔV_inb_waves25 = (V_inb_waves25 .- V_inb_waves0)[inirange]
ΔV_inb_waves58 = (V_inb_waves58 .- V_inb_waves0)[inirange]
ΔV_inb_waves811 = (V_inb_waves811 .- V_inb_waves0)[inirange]

ΔV_inb_waves25 = mean(ΔV_inb_waves[inirange, waves25],dims = (2,3))[:,1,1];
ΔV_inb_waves58 = mean(ΔV_inb_waves[inirange, waves58],dims = (2,3))[:,1,1];
ΔV_inb_waves811 = mean(ΔV_inb_waves[inirange, waves811], dims = (2,3))[:,1,1];

topo = curvedslope.(yb[1:ylength])

f1 = Figure(resolution = (1500, 1600), fontsize=30)
    ga = f1[1, 1] = GridLayout()    
    gb = f1[1, 2] = GridLayout()    

    axb = Axis(ga[2, 1], ylabel = "z [m]") 
    axb2 = Axis(ga[3, 1], ylabel = "z [m]") 
    axb3 = Axis(ga[4, 1], ylabel = "z [m]", xlabel = "y [m]") 

    gcb1 = ga[1, 1] = GridLayout()

    axv = Axis(gb[2,1], ylabel = "b [ms⁻²]")
    axv2 = Axis(gb[3,1], ylabel = "b [ms⁻²]",)
    axv3 = Axis(gb[4,1], ylabel = "b [ms⁻²]", xlabel = "ΔV = V(b,t) - V(b,0) [m³]")

    axb3.xticks = 500:500:2500
    axb.yticks = -500:250:0
    axb2.yticks = -500:250:-1
    axb3.yticks = -500:250:-1

    limits!(axb, 0, 2500, -500, 0)
    limits!(axb2, 0, 2500, -500, 0)
    limits!(axb3, 0, 2500, -500, 0)

    hidexdecorations!(axb)
    hidexdecorations!(axb2)
    b_bins_coarse = -5e-3:1e-3:-1e-3
    yb_cut = yb[1:ylength]
    Δblims = (-4.5e-4, 4.5e-4)

    hb = heatmap!(axb, yb_cut, zb, b_Wavg25_xavg .- b₀, colormap = :balance, colorrange = Δblims)
        contour!(axb, yb_cut, zb, b₀, color = :black, linewidth = 3, levels =b_bins_coarse)
        contour!(axb, yb_cut, zb, b_Wavg25_xavg, color = :gray35, linewidth = 5, levels =b_bins_coarse)
        lines!(axb, yb_cut, topo, color=:black, linewidth = 4)
        lines!(axb, yb_cut, linslope.(yb_cut) .+ pm.U₀/pm.Ñ, color=:black, linewidth = 3, linestyle = :dash)
        lines!(axb, yb_cut, linslope.(yb_cut) .+ hab, color=:black, linewidth = 7)
        #axislegend(axb, position = :lb, [LineElement(color = :gray41, linestyle = nothing, linewidth = 4)], ["Wave averaged: 2-5 Tσ"])
        text!(axb, Point.(100, -450), text = "2-5 Tσ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 30)

    heatmap!(axb2, yb_cut, zb, b_Wavg58_xavg .- b₀, colormap = :balance, colorrange = Δblims)
        contour!(axb2, yb_cut, zb, b₀, color = :black, linewidth = 3, levels =b_bins_coarse)
        contour!(axb2, yb_cut, zb, b_Wavg58_xavg, color = :gray35, linewidth = 5, levels =b_bins_coarse)
        lines!(axb2, yb_cut, linslope.(yb_cut) .+ pm.U₀/pm.Ñ, color=:black, linewidth = 3, linestyle = :dash)
        lines!(axb2, yb_cut, topo, color=:black, linewidth = 4)
        lines!(axb2, yb_cut, linslope.(yb_cut) .+ hab, color=:black, linewidth = 7)
        #axislegend(axb2, position = :lb, [LineElement(color = :gray41, linestyle = nothing, linewidth = 4)], ["Wave averaged: 5-8 Tσ"])
        text!(axb2, Point.(100, -450), text = "5-8 Tσ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 30)

    heatmap!(axb3, yb_cut, zb, b_Wavg811_xavg .- b₀, colormap = :balance, colorrange = Δblims)
        contour!(axb3, yb_cut, zb, b₀, color = :black, linewidth = 3, levels =b_bins_coarse)
        contour!(axb3, yb_cut, zb, b_Wavg811_xavg, color = :gray35, linewidth = 5, levels =b_bins_coarse) #bmin:0.0005:0, alpha = 0.5)
        lines!(axb3, yb_cut, linslope.(yb_cut) .+ pm.U₀/pm.Ñ, color=:black, linewidth = 3, linestyle = :dash)
        lines!(axb3, yb_cut, topo, color=:black, linewidth = 4)
        lines!(axb3, yb_cut, linslope.(yb_cut) .+ hab, color=:black, linewidth = 7)
        #axislegend(axb3, position = :lb, [LineElement(color = :gray41, linestyle = nothing, linewidth = 4)], ["Wave averaged: 8-11 Tσ"])
        text!(axb3, Point.(100, -450), text = "8-11 Tσ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 30)

    # create colorbars the size of the whole data set #label = "Δb = b - b₀ [ms⁻²]
    cb1 = Colorbar(gcb1[1,1], hb, ticks = (-4e-4:2e-4:4e-4), size =35, vertical = false) 
    cb1.alignmode = Mixed(top = 0)

    axv.yticks = -6e-3:1e-3:0
    axv.xticks = -2e6:2e6:4e6
    axv2.yticks = -6e-3:1e-3:0
    axv2.xticks = -2e6:2e6:4e6
    axv3.yticks = -6e-3:1e-3:0
    axv3.xticks = -2e6:2e6:4e6

    hidexdecorations!(axv, grid = false)
    hidexdecorations!(axv2, grid = false)

    limits!(axv, -2e6, 5e6, -5.5e-3, -5e-4)
    limits!(axv2, -2e6, 5e6, -5.5e-3, -5e-4)
    limits!(axv3, -2e6, 5e6, -5.5e-3, -5e-4)

    v2 = lines!(axv, ΔV_inb_waves25, b_bins, color = :black, linewidth = 6, label = "t = 2-5 Tσ")
    v5 = lines!(axv2, ΔV_inb_waves58, b_bins, color = :dodgerblue2, linewidth = 6, label = "t = 5-8 Tσ")
    v11 = lines!(axv3, ΔV_inb_waves811, b_bins, color = :firebrick2, linewidth = 6, label = "t = 8-11 Tσ")

    Legend( gb[1, 1], [v2, v5, v11], ["t = 2-5 Tσ", "t = 5-8 Tσ", "t = 8-11 Tσ"],
                    tellheight = false, tellwidth = false,# framevisible = false,
                    #margin = (10, 10, 10, 10), patchlabelgap = 7,
                    halign = :right, valign = :center, orientation = :horizontal)

    #leg = Legend(gb[2,1], axv, "Wave Average Range", margin = (100,100,0,0), framevisible = false, labelsize = 30)
    #axislegend(position=:rt, "Wave Average Range")
    #leg = Legend(gb[1,1],  [LineElement(color = :black, linewidth = 5), LineElement(color = :dodgerblue2, linewidth = 5), LineElement(color = :firebrick2, linewidth = 5)], 
    #        ["t = 2-5 Tσ", "t = 5-8 Tσ", "t = 8-11 Tσ"], "Wave Average Range",
    #        margin = (20,20,0,0), vertical = false)
    rowsize!(gb, 1, Relative(0.05))
    rowsize!(ga, 1, Relative(0.05))
    colsize!(f1.layout, 2, Relative(0.35))

    Label(ga[1, 1, TopLeft()], "a",
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :left)
    Label(gb[1, 1, TopLeft()], "b",
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :left)
    Label(ga[1, 1, Top()], "Δb = b - b₀ [ms⁻²]",
                    fontsize = 30,
                    padding = (5, 5, 5, 5),
                    halign = :center)
    Label(gb[1, 1, Top()], "Wave Average Range",
                    fontsize = 30,
                    padding = (5, 30, 5, 5),
                    halign = :center)

    colgap!(ga, 15)
    rowgap!(ga, 15)
    colgap!(gb, 15)
    rowgap!(gb, 15)

    #Box(f1[1, 1], color = (:red, 0.2), strokewidth = 0)
    #Box(f1[1, 2], color = (:blue, 0.2), strokewidth = 0)
    #Box(gb[1, 1], color = (:red, 0.2), strokewidth = 0)
    #Box(ga[1, 1], color = (:green, 0.2), strokewidth = 0)

save(path_name * "Analysis/b_Wavg_intime_Vol" * sn * ".png", f1, px_per_unit = 2)


topo = curvedslope.(yb[1:ylength])

f1 = Figure(resolution = (1200, 800), fontsize=30)
    ga = f1[1, 1] = GridLayout()    

    axb3 = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]") 

    gcb1 = ga[1, 1] = GridLayout()

    axb3.xticks = 500:500:2000
    axb3.yticks = -500:250:0

    limits!(axb3, 0, 2500, -500, 0)

    b_bins_coarse = -5e-3:1e-3:-1e-3
    yb_cut = yb[1:ylength]
    Δblims = (-4.5e-4, 4.5e-4)

    hb = heatmap!(axb3, yb_cut, zb, b_Wavg811_xavg .- b₀, colormap = :balance, colorrange = Δblims)
        contour!(axb3, yb_cut, zb, b₀, color = :black, linewidth = 3, levels =b_bins_coarse)
        contour!(axb3, yb_cut, zb, b_Wavg811_xavg, color = :gray35, linewidth = 5, levels =b_bins_coarse) #bmin:0.0005:0, alpha = 0.5)
        lines!(axb3, yb_cut, linslope.(yb_cut) .+ pm.U₀/pm.Ñ, color=:black, linewidth = 3, linestyle = :dash)
        band!(axb3, yb_cut, topo, -500, color=:black)

    # create colorbars the size of the whole data set #label = "Δb = b - b₀ [ms⁻²]
    cb1 = Colorbar(gcb1[1,1], hb, ticks = (-4e-4:2e-4:4e-4), size =35, vertical = false, label = "Δb = b - b₀ [ms⁻²]") 

    rowsize!(ga, 1, Relative(0.05))

save(path_name * "Analysis/b_Wavg_intime_Vol_end_" * sn * ".png", f1, px_per_unit = 2)
