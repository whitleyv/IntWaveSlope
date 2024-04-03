using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

apath = path_name * "Analysis/"

include("parameters.jl")

setname = "U250N100Lz100g100"

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

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

name1_prefix = "IntWave_mp_noS_" * setname # "IntWave_mp_U250N100Lz100g100"

filepath1 = path_name * name1_prefix * ".jld2"

b_timeseries = FieldTimeSeries(filepath1,"b");
v_timeseries = FieldTimeSeries(filepath1,"v");

ylength = 650

bi = interior(b_timeseries)[:,1:ylength,:,:];
vi = interior(v_timeseries)[:,1:ylength,:,:];

xb, yb, zb = nodes(b_timeseries) #CCC
(Δx, Δy, Δz) = (4.0, 4.0, 2.0)
lin_land = linslope.(yb[1:334])
xlength = length(xb)
tlength = length(b_timeseries.times)

### need to try different h values:
δ = pm.U₀/pm.Ñ
α = atan(pm.Tanα)
δp = δ/cos(α)
habrat = 5#1.1
hab = habrat* δ
δy = hab/pm.Tanα

@info "HAB = $hab"
# including all values out to horiozontal depth:
HABidxL = round(Int, δy/4)
yidx_start = 1
yidx_end = 334
yidx_length = length(yidx_start:yidx_end)

land = curvedslope.(yb[yidx_start:yidx_end])
zidx = sum( zb .< land', dims = 1)[1,:]
Idx_Array_z = reshape(repeat(zidx, HABidxL), yidx_length, HABidxL)
Idx_Array_y = reshape(repeat(yidx_start:yidx_end, HABidxL), yidx_length, HABidxL) .+ (0:HABidxL-1)'

b_RotatedArray_BelowHAB = zeros(xlength, yidx_length, HABidxL, tlength);
v_RotatedArray_BelowHAB = zeros(xlength, yidx_length, HABidxL, tlength);

for k = 1:yidx_length
    zdx = zidx[k] # each z index 
    for j = 1:HABidxL
        ydx = Idx_Array_y[k,j]
        b_RotatedArray_BelowHAB[:,k,j, :] = bi[:,ydx, zdx, :]
        v_RotatedArray_BelowHAB[:,k,j, :] = vi[:,ydx, zdx, :]
    end
    
end

########
# BIN FIRST THEN AVERAGE: binning buoyancy at every time step and calcualting the total volume
########
ini_Nisos = 75
Δbf = -500*pm.Ñ^2/ini_Nisos
Δb = round(1.04 * Δbf*1e7)*1e-7

V_inb = zeros(ini_Nisos, tlength)
vavg_inb = zeros(ini_Nisos, tlength)

for i = 1:tlength
    @info "Time $i of $tlength..."
    # at each time step get the data at (x,ycut,z)
    bt = b_RotatedArray_BelowHAB[:,:, :, i];
    vt = v_RotatedArray_BelowHAB[:,:,:, i];

    # for each buoyancy class:
    for n = 1:ini_Nisos
        # (1) CCC locations in the density class
        # starting with b < 0 should be almost the whole domain
        boolB = (bt.< Δb*(n-1)) .& (bt.>= Δb *n)
        # volume integrated dye concentration:
        V_inb[n,i] = sum(boolB) * Δx * Δy * Δz
        vavg_inb[n,i] = sum(vt[boolB])./sum(boolB)
    end
end

MeanV_inb = mean(V_inb[15:end - 15, 1])
ΔV_inb = V_inb .- MeanV_inb #V_inb[:, 1]

include("WaveValues.jl")

wave_info=get_wave_indices(b_timeseries, pm, tlength)
endwaves = wave_info.WavePeriods[:,8:11]

ΔV_inb_Wavg = mean(ΔV_inb[:,endwaves], dims = (2,3))[:,1,1];
V_inb_Wavg = mean(V_inb[:,endwaves], dims = (2,3))[:,1,1];
vavg_inb_Wavg = mean(vavg_inb[:,endwaves], dims = (2,3))[:,1,1];

b_Wavg = mean(b_RotatedArray_BelowHAB[:,:,:,endwaves], dims = (4,5))[:,:,:,1,1];
v_Wavg = mean(v_RotatedArray_BelowHAB[:,:,:,endwaves], dims = (4,5))[:,:,:,1,1];

V_Wavg_inb = zeros(ini_Nisos);
vavg_Wavg_inb = zeros(ini_Nisos);

# for each buoyancy class:
for n = 1:ini_Nisos
    # (1) CCC locations in the density class
    # starting with b < 0 should be almost the whole domain
    boolB = (b_Wavg.< Δb*(n-1)) .& (b_Wavg.>= Δb *n)
    # volume integrated dye concentration:
    V_Wavg_inb[n] = sum(boolB) * Δx * Δy * Δz
    vavg_Wavg_inb[n] = sum(v_Wavg[boolB])./sum(boolB)

end

ΔV_Wavg_inb = V_Wavg_inb .- MeanV_inb #V_inb[:, 1]

b_bins = (Δb .* (1:ini_Nisos))

habround = round(Int, hab)

f1 = Figure(resolution = (1500, 900), fontsize=26)
    ga = f1[1, 1] = GridLayout()   
    ax1 = Axis(ga[1, 1], xlabel = "ΔVol in ℛ̄(b) [m³]", ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[1, 2], xlabel = "Average Horizontal Velocity in ℛ̄(b) [ms⁻¹]", ylabel = "b [ms⁻²]") #vh

    ax1.yticks = -5e-3:1e-3:-1e-3
    ax2.yticks = -5e-3:1e-3:-1e-3
    limits!(ax1, -2.1e6, 2.1e6, Δb*ini_Nisos, Δb)
    limits!(ax2, -0.04, 0.04, Δb*ini_Nisos, Δb)

    lines!(ax1, ΔV_inb_Wavg, b_bins, color = :black, linewidth = 2, label = "binned then avg")
    lines!(ax1, ΔV_Wavg_inb, b_bins, color = :blue, linewidth = 2, label = "avg then binned")
    hspan!(ax1, -3e-3, (-3e-3+δ*pm.Ñ^2) , color = (:blue, 0.2), label = "δN²")

    lines!(ax2, vavg_inb_Wavg, b_bins, color = :black, linewidth = 2, label = "v̄, binned then avg")
    lines!(ax2, vavg_Wavg_inb, b_bins, color = :blue, linewidth = 2, label = "v̄, binned then avg")

   #axislegend(ax1, position = :lb)
   f1[1,2] = Legend(f1, ax1)

    Label(ga[1, 1:2, Top()], "Wave averaged over waves 7-11 Tσ, $(ini_Nisos) bins, h = $(habrat)δ",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)


    savename = apath * "TotalVolumeChange_nb_$(ini_Nisos)_$(habround)_Wavg_wv_"* setname 
    savename = apath * "TotalVolumeChange_nb_$(ini_Nisos)_$(habround)_Wavg_wv_constVol_"* setname 

save(savename * ".png", f1)

######
#   NORMALIZED VOLUME
#####
f1 = Figure(resolution = (1000, 900), fontsize=26)
    ga = f1[1, 1] = GridLayout()   
    ax1 = Axis(ga[1, 1], xlabel = "Normalized Total Changes in Volume into ℛ̄(b)", ylabel = "b [ms⁻²]") #vh

    ax1.yticks = -5e-3:1e-3:-1e-3

    lines!(ax1, ΔV_inb_Wavg[1:end-12]./ V_inb[1:end-12, 1], b_bins[1:end-12], color = :black, linewidth = 2, label = "ΔV, binned then avg")
    lines!(ax1, ΔV_Wavg_inb[1:end-12]./V_inb[1:end-12, 1], b_bins[1:end-12], color = :blue, linewidth = 2, label = "ΔV, avg then binned")
    hspan!(ax1, -3e-3, (-3e-3+δ*pm.Ñ^2) , color = (:blue, 0.2), label = "δN²")
   #axislegend(ax1, position = :lb)
   f1[1,2] = Legend(f1, ax1)

    Label(ga[1, 1, Top()], "Wave averaged over waves 7-11 Tσ, $(ini_Nisos) bins",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)


    savename = apath * "TotalVolumeChange_nb_$(ini_Nisos)_$(habround)_Wavg_norm_"* setname 

save(savename * ".png", f1)

#####
#  average Volume compared to initial
####
f1 = Figure(resolution = (1000, 900), fontsize=26)
    ga = f1[1, 1] = GridLayout()   
    ax1 = Axis(ga[1, 1], xlabel = "Total Volume into ℛ̄(b)", ylabel = "b [ms⁻²]") #vh

    ax1.yticks = -5e-3:1e-3:-1e-3

    lines!(ax1, V_inb_Wavg, b_bins, color = :black, linewidth = 2, label = "V̄, binned then avg")
    lines!(ax1, V_inb[:,1], b_bins, color = :red, linewidth = 2, label = "V₀, binned")
    hspan!(ax1, -3e-3, (-3e-3+δ*pm.Ñ^2) , color = (:blue, 0.2), label = "δN²")
   #axislegend(ax1, position = :lb)
   f1[1,2] = Legend(f1, ax1)

    Label(ga[1, 1, Top()], "Wave averaged over waves 7-11 Tσ, $(ini_Nisos) bins",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)


    savename = apath * "TotalVolume_nb_$(ini_Nisos)_$(habround)_Wavg_"* setname 

save(savename * ".png", f1)


######
#   TESTING VALUES
#####

yflat = vcat(Idx_Array_y...);
zflat = vcat(Idx_Array_z...);

btestrot = b_RotatedArray_BelowHAB[19,:,:,end]
btestrot = b_Wavg[19,:,:];
btestrotinit = b_RotatedArray_BelowHAB[19,:,:,1];

bflat = vcat(btestrot...);
btest = mean(bi[19,:,:,endwaves], dims = (3,4))[:,:,1,1]; #bi[19,:,:,end];
vtest = mean(vi[19,:,:,endwaves], dims = (3,4))[:,:,1,1]; #bi[19,:,:,end];

bini = bi[19,:,:,1]
Δb = btestrot .- btestrotinit;

Δbflat = vcat(Δb...);

# I need to start 8 grid poinrts in from the top and  
f = Figure(resolution = (1000,700)) 
    ga = f[1, 1] = GridLayout()    
    ax1 = Axis(ga[1, 1]) 
    limits!(ax1, 0, 2500, -500,0)
    scatter!(ax1, yb[yflat], zb[zflat], color = Δbflat, colormap = :balance, colorrange = (-maximum(abs.(Δbflat)), maximum(abs.(Δbflat))), markersize=4)    
    #scatter!(ax1, yb[yflat], zb[zflat], color = bflat, colormap = :thermal, markersize=4)
    lines!(ax1, yb[1:ylength], curvedslope.(yb[1:ylength]), color = :blue)
    lines!(ax1, yb[1:ylength], linslope.(yb[1:ylength]) .+ δ, color = :blue)
    contour!(ax1, yb[1:ylength], zb, bini, levels = b_bins, color = :green, linewidth = 1)
    contour!(ax1, yb[1:ylength], zb, btest, levels = b_bins, color = :black, linewidth = 1)

savename = apath * "btestrotationfullendavgdif" 

save(savename * ".png", f)

f = Figure(resolution = (1000,700)) 
    ga = f[1, 1] = GridLayout()    
    ax1 = Axis(ga[1, 1]) 
    limits!(ax1, 0, 2500, -500,0)
    heatmap!(ax1, yb[1:ylength], zb, vtest, colormap = :balance, colorrange = (-0.1, 0.1))    
    lines!(ax1, yb[1:ylength], curvedslope.(yb[1:ylength]), color = :blue)
    lines!(ax1, yb[1:ylength], linslope.(yb[1:ylength]) .+ δ, color = :blue)

savename = apath * "vavg_end" 

save(savename * ".png", f)