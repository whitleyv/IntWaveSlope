using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

apath = path_name * "Analysis/"

include("parameters.jl")

setname = "U350N100Lz100g100"

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

∇κ∇B_timeseries = FieldTimeSeries(filepath1, "SGS∇κ∇b");
b_timeseries = FieldTimeSeries(filepath1,"b");
v_timeseries = FieldTimeSeries(filepath1,"v");
w_timeseries = FieldTimeSeries(filepath1,"w");
N2_timeseries = FieldTimeSeries(filepath1, "N2");
e_timeseries = FieldTimeSeries(filepath1,"ϵ");

ylength = 650 #500

∇κ∇Bi = interior(∇κ∇B_timeseries)[:,1:ylength,:,:];
bi = interior(b_timeseries)[:,1:ylength+1,:,:];
vi = interior(v_timeseries)[:,1:ylength,:,:];
wi = interior(w_timeseries)[:,1:ylength,:,:];
N2i = interior(N2_timeseries)[:,1:ylength,:,:];
ei = interior(e_timeseries)[:,1:ylength,:,:];

xb, yb, zb = nodes(b_timeseries) #CCC

lin_land = linslope.(yb[1:334])
xlength = length(xb)
tlength = length(b_timeseries.times)

function deriv(A, h, dim_order)
    Ap = permutedims(A, dim_order)
    step = 1/h
    step2 = 1/(h*2)
    da = zeros(size(Ap));
    da[1, :] = (Ap[2,:] .- Ap[1,:]) .* step
    da[end, :] = (Ap[end, :] .- Ap[end-1,:]) .* step
    da[2:end-1,:] = (Ap[3:end,:] .- Ap[1:end-2,:]) .* step2
        
    dap = permutedims(da, dim_order)
    return dap
end

### need to try different h values:
δ = pm.U₀/pm.Ñ
α = atan(pm.Tanα)
δp = δ/cos(α)
habrat = 1.1
hab = habrat* δ
δy = hab/pm.Tanα
#HAB = (0.2:0.4:2.4).* δp

@info "HAB = $hab"
###################
# Creating Rotated Array near slope
####################

# including all values out to horiozontal depth:
HABidxL = round(Int, δy/4)
yidx_start = 2
yidx_end = 334
yidx_length = length(yidx_start:yidx_end)

land = curvedslope.(yb[yidx_start:yidx_end]) 
zidx = sum( zb .<= land', dims = 1)[1,:] .+ 1
Idx_Array_z = reshape(repeat(zidx, HABidxL), yidx_length, HABidxL)
Idx_Array_y = reshape(repeat(yidx_start:yidx_end, HABidxL), yidx_length, HABidxL) .+ (0:HABidxL-1)'

b_RotatedArray_BelowHAB = zeros(xlength, yidx_length, HABidxL, tlength);
∇κ∇B_RotatedArray_BelowHAB = zeros(xlength, yidx_length, HABidxL, tlength);
v_RotatedArray_HAB = zeros(xlength, yidx_length, tlength);
w_RotatedArray_HAB = zeros(xlength, yidx_length, tlength);
vfull_RotatedArray_HAB = zeros(xlength, yidx_length, HABidxL, tlength);
e_RotatedArray_HAB = zeros(xlength, yidx_length, HABidxL, tlength);
N2_RotatedArray_HAB = zeros(xlength, yidx_length, HABidxL, tlength);

for k = 1:yidx_length
    zdx = zidx[k] # each z index 

    for j = 1:HABidxL
        ydx = Idx_Array_y[k,j]

        N2_RotatedArray_HAB[:,k,j,:] = N2i[:,ydx, zdx, :]
        b_RotatedArray_BelowHAB[:,k,j, :] = bi[:,ydx, zdx, :]
        ∇κ∇B_RotatedArray_BelowHAB[:,k,j, :] = ∇κ∇Bi[:,ydx, zdx, :]
        vfull_RotatedArray_HAB[:,k,j, :] = vi[:,ydx, zdx, :]
        e_RotatedArray_HAB[:,k,j, :] = ei[:,ydx, zdx, :]

    end

    Hidx_y = Idx_Array_y[k, end]
    # create arrays of only the values at HAB:
    v_RotatedArray_HAB[:,k, :] = vi[:,Hidx_y, zdx, :]
    w_RotatedArray_HAB[:,k, :] = wi[:,Hidx_y, zdx, :]
    
end

##############
# Calculating BInning for Volume budget values
##############
# volume flux at each point, with normal vector normalized:
Hnorm = 1/sqrt(pm.Tanα^2 + 1);
M_integrand = (pm.Tanα .* v_RotatedArray_HAB .+ w_RotatedArray_HAB) .* Hnorm

ini_Nisos = 25
ini_Nisos_velocity = 50 # double resolution for velocity contours
Δb = -500*pm.Ñ^2/ini_Nisos
Δb_velocity = -500*pm.Ñ^2/ini_Nisos_velocity
∫B_inb = zeros(ini_Nisos, tlength);
∫M_inb = zeros(ini_Nisos, tlength);
V_inb = zeros(ini_Nisos, tlength);
b_inb = zeros(ini_Nisos, tlength);
vavg_inb = zeros(ini_Nisos_velocity, tlength);
eavg_inb = zeros(ini_Nisos_velocity, tlength);
N2avg_inb = zeros(ini_Nisos_velocity, tlength);
M2avg_inb = zeros(ini_Nisos_velocity, tlength);
∇bavg_inb = zeros(ini_Nisos_velocity, tlength);

for i = 1:tlength
    @info "Time $i of $tlength..."
    # at each time step get the data at (x,ycut,z)
    bt = b_RotatedArray_BelowHAB[:,:, :, i];
    ∇κ∇Bt = ∇κ∇B_RotatedArray_BelowHAB[:,:,:,i];
    M_integrandt = M_integrand[:,:,i];
    vt = vfull_RotatedArray_HAB[:,:,:,i];
    et = e_RotatedArray_HAB[:,:,:,i];
    N2t = N2_RotatedArray_HAB[:,:,:,i];

    # for each buoyancy class:
    for n = 1:ini_Nisos
        # (1) CCC locations in the density class
        # starting with b < 0 should be almost the whole domain
        boolB = (bt.< Δb*(n-1)) .& (bt.>= Δb *n)
        boolB_H =  boolB[:,:,end];
        # fluxes within the isopycnal layer
        ∇κ∇B_inb = ∇κ∇Bt[boolB]
        M_inb = M_integrandt[boolB_H]

        # volume integrated dye concentration:
        ∫B_inb[n,i] = sum(∇κ∇B_inb)*32
        V_inb[n,i] = sum(boolB) * 32
        ∫M_inb[n,i] = sum(M_inb)*16
    end

    for n = 1:ini_Nisos_velocity
        boolB = (bt.< Δb_velocity*(n-1)) .& (bt.>= Δb_velocity *n)
        v_inb = vt[boolB]
        e_inb = et[boolB]
        N2_inb = N2t[boolB]

        vol  = sum(boolB)
        if vol != 0
            vavg_inb[n,i] = sum(v_inb) ./ vol
            eavg_inb[n,i] = sum(e_inb) ./ vol
            N2avg_inb[n,i] = sum(N2_inb) ./ vol
        end
    end

end

#####
#  calculating Derivatives of Volume Budget terms
####

# now we need the buoyancy derivative deriv(A, h, dim_order)
∫A_inb = deriv(∫B_inb, Δb, [1, 2])

# n = 1 in the ini_Nisos dimension is the top of the domain,
# so for A(b₁,t) the first A(b₁,t)[1,:] = 0 since it's a wall,
# the rest of them come from A(b₁,t)[n,:] = A(b,t)[n-1,:]
∫B1_inb = vcat(zeros(tlength)', ∫B_inb[1:end-1,:])
∫A1_inb = deriv(∫B1_inb, Δb, [1, 2])

# now we need the volume changes in time
∂tV_inb = deriv(V_inb, 600.0, [2, 1])

######
# ONLY USING INDICES WHERE VOLUME IS SIMIMILAR
#####

volstart =  2 #10 #33 #findfirst(V_inb[:,1] .> 3e5)
volend = ini_Nisos - 2 #- 15 #40 #findfirst(reverse(V_inb[:,1]) .> 3e5) +1

ini_Nisos_cut = length(volstart:volend)
ini_Nisos_range = volstart:volend
ini_Nisos_range_velocity = volstart*2:volend*2
b_bins_top = (Δb .* (0:ini_Nisos-1))[ini_Nisos_range]
b_bins_bot = (Δb .* (1:ini_Nisos))[ini_Nisos_range]
b_bins_bot_velocity = (Δb_velocity .* (1:ini_Nisos_velocity))[ini_Nisos_range_velocity]

# now can plot a hovmoller of all the terms at specific h values:
waves = b_timeseries.times./pm.Tσ

# dV/dt      A1 - A - m
# A1 - A     - M
# A1         A

Adif = ∂tV_inb .+ ∫M_inb
Adif_true = ∫A1_inb .- ∫A_inb

habround = round(Int, hab)

f1 = Figure(resolution = (1500, 2000), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[1, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[1, 2], ) 

    ax3 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax4 = Axis(ga[2, 2], ) 

    ax5 = Axis(ga[3, 1], xlabel = "Tσ", ylabel = "b [ms⁻²]") #vh
    ax6 = Axis(ga[3, 2], xlabel = "Tσ" ) 

    gcb1 = ga[1, 3] = GridLayout()
    gcb2 = ga[2, 3] = GridLayout()
    gcb3 = ga[3, 3] = GridLayout()

    ax1.yticks = -5e-3:2e-3:-1e-3
    ax3.yticks = -5e-3:2e-3:-1e-3
    ax5.yticks = -5e-3:2e-3:-1e-3

    ax5.xticks = 2:2:10
    ax6.xticks = 2:2:10

    limits!(ax1, 0,11, -6e-3, 0)
    limits!(ax2, 0,11,-6e-3, 0)
    limits!(ax3, 0,11, -6e-3, 0)
    limits!(ax4, 0,11, -6e-3, 0)
    limits!(ax5, 0,11, -6e-3, 0)
    limits!(ax6, 0,11, -6e-3, 0)

    hidedecorations!(ax2)
    hidexdecorations!(ax3)
    hidexdecorations!(ax1)
    hidedecorations!(ax4)
    hideydecorations!(ax6)

    hV = heatmap!(ax1, waves, b_bins_bot, ∂tV_inb[ini_Nisos_range,:]', colormap = :balance, colorrange = (-900,900))
        text!(ax1, Point.(1, -5e-3), text = "∂V/∂t", align = (:left, :center), color = :black, font = :bold, fontsize = 30)

    heatmap!(ax2, waves, b_bins_bot, RHS[ini_Nisos_range,:]', colormap = :balance, colorrange = (-900,900))
        text!(ax2, Point.(1, -5e-3), text = "A(b-Δb) - A(b) - M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 30)

    heatmap!(ax3, waves, b_bins_bot, Adif[ini_Nisos_range,:]', colormap = :balance, colorrange = (-20,20))
        text!(ax3, Point.(1, -5e-3), text = "A(b-Δb) - A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 30)

    hM = heatmap!(ax4, waves, b_bins_bot, -1 .* ∫M_inb[ini_Nisos_range,:]', colormap = :balance, colorrange = (-500,500))
        text!(ax4, Point.(1, -5e-3), text = "- M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 30)
    
    hA = heatmap!(ax5, waves, b_bins_bot, ∫A1_inb[ini_Nisos_range,:]', colormap = :balance, colorrange = (-20,20))
        text!(ax5, Point.(1, -5e-3), text = "A(b-Δb)", align = (:left, :center), color = :black, font = :bold, fontsize = 30)

    heatmap!(ax6, waves, b_bins_bot, -1 .* ∫A_inb[ini_Nisos_range,:]', colormap = :balance, colorrange = (-20,20))
        text!(ax6, Point.(1, -5e-3), text = "- A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 30)

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-450:450:450), size =35, label = "∂V/∂t(b,t) [m³s⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hM, ticks = (-250:250:250), size =35, label = "M(b,t) [m³s⁻¹]")
    cb3 = Colorbar(gcb3[1,1], hA, ticks = (-10:10:10), size =35, label = "A(b,t) [m³s⁻¹]")

    Label(ga[1, 1:2, Top()], "HAB = $habrat δ m, no Corner Sponge",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)


#savename = apath * "VolumeTransBudget_Hovmoller_nb_$ini_Nisos" *"_"* setname 
savename = apath * "VolumeTransBudget_Hovmoller_nb_$(ini_Nisos)full_$(habround)_csp_"* setname 

save(savename * ".png", f1)

#####
#  AVERAGING IN TIME
#####

# try taking a phase average:
include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
endwaves = wave_info.WavePeriods[:,5:11]

∂tV_inb_Wavg = mean(∂tV_inb[ini_Nisos_range,endwaves], dims = (2,3))[:,1,1];
Adif_Wavg = mean(Adif[ini_Nisos_range,endwaves], dims = (2,3))[:,1,1];
∫M_inb_Wavg = mean(∫M_inb[ini_Nisos_range,endwaves], dims = (2,3))[:,1,1];
Adif_true_Wavg = mean(Adif_true[ini_Nisos_range,endwaves], dims = (2,3))[:,1,1];


f1 = Figure(resolution = (800, 900), fontsize=26)
    ga = f1[1, 1] = GridLayout()   
    ax1 = Axis(ga[1, 1], xlabel = "Volume Flux into ℛ̄(b)", ylabel = "b [ms⁻²]") #vh

    ax1.yticks = -5e-3:1e-3:-1e-3

    lines!(ax1, ∂tV_inb_Wavg, b_bins_bot, color = :black, linewidth = 2, label = "∂V/∂t")
    lines!(ax1, -1 .* ∫M_inb_Wavg, b_bins_bot, color = :blue, linewidth = 2, label = "- M(b)")
    lines!(ax1, Adif_Wavg, b_bins_bot, color = :red, linewidth = 2, label = "∂V/∂t + M(b)" )
    lines!(ax1, Adif_true_Wavg, b_bins_bot, color = :gray30, linewidth = 3, label = "A(b+Δb) - A(b)")
    hspan!(ax1, -3e-3, (-3e-3+δ*pm.Ñ^2) , color = (:blue, 0.2), label = "δN²")
    axislegend(ax1, position = :rt)

    Label(ga[1, 1, Top()], "HAB = $habrat δ m, Wave avg Tσ:4-11",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)

    savename = apath * "VolumeBudget_nb_$(ini_Nisos)_$(habround)_Wavg_fin_"* setname 

save(savename * ".png", f1)

#####
# Rolling avg
####
∂tV_inb_Wavg_rbavg = zeros(ini_Nisos_cut -2)
∫M_inb_Wavg_rbavg = zeros(ini_Nisos_cut -2)
Adif_Wavg_rbavg = zeros(ini_Nisos_cut -2)
Adif_true_Wavg_rbavg = zeros(ini_Nisos_cut -2)

for k = 2:ini_Nisos_cut-1
    ∂tV_inb_Wavg_rbavg[k-1] = mean(∂tV_inb_Wavg[k-1:k+1])
    ∫M_inb_Wavg_rbavg[k-1] = mean(∫M_inb_Wavg[k-1:k+1])
    Adif_Wavg_rbavg[k-1] = mean(Adif_Wavg[k-1:k+1])
    Adif_true_Wavg_rbavg[k-1] = mean(Adif_true_Wavg[k-1:k+1])
end

f1 = Figure(resolution = (800, 900), fontsize=26)
    ga = f1[1, 1] = GridLayout()   
    ax1 = Axis(ga[1, 1], xlabel = "Volume Flux into ℛ̄(b)", ylabel = "b [ms⁻²]") #vh

    ax1.yticks = -5e-3:1e-3:-1e-3

    lines!(ax1, ∂tV_inb_Wavg_rbavg, b_bins_bot[2:end-1], color = :black, linewidth = 2, label = "∂V/∂t")
    lines!(ax1, -1 .* ∫M_inb_Wavg_rbavg, b_bins_bot[2:end-1], color = :blue, linewidth = 2, label = "- M(b)")
    lines!(ax1, Adif_Wavg_rbavg, b_bins_bot[2:end-1], color = :red, linewidth = 2, label = "∂V/∂t + M(b)" )
    lines!(ax1, Adif_true_Wavg_rbavg, b_bins_bot[2:end-1], color = :gray30, linewidth = 3, label = "A(b+Δb) - A(b)")
    hspan!(ax1, -3e-3, (-3e-3+δ*pm.Ñ^2) , color = (:blue, 0.2), label = "δN²")
    axislegend(ax1, position = :rt)

    Label(ga[1, 1, Top()], "HAB = $habrat δ m, Wave avg Tσ:4-11",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)

    savename = apath * "VolumeBudget_nb_$(ini_Nisos)_$(habround)_Wavg_rbavg_"* setname 

save(savename * ".png", f1)

∂tV_inb_Phavg = mean(∂tV_inb[:,endwaves], dims = 3)[ini_Nisos_range,:,1];
Adif_Phavg = mean(Adif[:,endwaves], dims = 3)[ini_Nisos_range,:,1];
∫M_inb_Phavg = mean(∫M_inb[:,endwaves], dims = 3)[ini_Nisos_range,:,1];
vavg_inb_Phavg = mean(vavg_inb[:,endwaves], dims = 3)[ini_Nisos_range_velocity,:,1];
Adif_true_Phavg = mean(Adif_true[:,endwaves], dims = 3)[ini_Nisos_range,:,1];
eavg_inb_Phavg = mean(eavg_inb[:,endwaves], dims = 3)[ini_Nisos_range_velocity,:,1];
N2avg_inb_Phavg = mean(N2avg_inb[:,endwaves], dims = 3)[ini_Nisos_range_velocity,:,1] .- pm.Ñ.^ 2;
M2avg_inb_Phavg = mean(M2avg_inb[:,endwaves], dims = 3)[ini_Nisos_range_velocity,:,1] 
∇bavg_inb_Phavg = mean(∇bavg_inb[:,endwaves], dims = 3)[ini_Nisos_range_velocity,:,1] .- pm.Ñ.^ 2;

log_eavg_inb_Phavg = log10.(eavg_inb_Phavg)
phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

#(1500, 700)
f1 = Figure(resolution = (1500, 1400), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[2, 2]) 
    ax3 = Axis(ga[3, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax4 = Axis(ga[3, 2], xlabel = "Tσ" ) 

    gcb1 = ga[1, 1:2] = GridLayout()
    gcb2 = ga[4, 1:2] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3
    ax3.yticks = -5e-3:1e-3:-1e-3

    ax4.xticks = .2:.2:1
    ax3.xticks = .2:.2:1

    hidedecorations!(ax2)
    hidexdecorations!(ax1)
    hideydecorations!(ax4)

    scaling_levelsn = -.1:.1:.1
    scaling_levelsp = 0:.1:.1

    N2scaling_levelsn = -5e-6:5e-6:0
    N2scaling_levelsp = 0:1e-5:1e-5

    hV = heatmap!(ax1, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-600,600))
    vel = contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    text!(ax1, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax1, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax1, Point.(.9, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax1, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax1, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax1, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 15)

    heatmap!(ax2, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-600,600))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    text!(ax2, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax2, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax2, Point.(.9, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax2, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax2, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 15)
    text!(ax2, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 15)

    heatmap!(ax3, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-600,600))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    

    hE = heatmap!(ax4, phase_times, b_bins_bot_velocity, log_eavg_inb_Phavg', colormap = :thermal, colorrange = (-7,-5))
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    contour!(ax4, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg'; levels = N2scaling_levelsp, color = :white , linewidth = 3)    
    contour!(ax4, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg'; levels = N2scaling_levelsn, color = :white , linewidth = 3, linestyle = :dot)    

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-500:250:500), size =35, label = "∂V(b,t)/∂t [m³s⁻¹]", vertical = false)
    cb2 = Colorbar(gcb1[1,2], hV, ticks = (-500:250:500), size =35, label = "- M(b,t) [m³s⁻¹]", vertical = false)
    cb3 = Colorbar(gcb2[1,1], hV, ticks = (-500:250:500), size =35, label = "∂V(b,t)/∂t  + M(b,t)  = A(b+Δb,t) - A(b,t) [m³s⁻¹]", vertical = false, flipaxis = false)
    cb4 = Colorbar(gcb2[1,2], hE, ticks = (-7:1:-5, ["10⁻⁷", "10⁻⁸", "10⁻⁹"]), size =35, label = "ε [m³s⁻²]", vertical = false, flipaxis = false)

    rowsize!(ga, 1, Relative(0.05))
    rowsize!(ga, 4, Relative(0.05))

savename = apath * "VolumeMVelocityDissipN2_Hovmoller_Paper_nb_$(ini_Nisos)_$(habround)_Phavg_"* setname 

save(savename * ".png", f1)


f1 = Figure(resolution = (1500, 2000), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[1, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[1, 2], ) 

    ax3 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax4 = Axis(ga[2, 2], ) 

    ax5 = Axis(ga[3, 1], xlabel = "Tσ", ylabel = "b [ms⁻²]") #vh
    ax6 = Axis(ga[3, 2], xlabel = "Tσ" ) 

    gcb1 = ga[1, 3] = GridLayout()
    gcb2 = ga[2, 3] = GridLayout()
    gcb3 = ga[3, 3] = GridLayout()

    ax1.yticks = -5e-3:2e-3:-1e-3
    ax3.yticks = -5e-3:2e-3:-1e-3
    ax5.yticks = -5e-3:2e-3:-1e-3

    ax5.xticks = .2:.2:1
    ax6.xticks = .2:.2:1

    #limits!(ax1, 0,1, -5e-3, 0)
    #limits!(ax2, 0,1,-5e-3, 0)
    #limits!(ax3, 0,1, -5e-3, 0)
    #limits!(ax4, 0,1, -5e-3, 0)
    #limits!(ax5, 0,1, -6e-3, 0)
    #limits!(ax6, 0,1, -6e-3, 0)

    hidedecorations!(ax2)
    hidexdecorations!(ax3)
    hidexdecorations!(ax1)
    hidedecorations!(ax4)
    hideydecorations!(ax6)

    hV = heatmap!(ax1, phase_times, b_bins_bot, ∂tV_inb_Phavg[9:end-10,:]', colormap = :balance, colorrange = (-400,400))
        text!(ax1, Point.(.1, -5e-3), text = "∂V/∂t", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax2, phase_times, b_bins_bot, RHS_Phavg[9:end-10,:]', colormap = :balance, colorrange = (-400,400))
        text!(ax2, Point.(.1, -5e-3), text = "A(b-Δb) - A(b) - M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax3, phase_times, b_bins_bot, Adif_Phavg[9:end-10,:]', colormap = :balance, colorrange = (-8,8))
        text!(ax3, Point.(.1, -5e-3), text = "A(b-Δb) - A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    hM = heatmap!(ax4, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg[9:end-10,:]', colormap = :balance, colorrange = (-400,400))
        text!(ax4, Point.(.1, -5e-3), text = "- M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    
    hA = heatmap!(ax5, phase_times, b_bins_bot, ∫A1_inb_Phavg[9:end-10,:]', colormap = :balance, colorrange = (-8,8))
        text!(ax5, Point.(.1, -5e-3), text = "A(b-Δb)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax6, phase_times, b_bins_bot, -1 .* ∫A_inb_Phavg[9:end-10,:]', colormap = :balance, colorrange = (-8,8))
        text!(ax6, Point.(.1, -5e-3), text = "- A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-200:200:200), size =35, label = "∂V/∂t(b,t) [m³s⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hM, ticks = (-200:200:200), size =35, label = "M(b,t) [m³s⁻¹]")
    cb3 = Colorbar(gcb3[1,1], hA, ticks = (-4:4:4), size =35, label = "A(b,t) [m³s⁻¹]")

    Label(ga[1, 1:2, Top()], "Phase Averaged Tσ:6-11, HAB = $habrat δ m, w/ Corner Sponge",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)

#savename = apath * "VolumeTransBudget_Hovmoller_nb_$ini_Nisos" *"_Phavg_"* setname 
savename = apath * "VolumeTransBudget_Hovmoller_nb_$(ini_Nisos)full_$(habround)_Phavg_csp_"* setname 

save(savename * ".png", f1)

f1 = Figure(resolution = (1500, 2000), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[1, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[1, 2], ) 

    ax3 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax4 = Axis(ga[2, 2], ) 

    ax5 = Axis(ga[3, 1], xlabel = "Tσ", ylabel = "b [ms⁻²]") #vh
    ax6 = Axis(ga[3, 2], xlabel = "Tσ" ) 

    gcb1 = ga[1, 3] = GridLayout()
    gcb2 = ga[2, 3] = GridLayout()
    gcb3 = ga[3, 3] = GridLayout()

    ax1.yticks = -5e-3:2e-3:-1e-3
    ax3.yticks = -5e-3:2e-3:-1e-3
    ax5.yticks = -5e-3:2e-3:-1e-3

    ax5.xticks = .2:.2:1
    ax6.xticks = .2:.2:1

    #limits!(ax1, 0,1, -5e-3, 0)
    #limits!(ax2, 0,1,-5e-3, 0)
    #limits!(ax3, 0,1, -5e-3, 0)
    #limits!(ax4, 0,1, -5e-3, 0)
    #limits!(ax5, 0,1, -6e-3, 0)
    #limits!(ax6, 0,1, -6e-3, 0)

    hidedecorations!(ax2)
    hidexdecorations!(ax3)
    hidexdecorations!(ax1)
    hidedecorations!(ax4)
    hideydecorations!(ax6)

    hV = heatmap!(ax1, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-400,400))
        text!(ax1, Point.(.1, -5e-3), text = "∂V/∂t", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax2, phase_times, b_bins_bot, RHS_Phavg', colormap = :balance, colorrange = (-400,400))
        text!(ax2, Point.(.1, -5e-3), text = "A(b-Δb) - A(b) - M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax3, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-8,8))
        text!(ax3, Point.(.1, -5e-3), text = "A(b-Δb) - A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    hM = heatmap!(ax4, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-400,400))
        text!(ax4, Point.(.1, -5e-3), text = "- M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    
    hA = heatmap!(ax5, phase_times, b_bins_bot, -1 .* ∫A_inb_Phavg', colormap = :balance, colorrange = (-8,8))
        text!(ax5, Point.(.1, -5e-3), text = "-A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax6, phase_times, b_bins_bot,  -1 .* ∫Mwave_inb_Phavg', colormap = :balance, colorrange = (-400,400))
        text!(ax6, Point.(.1, -5e-3), text = "- Mw(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-200:200:200), size =35, label = "∂V/∂t(b,t) [m³s⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hM, ticks = (-200:200:200), size =35, label = "M(b,t) [m³s⁻¹]")
    cb3 = Colorbar(gcb3[1,1], hA, ticks = (-4:4:4), size =35, label = "A(b,t) [m³s⁻¹]")

    Label(ga[1, 1:2, Top()], "Phase Averaged Tσ:6-11, HAB = $habrat δ m, w/ Corner Sponge",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)


savename = apath * "VolumeTransBudget_Mwave_Hovmoller_nb_$(ini_Nisos)full_$(habround)_Phavg_csp_"* setname 

save(savename * ".png", f1)

begwaves = wave_info.WavePeriods[:,2:3]

∂tV_inb_Phavg = mean(∂tV_inb[:,begwaves], dims = 3)[:,:,1];
RHS_Phavg = mean(RHS[:,begwaves], dims = 3)[:,:,1];
Adif_Phavg = mean(Adif[:,begwaves], dims = 3)[:,:,1];
∫M_inb_Phavg = mean(∫M_inb[:,begwaves], dims = 3)[:,:,1];
∂b∫A1_inb_Phavg = mean(∂b∫A1_inb[:,begwaves], dims = 3)[:,:,1];
∂b∫A_inb_Phavg = mean(∂b∫A_inb[:,begwaves], dims = 3)[:,:,1];

phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

f1 = Figure(resolution = (1500, 2000), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[1, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[1, 2], ) 

    ax3 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax4 = Axis(ga[2, 2], ) 

    ax5 = Axis(ga[3, 1], xlabel = "Tσ", ylabel = "b [ms⁻²]") #vh
    ax6 = Axis(ga[3, 2], xlabel = "Tσ" ) 

    gcb1 = ga[1, 3] = GridLayout()
    gcb2 = ga[2, 3] = GridLayout()
    gcb3 = ga[3, 3] = GridLayout()

    ax1.yticks = -5e-3:2e-3:-1e-3
    ax3.yticks = -5e-3:2e-3:-1e-3
    ax5.yticks = -5e-3:2e-3:-1e-3

    ax5.xticks = .2:.2:1
    ax6.xticks = .2:.2:1

    limits!(ax1, 0,1, -6e-3, 0)
    limits!(ax2, 0,1,-6e-3, 0)
    limits!(ax3, 0,1, -6e-3, 0)
    limits!(ax4, 0,1, -6e-3, 0)
    limits!(ax5, 0,1, -6e-3, 0)
    limits!(ax6, 0,1, -6e-3, 0)

    hidedecorations!(ax2)
    hidexdecorations!(ax3)
    hidexdecorations!(ax1)
    hidedecorations!(ax4)
    hideydecorations!(ax6)

    hV = heatmap!(ax1, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-400,400))
        text!(ax1, Point.(.1, -5e-3), text = "∂V/∂t", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax2, phase_times, b_bins_bot, RHS_Phavg', colormap = :balance, colorrange = (-400,400))
        text!(ax2, Point.(.1, -5e-3), text = "A(b-Δb) - A(b) - M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax3, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-8,8))
        text!(ax3, Point.(.1, -5e-3), text = "A(b-Δb) - A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    hM = heatmap!(ax4, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-400,400))
        text!(ax4, Point.(.1, -5e-3), text = "- M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    
    hA = heatmap!(ax5, phase_times, b_bins_bot, ∂b∫A1_inb_Phavg', colormap = :balance, colorrange = (-8,8))
        text!(ax5, Point.(.1, -5e-3), text = "A(b-Δb)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax6, phase_times, b_bins_bot, -1 .* ∂b∫A_inb_Phavg', colormap = :balance, colorrange = (-8,8))
        text!(ax6, Point.(.1, -5e-3), text = "- A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-200:200:200), size =35, label = "∂V/∂t(b,t) [m³s⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hM, ticks = (-200:200:200), size =35, label = "M(b,t) [m³s⁻¹]")
    cb3 = Colorbar(gcb3[1,1], hA, ticks = (-4:4:4), size =35, label = "A(b,t) [m³s⁻¹]")

    Label(ga[1, 1:2, Top()], "Phase Averaged Tσ:1-3, HAB = $habrat δ m, no Corner Sponge",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)

#savename = apath * "VolumeTransBudget_Hovmoller_nb_$ini_Nisos" *"_Phavg_"* setname 
savename = apath * "VolumeTransBudget_Hovmoller_nb_$(ini_Nisos)_$(habround)_Phavg_early_"* setname 
save(savename * ".png", f1)

rWindow = 15
rWavg_length = tlength - 2*rWindow

# try taking a rolling average:
∂tV_inb_rWavg = zeros(ini_Nisos,rWavg_length);
RHS_rWavg = zeros(ini_Nisos,rWavg_length);
Adif_rWavg = zeros(ini_Nisos,rWavg_length);
∫M_inb_rWavg = zeros(ini_Nisos,rWavg_length);
∂b∫A1_inb_rWavg = zeros(ini_Nisos,rWavg_length);
∂b∫A_inb_rWavg = zeros(ini_Nisos,rWavg_length);

for (idx,i) in enumerate(rWindow+1:tlength-rWindow)
    ∂tV_inb_rWavg[:,idx] = mean(∂tV_inb[:,i-rWindow:i+rWindow],dims = 2)
    RHS_rWavg[:,idx] = mean(RHS[:,i-rWindow:i+rWindow],dims = 2)
    Adif_rWavg[:,idx] = mean(Adif[:,i-rWindow:i+rWindow],dims = 2)
    ∫M_inb_rWavg[:,idx] =mean(∫M_inb[:,i-rWindow:i+rWindow],dims = 2)
    ∂b∫A1_inb_rWavg[:,idx] = mean(∂b∫A1_inb[:,i-rWindow:i+rWindow],dims = 2)
    ∂b∫A_inb_rWavg[:,idx] = mean(∂b∫A_inb[:,i-rWindow:i+rWindow],dims = 2)

end

rollwaves = waves[rWindow+1:tlength-rWindow]

f1 = Figure(resolution = (1500, 2000), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[1, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[1, 2], ) 

    ax3 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax4 = Axis(ga[2, 2], ) 

    ax5 = Axis(ga[3, 1], xlabel = "Tσ", ylabel = "b [ms⁻²]") #vh
    ax6 = Axis(ga[3, 2], xlabel = "Tσ" ) 

    gcb1 = ga[1, 3] = GridLayout()
    gcb2 = ga[2, 3] = GridLayout()
    gcb3 = ga[3, 3] = GridLayout()

    ax1.yticks = -5e-3:2e-3:-1e-3
    ax3.yticks = -5e-3:2e-3:-1e-3
    ax5.yticks = -5e-3:2e-3:-1e-3

    ax5.xticks = 2:2:10
    ax6.xticks = 2:2:10

    #limits!(ax1, 1,10, -6e-3, 0)
    #limits!(ax2, 1,10,-6e-3, 0)
    #limits!(ax3, 1,10, -6e-3, 0)
    #limits!(ax4, 1,10, -6e-3, 0)
    #limits!(ax5, 1,10, -6e-3, 0)
    #limits!(ax6, 1,10, -6e-3, 0)

    hidexdecorations!(ax1)
    hidedecorations!(ax2)
    hidexdecorations!(ax3)
    hidedecorations!(ax4)
    hideydecorations!(ax6)

    hV = heatmap!(ax1, rollwaves, b_bins_bot, ∂tV_inb_rWavg[9:end-10,:]', colormap = :balance, colorrange = (-50,50))
        text!(ax1, Point.(2, -4e-3), text = "∂V/∂t", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax2, rollwaves, b_bins_bot, RHS_rWavg[9:end-10,:]', colormap = :balance, colorrange = (-50,50))
        text!(ax2, Point.(2, -4e-3), text = "A(b-Δb) - A(b) - M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax3, rollwaves, b_bins_bot, Adif_rWavg[9:end-10,:]', colormap = :balance, colorrange = (-4,4))
        text!(ax3, Point.(2, -4e-3), text = "A(b-Δb) - A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    hM = heatmap!(ax4, rollwaves, b_bins_bot, -1 .* ∫M_inb_rWavg[9:end-10,:]', colormap = :balance, colorrange = (-50,50))
        text!(ax4, Point.(2, -4e-3), text = "- M(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    
    hA = heatmap!(ax5, rollwaves, b_bins_bot, ∂b∫A1_inb_rWavg[9:end-10,:]', colormap = :balance, colorrange = (-4,4))
        text!(ax5, Point.(2, -4e-3), text = "A(b-Δb)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax6, rollwaves, b_bins_bot, -1 .* ∂b∫A_inb_rWavg[9:end-10,:]', colormap = :balance, colorrange = (-4,4))
        text!(ax6, Point.(2, -4e-3), text = "- A(b)", align = (:left, :center), color = :black, font = :bold, fontsize = 26)

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-25:25:25), size =35, label = "∂V/∂t(b,t) [m³s⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hM, ticks = (-25:25:25), size =35, label = "M(b,t) [m³s⁻¹]")
    cb3 = Colorbar(gcb3[1,1], hA, ticks = (-2:2:2), size =35, label = "A(b,t) [m³s⁻¹]")

    Label(ga[1, 1:2, Top()], "Rolling Wave Averaged Window: 2 Tσ, HAB = $habrat δ m, no Corner Sponge",
    fontsize = 30,
    font = :bold,
    padding = (5, 5, 5, 5),
    halign = :center)

savename = apath * "VolumeTransBudget_Hovmoller_nb_$(ini_Nisos)_$(habround)_rWavg_cut_"* setname 

save(savename * ".png", f1)
