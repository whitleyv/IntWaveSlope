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

b_timeseries = FieldTimeSeries(filepath1,"b");
v_timeseries = FieldTimeSeries(filepath1,"v");
w_timeseries = FieldTimeSeries(filepath1,"w");

ylength = 650 #500

bi = interior(b_timeseries)[:,1:ylength,:,:];
vi = interior(v_timeseries)[:,1:ylength,:,:];
wi = interior(w_timeseries)[:,1:ylength,:,:];

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
v_RotatedArray_HAB = zeros(xlength, yidx_length, tlength);
w_RotatedArray_HAB = zeros(xlength, yidx_length, tlength);

for k = 1:yidx_length
    zdx = zidx[k] # each z index 

    for j = 1:HABidxL
        ydx = Idx_Array_y[k,j]

        b_RotatedArray_BelowHAB[:,k,j, :] = bi[:,ydx, zdx, :]

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

ini_Nisos = 50
Δb = -500*pm.Ñ^2/ini_Nisos
∫M_inb = zeros(ini_Nisos, tlength);
V_inb = zeros(ini_Nisos, tlength);

for i = 1:tlength
    @info "Time $i of $tlength..."
    # at each time step get the data at (x,ycut,z)
    bt = b_RotatedArray_BelowHAB[:,:, :, i];
    M_integrandt = M_integrand[:,:,i];

    # for each buoyancy class:
    for n = 1:ini_Nisos
        # (1) CCC locations in the density class
        # starting with b < 0 should be almost the whole domain
        boolB = (bt.< Δb*(n-1)) .& (bt.>= Δb *n)
        boolB_H =  boolB[:,:,end];
        # fluxes within the isopycnal layer
        M_inb = M_integrandt[boolB_H]

        # volume integrated dye concentration:
        V_inb[n,i] = sum(boolB) * 32
        ∫M_inb[n,i] = sum(M_inb)*16
    end

end

#####
#  calculating Derivatives of Volume Budget terms
####

# now we need the volume changes in time
∂tV_inb = deriv(V_inb, 600.0, [2, 1])
Adif = ∂tV_inb .+ ∫M_inb
######
# ONLY USING INDICES WHERE VOLUME IS SIMIMILAR
#####

volstart =  4 #10 #33 #findfirst(V_inb[:,1] .> 3e5)
volend = ini_Nisos - 4 #- 15 #40 #findfirst(reverse(V_inb[:,1]) .> 3e5) +1

ini_Nisos_cut = length(volstart:volend)
ini_Nisos_range = volstart:volend
b_bins_top = (Δb .* (0:ini_Nisos-1))[ini_Nisos_range]
b_bins_bot = (Δb .* (1:ini_Nisos))[ini_Nisos_range]

habround = round(Int, hab)

#####
#  AVERAGING IN TIME
#####

# try taking a wave average:
include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
endwaves = wave_info.WavePeriods[:,5:11]

#####
# Rolling avg
####
∂tV_inb_Wavg_rbavg = zeros(ini_Nisos_cut -2)
∫M_inb_Wavg_rbavg = zeros(ini_Nisos_cut -2)
Adif_Wavg_rbavg = zeros(ini_Nisos_cut -2)

for k = 2:ini_Nisos_cut-1
    ∂tV_inb_Wavg_rbavg[k-1] = mean(∂tV_inb_Wavg[k-1:k+1])
    ∫M_inb_Wavg_rbavg[k-1] = mean(∫M_inb_Wavg[k-1:k+1])
    Adif_Wavg_rbavg[k-1] = mean(Adif_Wavg[k-1:k+1])
end

f1 = Figure(resolution = (700, 900), fontsize=26)
    ga = f1[1, 1] = GridLayout()   
    ax1 = Axis(ga[1, 1], xlabel = "Volume Flux into ℛ̄(b)", ylabel = "b [ms⁻²]") #vh

    ax1.yticks = -5e-3:1e-3:-1e-3

    lines!(ax1, ∂tV_inb_Wavg_rbavg, b_bins_bot[2:end-1], color = :black, linewidth = 4, label = "∂V/∂t")
    lines!(ax1, -1 .* ∫M_inb_Wavg_rbavg, b_bins_bot[2:end-1], color = :dodgerblue2, linewidth = 4, label = "- M(b)")
    lines!(ax1, Adif_Wavg_rbavg, b_bins_bot[2:end-1], color = :firebrick2, linewidth = 4, label = "∂V/∂t + M(b)" )
    hspan!(ax1, -4e-3, (-4e-3+δ*pm.Ñ^2) , color = (:blue, 0.2), label = "δN²")
    axislegend(ax1, position = :rt, labelsize = 30)

    savename = apath * "VolumeBudget_nb_$(ini_Nisos)_$(habround)_Wavg_rbavg_"* setname 

save(savename * ".png", f1, px_per_unit = 2)

