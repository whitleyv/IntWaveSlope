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
v_timeseries = FieldTimeSeries(filepath,"v");
w_timeseries = FieldTimeSeries(filepath,"w");

xb, yb, zb = nodes(b_timeseries) #CCC
 
ylength = 880 
tlength = 161
xlength = 38
land = curvedslope.(yb[1:ylength])

bi = interior(b_timeseries)[:,1:ylength,:,:];
vi = interior(v_timeseries)[:,1:ylength,:,:];
wi = interior(w_timeseries)[:,1:ylength,:,:];

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

function rotate_indices(habrat, δ)
    hab = habrat* δ
    δy = hab/pm.Tanα

    HABidxL = round(Int, δy/4)
    yidx_start = 2
    yidx_end = 334
    yidx_length = length(yidx_start:yidx_end)
    
    land = curvedslope.(yb[yidx_start:yidx_end])
    zidx = sum( zb .< land', dims = 1)[1,:] .+ 1
    Idx_Array_z = reshape(repeat(zidx, HABidxL), yidx_length, HABidxL)
    Idx_Array_y = reshape(repeat(yidx_start:yidx_end, HABidxL), yidx_length, HABidxL) .+ (0:HABidxL-1)'
    return (Idx_Array_y, Idx_Array_z, zidx, yidx_length, HABidxL)
end

function rotate_array(b, size_tuple, Idx_Array_y, zidx)
    _, yidx_length, HABidxL, _ = size_tuple

    b_RotatedArray_BelowHAB = zeros(size_tuple);

    for k = 1:yidx_length
        zdx = zidx[k] # each z index 
        for j = 1:HABidxL
            ydx = Idx_Array_y[k,j]
            b_RotatedArray_BelowHAB[:,k,j, :] = b[:,ydx, zdx, :]
        end
    end

    return b_RotatedArray_BelowHAB


end

function rotate_array(b, v, w, size_tuple, Idx_Array_y, zidx)
    xlength, yidx_length, HABidxL, tlength = size_tuple
    b_RotatedArray_BelowHAB = zeros(size_tuple);

    v_RotatedArray_HAB = zeros(xlength, yidx_length, tlength);
    w_RotatedArray_HAB = zeros(xlength, yidx_length, tlength);
    
    for k = 1:yidx_length
        zdx = zidx[k] # each z index 
        for j = 1:HABidxL
            ydx = Idx_Array_y[k,j]
            b_RotatedArray_BelowHAB[:,k,j, :] = b[:,ydx, zdx, :]
        end
        Hidx_y = Idx_Array_y[k, end]
        # create arrays of only the values at HAB:
        v_RotatedArray_HAB[:,k, :] = v[:,Hidx_y, zdx, :]
        w_RotatedArray_HAB[:,k, :] = w[:,Hidx_y, zdx, :]

    end

    return b_RotatedArray_BelowHAB, v_RotatedArray_HAB, w_RotatedArray_HAB


end
#=
function volume_buoyancy_binning(b_Rotated_Wavg, b0, ini_Nisos, Δb, inirange)
    V_inb_waves = zeros(ini_Nisos);
    V_inb_waves0 = zeros(ini_Nisos);

    # for each buoyancy class:
    for n = 1:ini_Nisos
        boolB = (b_Rotated_Wavg .< Δb*(n-1)) .& (b_Rotated_Wavg .>= Δb *n)
        boolB0 = (b0.< Δb*(n-1)) .& (b0.>= Δb *n)

        # volume integrated dye concentration:
        V_inb_waves[n] = sum(boolB) * 32
        V_inb_waves0[n] = sum(boolB0) * 32

    end
    b_bins = (Δb .* (1:ini_Nisos))[inirange]
    ΔV_inb_waves = (V_inb_waves .- V_inb_waves0)[inirange]


    return V_inb_waves, V_inb_waves0, ΔV_inb_waves, b_bins
end 
=#
function volume_buoyancy_binning(b_Rotated, ini_Nisos, Δb, inirange, tlength)
    V_inb = zeros(ini_Nisos, tlength);

    for i = 1:tlength
     #   @info "Time $i of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bt = b_Rotated[:,:, :, i];
    
        # for each buoyancy class:
        for n = 1:ini_Nisos
            # (1) CCC locations in the density class
            # starting with b < 0 should be almost the whole domain
            boolB = (bt.< Δb*(n-1)) .& (bt.>= Δb *n)

            # volume integrated dye concentration:
            V_inb[n,i] = sum(boolB) * 32
        end
    
    end

    b_bins = (Δb .* (1:ini_Nisos))[inirange]
    ΔV_inb = (V_inb .- V_inb[:,1])[inirange,:]


    return ΔV_inb, b_bins
end

function volume_buoyancy_binning(b_Rotated, M_integrand, ini_Nisos, Δb, inirange, tlength)
    V_inb = zeros(ini_Nisos, tlength);
    ∫M_inb = zeros(ini_Nisos, tlength);

    for i = 1:tlength
     #   @info "Time $i of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bt = b_Rotated[:,:, :, i];
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

    b_bins = (Δb .* (1:ini_Nisos))[inirange]
    

    return V_inb[inirange,:], ∫M_inb[inirange,:], b_bins
end

#########
#       FULL VOLUME OUT TO 5 DELTA and 1.1 DELTA
########
 
δ = pm.U₀/pm.Ñ
α = atan(pm.Tanα)
δp = δ/cos(α)

(Idx_Array_y5, Idx_Array_z5, zidx5, yidx_length5, HABidxL5) = rotate_indices(5, δ)
(Idx_Array_y1, Idx_Array_z1, zidx1, yidx_length1, HABidxL1) = rotate_indices(1.1, δ)

#########
#       ROTATE ARRAY
########

b_RotatedArray_BelowHAB5 = rotate_array(bi, (xlength, yidx_length5, HABidxL5, tlength), Idx_Array_y5, zidx5)
b_RotatedArray_BelowHAB1, v_RotatedArray_HAB1, w_RotatedArray_HAB1 = rotate_array(bi, vi, wi, (xlength, yidx_length1, HABidxL1, tlength), Idx_Array_y1, zidx1)

Hnorm = 1/sqrt(pm.Tanα^2 + 1);
M_integrand = (pm.Tanα .* v_RotatedArray_HAB1 .+ w_RotatedArray_HAB1) .* Hnorm

#######
#      WAVE AVERAGING
#######

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)

waves511 = wave_info.WavePeriods[:,6:8]#9:11]

b_RotatedArray_BelowHAB5_Wavg511 = mean(b_RotatedArray_BelowHAB5[:,:,:,waves511], dims = (4,5))[:,:,:,1,1];
b0 = b_RotatedArray_BelowHAB5[:,:,:,1]

#########
#       BINNING BUOYANCY
########

ini_Nisos = 50 #60#50
Δb = -500*pm.Ñ^2/ini_Nisos
inirange = 4:ini_Nisos - 4

_, _, ΔV_Wavg511_inb5, b_bins5 = volume_buoyancy_binning(b_RotatedArray_BelowHAB5_Wavg511, b0, ini_Nisos, Δb, inirange)

ΔV_inb5, b_bins5 = volume_buoyancy_binning(b_RotatedArray_BelowHAB5, ini_Nisos, Δb, inirange, tlength)

ΔV_inb5_dt = permutedims(permutedims(ΔV_inb5, [2,1]) ./ b_timeseries.times, [2,1])


#########
#       BINNING BUOYANCY
########

ini_Nisos = 50
Δb = -500*pm.Ñ^2/ini_Nisos
volstart =  4 #10 #33 #findfirst(V_inb[:,1] .> 3e5)
volend = ini_Nisos - 4 #- 15 #40 #findfirst(reverse(V_inb[:,1]) .> 3e5) +1
ini_Nisos_range = volstart:volend
ini_Nisos_cut = length(volstart:volend)

V_inb1, ∫M_inb1, b_bins1 = volume_buoyancy_binning(b_RotatedArray_BelowHAB1, M_integrand, ini_Nisos, Δb, ini_Nisos_range, tlength)

# now we need the volume changes in time
∂tV_inb1 = deriv(V_inb1, 600.0, [2, 1])
Adif1 = ∂tV_inb1 .+ ∫M_inb1

########
#      WAVE AVERAGING
#######
include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)

waves511 = wave_info.WavePeriods[:,5:11]#9:11]

∂tV_inb_Wavg = mean(∂tV_inb1[:,waves511], dims = (2,3))[:,1,1];
Adif_Wavg = mean(Adif1[:,waves511], dims = (2,3))[:,1,1];
∫M_inb_Wavg = mean(∫M_inb1[:,waves511], dims = (2,3))[:,1,1];

#waves511 = wave_info.WavePeriods[:,6:8]#9:11]
ΔV_inb5_Wavg = mean(ΔV_inb5[:,waves511], dims = (2,3))[:,1,1];
ΔV_inb5_dt_Wavg = mean(ΔV_inb5_dt[:,waves511], dims = (2,3))[:,1,1];

#####
# Rolling buoyancy bin avg
####
ΔV_inb5_Wavg_rbavg = zeros(ini_Nisos_cut -2)
ΔV_inb5_dt_Wavg_rbavg = zeros(ini_Nisos_cut -2)
∂tV_inb_Wavg_rbavg = zeros(ini_Nisos_cut -2)
∫M_inb_Wavg_rbavg = zeros(ini_Nisos_cut -2)
Adif_Wavg_rbavg = zeros(ini_Nisos_cut -2)

for k = 2:ini_Nisos_cut-1
    ∂tV_inb_Wavg_rbavg[k-1] = mean(∂tV_inb_Wavg[k-1:k+1])
    ∫M_inb_Wavg_rbavg[k-1] = mean(∫M_inb_Wavg[k-1:k+1])
    Adif_Wavg_rbavg[k-1] = mean(Adif_Wavg[k-1:k+1])
    ΔV_inb5_Wavg_rbavg[k-1] = mean(ΔV_inb5_Wavg[k-1:k+1])
    ΔV_inb5_dt_Wavg_rbavg[k-1] = mean(ΔV_inb5_dt_Wavg[k-1:k+1])

end

########
#      BUOYANCY CHANGE IN PHSYICAL SPACE
#######

# putting these in wave form
b_waves511 = bi[:,:,:,waves511];

# averaging over waves
b_Wavg511_xavg = mean(b_waves511, dims = (1,4,5))[1,:,:,1,1];
b₀ = bi[19,:,:, 1];

Δb_Wavg511_xavg = b_Wavg511_xavg .- b₀;

topo = curvedslope.(yb[1:ylength])

f1 = Figure(resolution = (2500, 1000), fontsize=30)
    ga = f1[1, 1] = GridLayout()    
    gb = f1[1, 2] = GridLayout()    

    axb = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]") 
    gcb1 = ga[1, 1] = GridLayout()

    axv = Axis(gb[2,1], ylabel = "b [ms⁻²]", xlabel = "ΔV/Δt = (V(b,t) - V(b,0)) / t [m³s⁻¹]")
    axvsm = Axis(gb[2,2], xlabel = "Volume Flux into R(b) [m³s⁻¹]")

    axb.xticks = 500:500:2500
    axb.yticks = -500:250:0
    limits!(axb, 0, 2500, -500, 0)

    b_bins_coarse = -5e-3:1e-3:-1e-3
    yb_cut = yb[1:ylength]
    Δblims = (-4.5e-4, 4.5e-4)

    hb = heatmap!(axb, yb_cut, zb, Δb_Wavg511_xavg, colormap = :balance, colorrange = Δblims)
        contour!(axb, yb_cut, zb, b₀, color = :black, linewidth = 3, levels =b_bins_coarse)
        contour!(axb, yb_cut, zb, b_Wavg511_xavg, color = :gray35, linewidth = 5, levels =b_bins_coarse) #bmin:0.0005:0, alpha = 0.5)
        lines!(axb, yb_cut, linslope.(yb_cut) .+ pm.U₀/pm.Ñ, color=:black, linewidth = 3, linestyle = :dash)
        lines!(axb, yb_cut, topo, color=:black, linewidth = 4)
        d5 = lines!(axb, yb_cut, linslope.(yb_cut) .+ 5* δ, color=:black, linewidth = 7)
        d1 = lines!(axb, yb_cut, linslope.(yb_cut) .+ 1.1* δ, color=:black, linewidth = 7, linestyle=:dot)

    # create colorbars the size of the whole data set #label = "Δb = b - b₀ [ms⁻²]
    cb1 = Colorbar(gcb1[1,1], hb, ticks = (-4e-4:2e-4:4e-4), size =35, vertical = false) 
    cb1.alignmode = Mixed(top = 0)

    axv.yticks = -6e-3:1e-3:0
    #axv.xticks = -1e6:1e6:1e6
    #axv.xticks = -5:5:15
    axv.xticks = -20:20:20
    axvsm.xticks = -20:20:20
    hideydecorations!(axvsm, grid = false)

    #limits!(axv, -1.2e6, 1.2e6, -5.5e-3, -5e-4)
    #limits!(axv, -6, 20, -5.5e-3, -5e-4)
    limits!(axv,-25, 40, -5.5e-3, -5e-4)
    limits!(axvsm, -25, 40, -5.5e-3, -5e-4)

    v11 = lines!(axv, ΔV_inb5_dt_Wavg_rbavg, b_bins5[2:end-1], color = :black, linewidth = 6, label = "t = 8-11 Tσ")

    vsm = lines!(axvsm, ∂tV_inb_Wavg_rbavg, b_bins1[2:end-1], color = :darkgreen, linewidth = 6, label = "∂V/∂t")
    lines!(axvsm, -1 .* ∫M_inb_Wavg_rbavg, b_bins1[2:end-1], color = :dodgerblue2, linewidth = 6, label = "- M(b)")
    lines!(axvsm, Adif_Wavg_rbavg, b_bins1[2:end-1], color = :firebrick2, linewidth = 6, label = "A(b+Δb) - A(b)" )
    
    axislegend(axvsm, position = :rt)

    rowsize!(gb, 1, Relative(0.1))
    rowsize!(ga, 1, Relative(0.1))
    #colsize!(f1.layout, 2, Relative(0.35))
    Legend( gb[1, 1], [d5], [rich("R(b) out to 5 h", subscript("w"))], linewidth=4,
                    tellheight = false, tellwidth = false, framevisible = false, 
                    "Boundary and Interior Region",
                    halign = :center, valign = :center, orientation = :horizontal)

    Legend( gb[1, 2], [d1], [rich("R(b) out to 1.1 h", subscript("w"))],linewidth=4,
                    "Near Boundary Region only",
                    tellheight = false, tellwidth = false, framevisible = false,
                    halign = :center, valign = :center, orientation = :horizontal)

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
    Label(gb[1, 2, TopLeft()], "c",
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :left)

    Label(ga[1, 1, Top()], "Δb = b - b₀ [ms⁻²]",
                    fontsize = 30, 
                    padding = (5, 5, 5, 10),
                    halign = :center)

    Label(gb[1, 1, Top()], "Total Volume Change",
                    fontsize = 30, font = :bold,
                    padding = (5, 5, 5, 10),
                    halign = :center)
    Label(gb[1, 2, Top()], "Instantaneous Volume Change",
                    fontsize = 30,font = :bold,
                    padding = (5, 5, 5, 10),
                    halign = :center)

    #colgap!(ga, 15)
    #rowgap!(ga, 15)
    colgap!(gb, 5)
    rowgap!(gb, 15)

save(path_name * "Analysis/b_Wavg_intime_Vol_wFluxes_dt_" * sn * ".png", f1, px_per_unit = 2)

