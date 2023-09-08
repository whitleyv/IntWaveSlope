using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using ArgParse

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
filepath2 = path_name * "IntWave_smdt_" * sn * ".jld2"

@info "getting data from: " * sn

b_timeseries = FieldTimeSeries(filepath, "b");
v_timeseries = FieldTimeSeries(filepath, "v");
w_timeseries = FieldTimeSeries(filepath, "w");
e_timeseries = FieldTimeSeries(filepath,"ϵ");
Cg_timeseries = FieldTimeSeries(filepath2,"Cg");

xb, yb, zb = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC
xw, yw, zw = nodes(w_timeseries) #CCF

ylength = 880 
tlength = 161

ci = interior(Cg_timeseries)[:, 1:ylength,:,1:2:end];
bi = interior(b_timeseries)[:,1:ylength,:,:];
vi = interior(v_timeseries)[:,1:ylength+1,:,:];
wi = interior(w_timeseries)[:,1:ylength,:,:];
ei = interior(e_timeseries)[:,1:ylength,:,:];

α = atan(pm.Tanα)

@info "interpolating..."

# interpolate to center
v_ccc = 0.5 .* (vi[:,1:end-1,:,:] .+ vi[:,2:end,:,:]);
w_ccc = 0.5 .* (wi[:,:,1:end-1,:] .+ wi[:,:,2:end,:]);

v̂ = v_ccc .* cos(α) .- w_ccc .* sin(α)
ŵ = v_ccc .* sin(α) .+ w_ccc .* cos(α)

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl


# array Wl x (7-11)nTσ
#__Period6 ____Period7_____Period8____...
#__ + 1/15
#__ + 2/15
#__ + 3/15
#    ...
#__ + 15/15
end5waves = wave_info.WavePeriods[:,7:end] # t = 7T : 11 T
beg3waves = wave_info.WavePeriods[:,2:4] # t = 1T : t = 4T
trans3waves = wave_info.WavePeriods[:,5:7] # t = 4T : t = 7T

# reorienting into phases
# x y z Ph W
function phase_orient(indexing, bi, v_ccc, w_ccc, v̂, ŵ, ei, ci)
    b_Ph = bi[:,:,:,indexing]
    v_Ph = v_ccc[:,:,:,indexing]
    w_Ph = w_ccc[:,:,:,indexing]
    v̂_Ph = v̂[:,:,:,indexing]
    ŵ_Ph = ŵ[:,:,:,indexing]
    e_Ph = ei[:,:,:, indexing]
    c_Ph = ci[:,:,:,indexing]
    return b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph, e_Ph, c_Ph
end

b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph, e_Ph, c_Ph  = phase_orient(end5waves, bi, v_ccc, w_ccc, v̂, ŵ, ei, ci)
b_Ph_b, v_Ph_b, w_Ph_b, v̂_Ph_b, ŵ_Ph_b, e_Ph_b, c_Ph_b  = phase_orient(beg3waves, bi, v_ccc, w_ccc, v̂, ŵ, ei, ci)
b_Ph_c, v_Ph_c, w_Ph_c, v̂_Ph_c, ŵ_Ph_c, e_Ph_c, c_Ph_c  = phase_orient(trans3waves, bi, v_ccc, w_ccc, v̂, ŵ, ei, ci)

@info " X averaging...."
function x_average(b_Ph, v_Ph, v̂_Ph, ŵ_Ph, e_Ph, c_Ph)

    # y, z, Ph, W
    b_xavg = mean(b_Ph, dims = 1)[1,:,:,:,:];
    v_xavg = mean(v_Ph, dims = 1)[1,:,:,:,:];
    v̂_xavg = mean(v̂_Ph, dims = 1)[1,:,:,:,:];
    ŵ_xavg = mean(ŵ_Ph, dims = 1)[1,:,:,:,:];
    e_xavg = mean(e_Ph, dims = 1)[1,:,:,:,:];
    c_xavg = mean(c_Ph, dims = 1)[1,:,:,:,:];
 
    xavg_values = (; b_xavg, v_xavg, v̂_xavg, ŵ_xavg,e_xavg, c_xavg)
    return xavg_values
end 

xavg_values = x_average(b_Ph, v_Ph, v̂_Ph, ŵ_Ph, e_Ph, c_Ph)
xavg_values_b = x_average(b_Ph_b, v_Ph_b, v̂_Ph_b, ŵ_Ph_b, e_Ph_b, c_Ph_b )
xavg_values_c = x_average(b_Ph_c, v_Ph_c, v̂_Ph_c, ŵ_Ph_c, e_Ph_c, c_Ph_c)

function fluxes(b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph, xavg_values)
    @info "Finding Perturbations..."
    # (x,y,z,Ph,W) - (y,z,Ph,W)
    b_xpert = b_Ph .- xavg_values.b_xavg;
    v_xpert = v_Ph .- xavg_values.v_xavg;
    w_xpert = w_Ph .- xavg_values.w_xavg;
    v̂_xpert = v̂_Ph .- xavg_values.v̂_xavg;
    ŵ_xpert = ŵ_Ph .- xavg_values.ŵ_xavg;

    # (x,y,z,Ph,W)
    vb_xpert = b_xpert .* v_xpert;
    wb_xpert = b_xpert .* w_xpert;
    v̂b_xpert = b_xpert .* v̂_xpert;
    ŵb_xpert = b_xpert .* ŵ_xpert;

    @info "Averaging fluxes..."
    # (y, z, Ph, W)
    vb_xpertxavg = mean(vb_xpert, dims = (1))[1,:,:,:,:];
    wb_xpertxavg = mean(wb_xpert, dims = (1))[1,:,:,:,:];
    v̂b_xpertxavg = mean(v̂b_xpert, dims = (1))[1,:,:,:,:];
    ŵb_xpertxavg = mean(ŵb_xpert, dims = (1))[1,:,:,:,:];

    flux_values = (; vb_xpertxavg, wb_xpertxavg, v̂b_xpertxavg,ŵb_xpertxavg)

    return flux_values
end

flux_values = fluxes(b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph, xavg_values)
flux_values_b = fluxes(b_Ph_b, v_Ph_b, w_Ph_b, v̂_Ph_b, ŵ_Ph_b, xavg_values_b)
flux_values_c = fluxes(b_Ph_c, v_Ph_c, w_Ph_c, v̂_Ph_c, ŵ_Ph_c, xavg_values_c)

function phase_average(xavg_values, flux_values)
    # x, y, z, Ph
    b_xavgPhavg = mean(xavg_values.b_xavg, dims = 4)[:,:,:,1];
    v_xavgPhavg = mean(xavg_values.v_xavg, dims = 4)[:,:,:,1];
    w_xavgPhavg = mean(xavg_values.w_xavg, dims = 4)[:,:,:,1];
    v̂_xavgPhavg = mean(xavg_values.v̂_xavg, dims = 4)[:,:,:,1];
    ŵ_xavgPhavg = mean(xavg_values.ŵ_xavg, dims = 4)[:,:,:,1];
    e_xavgPhavg = mean(xavg_values.e_xavg, dims = 4)[:,:,:,1];
    c_xavgPhavg = mean(xavg_values.c_xavg, dims = 4)[:,:,:,1];
    vb_xpertxavgPhavg = mean(flux_values.vb_xpertxavg, dims = 4)[:,:,:,1];
    wb_xpertxavgPhavg = mean(flux_values.wb_xpertxavg, dims = 4)[:,:,:,1];
    v̂b_xpertxavgPhavg = mean(flux_values.v̂b_xpertxavg, dims = 4)[:,:,:,1];
    ŵb_xpertxavgPhavg = mean(flux_values.ŵb_xpertxavg, dims = 4)[:,:,:,1];

    PhaseAveragedVals = (; b_xavgPhavg, v_xavgPhavg, w_xavgPhavg, 
    v̂_xavgPhavg, ŵ_xavgPhavg, e_xavgPhavg, c_xavgPhavg,
    vb_xpertxavgPhavg, wb_xpertxavgPhavg, v̂b_xpertxavgPhavg,
    ŵb_xpertxavgPhavg)

    return PhaseAveragedVals
end

function wave_average(xavg_values, flux_values)
    # x, y, z
    v_xavgWavg = mean(xavg_values.v_xavg, dims = (3,4))[:,:,1,1];
    w_xavgWavg = mean(xavg_values.w_xavg, dims = (3,4))[:,:,1,1];
    v̂_xavgWavg = mean(xavg_values.v̂_xavg, dims = (3,4))[:,:,1,1];
    ŵ_xavgWavg = mean(xavg_values.ŵ_xavg, dims = (3,4))[:,:,1,1];
    vb_xpertxavgWavg = mean(flux_values.vb_xpertxavg, dims = (3,4))[:,:,1,1];
    wb_xpertxavgWavg = mean(flux_values.wb_xpertxavg, dims = (3,4))[:,:,1,1];
    v̂b_xpertxavgWavg = mean(flux_values.v̂b_xpertxavg, dims = (3,4))[:,:,1,1];
    ŵb_xpertxavgWavg = mean(flux_values.ŵb_xpertxavg, dims = (3,4))[:,:,1,1];

    WaveAveragedVals = (; b_xavgWavg, v_xavgWavg, w_xavgWavg, 
    v̂_xavgWavg, ŵ_xavgWavg, vb_xpertxavgWavg, wb_xpertxavgWavg, v̂b_xpertxavgWavg,
    ŵb_xpertxavgWavg)

    return WaveAveragedVals
end

PhaseAveragedVals = phase_average(xavg_values, flux_values)
PhaseAveragedVals_b = phase_average(xavg_values_b, flux_values_b)
PhaseAveragedVals_c = phase_average(xavg_values_c, flux_values_c)

WaveAveragedVals = wave_average(xavg_values, flux_values)
WaveAveragedVals_b = wave_average(xavg_values_b, flux_values_b)
WaveAveragedVals_c = wave_average(xavg_values_c, flux_values_c)

zlength = length(zb)
ylength = 880

Ygrid = reshape(repeat(yb[2:ylength], zlength-1), ylength-1, zlength-1)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

Δy = yb[1]-yb[2]
Δz = zb[1]-zb[2]

function Flux_Divergence(Δy, Δz, boolZY, WaveAveragedVals, PhaseAveragedVals)

    wb_xpertxavgWavg_dz = ((WaveAveragedVals.wb_xpertxavgWavg[2:end,1:end-1] .- WaveAveragedVals.wb_xpertxavgWavg[2:end,2:end])./Δz) .* boolZY;
    vb_xpertxavgWavg_dy = ((WaveAveragedVals.vb_xpertxavgWavg[1:end-1,2:end] .- WaveAveragedVals.vb_xpertxavgWavg[2:end,2:end])./Δy) .* boolZY;
    ∇_xpertxavgWavg = wb_xpertxavgWavg_dz .+ vb_xpertxavgWavg_dy;

    FluxDiv_Wavg = (; wb_xpertxavgWavg_dz, vb_xpertxavgWavg_dy, ∇_xpertxavgWavg)

    wb_xpertxavgPhavg_dz = ((PhaseAveragedVals.wb_xpertxavgPhavg[2:end,1:end-1,:] .- PhaseAveragedVals.wb_xpertxavgPhavg[2:end,2:end,:])./Δz) .* boolZY;
    vb_xpertxavgPhavg_dy = ((PhaseAveragedVals.vb_xpertxavgPhavg[1:end-1,2:end,:] .- PhaseAveragedVals.vb_xpertxavgPhavg[2:end,2:end,:])./Δy) .* boolZY;
    ∇_xpertxavgPhavg = wb_xpertxavgPhavg_dz .+ vb_xpertxavgPhavg_dy;

    FluxDiv_Phavg = (; wb_xpertxavgPhavg_dz, vb_xpertxavgPhavg_dy, ∇_xpertxavgPhavg)

    return FluxDiv_Wavg, FluxDiv_Phavg

end

FluxDiv_Wavg, FluxDiv_Phavg = Flux_Divergence(Δy, Δz, boolZY, WaveAveragedVals, PhaseAveragedVals)
FluxDiv_Wavg_b, FluxDiv_Phavg_b = Flux_Divergence(Δy, Δz, boolZY, WaveAveragedVals_b, PhaseAveragedVals_b)
FluxDiv_Wavg_c, FluxDiv_Phavg_c = Flux_Divergence(Δy, Δz, boolZY, WaveAveragedVals_c, PhaseAveragedVals_c)

phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

@info "Saving Files..."

function savefluxes(savename, apath, FluxDiv_Wavg, FluxDiv_Phavg, PhaseAveragedVals, WaveAveragedVals, yb, zb, phase_times)
    
    filescalename = apath * savename * ".jld2"

    jldsave(filescalename; 
        FluxDiv_Wavg,
        FluxDiv_Phavg,
        PhaseAveragedVals,
        WaveAveragedVals,
        yb, zb, phase_times)
end

savefluxes("BFluxTurb_end_" * sn, path_name * "Analysis/", FluxDiv_Wavg, FluxDiv_Phavg,
 PhaseAveragedVals, WaveAveragedVals, yb, zb, phase_times)

savefluxes("BFluxTurb_beg_" * sn, path_name * "Analysis/", FluxDiv_Wavg_b, FluxDiv_Phavg_b, 
PhaseAveragedVals_b, WaveAveragedVals_b, yb, zb, phase_times)

savefluxes("BFluxTurb_mid_" * sn, path_name * "Analysis/", FluxDiv_Wavg_c, FluxDiv_Phavg_c, 
PhaseAveragedVals_c, WaveAveragedVals_c, yb, zb, phase_times)