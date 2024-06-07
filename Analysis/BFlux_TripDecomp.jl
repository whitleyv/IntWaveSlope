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

#name_prefix = "IntWave_" * sn
name_prefix = "IntWave_mp_" * setname
#name_prefix = "IntWave_mp_nExD_" * setname

filepath = path_name * name_prefix * ".jld2"
#filepath2 = path_name * "cgl_" * sn * ".jld2"

@info "getting data from: " * sn

b_timeseries = FieldTimeSeries(filepath, "b");
v_timeseries = FieldTimeSeries(filepath, "v");
w_timeseries = FieldTimeSeries(filepath, "w");
Cg_timeseries = FieldTimeSeries(filepath,"Cg");
SGS_timeseries = FieldTimeSeries(filepath, "SGS∇κ∇b");

xb, yb, zb = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC
xw, yw, zw = nodes(w_timeseries) #CCF

ylength = 880 
tlength = 161

ci = interior(Cg_timeseries)[:, 1:ylength,:, :];
bi = interior(b_timeseries)[:,1:ylength,:,:];
vi = interior(v_timeseries)[:,1:ylength+1,:,:];
wi = interior(w_timeseries)[:,1:ylength,:,:];
SGSi = interior(SGS_timeseries)[:,1:ylength,:,:];

#f3 = jldopen(filepath2);
#CGli = f3["Cgli"];
#spinlength = size(CGli)[4]
#spinuplength = tlength - spinlength;

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl

v_ccc = 0.5 .* (vi[:,1:end-1,:,:] .+ vi[:,2:end,:,:]);
w_ccc = 0.5 .* (wi[:,:,1:end-1,:] .+ wi[:,:,2:end,:]);

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

fullwaves = wave_info.WavePeriods[:,3:end] # t = 4T : t = 7T
fullwaves = wave_info.WavePeriods[:,4:end] # t = 4T : t = 7T

# reorienting into phases
# x y z Ph W
function phase_orient(indexing, bi, v_ccc, w_ccc)
    b_Ph = bi[:,:,:,indexing]
    v_Ph = v_ccc[:,:,:,indexing]
    w_Ph = w_ccc[:,:,:,indexing]

    PhaseOrientVals = (; b_Ph, v_Ph, w_Ph)

    return PhaseOrientVals
end

function phase_orient(indexing, bi, v_ccc, w_ccc, ci, SGSi)
    b_Ph = bi[:,:,:,indexing]
    v_Ph = v_ccc[:,:,:,indexing]
    w_Ph = w_ccc[:,:,:,indexing]
    c_Ph = ci[:,:,:,indexing]
    SGS_Ph = SGSi[:,:,:,indexing]

    PhaseOrientVals = (; b_Ph, v_Ph, w_Ph, c_Ph, SGS_Ph)

    return PhaseOrientVals
end

function phase_orient(indexing, bi, v_ccc, w_ccc, ci, SGSi, CGli, minlength)
    b_Ph = bi[:,:,:,indexing]
    v_Ph = v_ccc[:,:,:,indexing]
    w_Ph = w_ccc[:,:,:,indexing]
    c_Ph = ci[:,:,:,indexing]
    SGS_Ph = SGSi[:,:,:,indexing]
    CGl_Ph = CGli[:,:,:, indexing .- minlength]

    PhaseOrientVals = (; b_Ph, v_Ph, w_Ph, c_Ph, SGS_Ph, CGl_Ph)

    return PhaseOrientVals
end

PhaseOrientVals_mid  = phase_orient(trans3waves, bi, v_ccc, w_ccc, ci, SGSi, CGli, spinuplength)
PhaseOrientVals_end  = phase_orient(end5waves, bi, v_ccc, w_ccc, ci, SGSi, CGli, spinuplength)

PhaseOrientVals_beg  = phase_orient(beg3waves, bi, v_ccc, w_ccc, ci, SGSi)
PhaseOrientVals_mid  = phase_orient(trans3waves, bi, v_ccc, w_ccc, ci, SGSi)
PhaseOrientVals_end  = phase_orient(end5waves, bi, v_ccc, w_ccc, ci, SGSi)

PhaseOrientVals_full  = phase_orient(fullwaves, bi, v_ccc, w_ccc, ci, SGSi)

@info "Wave averaging..."
# x y z
function waveaveraging(PhaseOrientVals)
    
    v_Wavg = mean(PhaseOrientVals.v_Ph, dims = (4,5))[:,:,:,1,1];
    w_Wavg = mean(PhaseOrientVals.w_Ph, dims = (4,5))[:,:,:,1,1];
    b_Wavg = mean(PhaseOrientVals.b_Ph, dims = (4,5))[:,:,:,1,1];
    c_Wavg = mean(PhaseOrientVals.c_Ph, dims = (4,5))[:,:,:,1,1];
    SGS_Wavg = mean(PhaseOrientVals.SGS_Ph, dims = (4,5))[:,:,:,1,1];

    if keys(PhaseOrientVals)[end] == :CGl_Ph
        CGl_Wavg = mean(PhaseOrientVals.CGl_Ph, dims = (4,5))[:,:,:,1,1];
        WaveAveragedVals = (; v_Wavg, w_Wavg, c_Wavg, b_Wavg, SGS_Wavg, CGl_Wavg)
    else
        WaveAveragedVals = (; v_Wavg, w_Wavg, c_Wavg, b_Wavg, SGS_Wavg)
    end

    return WaveAveragedVals
end

WaveAveragedVals_end = waveaveraging(PhaseOrientVals_end)
WaveAveragedVals_beg = waveaveraging(PhaseOrientVals_beg)
WaveAveragedVals_mid = waveaveraging(PhaseOrientVals_mid)

WaveAveragedVals_full = waveaveraging(PhaseOrientVals_full)

@info "Phase averaging..."

# x y z Ph
function phaseaveraging(PhaseOrientVals, WaveAveragedVals)
    v_Phavg = mean(PhaseOrientVals.v_Ph, dims = (5))[:,:,:,:,1];
    w_Phavg = mean(PhaseOrientVals.w_Ph, dims = (5))[:,:,:,:,1];
    b_Phavg = mean(PhaseOrientVals.b_Ph, dims = (5))[:,:,:,:,1];
    c_Phavg = mean(PhaseOrientVals.c_Ph, dims = (5))[:,:,:,:,1];

    # x y z Ph .- x y z = x y z Ph
    v_Phdep = v_Phavg .- WaveAveragedVals.v_Wavg
    w_Phdep = w_Phavg .- WaveAveragedVals.w_Wavg
    b_Phdep = b_Phavg .- WaveAveragedVals.b_Wavg
    c_Phdep = c_Phavg .- WaveAveragedVals.c_Wavg
    PhaseDependentVals = (; v_Phdep, w_Phdep, b_Phdep, c_Phdep)

    if keys(PhaseOrientVals)[end] == :CGl_Ph
        CGl_Phavg = mean(PhaseOrientVals.CGl_Ph, dims = (5))[:,:,:,:,1];
        CGl_Phdep = CGl_Phavg .- WaveAveragedVals.CGl_Wavg

        PhaseAveragedVals = (; v_Phavg, w_Phavg, c_Phavg, b_Phavg, CGl_Phavg)
        PhaseDependentVals = (; v_Phdep, w_Phdep, b_Phdep, c_Phdep, CGl_Phdep)
    else
        PhaseAveragedVals = (;v_Phavg, w_Phavg, c_Phavg, b_Phavg)
        PhaseDependentVals = (; v_Phdep, w_Phdep, b_Phdep, c_Phdep)
    end
    
    return PhaseAveragedVals, PhaseDependentVals
end

function phaseaveraging(PhaseOrientVals, WaveAveragedVals)
    v_Phavg = mean(PhaseOrientVals.v_Ph, dims = (5))[:,:,:,:,1];
    w_Phavg = mean(PhaseOrientVals.w_Ph, dims = (5))[:,:,:,:,1];
    b_Phavg = mean(PhaseOrientVals.b_Ph, dims = (5))[:,:,:,:,1];

    # x y z Ph .- x y z = x y z Ph
    v_Phdep = v_Phavg .- WaveAveragedVals.v_Wavg
    w_Phdep = w_Phavg .- WaveAveragedVals.w_Wavg
    b_Phdep = b_Phavg .- WaveAveragedVals.b_Wavg
    PhaseDependentVals = (; v_Phdep, w_Phdep, b_Phdep)

    PhaseAveragedVals = (;v_Phavg, w_Phavg, b_Phavg)
    PhaseDependentVals = (; v_Phdep, w_Phdep, b_Phdep)
    
    return PhaseAveragedVals, PhaseDependentVals
end

PhaseAveragedVals_end, PhaseDependentVals_end = phaseaveraging(PhaseOrientVals_end, WaveAveragedVals_end)
PhaseAveragedVals_beg, PhaseDependentVals_beg = phaseaveraging(PhaseOrientVals_beg, WaveAveragedVals_beg)
PhaseAveragedVals_mid, PhaseDependentVals_mid = phaseaveraging(PhaseOrientVals_mid, WaveAveragedVals_mid)

PhaseAveragedVals_full, PhaseDependentVals_full = phaseaveraging(PhaseOrientVals_full, WaveAveragedVals_full)

@info "turbulent qunatities..."
# x y z Ph W .-  x y z Ph =  x y z Ph W
function turbulentquants(PhaseOrientVals, PhaseAveragedVals)
    v_turb = PhaseOrientVals.v_Ph .- PhaseAveragedVals.v_Phavg
    w_turb = PhaseOrientVals.w_Ph .- PhaseAveragedVals.w_Phavg
    b_turb = PhaseOrientVals.b_Ph .- PhaseAveragedVals.b_Phavg
    c_turb = PhaseOrientVals.c_Ph .- PhaseAveragedVals.c_Phavg

    TurbVals = (;v_turb, w_turb, c_turb, b_turb)

    return TurbVals
end

function turbulentquants(PhaseOrientVals, PhaseAveragedVals)
    v_turb = PhaseOrientVals.v_Ph .- PhaseAveragedVals.v_Phavg
    w_turb = PhaseOrientVals.w_Ph .- PhaseAveragedVals.w_Phavg
    b_turb = PhaseOrientVals.b_Ph .- PhaseAveragedVals.b_Phavg

    TurbVals = (;v_turb, w_turb, b_turb)

    return TurbVals
end

TurbVals_end = turbulentquants(PhaseOrientVals_end, PhaseAveragedVals_end)
TurbVals_beg = turbulentquants(PhaseOrientVals_beg, PhaseAveragedVals_beg)
TurbVals_mid = turbulentquants(PhaseOrientVals_mid, PhaseAveragedVals_mid)

TurbVals_full = turbulentquants(PhaseOrientVals_full, PhaseAveragedVals_full)

function fluxes(TurbVals, PhaseDependentVals)
   
    @info "Multiplying fluxes"
    # (x,y,z,Ph,W)
    vb_turb = TurbVals.b_turb .* TurbVals.v_turb;
    wb_turb = TurbVals.b_turb .* TurbVals.w_turb;
    # (x,y,z,Ph)
    vb_phasedep = PhaseDependentVals.b_Phdep .* PhaseDependentVals.v_Phdep;
    wb_phasedep = PhaseDependentVals.b_Phdep .* PhaseDependentVals.w_Phdep;

    @info "Averaging fluxes..."
    # (x,y,z)
    vb_turbWavg = mean(vb_turb, dims = (4,5))[:,:,:,1,1];
    wb_turbWavg = mean(wb_turb, dims = (4,5))[:,:,:,1,1];
    # (x,y,z)
    vb_phasedepWavg = mean(vb_phasedep, dims = (4))[:,:,:,1];
    wb_phasedepWavg = mean(wb_phasedep, dims = (4))[:,:,:,1];

    flux_values = (; vb_turbWavg, wb_turbWavg, vb_phasedepWavg, wb_phasedepWavg)

    return flux_values
end

FluxVals_end = fluxes(TurbVals_end, PhaseDependentVals_end)
FluxVals_beg = fluxes(TurbVals_beg, PhaseDependentVals_beg)
FluxVals_mid = fluxes(TurbVals_mid, PhaseDependentVals_mid)

FluxVals_full = fluxes(TurbVals_full, PhaseDependentVals_full)

zlength = length(zb)
ylength = 880
xlength = length(xb)

Ygrid = reshape(repeat(yb[2:ylength], xlength*(zlength-1)), ylength-1, zlength-1, xlength)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

boolZYX = permutedims(boolZY, [3, 1, 2])

Δy = yb[1]-yb[2]
Δz = zb[1]-zb[2]

function Flux_Divergence(Δy, Δz, boolZY, FluxVals)

    # (x,y,z)
    wb_turbWavg_dz = ((FluxVals.wb_turbWavg[:,2:end,1:end-1] .- FluxVals.wb_turbWavg[:,2:end,2:end])./Δz) .* boolZY;
    vb_turbWavg_dy = ((FluxVals.vb_turbWavg[:,1:end-1,2:end] .- FluxVals.vb_turbWavg[:,2:end,2:end])./Δy) .* boolZY;
    ∇_turbWavg = wb_turbWavg_dz .+ vb_turbWavg_dy;
    
    FluxDiv_turbWavg = (; wb_turbWavg_dz, vb_turbWavg_dy, ∇_turbWavg)

    wb_phasedepWavg_dz = ((FluxVals.wb_phasedepWavg[:,2:end,1:end-1] .- FluxVals.wb_phasedepWavg[:,2:end,2:end])./Δz) .* boolZY;
    vb_phasedepWavg_dy = ((FluxVals.vb_phasedepWavg[:,1:end-1,2:end] .- FluxVals.vb_phasedepWavg[:,2:end,2:end])./Δy) .* boolZY;
    ∇_phasedepWavg = wb_phasedepWavg_dz .+ vb_phasedepWavg_dy;

    FluxDiv_phasedepWavg = (; wb_phasedepWavg_dz, vb_phasedepWavg_dy, ∇_phasedepWavg)

    return FluxDiv_turbWavg, FluxDiv_phasedepWavg

end

FluxDiv_turbWavg_end, FluxDiv_phasedepWavg_end = Flux_Divergence(Δy, Δz, boolZYX, FluxVals_end)
FluxDiv_turbWavg_beg, FluxDiv_phasedepWavg_beg = Flux_Divergence(Δy, Δz, boolZYX, FluxVals_beg)
FluxDiv_turbWavg_mid, FluxDiv_phasedepWavg_mid = Flux_Divergence(Δy, Δz, boolZYX, FluxVals_mid)

FluxDiv_turbWavg_full, FluxDiv_phasedepWavg_full = Flux_Divergence(Δy, Δz, boolZYX, FluxVals_full)

phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

function Advective_Flux(Δy, Δz, boolZY, WaveAveragedVals)

    # (x,y,z)
    b_Wavg_dz = ((WaveAveragedVals.b_Wavg[:,2:end,1:end-1] .- WaveAveragedVals.b_Wavg[:,2:end,2:end])./Δz) .* boolZY;
    b_Wavg_dy = ((WaveAveragedVals.b_Wavg[:,1:end-1,2:end] .- WaveAveragedVals.b_Wavg[:,2:end,2:end])./Δy) .* boolZY;
    
    v_b_Wavg_dy = WaveAveragedVals.v_Wavg[:,2:end,2:end] .* b_Wavg_dy
    w_b_Wavg_dz = WaveAveragedVals.w_Wavg[:,2:end,2:end] .* b_Wavg_dz
    u∇b_Wavg = v_b_Wavg_dy .+ w_b_Wavg_dz;
    
    return u∇b_Wavg
end

u∇b_Wavg_full = Advective_Flux(Δy, Δz, boolZYX, WaveAveragedVals_full)
@info "Saving Files..."

function savefluxes(savename, apath, FluxDiv_turbWavg, FluxDiv_phasedepWavg, PhaseDependentVals, TurbVals, WaveAveragedVals, FluxVals, yb, zb, phase_times)
    
    filescalename = apath * savename * ".jld2"

    jldsave(filescalename; 
    FluxDiv_turbWavg,
    FluxDiv_phasedepWavg,
    PhaseDependentVals,
    TurbVals,
    WaveAveragedVals,
    FluxVals,
    yb, zb, phase_times)
end

savefluxes("BFluxTrip_end_mp_" * sn, path_name * "Analysis/", FluxDiv_turbWavg_end, FluxDiv_phasedepWavg_end, 
    PhaseDependentVals_end, TurbVals_end, WaveAveragedVals_end, FluxVals_end, yb, zb, phase_times)

savefluxes("BFluxTrip_beg_mp_" * sn, path_name * "Analysis/", FluxDiv_turbWavg_beg, FluxDiv_phasedepWavg_beg,
    PhaseDependentVals_beg, TurbVals_beg, WaveAveragedVals_beg, FluxVals_beg, yb, zb, phase_times)

savefluxes("BFluxTrip_mid_mp_" * sn, path_name * "Analysis/", FluxDiv_turbWavg_mid, FluxDiv_phasedepWavg_mid, 
    PhaseDependentVals_mid, TurbVals_mid, WaveAveragedVals_mid, FluxVals_mid, yb, zb, phase_times)

savefluxes("BFluxTrip_full_mp_noS_" * sn, path_name * "Analysis/", FluxDiv_turbWavg_full, FluxDiv_phasedepWavg_full, 
    PhaseDependentVals_full, TurbVals_full, WaveAveragedVals_full, FluxVals_full, u∇b_Wavg_full, yb, zb, phase_times)

function savefluxes(savename, apath, FluxDiv_turbWavg, FluxDiv_phasedepWavg, PhaseDependentVals, TurbVals, WaveAveragedVals, FluxVals, u∇b_Wavg, yb, zb, phase_times)
    
        filescalename = apath * savename * ".jld2"
    
        jldsave(filescalename; 
        FluxDiv_turbWavg,
        FluxDiv_phasedepWavg,
        PhaseDependentVals,
        TurbVals,
        WaveAveragedVals,
        FluxVals, u∇b_Wavg,
        yb, zb, phase_times)
end

jldsave(path_name * "Analysis/BFluxTrip_full_mp_noS_xavg_" * sn; 
∇_phasedepWavg_xavg,
∇_turbWavg_xavg,
vb_phasedepWavg_dy_xavg,
wb_phasedepWavg_dz_xavg,
vb_turbWavg_dy_xavg,
wb_turbWavg_dz_xavg,
SGS_Wavg_xavg,
u∇b_Wavg_xavg,
yb, zb)

