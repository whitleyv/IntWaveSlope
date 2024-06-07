using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "paramset"
            help = "sets which parameters to use"
            default = "U250N100Lz100g100"
    end
    return parse_args(s)
end

args=parse_commandline()

sn = args["paramset"]

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

#sn1 = "U150N100Lz100g100"
sn1 = sn

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

name_prefix = "IntWave_mp_" * sn1
filepath = path_name * name_prefix * ".jld2"

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

v_ccc = 0.5 .* (vi[:,1:end-1,:,:] .+ vi[:,2:end,:,:]);
w_ccc = 0.5 .* (wi[:,:,1:end-1,:] .+ wi[:,:,2:end,:]);

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl

lin_land = linslope.(yb) # the actual z value the land hits for each value of y.

# array Wl x (7-11)nTσ
#__Period6 ____Period7_____Period8____...
#__ + 1/15
#__ + 2/15
#__ + 3/15
#    ...
#__ + 15/15

# reorienting into phases
# x y z Ph W
function phase_orient(indexing, bi, v_ccc, w_ccc, ci, SGSi)
    b_Ph = bi[:,:,:,indexing]
    v_Ph = v_ccc[:,:,:,indexing]
    w_Ph = w_ccc[:,:,:,indexing]
    c_Ph = ci[:,:,:,indexing]
    SGS_Ph = SGSi[:,:,:,indexing]

    PhaseOrientVals = (; b_Ph, v_Ph, w_Ph, c_Ph, SGS_Ph)

    return PhaseOrientVals
end

@info "Wave averaging..."
# x y z
# x y z Ph W[2:end-1] if rolling

function rollingwaveaveraging(bi, v_ccc, w_ccc, ci, SGSi, xL, yL, zL, WL, nWl, WavePeriods)
    
    v_rWavg = zeros(xL, yL, zL, WL*nWl - 2*Wl);
    b_rWavg = zeros(xL, yL, zL,WL*nWl - 2*Wl);
    c_rWavg = zeros(xL, yL, zL, WL*nWl - 2*Wl);
    w_rWavg = zeros(xL, yL, zL, WL*nWl - 2*Wl);
    SGS_rWavg = zeros(xL, yL, zL, WL*nWl - 2*Wl);

    # xL, yL, zL, tL[Wl:165-Wl]
    for i = WL+1:(nWl*Wl - Wl)
        @info "$i/135"
        Windexs = WavePeriods[i-Wl:i+Wl]
        v_rWavg[:,:,:,i - Wl] = mean(v_ccc[:,:,:,Windexs], dims = 4)
        w_rWavg[:,:,:,i - Wl] = mean(w_ccc[:,:,:,Windexs], dims = 4)
        b_rWavg[:,:,:,i - Wl] = mean(bi[:,:,:,Windexs], dims = 4)
        c_rWavg[:,:,:,i - Wl] = mean(ci[:,:,:,Windexs], dims = 4)
        SGS_rWavg[:,:,:,i - Wl] = mean(SGSi[:,:,:,Windexs], dims = 4)
    end

    WaveAveragedVals = (; v_rWavg, w_rWavg, c_rWavg, b_rWavg, SGS_rWavg)

    return WaveAveragedVals
end

PhaseOrientVals = phase_orient(wave_info.WavePeriods, bi, v_ccc, w_ccc, ci, SGSi)
(xL, yL, zL, WL, nWl) = size(PhaseOrientVals.v_Ph)

RollingWaveAveragedVals = rollingwaveaveraging(bi, v_ccc, w_ccc, ci, SGSi, xL, yL, zL, WL, nWl, wave_info.WavePeriods)

@info "Phase averaging..."
# x y z Ph

function rollingphaseaveraging(PhaseOrientVals, RollingWaveAveragedVals, xL, yL, zL, WL, nWl)

    v_rPhavg = zeros(xL, yL, zL, WL, nWl -2);
    b_rPhavg = zeros(xL, yL, zL, WL, nWl -2);
    c_rPhavg = zeros(xL, yL, zL, WL, nWl -2);
    w_rPhavg = zeros(xL, yL, zL, WL, nWl -2);

    for i = 2:nWl-1
        @info "$i/10"
        v_rPhavg[:,:,:,:,i-1] = mean(PhaseOrientVals.v_Ph[:,:,:,:,i-1:i+1], dims = (5))
        w_rPhavg[:,:,:,:,i-1] = mean(PhaseOrientVals.w_Ph[:,:,:,:,i-1:i+1], dims = (5))
        b_rPhavg[:,:,:,:,i-1] = mean(PhaseOrientVals.b_Ph[:,:,:,:,i-1:i+1], dims = (5))
        c_rPhavg[:,:,:,:,i-1] = mean(PhaseOrientVals.c_Ph[:,:,:,:,i-1:i+1], dims = (5))
    end

    v_rPhavg_tl = reshape(v_rPhavg, xL, yL, zL, WL*nWl - 2*Wl)
    w_rPhavg_tl = reshape(w_rPhavg, xL, yL, zL, WL*nWl - 2*Wl)
    b_rPhavg_tl = reshape(b_rPhavg, xL, yL, zL, WL*nWl - 2*Wl)
    c_rPhavg_tl = reshape(c_rPhavg, xL, yL, zL, WL*nWl - 2*Wl)

    # x y z tL[Wl:end-Wl] .- x y z = x y z Ph
    v_Phdep = v_rPhavg_tl .- RollingWaveAveragedVals.v_rWavg
    w_Phdep = w_rPhavg_tl .- RollingWaveAveragedVals.w_rWavg
    b_Phdep = b_rPhavg_tl .- RollingWaveAveragedVals.b_rWavg
    c_Phdep = c_rPhavg_tl .- RollingWaveAveragedVals.c_rWavg

    RollingPhaseDependentVals = (; v_Phdep, w_Phdep, b_Phdep, c_Phdep)

    RollingPhaseAveragedVals = (; v_rPhavg, w_rPhavg, c_rPhavg, b_rPhavg)

    return RollingPhaseAveragedVals, RollingPhaseDependentVals
end

RollingPhaseAveragedVals, RollingPhaseDependentVals = rollingphaseaveraging(PhaseOrientVals, RollingWaveAveragedVals, xL, yL, zL, WL, nWl)

@info "turbulent qunatities..."
# x y z Ph W .-  x y z Ph =  x y z Ph W

function rollingturbulentquants(PhaseOrientVals, RollingPhaseAveragedVals)
    v_turb = PhaseOrientVals.v_Ph[:,:,:,:, 2:end-1] .- RollingPhaseAveragedVals.v_rPhavg;
    w_turb = PhaseOrientVals.w_Ph[:,:,:,:, 2:end-1] .- RollingPhaseAveragedVals.w_rPhavg;
    b_turb = PhaseOrientVals.b_Ph[:,:,:,:, 2:end-1] .- RollingPhaseAveragedVals.b_rPhavg;
    c_turb = PhaseOrientVals.c_Ph[:,:,:,:, 2:end-1] .- RollingPhaseAveragedVals.c_rPhavg;

    TurbVals = (;v_turb, w_turb, c_turb, b_turb)

    return TurbVals
end

RollingTurbVals = rollingturbulentquants(PhaseOrientVals, RollingPhaseAveragedVals)


function rollingfluxes(RollingTurbVals, RollingPhaseDependentVals, xL, yL, zL, WL, nWl)
   
    @info "Multiplying fluxes"
    # (x,y,z,tL- 2Wl)
    vb_turb = reshape(RollingTurbVals.b_turb .* RollingTurbVals.v_turb, xL, yL, zL, WL*nWl - 2*Wl);
    wb_turb = reshape(RollingTurbVals.b_turb .* RollingTurbVals.w_turb, xL, yL, zL, WL*nWl - 2*Wl);
    # (x,y,z,tL - 2Wl)
    vb_phasedep = RollingPhaseDependentVals.b_Phdep .* RollingPhaseDependentVals.v_Phdep;
    wb_phasedep = RollingPhaseDependentVals.b_Phdep .* RollingPhaseDependentVals.w_Phdep;

    
    @info "Averaging fluxes..."
    vb_turb_rWavg = zeros(xL, yL, zL, WL*nWl - 4*Wl)
    wb_turb_rWavg = zeros(xL, yL, zL, WL*nWl - 4*Wl)
    vb_phasedep_rWavg = zeros(xL, yL, zL, WL*nWl - 4*Wl)
    wb_phasedep_rWavg = zeros(xL, yL, zL, WL*nWl - 4*Wl)

    # xL, yL, zL, tL[Wl:165-Wl]
    for i = WL+1:(nWl*Wl - 3*Wl)
        @info "$i/120"
        vb_turb_rWavg[:,:,:,i - Wl] = mean(vb_turb[:,:,:,i-Wl:i+Wl], dims = 4)
        wb_turb_rWavg[:,:,:,i - Wl] = mean(wb_turb[:,:,:,i-Wl:i+Wl], dims = 4)
        vb_phasedep_rWavg[:,:,:,i - Wl] = mean(vb_phasedep[:,:,:,i-Wl:i+Wl], dims = 4)
        wb_phasedep_rWavg[:,:,:,i - Wl] = mean(wb_phasedep[:,:,:,i-Wl:i+Wl], dims = 4)
    end


    flux_values = (; vb_turb_rWavg, wb_turb_rWavg, vb_phasedep_rWavg, wb_phasedep_rWavg)

    return flux_values
end

RollingFluxVals = rollingfluxes(RollingTurbVals, RollingPhaseDependentVals, xL, yL, zL, WL, nWl)

zlength = length(zb)
ylength = 880
xlength = length(xb)

Ygrid = reshape(repeat(yb[2:ylength], xlength*(zlength-1)), ylength-1, zlength-1, xlength)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

boolZYX = permutedims(boolZY, [3, 1, 2])

Δy = yb[1]-yb[2]
Δz = zb[1]-zb[2]


function RollingFlux_Divergence(Δy, Δz, boolZY, RollingFluxVals)

    # (x,y,z)
    wb_turb_rWavg_dz = ((RollingFluxVals.wb_turb_rWavg[:,2:end,1:end-1,:] .- RollingFluxVals.wb_turb_rWavg[:,2:end,2:end, :])./Δz) .* boolZY;
    vb_turb_rWavg_dy = ((RollingFluxVals.vb_turb_rWavg[:,1:end-1,2:end, :] .- RollingFluxVals.vb_turb_rWavg[:,2:end,2:end,:])./Δy) .* boolZY;
    ∇_turb_rWavg = wb_turb_rWavg_dz .+ vb_turb_rWavg_dy;
    
    FluxDiv_turb_rWavg = (; wb_turb_rWavg_dz, vb_turb_rWavg_dy, ∇_turb_rWavg)

    wb_phasedep_rWavg_dz = ((RollingFluxVals.wb_phasedep_rWavg[:,2:end,1:end-1,:] .- RollingFluxVals.wb_phasedep_rWavg[:,2:end,2:end, :])./Δz) .* boolZY;
    vb_phasedep_rWavg_dy = ((RollingFluxVals.vb_phasedep_rWavg[:,1:end-1,2:end, :] .- RollingFluxVals.vb_phasedep_rWavg[:,2:end,2:end, :])./Δy) .* boolZY;
    ∇_phasedep_rWavg = wb_phasedep_rWavg_dz .+ vb_phasedep_rWavg_dy;

    FluxDiv_phasedep_rWavg = (; wb_phasedep_rWavg_dz, vb_phasedep_rWavg_dy, ∇_phasedep_rWavg)

    return FluxDiv_turb_rWavg, FluxDiv_phasedep_rWavg

end

FluxDiv_turb_rWavg, FluxDiv_phasedep_rWavg = RollingFlux_Divergence(Δy, Δz, boolZYX, RollingFluxVals)


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


rolling_phase_times = b_timeseries.times[wave_info.WavePeriods[2*Wl+1:end-2*Wl]]/pm.Tσ

savefluxes("BFluxTrip_fullrAvg_mp_" * sn1, path_name * "Analysis/", FluxDiv_turb_rWavg, FluxDiv_phasedep_rWavg, 
    RollingPhaseDependentVals, RollingTurbVals, RollingWaveAveragedVals, RollingFluxVals, yb, zb, rolling_phase_times)

savepath = path_name * "Analysis/" * "BFluxTrip_fullrAvg_mp_" * sn1 * ".jld2"
scale_file = jldopen(savepath, "r+")
skeys = keys(scale_file)

RollingWaveAveragedVals = scale_file[skeys[5]];
FluxDiv_phasedep_rWavg = scale_file[skeys[2]];
FluxDiv_turb_rWavg = scale_file[skeys[1]];

ΔWt = pm.Tσ/ wave_info.Wl
dt_b_rWavg = (RollingWaveAveragedVals.b_rWavg[:,:,:,1:end-1] .- RollingWaveAveragedVals.b_rWavg[:,:,:,2:end]) ./ ΔWt;

b_rWavg_xavg = mean(RollingWaveAveragedVals.b_rWavg, dims = 1)[1,:,:,:];
v_rWavg_xavg = mean(RollingWaveAveragedVals.v_rWavg, dims = 1)[1,:,:,:];
∇_phasedep_rWavg_xavg = mean(FluxDiv_phasedep_rWavg.∇_phasedep_rWavg, dims = 1)[1,:,:, :];
wb_phasedep_rWavg_dz_xavg = mean(FluxDiv_phasedep_rWavg.wb_phasedep_rWavg_dz, dims = 1)[1,:,:, :];
vb_phasedep_rWavg_dy_xavg = mean(FluxDiv_phasedep_rWavg.vb_phasedep_rWavg_dy, dims = 1)[1,:,:, :];
wb_turb_rWavg_dz_xavg = mean(FluxDiv_turb_rWavg.wb_turb_rWavg_dz, dims = 1)[1,:,:, :];
vb_turb_rWavg_dy_xavg = mean(FluxDiv_turb_rWavg.vb_turb_rWavg_dy, dims = 1)[1,:,:, :];
∇_turb_rWavg_xavg = mean(FluxDiv_turb_rWavg.∇_turb_rWavg, dims = 1)[1,:,:, :];
c_rWavg_xavg = mean(RollingWaveAveragedVals.c_rWavg, dims = 1)[1,:,:, :];
SGS_rWavg_xavg = mean(RollingWaveAveragedVals.SGS_rWavg, dims = 1)[1,:,:, :];
dt_b_rWavg_xavg = mean(dt_b_rWavg, dims = 1)[1,:,:,:];

function normal_profile(zht, lin_land, zc, bi, start_index, ylength, tlength)

    b_profile = zeros(zht, tlength)

    for hdx = 1:zht
        hab = hdx*2
        # where we want the hovmoller data to come from
        land_pz = lin_land .+ hab # where we want the hovmoller data to come from

        # for each y value, find the first z value above that number
        indexLength = length(start_index:ylength)
        depth_indices = zeros(Int, indexLength);

        for (zdx, i) in enumerate(start_index:ylength)
            depth_indices[zdx] = sum(zc .< land_pz[i]) +1
        end

        # the line at this height above bottom
        b_slopelines = zeros(indexLength, tlength);

        # at each y value
        for (zdx, j) in enumerate(start_index:ylength)
            b_slopelines[zdx, :] = bi[j, depth_indices[zdx], :]
        end

        b_profile[hdx, :] = mean(b_slopelines, dims = 1)
    end

    return b_profile

end

st_y = 140 # 107
en_y = 212 # 250
zL = 75

b_rWavg_prof =  normal_profile(zL, lin_land, zb, b_rWavg_xavg, st_y, en_y, size(b_rWavg_xavg)[end])

SGS_rWavg_prof=  normal_profile(zL, lin_land, zb, SGS_rWavg_xavg, st_y, en_y, size(SGS_rWavg_xavg)[end])
c_rWavg_prof =  normal_profile(zL, lin_land, zb, c_rWavg_xavg, st_y, en_y, size(c_rWavg_xavg)[end])
v_rWavg_prof =  normal_profile(zL, lin_land, zb, v_rWavg_xavg, st_y, en_y, size(v_rWavg_xavg)[end])
∇_phasedep_rWavg_prof =  normal_profile(zL-1, lin_land, zb[2:end], ∇_phasedep_rWavg_xavg, st_y, en_y, size(∇_phasedep_rWavg_xavg)[end])
∇_turb_rWavg_prof =  normal_profile(zL-1, lin_land, zb[2:end], ∇_turb_rWavg_xavg, st_y, en_y, size(∇_turb_rWavg_xavg)[end])
wb_phasedep_rWavg_dz_prof =  normal_profile(zL-1, lin_land, zb[2:end], wb_phasedep_rWavg_dz_xavg, st_y, en_y, size(wb_phasedep_rWavg_dz_xavg)[end])
vb_phasedep_rWavg_dy_prof =  normal_profile(zL-1, lin_land, zb[2:end], vb_phasedep_rWavg_dy_xavg, st_y, en_y, size(vb_phasedep_rWavg_dy_xavg)[end])
wb_turb_rWavg_dz_prof =  normal_profile(zL-1, lin_land, zb[2:end], wb_turb_rWavg_dz_xavg, st_y, en_y, size(wb_turb_rWavg_dz_xavg)[end])
vb_turb_rWavg_dy_prof =  normal_profile(zL-1, lin_land, zb[2:end], vb_turb_rWavg_dy_xavg, st_y, en_y, size(vb_turb_rWavg_dy_xavg)[end])
dt_b_rWavg_prof =  normal_profile(zL, lin_land, zb, dt_b_rWavg_xavg, st_y, en_y, size(dt_b_rWavg_xavg)[end])

dt_b_rWavg_prof250 = normal_profile(zL, lin_land, zb, dt_b_rWavg_xavg, st_y, en_y, size(dt_b_rWavg_xavg)[end])
dt_b_rWavg_prof150 = normal_profile(zL, lin_land, zb, dt_b_rWavg_xavg, st_y, en_y, size(dt_b_rWavg_xavg)[end])
dt_b_rWavg_prof300 = normal_profile(zL, lin_land, zb, dt_b_rWavg_xavg, st_y, en_y, size(dt_b_rWavg_xavg)[end])

filescalename = path_name * "Analysis/" * "BFluxTrip_rAvg_prof_mp_" * sn1  * ".jld2"
filescalename = path_name * "Analysis/" * "BFluxTripdbdt_rAvg_prof_mp"  * ".jld2"
filescalename = path_name * "Analysis/" * "BFluxTripb_rAvg_prof_mp"  * ".jld2"

rolling_phase_times1 = b_timeseries.times[wave_info.WavePeriods[Wl+1:end-Wl]]/pm.Tσ
rolling_phase_times2 = b_timeseries.times[wave_info.WavePeriods[2*Wl+1:end-2*Wl]]/pm.Tσ

jldsave(filescalename; 
SGS_rWavg_prof,
c_rWavg_prof,
v_rWavg_prof,
∇_phasedep_rWavg_prof,
∇_turb_rWavg_prof,
wb_phasedep_rWavg_dz_prof,
vb_phasedep_rWavg_dy_prof,
wb_turb_rWavg_dz_prof,
vb_turb_rWavg_dy_prof,
yb, zb, rolling_phase_times1,
rolling_phase_times2)

jldsave(filescalename; 
dt_b_rWavg_prof150,
dt_b_rWavg_prof250,
dt_b_rWavg_prof300)


jldsave(filescalename; 
b_rWavg_prof150,
b_rWavg_prof250,
b_rWavg_prof300)