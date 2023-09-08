using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using ArgParse

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "paramset"
            help = "sets which parameters to use"
            default = "U100N100Lz100g100"
    end
    return parse_args(s)
end

args=parse_commandline()

sn = args["paramset"]

ENV["GKSwstype"] = "nul" # if on remote HPC

include("../parameters.jl")
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

b_timeseries = FieldTimeSeries(filepath, "b");
v_timeseries = FieldTimeSeries(filepath, "v");
w_timeseries = FieldTimeSeries(filepath, "w");
e_timeseries = FieldTimeSeries(filepath,"ϵ");

xb, yb, zb = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC
xw, yw, zw = nodes(w_timeseries) #CCF

ylength = 880 
tlength = 161
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

include("../WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl

# array Wl x (7-11)nTσ
#__Period6 ____Period7_____Period8____...
#__ + 1/15
#__ + 2/15
#__ + 3/15
#    ...
#__ + 15/15
end5waves = wave_info.WavePeriods[:,7:end] # t = 7 T :11 T
beg4waves = wave_info.WavePeriods[:,1:4] # t = 0 : t = 4T
#beg4waves = wave_info.WavePeriods[:,1:3] # t = 0 : t = 3T
trans2waves = wave_info.WavePeriods[:,5:6] # t = 4T : t = 6T

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
b_Ph_b, v_Ph_b, w_Ph_b, v̂_Ph_b, ŵ_Ph_b, e_Ph_b, c_Ph_b  = phase_orient(beg4waves, bi, v_ccc, w_ccc, v̂, ŵ, ei, ci)
b_Ph_c, v_Ph_c, w_Ph_c, v̂_Ph_c, ŵ_Ph_c, e_Ph_c, c_Ph_c  = phase_orient(trans2waves, bi, v_ccc, w_ccc, v̂, ŵ, ei, ci)

function phase_average(b_Ph, v_Ph, v̂_Ph, ŵ_Ph, e_Ph, c_Ph)
    # x, y, z, Ph
    b_Phavg = mean(b_Ph, dims = 5)[:,:,:,:,1];
    v_Phavg = mean(v_Ph, dims = 5)[:,:,:,:,1];
    v̂_Phavg = mean(v̂_Ph, dims = 5)[:,:,:,:,1];
    ŵ_Phavg = mean(ŵ_Ph, dims = 5)[:,:,:,:,1];
    e_Phavg = mean(e_Ph, dims = 5)[:,:,:,:,1];
    c_Phavg = mean(c_Ph, dims = 5)[:,:,:,:,1];

    return b_Phavg, v_Phavg, v̂_Phavg, ŵ_Phavg, e_Phavg, c_Phavg
end

b_Phavg, v_Phavg, v̂_Phavg, ŵ_Phavg, e_Phavg, c_Phavg = phase_average(b_Ph, v_Ph, v̂_Ph, ŵ_Ph, e_Ph, c_Ph)
b_Phavg_b, v_Phavg_b, v̂_Phav_b, ŵ_Phavg_b, e_Phavg_b, c_Phavg_b = phase_average(b_Ph_b, v_Ph_b, v̂_Ph_b, ŵ_Ph_b, e_Ph_b, c_Ph_b )
b_Phavg_c, v_Phavg_c, v̂_Phavg_c, ŵ_Phavg_c, e_Phavg_c, c_Phavg_c = phase_average(b_Ph_c, v_Ph_c, v̂_Ph_c, ŵ_Ph_c, e_Ph_c,c_Ph_c)

@info "Wave Averaging..."

function wave_average(b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph)
    # x, y, z
    b_Wavg = mean(b_Ph, dims = (4,5))[:,:,:,1,1];
    v_Wavg = mean(v_Ph, dims = (4,5))[:,:,:,1,1];
    w_Wavg = mean(w_Ph, dims = (4,5))[:,:,:,1,1];
    v̂_Wavg = mean(v̂_Ph, dims = (4,5))[:,:,:,1,1];
    ŵ_Wavg = mean(ŵ_Ph, dims = (4,5))[:,:,:,1,1];

    return b_Wavg, v_Wavg, w_Wavg, v̂_Wavg, ŵ_Wavg
end

b_Wavg, v_Wavg, w_Wavg, v̂_Wavg, ŵ_Wavg = wave_average(b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph)
b_Wavg_b, v_Wavg_b, w_Wavg_b, v̂_Wavg_b, ŵ_Wavg_b = wave_average(b_Ph_b, v_Ph_b, w_Ph_b, v̂_Ph_b, ŵ_Ph_b)
b_Wavg_c, v_Wavg_c, w_Wavg_c, v̂_Wavg_c, ŵ_Wavg_c = wave_average(b_Ph_c, v_Ph_c, w_Ph_c, v̂_Ph_c, ŵ_Ph_c)

@info "Finding Perturbations..."
function fluxes(b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph, b_Wavg, v_Wavg, w_Wavg, v̂_Wavg, ŵ_Wavg)
    @info "Finding Perturbations..."
    # (x,y,z,Ph,W) - (x, y,z)
    b_Wpert = b_Ph .- b_Wavg;
    v_Wpert = v_Ph .- v_Wavg;
    w_Wpert = w_Ph .- w_Wavg;
    v̂_Wpert = v̂_Ph .- v̂_Wavg;
    ŵ_Wpert = ŵ_Ph .- ŵ_Wavg;

    # (x,y,z,Ph,W)
    vb_Wpert = b_Wpert .* v_Wpert;
    wb_Wpert = b_Wpert .* w_Wpert;
    v̂b_Wpert = b_Wpert .* v̂_Wpert;
    ŵb_Wpert = b_Wpert .* ŵ_Wpert;

    @info "Averaging fluxes..."
    # (x,y,z)
    vb_WpertWavg = mean(vb_Wpert, dims = (4,5))[:,:,:,1,1];
    wb_WpertWavg = mean(wb_Wpert, dims = (4,5))[:,:,:,1,1];
    v̂b_WpertWavg = mean(v̂b_Wpert, dims = (4,5))[:,:,:,1,1];
    ŵb_WpertWavg = mean(ŵb_Wpert, dims = (4,5))[:,:,:,1,1];

    # (x,y,z, W)
    vb_WpertPhavg = mean(vb_Wpert, dims = 5)[:,:,:,:,1];
    wb_WpertPhavg = mean(wb_Wpert, dims = 5)[:,:,:,:,1];
    v̂b_WpertPhavg = mean(v̂b_Wpert, dims = 5)[:,:,:,:,1];
    ŵb_WpertPhavg = mean(ŵb_Wpert, dims = 5)[:,:,:,:,1];

    return vb_WpertPhavg, wb_WpertPhavg, v̂b_WpertPhavg,ŵb_WpertPhavg,vb_WpertWavg, wb_WpertWavg, v̂b_WpertWavg, ŵb_WpertWavg
end

vb_WpertPhavg, wb_WpertPhavg, v̂b_WpertPhavg, ŵb_WpertPhavg, vb_WpertWavg, wb_WpertWavg, v̂b_WpertWavg, ŵb_WpertWavg = fluxes(b_Ph, v_Ph, w_Ph, v̂_Ph, ŵ_Ph, b_Wavg, v_Wavg, w_Wavg, v̂_Wavg, ŵ_Wavg)
vb_WpertPhavg_b, wb_WpertPhavg_b, v̂b_WpertPhavg_b, ŵb_WpertPhavg_b, vb_WpertWavg_b, wb_WpertWavg_b, v̂b_WpertWavg_b, ŵb_WpertWavg_b = fluxes(b_Ph_b, v_Ph_b, w_Ph_b, v̂_Ph_b, ŵ_Ph_b, b_Wavg_b, v_Wavg_b, w_Wavg_b, v̂_Wavg_b, ŵ_Wavg_b)
vb_WpertPhavg_c, wb_WpertPhavg_c, v̂b_WpertPhavg_c, ŵb_WpertPhavg_c, vb_WpertWavg_c, wb_WpertWavg_c, v̂b_WpertWavg_c, ŵb_WpertWavg_c = fluxes(b_Ph_c, v_Ph_c, w_Ph_c, v̂_Ph_c, ŵ_Ph_c, b_Wavg_c, v_Wavg_c, w_Wavg_c, v̂_Wavg_c, ŵ_Wavg_c)

phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

@info "Saving Files..."

function savefluxes(savename, apath, vb_WpertPhavg, wb_WpertPhavg, v̂b_WpertPhavg, ŵb_WpertPhavg,
    vb_WpertWavg, wb_WpertWavg, v̂b_WpertWavg, ŵb_WpertWavg, v_Phavg, v̂_Phavg, b_Phavg, e_Phavg,
    v_Wavg, v̂_Wavg, yb, zb, phase_times)

    filescalename = apath * savename * ".jld2"

    jldsave(filescalename; 
        vb_WpertPhavg, # perturbations via wave average, but still split in phases for <>
        wb_WpertPhavg,
        v̂b_WpertPhavg,
        ŵb_WpertPhavg,
        vb_WpertWavg, # perturbations via wave average and <> wave averaged
        wb_WpertWavg,
        v̂b_WpertWavg,
        ŵb_WpertWavg,
        v_Phavg, 
        v_Wavg,
        v̂_Phavg,
        ŵ_Phavg,
        v̂_Wavg,
        b_Phavg,
        e_Phavg,
        yb, zb, phase_times)
end

savefluxes("BFlux_end_" * sn, path_name * "Analysis/", vb_WpertPhavg, wb_WpertPhavg, v̂b_WpertPhavg, ŵb_WpertPhavg,
    vb_WpertWavg, wb_WpertWavg, v̂b_WpertWavg, ŵb_WpertWavg, v_Phavg, v̂_Phavg, b_Phavg, e_Phavg,
    v_Wavg, v̂_Wavg, yb, zb, phase_times)

savefluxes("BFlux_beg_" * sn, path_name * "Analysis/", vb_WpertPhavg_b, wb_WpertPhavg_b, v̂b_WpertPhavg_b, ŵb_WpertPhavg_b,
    vb_WpertWavg_b, wb_WpertWavg_b, v̂b_WpertWavg_b, ŵb_WpertWavg_b, v_Phavg_b, v̂_Phav_b, b_Phavg_b, e_Phavg_b,
    v_Wavg_b, v̂_Wavg_b, yb, zb, phase_times)

savefluxes("BFlux_mid_" * sn, path_name * "Analysis/", vb_WpertPhavg_c, wb_WpertPhavg_c, v̂b_WpertPhavg_c, ŵb_WpertPhavg_c,
    vb_WpertWavg_c, wb_WpertWavg_c, v̂b_WpertWavg_c, ŵb_WpertWavg_c, v_Phavg_c, v̂_Phavg_c, b_Phavg_c, e_Phavg_c,
    v_Wavg_c, v̂_Wavg_c, yb, zb, phase_times)