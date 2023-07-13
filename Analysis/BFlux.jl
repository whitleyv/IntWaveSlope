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

xb, yb, zb = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC
xw, yw, zw = nodes(w_timeseries) #CCF

ylength = 880 
tlength = 161
bi = interior(b_timeseries)[:,1:ylength,:,:];
vi = interior(v_timeseries)[:,1:ylength+1,:,:];
wi = interior(w_timeseries)[:,1:ylength,:,:];

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
end5waves = wave_info.WavePeriods[:,7:end]

# reorienting into phases
# x y z Ph W
b_Ph = bi[:,:,:,end5waves];
v_Ph = v_ccc[:,:,:,end5waves];
w_Ph = w_ccc[:,:,:,end5waves];
v̂_Ph = v̂[:,:,:,end5waves];
ŵ_Ph = ŵ[:,:,:,end5waves];

# x, y, z, Ph
b_Phavg = mean(b_Ph, dims = 5)[:,:,:,:,1];
v_Phavg = mean(v_Ph, dims = 5)[:,:,:,:,1];
v̂_Phavg = mean(v̂_Ph, dims = 5)[:,:,:,:,1];
w_Phavg = mean(ŵ_Ph, dims = 5)[:,:,:,:,1];
@info "Wave Averaging..."

# x, y, z
b_Wavg = mean(b_Ph, dims = (4,5))[:,:,:,1,1];
v_Wavg = mean(v_Ph, dims = (4,5))[:,:,:,1,1];
w_Wavg = mean(w_Ph, dims = (4,5))[:,:,:,1,1];
v̂_Wavg = mean(v̂_Ph, dims = (4,5))[:,:,:,1,1];
ŵ_Wavg = mean(ŵ_Ph, dims = (4,5))[:,:,:,1,1];

# (x, y,z,Ph,W) - (x,y,z)
#b_Phpert = b_Ph .- b_Phavg;
#v_Phpert = v_Ph .- v_Phavg;
#w_Phpert = w_Ph .- w_Phavg;
#v̂_Phpert = v̂_Ph .- v̂_Phavg;
#ŵ_Phpert = ŵ_Ph .- ŵ_Phavg;

@info "Finding Perturbations..."

# (x,y,z,Ph,W) - (x, y,z)
b_Wpert = b_Ph .- b_Wavg;
v_Wpert = v_Ph .- v_Wavg;
w_Wpert = w_Ph .- w_Wavg;
v̂_Wpert = v̂_Ph .- v̂_Wavg;
ŵ_Wpert = ŵ_Ph .- ŵ_Wavg;

# (y,z,Ph,W) * (y,z,Ph,W)
#vb_Phpert = b_Phpert .* v_Phpert;
#wb_Phpert = b_Phpert .* w_Phpert;
#v̂b_Phpert = b_Phpert .* v̂_Phpert;
#ŵb_Phpert = b_Phpert .* ŵ_Phpert;

# (x,y,z,Ph,W)
vb_Wpert = b_Wpert .* v_Wpert;
wb_Wpert = b_Wpert .* w_Wpert;
v̂b_Wpert = b_Wpert .* v̂_Wpert;
ŵb_Wpert = b_Wpert .* ŵ_Wpert;

# (x,y,z,Ph)
#vb_PhpertPhavg = mean(vb_Phpert, dims = 4)[:,:,:,1];
#wb_PhpertPhavg = mean(wb_Phpert, dims = 4)[:,:,:,1];
#v̂b_PhpertPhavg = mean(v̂b_Phpert, dims = 4)[:,:,:,1];
#ŵb_PhpertPhavg = mean(ŵb_Phpert, dims = 4)[:,:,:,1];

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

phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

@info "Saving Files..."

savename = "BFluxFull_" * sn
apath  = path_name * "Analysis/"

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
v̂_Wavg,
b_Phavg,
yb, zb, phase_times)

