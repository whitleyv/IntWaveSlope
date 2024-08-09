using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using CairoMakie

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

filesetnames = "SetList_mp.jld2"
scale_file = jldopen(filesetnames, "r+")
sns = scale_file["setnames"]
sfiles = scale_file["setfilenames"]

name_prefix =  sfiles[7] * sns[7]

filepath = path_name * name_prefix * ".jld2"
apath = path_name * "Analysis/"

b_timeseries = FieldTimeSeries(filepath, "b");
v_timeseries = FieldTimeSeries(filepath, "v");
w_timeseries = FieldTimeSeries(filepath, "w");

xb, yb, zb = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC
xw, yw, zw = nodes(w_timeseries) #CCF

ylength = 880 
tlength = 161
zlength = length(zb)
xlength = length(xb)

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

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl
fullwaves = wave_info.WavePeriods[:,4:end] # t = 4T : t = 7T

v̂_Ph = v̂[:,:,:,fullwaves];
ŵ_Ph = ŵ[:,:,:,fullwaves];
b_Ph = bi[:,:,:,fullwaves];

v̂_Phavg = mean(v̂_Ph, dims = 5)[:,:,:,:,1];
ŵ_Phavg = mean(ŵ_Ph, dims = 5)[:,:,:,:,1];
b_Phavg = mean(b_Ph, dims = 5)[:,:,:,:,1];

v̂_Wavg = mean(v̂_Ph, dims = (4,5))[:,:,:,1,1];
ŵ_Wavg = mean(ŵ_Ph, dims = (4,5))[:,:,:,1,1];
b_Wavg = mean(b_Ph, dims = (4,5))[:,:,:,1,1];

v̂_Phdep = v̂_Phavg .- v̂_Wavg;
ŵ_Phdep = ŵ_Phavg .- ŵ_Wavg;
b_Phdep = b_Phavg .- b_Wavg;

vb_phasedep = b_Phdep .* v̂_Phdep;
wb_phasedep = b_Phdep .* ŵ_Phdep;

# (x,y,z)
vb_phasedepWavg = mean(vb_phasedep, dims = (4))[:,:,:,1];
wb_phasedepWavg = mean(wb_phasedep, dims = (4))[:,:,:,1];

zlength = length(zb)
ylength = 880
xlength = length(xb)

Ygrid = reshape(repeat(yb[2:ylength], xlength*(zlength-1)), ylength-1, zlength-1, xlength)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

boolZYX = permutedims(boolZY, [3, 1, 2])

filescalename = apath * "BFlux_Rotated_" * sn * ".jld2"

phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

jldsave(filescalename; 
    vb_phasedepWavg,
    wb_phasedepWavg,
    ŵ_Wavg,
    v̂_Wavg,
    phase_times,
    boolZYX,
    xb,yb,zb) 
    