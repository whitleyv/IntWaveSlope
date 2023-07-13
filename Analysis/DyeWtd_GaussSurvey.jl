using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using JLD2

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
apath = path_name * "Analysis/"

include("../parameters.jl")

setname = "U300N100Lz100g100"

# how many times were saved?
zlength = 250
ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)
xlength = 38

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

function calculate_tracer_weight(CGi, Bi, ∇b_mg², sqt_∇b_mg², tlength)
    
    @info "Calculating Denominator of tracer wtd average..."
    # ∫ c dV
    csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

    @info "Calculating first moment..."
    # ∫ cb dV
    cb = CGi .* Bi;
    cbsum = sum(cb, dims=(1,2,3))[1,1,1,:];

    # b̄ =  ∫ cb dV / ∫ c dV
    c_weighted_bavg = cbsum ./ csum;

    @info "Calculating second moment..."
    # second moment:
    b2 = zeros(size(Bi));
    for i = 1:tlength
        b2[:,:,:, i] = (Bi[:,:,:,i] .- c_weighted_bavg[i])
    end

    cb2sum = sum(b2.^2 .* CGi, dims = (1,2,3))[1,1,1,:];

    # b̄² = ∫ (b - b̄)² c dV / ∫ c dV 
    c_weighted_b2avg = cb2sum ./ csum;
    
    @info "Calculating Denoms... "
    # 2<|∇b|²>  = 2 ∫|∇b|²c dV / ∫ c dV
    KTracer_denom = 2 .* sum(∇b_mg² .* CGi, dims = (1,2,3))[1,1,1,:] ./ csum;
    # 2<|∇b|>  = 2 ∫|∇b|c dV / ∫ c dV
    WTracer_denom = 2 .* sum(sqt_∇b_mg² .* CGi, dims = (1,2,3))[1,1,1,:] ./ csum;

    return c_weighted_b2avg, c_weighted_bavg, KTracer_denom, WTracer_denom
end

function rollingwaveavg(wave_info, Wt, Kt, dbt, Wdenom)

    # rolling wave average:
    Wl = wave_info.Wl

    wav_arr_length = Wl*wave_info.nTσ
    wav_Avgarr_length = wav_arr_length - 2*Wl - 1

    # global value in time
    Kt_Wavg = zeros(wav_Avgarr_length)
    Wt_Wavg = zeros(wav_Avgarr_length)
    dbt_Wavg = zeros(wav_Avgarr_length)
    Wdenom_Wavg = zeros(wav_Avgarr_length+1)

    for (ki, _) in enumerate(wave_info.WavePeriods[Wl+2:end-Wl])
        Widxs = wave_info.WavePeriods[ki:ki+2*Wl] 
        Kt_Wavg[ki] = mean(Kt[Widxs])
        Wt_Wavg[ki] = mean(Wt[Widxs])
        dbt_Wavg[ki] = mean(dbt[Widxs])
        Wdenom_Wavg[ki] = mean(Wdenom[Widxs])
    end

    Widxf = wave_info.WavePeriods[wav_Avgarr_length+1:wav_Avgarr_length+1+2*Wl] 
    Wdenom_Wavg[end] = mean(Wdenom[Widxf])

    return Wt_Wavg, Kt_Wavg, dbt_Wavg, Wdenom_Wavg

end

pm2 = getproperty(SimParams(), Symbol(setname))

pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
m = -π/pm2.Lz,
l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
Tf = 2*π/pm2.f, 
Tσ = 2*π/pm2.σ))

zSlopeSameˢ = -pm2.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm2.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame
    
@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm2.Lzˢ + pm2.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm2.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

name1_prefix = "IntWave_" * setname
name2_prefix = "IntWave_smdt_" * setname

filepath1 = path_name * name1_prefix * ".jld2"
filepath2 = path_name * name2_prefix * ".jld2"

@info "getting data from: " * setname

b_timeseries = FieldTimeSeries(filepath1,"b");
xb, yb, zb = nodes(b_timeseries) #CCC

Cg_timeseries = FieldTimeSeries(filepath2,"Cg");

Cgr_timeseries = FieldTimeSeries(filepath1,"Cgr");
Cgs_timeseries = FieldTimeSeries(filepath1,"Cgs");
∇b_mg²_timeseries = FieldTimeSeries(filepath1, "∇b²");

tlength1  = length(b_timeseries.times)
CGi = interior(Cg_timeseries)[1:xlength, 1:ylength, 1:zlength,1:2:end];
CGri = interior(Cgr_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
CGsi = interior(Cgs_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
Bi = interior(b_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
∇b_mg²i = interior(∇b_mg²_timeseries)[1:xlength, 1:ylength, 1:zlength,:];

sqt_∇b_mg² = sqrt.(∇b_mg²i);

(b̄2_Cg, b̄_Cg, ∇b̄²_Cg, ∇b̄_Cg) = calculate_tracer_weight(CGi, Bi, 
                                    ∇b_mg²i, sqt_∇b_mg², tlength1)
(b̄2_Cgr, b̄_Cgr, ∇b̄²_Cgr, ∇b̄_Cgr) = calculate_tracer_weight(CGri, Bi, 
                                    ∇b_mg²i, sqt_∇b_mg², tlength1)
(b̄2_Cgs, b̄_Cgs, ∇b̄²_Cgs, ∇b̄_Cgs) = calculate_tracer_weight(CGsi, Bi, 
                                    ∇b_mg²i, sqt_∇b_mg², tlength1)


@info "Calculating instantaneous K/W..."
Δt = b_timeseries.times[2:end] - b_timeseries.times[1:end-1]

# ∂ₜb̄
∂ₜb̄_Cg = (b̄_Cg[2:end] .- b̄_Cg[1:end-1])./Δt
∂ₜb̄_Cgr = (b̄_Cgr[2:end] .- b̄_Cgr[1:end-1])./Δt
∂ₜb̄_Cgs = (b̄_Cgs[2:end] .- b̄_Cgs[1:end-1])./Δt

# ∂ₜ<(b-b̄)²> 
∂ₜb̄2_Cg = (b̄2_Cg[2:end] .- b̄2_Cg[1:end-1])./Δt
∂ₜb̄2_Cgr = (b̄2_Cgr[2:end] .- b̄2_Cgr[1:end-1])./Δt
∂ₜb̄2_Cgs = (b̄2_Cgs[2:end] .- b̄2_Cgs[1:end-1])./Δt

# Wₜ = 0.5 ∂ₜb̄ / <|∇b|> = ∂ₜ(c_weighted_bavg) / WTracer_denom
Wₜ_Cg = ∂ₜb̄_Cg ./ ∇b̄_Cg[2:end] 
Wₜ_Cgr = ∂ₜb̄_Cgr ./ ∇b̄_Cgr[2:end] 
Wₜ_Cgs = ∂ₜb̄_Cgs ./ ∇b̄_Cgs[2:end] 

# Κₜ = 0.5 ∂ₜ<(b-b̄)²> / <|∇b|²> = ∂ₜ(c_weighted_b2avg) / KTracer_denom
Kₜ_Cg = ∂ₜb̄2_Cg ./ ∇b̄²_Cg[2:end]
Kₜ_Cgr = ∂ₜb̄2_Cgr ./ ∇b̄²_Cgr[2:end]
Kₜ_Cgs = ∂ₜb̄2_Cgs ./ ∇b̄²_Cgs[2:end]

@info "Wave Averaging..."

include("../WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm2, tlength1)
Wl = wave_info.Wl

(Wₜ_Wavg_Cg, Kₜ_Wavg_Cg, ∂ₜb̄_Wavg_Cg, ∇b̄_Wavg_Cg) = rollingwaveavg(wave_info, Wₜ_Cg, Kₜ_Cg, ∂ₜb̄_Cg, ∇b̄_Cg)
(Wₜ_Wavg_Cgr, Kₜ_Wavg_Cgr, ∂ₜb̄_Wavg_Cgr, ∇b̄_Wavg_Cgr) = rollingwaveavg(wave_info, Wₜ_Cgr, Kₜ_Cgr, ∂ₜb̄_Cgr, ∇b̄_Cgr)
(Wₜ_Wavg_Cgs, Kₜ_Wavg_Cgs, ∂ₜb̄_Wavg_Cgs, ∇b̄_Wavg_Cgs) = rollingwaveavg(wave_info, Wₜ_Cgs, Kₜ_Cgs, ∂ₜb̄_Cgs, ∇b̄_Cgs)

filescalename = apath * "DCwtd_Gauss.jld2"

jldsave(filescalename; 
Wₜ_Wavg_Cgr, Kₜ_Wavg_Cgr, ∂ₜb̄_Wavg_Cgr, ∇b̄_Wavg_Cgr,
Wₜ_Wavg_Cg, Kₜ_Wavg_Cg, ∂ₜb̄_Wavg_Cg, ∇b̄_Wavg_Cg,
Wₜ_Wavg_Cgs, Kₜ_Wavg_Cgs, ∂ₜb̄_Wavg_Cgs, ∇b̄_Wavg_Cgs,
Wₜ_Cg, ∂ₜb̄_Cg, ∇b̄_Cg,
Wₜ_Cgr, ∂ₜb̄_Cgr, ∇b̄_Cgr,
Wₜ_Cgs, ∂ₜb̄_Cgs, ∇b̄_Cgs)