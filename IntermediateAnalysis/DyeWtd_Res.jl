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

filescalename = apath * "DCwtdB_Res.jld2"

f = jldopen(filescalename, "a+")

setname = "U300N100Lz100g100"
 
pm2 = getproperty(SimParams(), Symbol(setname))

pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
m = -π/pm2.Lz,
l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
Tf = 2*π/pm2.f, 
Tσ = 2*π/pm2.σ))

δ = pm2.U₀/pm2.Ñ

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

name_prefix = "IntWave_smdt_" * setname
filepath = path_name * name_prefix * ".jld2"
@info "getting data from: " * setname

b_timeseries = FieldTimeSeries(filepath,"b");
xb, yb, zb = nodes(b_timeseries) #CCC
Cg_timeseries = FieldTimeSeries(filepath,"Cg");
∇b_mg²_timeseries = FieldTimeSeries(filepath, "∇b²");

tlength  = length(b_timeseries.times)
CGi = interior(Cg_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
Bi = interior(b_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
∇b_mg² = interior(∇b_mg²_timeseries)[1:xlength, 1:ylength, 1:zlength,:];

@info "Calculating buoyancy gradient averaged..."

sqt_∇b_mg² = sqrt.(∇b_mg²);

(c_weighted_b2avg, c_weighted_bavg, KTracer_denom, WTracer_denom) = calculate_tracer_weight(CGi, Bi, 
                                        ∇b_mg², sqt_∇b_mg², tlength)

@info "Calculating instantaneous K/W..."
Δt = b_timeseries.times[2:end] - b_timeseries.times[1:end-1]
# ∂ₜb̄
∂ₜb̄ = (c_weighted_bavg[2:end] .- c_weighted_bavg[1:end-1])./Δt
# ∂ₜ<(b-b̄)²> 
∂ₜb̄2 = (c_weighted_b2avg[2:end] .- c_weighted_b2avg[1:end-1])./Δt
# Wₜ = 0.5 ∂ₜb̄ / <|∇b|> = ∂ₜ(c_weighted_bavg) / WTracer_denom
Wₜ = ∂ₜb̄ ./ WTracer_denom[2:end] 
# Κₜ = 0.5 ∂ₜ<(b-b̄)²> / <|∇b|²> = ∂ₜ(c_weighted_b2avg) / KTracer_denom
Kₜ = ∂ₜb̄2 ./ KTracer_denom[2:end]

@info "Wave Averaging..."

include("../WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm2, tlength)
Wl = wave_info.Wl

(Wₜ_Wavg, Kₜ_Wavg, ∂ₜb̄_Wavg, WTracer_denom_Wavg) =rollingwaveavg(wave_info, Wₜ, Kₜ, ∂ₜb̄, WTracer_denom)

@info "Plots..."

big_title = @sprintf("Tracer Weighted Buoyancy, U₀=%0.2f, N=%0.2f×10⁻³, δ=%0.1f", pm2.U₀, 10^3*pm2.Ñ, pm2.U₀/pm2.Ñ)

bp1 = plot(b_timeseries.times[wave_info.WavePeriods[Wl+2:end-Wl]]/pm2.Tσ, Wₜ_Wavg, lw = 5, color = :dodgerblue2,
        xticks = false, ylabel="K [m²s⁻³]", 
        guidefontsize = 20, titlefont=20, tickfont = 14, left_margin=10.0mm,
        legend = false, size = (1500,1000), title = "Centered 2nd Moment Diffusivity")

bp2 = plot(b_timeseries.times[wave_info.WavePeriods[Wl+2:end-Wl]]/pm2.Tσ, Kₜ_Wavg, lw = 5, color = :firebrick2,
        xticks = false, ylabel="W [ms⁻¹]", 
        guidefontsize = 20, titlefont=20, tickfont = 14, left_margin=10.0mm,
        legend = false, size = (1500,1000), title = "1st Moment Diapycnal Velocity")

bp3 = plot(b_timeseries.times[wave_info.WavePeriods[Wl+2:end-Wl]]/pm2.Tσ, ∂ₜb̄_Wavg, lw = 5, color = :firebrick2,
xlabel = "Tσ", ylabel="∂ₜb̄ [ms⁻³]", 
guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=5.0mm, left_margin=10.0mm,
legend = false, size = (1500,1000), title = "Evolution of 1st moment")
b1plots = plot(bp1, bp2, layout= (2,1))

bp4 = plot(b_timeseries.times[wave_info.WavePeriods[Wl+2:end-Wl]]/pm2.Tσ, WTracer_denom_Wavg[2:end], lw = 5, color = :firebrick2,
xlabel = "Tσ", ylabel="|∇b̄| [s⁻²]", 
guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=5.0mm, left_margin=10.0mm,
legend = false, size = (1500,1000), title = "1st Moment Diapycnal Velocity")

b1plots = plot(bp1, bp2, bp3, bp4, layout= (2,2))

yy = 1:3
BigT = Plots.scatter(yy, marker=0, markeralpha=0, annotations=(2, yy[2],
        Plots.text(big_title, 20)), axis=nothing, grid=false, leg=false,
        foreground_color_subplot=colorant"white")

fin = Plots.plot(BigT, b1plots, layout=grid(2, 1, heights=[0.05,0.95]))

savename = "cWtdB_KWinstWavg_Res_" * setname
savefig(apath * savename * ".png")

write(f, "Kₜ_" * setname, Kₜ_Wavg)
write(f, "Wₜ_" * setname, Wₜ_Wavg)
write(f, "Wtdenom_" * setname, WTracer_denom_Wavg)
write(f, "dbt_" * setname, ∂ₜb̄_Wavg)




