using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using CurveFit
using JLD2

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/"
#apath = path_name * "Analysis/"

include("parameters.jl")

# how many times were saved?
zlength = 250
ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)
xlength = 38

setname =  "U250N100Lz100g100"
pm2 = getproperty(SimParams(), Symbol(setname))

pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
m = -π/pm2.Lz,
l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
Tf = 2*π/pm2.f, 
Tσ = 2*π/pm2.σ))

δ = pm2.U₀/pm2.Ñ
Ñn = pm2.Ñ

name_prefix1 = "IntWave_spinup_" * setname
filepath1 = path_name * name_prefix1 * ".jld2"
b_timeseries1 = FieldTimeSeries(filepath1,"b");

name_prefix = "IntWave_postspin_" * setname
filepath = path_name * name_prefix * ".jld2"
@info "getting data from: " * setname

b_timeseries = FieldTimeSeries(filepath,"b");
xb, yb, zb = nodes(b_timeseries) #CCC

Cg_timeseries = FieldTimeSeries(filepath,"Cg");
bgrad_timeseries = FieldTimeSeries(filepath, "∇b²");

tlength  = length(b_timeseries.times)

CGi = interior(Cg_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
Bi = interior(b_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
∇b_mg² = interior(bgrad_timeseries)[1:xlength, 1:ylength, 1:zlength,:];

sqt_∇b_mg² = sqrt.(∇b_mg²)

# ∫ c dV
@info "Calculating Denominator of tracer wtd average..."
csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

@info "Calculating first moment..."
# ∫ cb dV
cb = CGi .* Bi;
cbsum = sum(cb, dims=(1,2,3))[1,1,1,:];

# b̄ =  ∫ cb dV / ∫ c dV

c_weighted_bavg = cbsum ./ csum;

@info "Calculating second moment..."

# ∫ (b - b̄)² c dV
b2 = zeros(size(Bi));

for i = 1:tlength
    b2[:,:,:, i] = (Bi[:,:,:,i] .- c_weighted_bavg[i])
end

cb2sum = sum(b2.^2 .* CGi, dims = (1,2,3))[1,1,1,:];

# b̄² = ∫ (b - b̄)² c dV / ∫ c dV 
c_weighted_b2avg = cb2sum ./ csum;

# 2<|∇b|²>  = 2 ∫|∇b|²c dV / ∫ c dV
KTracer_denom = 2 .* sum(∇b_mg² .* CGi, dims = (1,2,3))[1,1,1,:] ./ csum;

# 2<|∇b|>  = 2 ∫|∇b|c dV / ∫ c dV
WTracer_denom = 2 .* sum(sqt_∇b_mg² .* CGi, dims = (1,2,3))[1,1,1,:] ./ csum;

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

# averaging over last wave lengths

spinlength = length(b_timeseries1.times)
fulltime = (; times = vcat(b_timeseries1.times, b_timeseries.times))
fulltimelength = length(fulltime.times)

include("WaveValues.jl")
wave_info=get_wave_indices(fulltime, pm2, fulltimelength)

Tσ6_idx = wave_info.T_Tσs[7] .- spinlength
Tσ10_idx = wave_info.T_Tσs[11] .- spinlength

Wₜ_endAvg = mean(Wₜ[Tσ6_idx:Tσ10_idx]) 
Kₜ_endAvg = mean(Kₜ[Tσ6_idx:Tσ10_idx])

@info "Wave Averaging tracer weighted moments first..."

# rolling wave average:
Wl = wave_info.Wl
# take out spinup period of 4 wave lengths
wav_arr_length = Wl*(wave_info.nTσ - 4)
# take out 2 wave lengths at either end for rolling wave average
wav_Avgarr_length = wav_arr_length - 2*Wl

# global value in time
c_weighted_bavg_Wavg = zeros(wav_Avgarr_length)
c_weighted_b2avg_Wavg = zeros(wav_Avgarr_length)
KTracer_denom_Wavg = zeros(wav_Avgarr_length)
WTracer_denom_Wavg = zeros(wav_Avgarr_length)

for (ki, _) in enumerate(wave_info.WavePeriods[5*Wl+1:end-Wl])
    Widxs = wave_info.WavePeriods[5*Wl + ki: 6*Wl-1 + ki] .- spinlength
    KTracer_denom_Wavg[ki] = mean(KTracer_denom[Widxs])
    WTracer_denom_Wavg[ki] = mean(WTracer_denom[Widxs])
    c_weighted_bavg_Wavg[ki] = mean(c_weighted_bavg[Widxs])
    c_weighted_b2avg_Wavg[ki] = mean(c_weighted_b2avg[Widxs])

end

# wave value index in terms of rolling array
Tσ6_idxW = 5*Wl-1 .- spinlength
Tσ10_idxW = 9*Wl-1 .- spinlength

tim_Wavg = fulltime.times[wave_info.WavePeriods[5*Wl+1:end-Wl]]
ΔtW = 600
∂ₜb̄_Wavg = (c_weighted_bavg_Wavg[2:end] .- c_weighted_bavg_Wavg[1:end-1])./ΔtW
∂ₜb̄2_Wavg = (c_weighted_b2avg_Wavg[2:end] .- c_weighted_b2avg_Wavg[1:end-1])./ΔtW

Wₜ_Wavg = ∂ₜb̄_Wavg ./ WTracer_denom_Wavg[2:end] 
Kₜ_Wavg = ∂ₜb̄2_Wavg ./ KTracer_denom_Wavg[2:end]

@info "Average Values at end..."
Wₜ_Wavg_endAvg = mean(Wₜ_Wavg[Tσ6_idxW:Tσ10_idxW-1]) 
Kₜ_Wavg_endAvg = mean(Kₜ_Wavg[Tσ6_idxW:Tσ10_idxW-1])

big_title = @sprintf("Tracer Weighted Buoyancy Post Spinup, U₀=%0.2f, N=%0.2f×10⁻³, δ=%0.1f", pm2.U₀, 10^3*pm2.Ñ, pm2.U₀/pm2.Ñ)

bp6 = plot(b_timeseries.times[2:end]/pm2.Tσ, Kₜ, lw = 5, 
                label=@sprintf("<K>₆_₁₀ = %0.3e", Kₜ_endAvg),
                xticks = false, ylabel="K [m²s⁻³]", legend_font = font(16),
                guidefontsize = 20, titlefont=20, tickfont = 14, left_margin=10.0mm,
                legend = :topright, size = (1100,800), title = "Centered 2nd Moment Diffusivity")
    plot!(fulltime.times[wave_info.WavePeriods[5*Wl+2:end-Wl]]/pm2.Tσ, Kₜ_Wavg, lw = 5, 
    label=@sprintf("<K>₆_₁₀ = %0.3e", Kₜ_Wavg_endAvg),color = :gray)

bp5 = plot(b_timeseries.times[2:end]/pm2.Tσ, Wₜ, lw = 5, 
                label=@sprintf("<W>₆_₁₀ =  %0.3e", Wₜ_endAvg),
                xlabel = "Tσ", ylabel="W [ms⁻³]", legend_font = font(16),
                guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=5.0mm, left_margin=10.0mm,
                legend = :topright, size = (1100,800), title = "1st Moment Diapycnal Velocity")
    plot!(fulltime.times[wave_info.WavePeriods[5*Wl+2:end-Wl]]/pm2.Tσ, Wₜ_Wavg, lw = 5, 
                label=@sprintf("<W>₆_₁₀ =  %0.3e", Wₜ_Wavg_endAvg), color = :gray)
bplots = plot(bp6, bp5, layout= (2,1))

yy = 1:3
BigT = Plots.scatter(yy, marker=0, markeralpha=0, annotations=(2, yy[2],
                Plots.text(big_title, 20)), axis=nothing, grid=false, leg=false,
                foreground_color_subplot=colorant"white")
        
fin = Plots.plot(BigT, bplots, layout=grid(2, 1, heights=[0.05,0.95]))

savename = "cWtdB_KW_" * setname * "spin"
savefig(path_name * savename * ".png")
    

bp6 = plot(fulltime.times[wave_info.WavePeriods[5*Wl+2:end-Wl]]/pm2.Tσ, Kₜ_Wavg, lw = 5, 
                label=@sprintf("<K>₆_₁₀ = %0.3e", Kₜ_Wavg_endAvg),
                xticks = false, ylabel="K [m²s⁻³]", legend_font = font(16),
                guidefontsize = 20, titlefont=20, tickfont = 14, left_margin=10.0mm,
                legend = :topright, size = (900,800), title = "Centered 2nd Moment Diffusivity")


bp5 = plot(fulltime.times[wave_info.WavePeriods[5*Wl+2:end-Wl]]/pm2.Tσ, Wₜ_Wavg, lw = 5, 
                label=@sprintf("<W>₆_₁₀ =  %0.3e", Wₜ_Wavg_endAvg),
                xlabel = "Tσ", ylabel="W [ms⁻³]", legend_font = font(16),
                guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=5.0mm, left_margin=10.0mm,
                legend = :topright, size = (900,800), title = "1st Moment Diapycnal Velocity")

bplots = plot(bp6, bp5, layout= (2,1))

yy = 1:3
BigT = Plots.scatter(yy, marker=0, markeralpha=0, annotations=(2, yy[2],
                Plots.text(big_title, 20)), axis=nothing, grid=false, leg=false,
                foreground_color_subplot=colorant"white")
        
fin = Plots.plot(BigT, bplots, layout=grid(2, 1, heights=[0.05,0.95]))

savename = "cWtdB_KW_wavg_" * setname * "spin"
savefig(path_name * savename * ".png")
    
