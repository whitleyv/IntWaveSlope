using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using CurveFit
ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/"
apath = path_name * "Analysis/"

include("parameters.jl")

# how many times were saved?
zlength = 250
ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)
xlength = 38

global Us = .05:.05:.55
global setnames=[]
for u = 50:50:550
    global setnames = [setnames ; @sprintf("U%dN100Lz100g100", u)]
end
for u = [150,350, 450]
    global setnames = [setnames ; @sprintf("U250Nfd%dLz100g100", u)]
end

Lvals = length(setnames)
Kₜ_endAvg = zeros(Lvals)
Wₜ_endAvg = zeros(Lvals)
δ = zeros(Lvals)
Ñn = zeros(Lvals)

function calculate_tracer_weight(CGi, Bi, ∇b_mg², sqt_∇b_mg², tlength)
    @info "Calculating Denominator of tracer wtd average..."
    csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

    @info "Calculating first moment..."
    # volume integral of dye * buoyancy
    cb = CGi .* Bi;
    cbsum = sum(cb, dims=(1,2,3))[1,1,1,:];

    # ratio of int cb dV / int c dV
    c_weighted_bavg = cbsum ./ csum;

    @info "Calculating second moment..."
    # second moment:
    b2 = zeros(size(Bi));
    for i = 1:tlength
        b2[:,:,:, i] = (Bi[:,:,:,i] .- c_weighted_bavg[i])
    end

    cb2sum = sum(b2.^2 .* CGi, dims = (1,2,3))[1,1,1,:];

    c_weighted_b2avg = cb2sum ./ csum;

    KTracer_denom = 2 .* sum(∇b_mg² .* CGi[2:end, 2:end, 2:end, :], dims = (1,2,3))[1,1,1,:] ./ csum;
    WTracer_denom = 2 .* sum(sqt_∇b_mg² .* CGi[2:end, 2:end, 2:end, :], dims = (1,2,3))[1,1,1,:] ./ csum;

    return c_weighted_b2avg, c_weighted_bavg, KTracer_denom, WTracer_denom

end

function rollingwaveavg(wave_info, KTracer_denom, WTracer_denom, c_weighted_bavg, c_weighted_b2avg)
    # rolling wave average:
    Wl = wave_info.Wl

    wav_arr_length = Wl*wave_info.nTσ
    wav_Avgarr_length = wav_arr_length - 2*Wl

    # global value in time
    c_weighted_bavg_Wavg = zeros(wav_Avgarr_length)
    c_weighted_b2avg_Wavg = zeros(wav_Avgarr_length)
    KTracer_denom_Wavg = zeros(wav_Avgarr_length)
    WTracer_denom_Wavg = zeros(wav_Avgarr_length)

    for (ki, _) in enumerate(wave_info.WavePeriods[Wl+1:end-Wl])
        Widxs = wave_info.WavePeriods[ki:ki+2*Wl] 
        KTracer_denom_Wavg[ki] = mean(KTracer_denom[Widxs])
        WTracer_denom_Wavg[ki] = mean(WTracer_denom[Widxs])
        c_weighted_bavg_Wavg[ki] = mean(c_weighted_bavg[Widxs])
        c_weighted_b2avg_Wavg[ki] = mean(c_weighted_b2avg[Widxs])

    end
    return KTracer_denom_Wavg, WTracer_denom_Wavg, c_weighted_bavg_Wavg, c_weighted_b2avg_Wavg
end


for (m, setname) in enumerate(setnames)
 
    pm2 = getproperty(SimParams(), Symbol(setname))

    pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
    Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
    m = -π/pm2.Lz,
    l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
    Tf = 2*π/pm2.f, 
    Tσ = 2*π/pm2.σ))

    δ[m] = pm2.U₀/pm2.Ñ
    Ñn[m] = pm2.Ñ

    name_prefix = "IntWave_" * setname
    filepath = path_name * name_prefix * ".jld2"
    @info "getting data from: " * setname

    b_timeseries = FieldTimeSeries(filepath,"b");
    xb, yb, zb = nodes(b_timeseries) #CCC
    Cg_timeseries = FieldTimeSeries(filepath,"Cg");

    tlength  = length(b_timeseries.times)
    CGi = interior(Cg_timeseries)[1:xlength, 1:ylength, 1:zlength,:];
    Bi = interior(b_timeseries)[1:xlength, 1:ylength, 1:zlength,:];

    @info "Calculating buoyancy gradient averaged..."
    ∂xb = Bi[1:end-1,2:end,2:end,:] .- Bi[2:end,2:end,2:end,:];
    ∂yb = Bi[2:end,1:end-1,2:end,:] .- Bi[2:end,2:end,2:end,:];
    ∂zb = Bi[2:end,2:end,1:end-1,:] .- Bi[2:end,2:end,2:end,:];

    ∇b_mg² = ∂xb.^2 .+ ∂yb.^2 .+ ∂zb.^2
    sqt_∇b_mg² = sqrt.(∇b_mg²)

    (c_weighted_b2avg, c_weighted_bavg, KTracer_denom, WTracer_denom) = calculate_tracer_weight(CGi, Bi, ∇b_mg², sqt_∇b_mg², tlength)

    @info "Wave Averaging..."

    include("WaveValues.jl")
    wave_info=get_wave_indices(b_timeseries, pm2, tlength)
    Wl = wave_info.Wl
    (KTracer_denom_Wavg, WTracer_denom_Wavg, c_weighted_bavg_Wavg, c_weighted_b2avg_Wavg) =rollingwaveavg(wave_info, KTracer_denom, WTracer_denom, c_weighted_bavg, c_weighted_b2avg)

    Tσ6_idx = wave_info.T_Tσs[7]
    Tσ10_idx = wave_info.T_Tσs[11]

    @info "Finding the rate of change..."
    Δt = b_timeseries.times[2] - b_timeseries.times[1]

    ∂ₜb̄_Wavg = (c_weighted_bavg_Wavg[2:end] .- c_weighted_bavg_Wavg[1:end-1])./Δt
    ∂ₜb̄2_Wavg = (c_weighted_b2avg_Wavg[2:end] .- c_weighted_b2avg_Wavg[1:end-1])./Δt

    @info "Determining the diffusivity/ velocity..."
    Wₜ = ∂ₜb̄_Wavg ./ WTracer_denom_Wavg[2:end] # average in time first to get one value like linear scaling in slope?
    Kₜ = ∂ₜb̄2_Wavg ./ KTracer_denom_Wavg[2:end]

    @info "Average Values at end..."
    Wₜ_endAvg[m] = mean(Wₜ[Tσ6_idx-Wl:Tσ10_idx-Wl-1]) 
    Kₜ_endAvg[m] = mean(Kₜ[Tσ6_idx-Wl:Tσ10_idx-Wl-1])

    big_title = @sprintf("Tracer Weighted Buoyancy, U₀=%0.2f, N=%0.2f×10⁻³, δ=%0.1f", pm2.U₀, 10^3*pm2.Ñ, pm2.U₀/pm2.Ñ)

    bp6 = plot(b_timeseries.times[Wl+2:tlength-Wl]/pm2.Tσ, Kₜ, lw = 5, 
            label=@sprintf("<K>₆_₁₀ = %0.3e", Kₜ_endAvg[m]),
            xticks = false, ylabel="K [m²s⁻³]", legend_font = font(16),
            guidefontsize = 20, titlefont=20, tickfont = 14, left_margin=10.0mm,
            legend = :topleft, size = (1100,800), title = "Centered 2nd Moment Diffusivity")

    bp5 = plot(b_timeseries.times[Wl+2:tlength-Wl]/pm2.Tσ, Wₜ, lw = 5, 
            label=@sprintf("<W>₆_₁₀ =  %0.3e", Wₜ_endAvg[m]),
            xlabel = "Tσ", ylabel="W [ms⁻³]", legend_font = font(16),
            guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=5.0mm, left_margin=10.0mm,
            legend = :topleft, size = (1100,800), title = "1st Moment Diapycnal Velocity")

    bplots = plot(bp5, bp6, layout= (2,1))

    yy = 1:3
    BigT = Plots.scatter(yy, marker=0, markeralpha=0, annotations=(2, yy[2],
            Plots.text(big_title, 20)), axis=nothing, grid=false, leg=false,
            foreground_color_subplot=colorant"white")
    
    fin = Plots.plot(BigT, bplots, layout=grid(2, 1, heights=[0.05,0.95]))

    savename = "cWtdB_KWonly_" * setname
    savefig(apath * savename * ".png")

end

filescalename = apath * "DeltavCwtdB.jld2"

jldsave(filescalename; setnames, 
            Kₜ_endAvg, Wₜ_endAvg,
            δ, Ñn)

