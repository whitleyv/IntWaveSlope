using Plots
using Measures
using Statistics
using Printf
using JLD2

sn = "U250N100Lz100g100"

apath = "Analysis/Plots/"
dpath = "Data/"
dpath = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/Analysis/"
apath = dpath
filescalename = dpath * "PE_mixingrate.jld2"

scale_file = jldopen(filescalename, "r+")

ϕd = scale_file[sn*"/ϕd"]
ϕi = scale_file[sn*"/ϕi"]
ϕe = scale_file[sn*"/ϕe"]
N_b = scale_file[sn*"/N_b"]
E_b = scale_file[sn*"/E_b"]

tlength = length(ϕd)
timeseries = (;times=0:600:(tlength-1)*600)

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

include("WaveValues.jl")
wave_info=get_wave_indices(timeseries, pm, tlength)

Tσ6_idx = wave_info.T_Tσs[7] # completion of 6 waves
Tσ10_idx = wave_info.T_Tσs[11] # completion of 6 waves

ϕe_610 = mean(ϕe[Tσ6_idx:Tσ10_idx])

@info "Plotting..."

big_title = @sprintf("Irreversible Mixing Rate, U₀=%0.2f, N=%0.2f×10⁻³, δ=%0.1f", pm.U₀, 10^3*pm.Ñ, pm.U₀/pm.Ñ)

ep2 = plot(wavetimes, ϕe, lw = 5, color = :green,
        label=@sprintf("<ϕe>₆_₁₀= %0.3f", ϕe_610),
        xlabel = "Tσ", ylabel="ϕe [m²s⁻³]", legend_font = font(16), 
        guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=10.0mm, left_margin=10.0mm, right_margin=20.0mm,
        legend = :bottomleft, size = (1000,800), title = big_title)
    plot!(wavetimes, ϕd, lw = 5, color = :gray50, label="ϕd")
    plot!(wavetimes, ϕi, lw = 5, color = :gray30, label="ϕi")


savefig(ep2, apath * "PE_rate" * setname * ".png")

# rolling wave average?
@info "Computing Rolling Wave Averages..."
# averages based on true wave times not consistent across N varying sims

Wl = wave_info.Wl
wav_arr_length = Wl*wave_info.nTσ
wav_Avgarr_length = wav_arr_length - 2*Wl

# global value in time
ϕe_Wavg = zeros(wav_Avgarr_length);
ϕd_Wavg = zeros(wav_Avgarr_length);
ϕi_Wavg = zeros(wav_Avgarr_length);

for (ki, tk) in enumerate(wave_info.WavePeriods[Wl+1:end-Wl])
    # glo2Wlbal value in time
    Widxs = wave_info.WavePeriods[ki:ki+2*Wl] 
    ϕe_Wavg[ki] = mean(ϕe[Widxs])
    ϕd_Wavg[ki] = mean(ϕd[Widxs])
    ϕi_Wavg[ki] = mean(ϕi[Widxs])

end

wavelength = 1+1/Wl:1/Wl:wave_info.nTσ-1


ep2 = plot(wavelength, ϕe_Wavg, lw = 5, color = :green,
        label=@sprintf("<ϕe>₆_₁₀= %0.3f", ϕe_610),
        xlabel = "Tσ", ylabel="ϕe [m²s⁻³]", legend_font = font(20), 
        guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=10.0mm, left_margin=10.0mm, right_margin=20.0mm,
        legend = :topleft, size = (1000,800), title = big_title)
    plot!(wavelength, ϕd_Wavg, lw = 5, color = :gray50, label="ϕd")
    plot!(wavelength, ϕi_Wavg, lw = 5, color = :gray30, label="ϕi")


savefig(ep2, apath * "PE_rateW" * setname * ".png")
