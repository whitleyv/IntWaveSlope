using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using CurveFit
using ArgParse
using JLD2
 
ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############
 
#apath = "/glade/scratch/whitleyv/NewAdvection/Parameters/Dissipation/"
path = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
apath = path * "Analysis/"

sn = "U250N100Lz100g100"

include("../parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

        
filesetnames =  "../SetnamesList.jld2"

scale_file = jldopen(filesetnames, "r+")

sns = scale_file["setnames"]
setnames = sns[1:25]

Lvals = length(setnames)

zlength = pm.nz
ylength = 880

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

# average dissipation over waves 2-5
eps_beginAvg = zeros(Lvals)
# average dissipation over waves 6-10
eps_endAvg = zeros(Lvals)
# maximum average dissipation over waves 2-5
eps_beginMaxAvg = zeros(Lvals)
# maximum average dissipation over waves 6-10
eps_endMaxAvg = zeros(Lvals)

for (m, setname) in enumerate(setnames)

    @info "getting data from: " * setname
   path_name = path

    # need to recalculate parameters each time because N, and wave length will change!
    pm2 = getproperty(SimParams(), Symbol(setname))

    pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
                Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
                nz = round(Int,pm2.Lz/2),
                m = -π/pm2.Lz,
                l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
                Tf = 2*π/pm2.f, 
                Tσ = 2*π/pm2.σ))
    
    name_prefix = "IntWave_" * setname
    filepath = path_name * name_prefix * ".jld2"

    e_timeseries = FieldTimeSeries(filepath, "ϵ");

    # only include simulation in count if it ran long enough

    tfin = e_timeseries.times[end]
    wave10 = pm2.Tσ*10
    if tfin >= wave10
            
        xe, ye, ze = nodes(e_timeseries) #CCC

        ei = interior(e_timeseries)[:,1:ylength,1:zlength,:];

        @info "Calculating Global Averages..."
        e_xyzmax = maximum(ei,dims =(1,2,3))[1,1,1,:]
        e_xavg = mean(ei, dims=1)[1,:,:,:];

        Ygrid = reshape(repeat(ye[1:ylength], zlength), ylength, zlength)
        SlopeGrid = curvedslope.(Ygrid)
        boolZ = (SlopeGrid .+ 4) .> ze[1:zlength]'
        # above slope values at each time step
        e_fvals = e_xavg[boolZ,:]

        e_xyzavg = mean(e_fvals, dims = 1)[1,:]

        tlength = length(e_timeseries.times)

        e_xyzavg_cutoff = zeros(tlength)

        for i = 1:tlength
            tbool = e_fvals[:,i] .> 1e-6 # changing up a order to -4 orders less than expected avg for smallest avg
            if sum(tbool) > 0
                e_xyzavg_cutoff[i] = mean(e_fvals[tbool, i])
            end
        end

        @info "Computing Rolling Wave Averages..."
        # averages based on true wave times not consistent across N varying sims
        include("../WaveValues.jl")
        wave_info=get_wave_indices(e_timeseries, pm2, tlength)

        Wl = wave_info.Wl
        wav_arr_length = Wl*wave_info.nTσ
        wav_Avgarr_length = wav_arr_length - 2*Wl

        # global value in time
        eps_glob_Wavg_cut = zeros(wav_Avgarr_length);
        eps_glob_max_Wavg = zeros(wav_Avgarr_length);

        for (ki, tk) in enumerate(wave_info.WavePeriods[Wl+1:end-Wl])
            # glo2Wlbal value in time
            Widxs = wave_info.WavePeriods[ki:ki+2*Wl] 
            eps_glob_Wavg_cut[ki] = mean(e_xyzavg_cutoff[Widxs])
            eps_glob_max_Wavg[ki] = mean(e_xyzmax[Widxs])

        end

        @info "Averaging over waves 6-10 and 2-5..."
        # you have completed 1 wave by T_Tσs[2]
        # this is the time index for full array
        Tσ6_idx = wave_info.T_Tσs[7] # completion of 6 waves
        Tσ10_idx = wave_info.T_Tσs[11] # completion of 10 waves

        Tσ2_idx = wave_info.T_Tσs[3] # completion of 2 waves
        Tσ5_idx = wave_info.T_Tσs[6]

        eps_avg_Wavg_25Avg_cut = mean(e_xyzavg_cutoff[Tσ2_idx:Tσ5_idx])
        eps_avg_Wavg_610Avg_cut = mean(e_xyzavg_cutoff[Tσ6_idx:Tσ10_idx])

        eps_max_Wavg_25Avg = mean(e_xyzmax[Tσ2_idx:Tσ5_idx])
        eps_max_Wavg_610Avg = mean(e_xyzmax[Tσ6_idx:Tσ10_idx])

        wavetimes = 1:1/Wl:((wav_Avgarr_length+(Wl-1))/Wl)
        
        @info "Plotting..."

        big_title = @sprintf("Mean Dissipation>10⁻⁹, U₀=%0.2f, N=%0.2f×10⁻³, δ=%0.1f", pm2.U₀, 10^3*pm2.Ñ, pm2.U₀/pm2.Ñ)
        
        ep2 = plot(wavetimes, eps_glob_Wavg_cut, lw = 5, color = :green,
                label=@sprintf("<ϵ̄>₆_₁₀= %0.4f × 10⁻⁶", eps_avg_Wavg_610Avg_cut*1e6),
                xlabel = "Tσ", ylabel="ϵ̄ [m²s⁻³]", legend_font = font(16), 
                guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=10.0mm, left_margin=10.0mm, right_margin=20.0mm,
                legend = :bottomleft, size = (1000,800), title = big_title)
            plot!(wavetimes, eps_glob_Wavg_cut, lw = 5, color = :green,
                label=@sprintf("<ϵ̄>₂₋₅ = %0.4f × 10⁻⁶", eps_avg_Wavg_25Avg_cut*1e6))

        savefig(ep2, apath * "Dissip_Mean_" * setname * ".png")
     
        ep3 = plot(wavetimes, eps_glob_max_Wavg, lw = 5, color = :green,
                label=@sprintf("<ϵₘₐₓ>₆_₁₀= %0.4f × 10⁻⁶", eps_max_Wavg_610Avg*1e6),
                xlabel = "Tσ", ylabel="ϵ̄ [m²s⁻³]", legend_font = font(16), 
                guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=10.0mm, 
                left_margin=10.0mm, right_margin=20.0mm,
                legend = :bottomleft, size = (1000,800), title = big_title)
            plot!(wavetimes, eps_glob_max_Wavg, lw = 5, color = :green,
                label=@sprintf("<ϵₘₐₓ>₂₋₅ = %0.4f × 10⁻⁶", eps_max_Wavg_25Avg*1e6))

        savefig(ep3, apath * "DissipMax_" * setname * ".png")

        # average dissipation over waves 2-5
        eps_beginAvg[m] = eps_avg_Wavg_25Avg_cut
        # average dissipation over waves 6-10
        eps_endAvg[m] = eps_avg_Wavg_610Avg_cut

        # average dissipation over waves 2-5
        eps_beginMaxAvg[m] = eps_max_Wavg_25Avg
        # average dissipation over waves 6-10
        eps_endMaxAvg[m] = eps_max_Wavg_610Avg
    end
end

filescalename = apath * "DeltavDissip.jld2"

jldsave(filescalename; setnames, 
    eps_beginAvg, eps_endAvg,
    eps_beginMaxAvg, eps_endMaxAvg)