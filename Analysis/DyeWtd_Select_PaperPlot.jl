using Plots
using Measures
using Statistics
using Printf
using JLD2

dpath = "Data/"
apath = "Analysis/Plots/"

fname= "DCwtdB_Select.jld2"
filepath = dpath * fname

wtd_file = jldopen(filepath, "r+")

fkeys = keys(wtd_file)

setnames = wtd_file["setnames"]
δs = wtd_file["δs"]
Ñs = wtd_file["Ñs"]

Kt_VaryU = zeros(4,134)
Wt_VaryU = zeros(4,134)
for (m, setname) in enumerate(setnames[1:4])
    Kt_VaryU[m,:] = wtd_file["Kₜ_" * setname]
    Wt_VaryU[m,:] = wtd_file["Wₜ_" * setname]
end

Kt_varyU2 = Kt_VaryU[1,:]
Kt_varyU3 = Kt_VaryU[2,:]
Kt_varyU4 = Kt_VaryU[3,:]
Kt_varyU5 = Kt_VaryU[4,:]

Wt_varyU2 = Wt_VaryU[1,:]
Wt_varyU3 = Wt_VaryU[2,:]
Wt_varyU4 = Wt_VaryU[3,:]
Wt_varyU5 = Wt_VaryU[4,:]

Kt_varyN2 = wtd_file["Kₜ_U250Nfd200Lz100g100"]
Wt_varyN2 = wtd_file["Wₜ_U250Nfd200Lz100g100"]
Kt_varyN3 = wtd_file["Kₜ_U250Nfd300Lz100g100"]
Wt_varyN3 = wtd_file["Wₜ_U250Nfd300Lz100g100"]
Kt_varyN4 = wtd_file["Kₜ_U250Nfd400Lz130g100"]
Wt_varyN4 = wtd_file["Wₜ_U250Nfd400Lz130g100"]
Kt_varyN5 = wtd_file["Kₜ_U250Nfd500Lz130g100"]
Wt_varyN5 = wtd_file["Wₜ_U250Nfd500Lz130g100"]

Ktvary = [Kt_varyU2,Kt_varyU3, Kt_varyU4, Kt_varyU5, Kt_varyN2,Kt_varyN3, Kt_varyN4, Kt_varyN5]
Wtvary = [Wt_varyU2,Wt_varyU3, Wt_varyU4, Wt_varyU5, Wt_varyN2,Wt_varyN3, Wt_varyN4, Wt_varyN5]

lengthWavs = zeros(Int64,8)
for (m, setname) in enumerate(setnames)
    lengthWavs[m] = length(wtd_file["Kₜ_" * setname])
end
Wls = (lengthWavs.+1)/9

halfwave = round.(Int,lengthWavs/2)

Ktheory = δs.^2 .* Ñs

Ktvary_sorted = Ktvary[sortperm(Ktheory, rev=true)]
Wls_sorted = Wls[sortperm(Ktheory, rev=true)]
halfwave_sorted = halfwave[sortperm(Ktheory, rev=true)]
Ktheory_sorted = sort(Ktheory, rev=true)
Wtheory_sorted = sort(δs, rev=true)
Wtvary_sorted = Wtvary[sortperm(δs, rev=true)]
Wls_sortedW = Wls[sortperm(δs, rev=true)]
halfwave_sortedW = halfwave[sortperm(δs, rev=true)]


#######################

f = Figure(resolution = (900, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat
gb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "K [m²s⁻¹]",
 title = "Tracer Weighted Buoyancy Diffusivity") 
ax1.yticks = (0:4e-4:1.2e-3 , ["0", "4×10⁻⁴", "8×10⁻⁴", "1.2×10⁻³"])
hidexdecorations!(ax1)
limits!(ax1, 1, 10, -3e-4, 1.3e-3)

ax2 = Axis(ga[2, 1], ylabel = "W [ms⁻¹]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Velocity") 
ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
ax2.xticks = 2:2:10
limits!(ax2, 1, 10, -1.2e-4, 1.2e-4)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace

rowgap!(ga, 10)

Kt_avg = mean(Kt_VaryU, dims = 1)[1,:]
Wt_avg = mean(Wt_VaryU, dims = 1)[1,:]

wavavg = (Wls[1] .+ 1)./Wls[1] : 1 ./Wls[1]: 10 .- 1 ./Wls[1]

colorgrad = :grayC #cgrad([:dodgerblue, :slateblue4, :gray20], [0.3, 0.8,1.0])

wavs1 = (Wls_sorted[1] .+ 1)./Wls_sorted[1] : 1 ./Wls_sorted[1]: 10 .- 1 ./Wls_sorted[1]

kp = lines!(ax1, wavs1, Ktvary_sorted[1], linewidth = 3, color = Ktheory_sorted[1].*ones(length(wavs1)), 
        colormap = colorgrad, colorrange=(0,70))
        for m = 2:8
            wavs = (Wls_sorted[m] .+ 1)./Wls_sorted[m] : 1 ./Wls_sorted[m]: 10 .- 1 ./Wls_sorted[m]
            lines!(ax1, wavs, Ktvary_sorted[m], linewidth = 3, color = Ktheory_sorted[m].*ones(length(wavs)),
            colormap = colorgrad,colorrange=(0,70))
        end
        lines!(ax1, wavavg, Kt_avg, linewidth = 5, color = :dodgerblue2)

wavs2 = (Wls_sortedW[1] .+ 1)./Wls_sortedW[1] : 1 ./Wls_sortedW[1]: 10 .- 1 ./Wls_sortedW[1]

wp = lines!(ax2, wavs2, Wtvary_sorted[1], linewidth = 3, color = Wtheory_sorted[1].*ones(length(wavs2)), 
        colormap = colorgrad, colorrange = (30,145))
        for m = 2:8
            wavs = (Wls_sortedW[m] .+ 1)./Wls_sortedW[m] : 1 ./Wls_sortedW[m]: 10 .- 1 ./Wls_sortedW[m]
            lines!(ax2, wavs, Wtvary_sorted[m], linewidth = 3, color = Wtheory_sorted[m].*ones(length(wavs)), 
                        colormap = colorgrad, colorrange = (30,145))
        end
        lines!(ax2, wavavg, Wt_avg, linewidth = 5, color = :firebrick2)

cb1 = Colorbar(gb[1,1], kp, ticks = 20:20:70, labelpadding = 20.0,
    size =25, label = "δ²N [m²s⁻¹]", labelsize = 30)
cb2 = Colorbar(gb[2,1], wp, ticks = 60:30:140, 
    size =25, label = "δ [m]", labelsize = 30)

colsize!(f.layout, 2, Relative(0.05))

savename = "Paper_Cwtd_Select" 

save(apath * savename * ".png", f)


