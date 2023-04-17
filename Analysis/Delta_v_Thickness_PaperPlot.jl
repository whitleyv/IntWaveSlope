using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using CurveFit
using ArgParse
using JLD2

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename = dpath * "DeltavAllScale.jld2"

scale_file = jldopen(filescalename, "r+")

setnames = vcat(scale_file["setnames"],scale_file["setnames_oth"])
δ = vcat(scale_file["δ"],scale_file["δ_oth"])
Cheight_havg_tavg =vcat(scale_file["Cheight_havg_tavg"],scale_file["Cheight_havg_tavg_oth"])
Nheight_havg_tavg =vcat(scale_file["Nheight_havg_tavg"],scale_file["Nheight_havg_tavg_oth"])

@info "Find Bestfit Lines"

############################### INTRUSIONS v Delta
(b_Chaδ,m_Chaδ) = linear_fit(δ, Cheight_havg_tavg)
################################ Isopycnals v Delta
(b_Nhaδ,m_Nhaδ) = linear_fit(δ, Nheight_havg_tavg)

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryN = 1:11
VaryU = 12:22
# varying sigma to change criticality, holding N and f constant
Varysg = 23:26
# varying sigma but not changing criticaility, holding N and f const but changing topography
Varys = 27:30
# ones that are subcritical out of sigma varying set
Varysg_sub = 25:26


# simulations where N is varied to change delta
Cp = plot(δ, m_Chaδ .* δ .+ b_Chaδ, label =  "", color=:gray50, lw = 4)
scatter!(δ[Varysg], Cheight_havg_tavg[Varysg], label ="", markersize = 8, marker=:pentagon, color= :lightgreen)
annotate!(δ[Varysg_sub] .+ 5, Cheight_havg_tavg[Varysg_sub] .- 5, text("Subcritical γ", 10), aspect_ratio=:true)
# varying sigma without changing criticality, changing slope instead
scatter!(δ[Varys], Cheight_havg_tavg[Varys], label ="", markersize = 8, marker=:pentagon, color= :lightgreen)
scatter!(δ[VaryN], Cheight_havg_tavg[VaryN], label ="Vary N", markersize = 8, color =:green, xlabel = "δ [m]", ylabel = "Lₕ [m]",
tickfont = 15, guidefontsize = 20, titlefont = 22, title = "Dye Intrusion Thickness, Lₕ", 
bottom_margin=5mm, left_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c,
xlims = (0, 160), ylims = (0, 160), size = (1400, 700), yticks = 0:40:160, xticks = 0:40:160)
# simulations where U is varied to change delta
scatter!(δ[VaryU], Cheight_havg_tavg[VaryU], label ="Vary V₀", markersize = 8, marker=:d, color = :green)
scatter!(-1 .* δ[Varysg], Cheight_havg_tavg[Varysg], label ="Vary σ", markersize = 8, marker=:pentagon, color= :lightgreen)
# varying sigma to change the criticality 
plot!(δ, -1 .* m_Chaδ .* δ .- b_Chaδ, label =  @sprintf("%0.1f δ + %0.1f", m_Chaδ, b_Chaδ),
        color=:gray50, lw = 4)

Np = plot(δ, m_Nhaδ .* δ .+ b_Nhaδ, label =  "", color=:gray50, lw = 4)
scatter!(δ[Varysg], Nheight_havg_tavg[Varysg], label ="", markersize = 8, marker=:pentagon, color =:lightgreen,)
annotate!(δ[Varysg_sub].+ 18, Nheight_havg_tavg[Varysg_sub], text("Subcritical γ", 10), aspect_ratio=:true)
# varying sigma without changing criticality, changing slope instead
scatter!(δ[Varys], Nheight_havg_tavg[Varys], label ="", markersize = 8, marker=:pentagon, color =:lightgreen,)
scatter!(δ[VaryN], Nheight_havg_tavg[VaryN], label ="Vary N", markersize = 8, color =:green, xlabel = "δ [m]", ylabel = "Lₙ₂ [m]",
tickfont = 15, guidefontsize = 20, titlefont = 22, title = "N'²<0 Intrusion Thickness, Lₙ₂", 
bottom_margin=5mm, left_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c,
xlims = (0, 160), ylims = (0, 160), size = (1400, 700), yticks = 0:40:160, xticks = 0:40:160)
# simulations where U is varied to change delta
scatter!(δ[VaryU], Nheight_havg_tavg[VaryU], label ="Vary V₀", markersize = 8, marker=:d, color = :green)
scatter!(-1 .* δ[Varysg], Nheight_havg_tavg[Varysg], label ="Vary σ", markersize = 8, marker=:pentagon, color =:lightgreen)
# varying sigma to change the criticality 
plot!(δ, -1 .* m_Nhaδ .* δ .- b_Nhaδ, label =  @sprintf("%0.1f δ + %0.1f", m_Nhaδ, b_Nhaδ),
        color=:gray50, lw = 4)


bothplots = plot(Cp, Np,layout = (1,2))
y = 1:5
big_title = "Intrusion Thickness Scales Linearly with δ"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_Intrusion_v_Delta.png")