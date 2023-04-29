using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using CurveFit
using ArgParse
using JLD2

#ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename = dpath * "DeltavAllScale.jld2"

scale_file = jldopen(filescalename, "r+")

setnames = vcat(scale_file["setnames"],scale_file["setnames_oth"])
δ = vcat(scale_file["δ"],scale_file["δ_oth"])
thorpe_hmax_tavg = vcat(scale_file["thorpe_hmax_tavg"],scale_file["thorpe_hmax_tavg_oth"])

@info "Find Bestfit Lines"

(b_Ltmδ,m_Ltmδ) = linear_fit(δ, thorpe_hmax_tavg)


@info "Plotting!"

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
plot(δ, m_Ltmδ .* δ .+ b_Ltmδ, label =  "", color=:gray50, lw = 4)
scatter!(δ[Varysg], thorpe_hmax_tavg[Varysg], label ="", markersize = 8, marker=:star4, color =:darkgreen)
scatter!(δ[Varys], thorpe_hmax_tavg[Varys], label ="", markersize = 8, marker=:star4, color= :darkgreen)
scatter!(δ[Varysg_sub], thorpe_hmax_tavg[Varysg_sub],label ="", markersize = 8, marker=:star4, color= :darkgreen, markerstrokecolor=:darkgoldenrod, markerstrokewidth=1.5)
annotate!(δ[Varysg_sub], thorpe_hmax_tavg[Varysg_sub] .+ 2, text("Subcritical γ", 10, color=:darkgoldenrod))
# varying sigma without changing criticality, changing slope instead
scatter!(δ[VaryN], thorpe_hmax_tavg[VaryN], label ="Vary N", markersize = 8, color =:firebrick2, xlabel = "δ [m]", ylabel = "Lₜ [m]",
tickfont = 15, guidefontsize = 20, titlefont = 22, title = "Maximum RMS Thorpe Scale (Lₜ) Scales Linearly with δ", 
bottom_margin=5mm, left_margin=10mm, right_margin = 10mm, legendfont=15, legend = :topleft, marker=:c,
xlims = (0, 160), ylims = (0, 70), size = (1100, 700), yticks = 20:20:60, xticks = 0:40:160)
# simulations where U is varied to change delta
scatter!(δ[VaryU], thorpe_hmax_tavg[VaryU], label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
scatter!(-1 .* δ[Varysg], thorpe_hmax_tavg[Varysg], label ="Vary σ", markersize = 8, marker=:star4, color =:darkgreen)
# varying sigma to change the criticality 
plot!(δ, -1 .* m_Ltmδ .* δ .- b_Ltmδ, label =  @sprintf("%0.1f δ + %0.1f", m_Ltmδ, b_Ltmδ),
        color=:gray50, lw = 4)

savefig(apath * "Paper_Thorpe_v_Delta_All.png")
        