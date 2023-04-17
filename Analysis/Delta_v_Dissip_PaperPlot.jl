using Plots
using Measures
using Statistics
using Printf
#using Oceananigans
using CurveFit
#using ArgParse
using JLD2

#ENV["GKSwstype"] = "nul" # if on remote HPC

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename_ε = dpath * "DeltavDissip.jld2"
filescalename_Lt = dpath * "DeltavAllScale.jld2"

scale_file_ε = jldopen(filescalename_ε, "r+")
scale_file_Lt = jldopen(filescalename_Lt, "r+")

setnames = scale_file_ε["setnames"]
δ = scale_file_ε["δ"]
eps_beginAvg = scale_file_ε["eps_beginAvg"]
eps_endAvg = scale_file_ε["eps_endAvg"]
thorpe_hmax_tavg = scale_file_Lt["thorpe_hmax_tavg"]
# Thorpe vals did vary N first, switching so same as dissipation
Lt_corrected = vcat(thorpe_hmax_tavg[12:end], thorpe_hmax_tavg[1:11])

Ñ_U = ones(11) .* (3.5*1e-3)
Ñ_N = 0.25 ./ δ[12:22]
Ñ = vcat(Ñ_U, Ñ_N)
δN = δ.^2 .* Ñ.^3
LtN = Lt_corrected.^2 .* Ñ.^3

deleteat!(LtN, 12:13)
deleteat!(δN, 12:13)
deleteat!(eps_endAvg, 12:13)
@info "Find Bestfit Lines"
(b_εδ,m_εδ) = linear_fit(δN, eps_endAvg)
(b_εLt,m_εLt) = linear_fit(LtN, eps_endAvg)

@info "Plotting!"

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryU = 1:11
VaryN = 12:20


# simulations where N is varied to change delta
dp = plot(δN, m_εδ .* δN .+ b_εδ, label =  "", color=:gray50, lw = 4)
scatter!(δN[VaryN], eps_endAvg[VaryN], label ="Vary N", markersize = 8, color =:green, xlabel = "δ²N³ [m²s⁻³]", yticks=false,
    tickfont = 15, guidefontsize = 20, titlefont = 20, legendfont = 15, title = "", 
    bottom_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-5, 1.1*1e-3), ylims = (-1e-6, 1.7*1e-5), size = (1500, 600), xformatter = :scientific,
    xticks = 0:5*1e-4:1e-3)
# simulations where U is varied to change delta
scatter!(δN[VaryU], eps_endAvg[VaryU], label ="Vary V₀", markersize = 8, marker=:d, color = :green)
plot!(LtN, -1 .* m_εLt .* LtN .- 5, label =  @sprintf("%0.3f Lₜ²N³", m_εLt),
        color=:gray50, lw = 4)
plot!(δ, -1 .* m_εδ .* δ .- b_εδ, label =  @sprintf("%0.3f δ²N³", m_εδ),
        color=:gray50, lw = 4)

lp = plot(LtN, m_εLt .* LtN .+ b_εLt, color=:gray50, lw = 4)
scatter!(LtN[VaryN], eps_endAvg[VaryN], markersize = 8, color =:green, xlabel = "Lₜ²N³ [m²s⁻³]", ylabel = "ε [m²s⁻³]",
    tickfont = 15, guidefontsize = 20, titlefont = 20, title = "", 
    bottom_margin=10mm, left_margin=10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-6, 1.4*1e-4), ylims = (-1e-6, 1.7*1e-5), size = (1500, 600), xformatter = :scientific,
    xticks = 0:5*1e-5:1e-4)
# simulations where U is varied to change delta
scatter!(LtN[VaryU], eps_endAvg[VaryU],  markersize = 8, marker=:d, color = :green, legend = false)



bothplots = plot(lp, dp,)
y = 1:5
big_title = "Mean Dissipation Scales With Energetically Motivated Terms"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_DissipLt_v_Delta.png")