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
filesetnames =  "SetnamesList.jld2"

scale_file_sn = jldopen(filesetnames, "r+")
δ = scale_file_sn["δs"][1:25]

scale_file_ε = jldopen(filescalename_ε, "r+")
scale_file_Lt = jldopen(filescalename_Lt, "r+")

setnames = scale_file_ε["setnames"]

include("parameters.jl")

Ñ = zeros(length(setnames))
for (m,setname) in enumerate(setnames)
    pm = getproperty(SimParams(), Symbol(setname))
    Ñ[m] = pm.Ñ
end

eps_beginAvg = scale_file_ε["eps_beginAvg"]
eps_endAvg = scale_file_ε["eps_endAvg"]
eps_beginMax = scale_file_ε["eps_beginMaxAvg"]
eps_endMax = scale_file_ε["eps_endMaxAvg"]
thorpe_hmax_tavg = scale_file_Lt["thorpe_hmax_tavg"]

δN = δ.^2 .* Ñ.^3
LtN = thorpe_hmax_tavg.^2 .* Ñ.^3

@info "Find Bestfit Lines for full data"
(b_εδ,m_εδ) = linear_fit(δN, eps_endAvg)
(b_εLt,m_εLt) = linear_fit(LtN, eps_endAvg)
(b_Mεδ,m_Mεδ) = linear_fit(δN, eps_endMax)
(b_MεLt,m_MεLt) = linear_fit(LtN, eps_endMax)

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryU = 1:11
VaryN = 12:20

VaryU = 1:11
VaryN = 12:25

# removing last N and first U val from data
VaryU2 = 1:10
VaryN2 = 13:22

# indices now avai;able for two data sets after cutting out above points
VaryU_cut = 1:10
VaryN_cut = 11:20

@info "Find Bestfit Lines for data cut"

eps_endAvg_cut = vcat(eps_endAvg[VaryU2], eps_endAvg[VaryN2])
eps_endMax_cut = vcat(eps_endMax[VaryU2], eps_endMax[VaryN2])
δN_cut = vcat(δN[VaryU2], δN[VaryN2])
LtN_cut = vcat(LtN[VaryU2], LtN[VaryN2])

(b_εδ_cut,m_εδ_cut) = linear_fit(δN_cut, eps_endAvg_cut)
(b_εLt_cut,m_εLt_cut) = linear_fit(LtN_cut, eps_endAvg_cut)
(b_Mεδ_cut,m_Mεδ_cut) = linear_fit(δN_cut, eps_endMax_cut)
(b_MεLt_cut,m_MεLt_cut) = linear_fit(LtN_cut, eps_endMax_cut)

@info "Plotting..."

""" Plotting Full Plots for Avg Dissipation"""

xlimLt =maximum(LtN)
xlimd = maximum(δN)
xmind = minimum(δN)
xminLt = minimum(LtN)
ylime = maximum(eps_endAvg)
ymine = minimum(eps_endAvg)

# simulations where N is varied to change delta
dp = plot(δN, m_εδ .* δN .+ b_εδ, label =  "", color=:gray30, lw = 4)
scatter!(δN[VaryN], eps_endAvg[VaryN], label ="Vary N", markersize = 8, 
    color =:firebrick2, xlabel = "δ²N³ [m²s⁻³]", yticks=false,
    tickfont = 15, guidefontsize = 20, titlefont = 20, legendfont = 15, title = "", 
    bottom_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-4, 1.2*1e-3), ylims = (-1e-6, 3.8e-5), size = (1500, 600), xformatter = :scientific,
    xticks = 0:5*1e-4:1e-3)
# simulations where U is varied to change delta
scatter!(δN[VaryU], eps_endAvg[VaryU], label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
plot!(LtN, -1 .* m_εLt .* LtN .- 5, label =  @sprintf("%0.3f Lₜ²N³", m_εLt),
        color=:gray50, lw = 4)
plot!(δ, -1 .* m_εδ .* δ .- b_εδ, label =  @sprintf("%0.3f δ²N³", m_εδ),
        color=:gray30, lw = 4)

lp = plot(LtN, m_εLt .* LtN .+ b_εLt, color=:gray50, lw = 4)
scatter!(LtN[VaryN], eps_endAvg[VaryN], markersize = 8, color =:firebrick2, xlabel = "Lₜ²N³ [m²s⁻³]", ylabel = "⟨ε⟩ [m²s⁻³]",
    tickfont = 15, guidefontsize = 20, titlefont = 20, title = "", 
    bottom_margin=10mm, left_margin=10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-4, 1.2*1e-3), ylims = (-1e-6, 3.8*1e-5), size = (1500, 600), xformatter = :scientific,
    xticks = 0:5*1e-4:1e-3)
# simulations where U is varied to change delta
scatter!(LtN[VaryU], eps_endAvg[VaryU],  markersize = 8, marker=:utriangle, color = :dodgerblue2, legend = false)

bothplots = plot(lp, dp,)
y = 1:5
big_title = "Mean Dissipation Scales With Energetically Motivated Terms"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_AvgDissip_v_LtDelta_Full.png")

""" Plotting Cut Plots for Avg Dissipation"""

xlimLt_cut =maximum(LtN_cut)
xlimd_cut = maximum(δN_cut)
ylime_cut = maximum(eps_endAvg_cut)

# simulations where N is varied to change delta
dp = plot(δN_cut, m_εδ_cut .* δN_cut .+ b_εδ_cut, label =  "", color=:gray30, lw = 4)
scatter!(δN_cut[VaryN_cut], eps_endAvg_cut[VaryN_cut], label ="Vary N", markersize = 8, color =:firebrick2, xlabel = "δ²N³ [m²s⁻³]", yticks=false,
    tickfont = 15, guidefontsize = 20, titlefont = 20, legendfont = 15, title = "", 
    bottom_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-6, 1.0*1e-3), ylims = (-1e-6, 1.2*1e-5), size = (1500, 600), xformatter = :scientific,
    xticks = 0:5*1e-4:1e-3)
# simulations where U is varied to change delta
scatter!(δN_cut[VaryU_cut], eps_endAvg_cut[VaryU_cut], label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
plot!(LtN_cut, -1 .* m_εLt_cut .* LtN_cut .- 5, label =  @sprintf("%0.3f Lₜ²N³", m_εLt_cut),
        color=:gray50, lw = 4)
plot!(δN_cut, -1 .* m_εδ_cut .* δN_cut .- 5, label =  @sprintf("%0.3f δ²N³", m_εδ_cut),
        color=:gray30, lw = 4)

lp = plot(LtN_cut, m_εLt_cut .* LtN_cut .+ b_εLt_cut, color=:gray50, lw = 4)
scatter!(LtN_cut[VaryN_cut], eps_endAvg_cut[VaryN_cut], markersize = 8, color =:firebrick2, xlabel = "Lₜ²N³ [m²s⁻³]", ylabel = "⟨ε⟩ [m²s⁻³]",
    tickfont = 15, guidefontsize = 20, titlefont = 20, title = "", 
    bottom_margin=10mm, left_margin=10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-6, 3.0*1e-4), ylims = (-1e-6, 1.2*1e-5), size = (1500, 600), xformatter = :scientific,
    xticks = 0:1*1e-4:2e-4)
# simulations where U is varied to change delta
scatter!(LtN_cut[VaryU_cut], eps_endAvg_cut[VaryU_cut],  markersize = 8, marker=:utriangle, color = :dodgerblue2, legend = false)

bothplots = plot(lp, dp,)
y = 1:5
big_title = "Mean Dissipation Scales With Energetically Motivated Terms"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_AvgDissip_v_LtDelta_Cut.png")


""" Plotting Full Plots for Max Dissipation"""

ylime = maximum(eps_endMax)

# simulations where N is varied to change delta
dp = plot(δN, m_Mεδ .* δN .+ b_Mεδ, label =  "", color=:gray30, lw = 4)
scatter!(δN[VaryN], eps_endMax[VaryN], label ="Vary N", markersize = 8, color =:firebrick2, xlabel = "δ²N³ [m²s⁻³]", yticks=false,
    tickfont = 15, guidefontsize = 20, titlefont = 20, legendfont = 15, title = "", 
    bottom_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-5, 1.2*1e-3), ylims = (-1e-5, 1.2*1e-3), size = (1500, 600), xformatter = :scientific, yformatter = :scientific,
    xticks = 0:5*1e-4:1e-3)
# simulations where U is varied to change delta
scatter!(δN[VaryU], eps_endMax[VaryU], label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
plot!(LtN, -1 .* m_MεLt .* LtN .- 5, label =  @sprintf("%0.3f Lₜ²N³", m_MεLt), 
        color=:gray50, lw = 4)
plot!(δ, -1 .* m_Mεδ .* δ .- 5, label =  @sprintf("%0.3f δ²N³", m_Mεδ),
        color=:gray30, lw = 4)

lp = plot(LtN, m_MεLt .* LtN .+ b_MεLt, color=:gray50, lw = 4)
scatter!(LtN[VaryN], eps_endMax[VaryN], markersize = 8, color =:firebrick2, xlabel = "Lₜ²N³ [m²s⁻³]", ylabel = "⟨εₘₐₓ⟩ [m²s⁻³]",
    tickfont = 15, guidefontsize = 20, titlefont = 20, title = "", 
    bottom_margin=10mm, left_margin=10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-5, 1.2*1e-3), ylims = (-1e-5, 1.2*1e-3), size = (1500, 600), xformatter = :scientific, yformatter = :scientific,
    xticks = 0:5*1e-4:1e-3)
# simulations where U is varied to change delta
scatter!(LtN[VaryU], eps_endMax[VaryU],  markersize = 8, marker=:utriangle, color = :dodgerblue2, legend = false)

bothplots = plot(lp, dp,)
y = 1:5
big_title = "Max Dissipation Scales With Energetically Motivated Terms"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_MaxDissip_v_LtDelta_Full.png")

ylime = maximum(eps_endMax_cut)

# simulations where N is varied to change delta
dp = plot(δN_cut, m_Mεδ_cut .* δN_cut .+ b_Mεδ_cut, label =  "", color=:gray30, lw = 4)
scatter!(δN_cut[VaryN_cut], eps_endMax_cut[VaryN_cut], label ="Vary N", markersize = 8, color =:firebrick2, xlabel = "δ²N³ [m²s⁻³]", yticks=false,
    tickfont = 15, guidefontsize = 20, titlefont = 20, legendfont = 15, title = "", 
    bottom_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-6, 1.0*1e-3), ylims = (-1e-5, 5*1e-4), size = (1500, 600), xformatter = :scientific,
    xticks = 0:3e-4:6e-4)
# simulations where U is varied to change delta
scatter!(δN_cut[VaryU_cut], eps_endMax_cut[VaryU_cut], label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
plot!(LtN_cut, -1 .* m_MεLt_cut .* LtN_cut .- 5, label =  @sprintf("%0.3f Lₜ²N³", m_MεLt_cut),
        color=:gray50, lw = 4)
plot!(δN_cut, -1 .* m_Mεδ_cut .* δN_cut .- 5, label =  @sprintf("%0.3f δ²N³", m_Mεδ_cut),
        color=:gray30, lw = 4)

lp = plot(LtN_cut, m_MεLt_cut .* LtN_cut .+ b_MεLt_cut, color=:gray50, lw = 4)
scatter!(LtN_cut[VaryN_cut], eps_endMax_cut[VaryN_cut], markersize = 8, color =:firebrick2, xlabel = "Lₜ²N³ [m²s⁻³]", ylabel = "ε [m²s⁻³]",
    tickfont = 15, guidefontsize = 20, titlefont = 20, title = "", 
    bottom_margin=10mm, left_margin=10mm, legend = :topleft, marker=:c, 
    xlims = (-1e-6, 3.0*1e-4), ylims = (-1e-5, 5*1e-4), size = (1500, 600), xformatter = :scientific,
    xticks = 0:1*1e-4:2e-4)
# simulations where U is varied to change delta
scatter!(LtN_cut[VaryU_cut], eps_endMax_cut[VaryU_cut],  markersize = 8, marker=:utriangle, color = :dodgerblue2, legend = false)

bothplots = plot(lp, dp,)
y = 1:5
big_title = "Max Dissipation Scales With Energetically Motivated Terms"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_MaxDissip_v_LtDelta_Cut.png")