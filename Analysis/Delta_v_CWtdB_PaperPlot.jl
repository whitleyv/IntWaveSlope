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
filescalename = dpath * "DeltavCwtdB.jld2"

scale_file = jldopen(filescalename, "r+")

setnames = scale_file["setnames"]
δ = scale_file["δ"]
Ñ = scale_file["Ñn"]
Kₜ = scale_file["Kₜ_endAvg"]
Wₜ = scale_file["Wₜ_endAvg"]

nZ(q) = q!=0
δN = δ.^2 .* Ñ
nZidx = findall(nZ, Kₜ)

δ_cut = δ[nZidx]
Ñ_cut = Ñ[nZidx]
δN_cut = δN[nZidx]
Kₜ_cut = Kₜ[nZidx]
Wₜ_cut = Wₜ[nZidx]

@info "Find Bestfit Lines"
(b_Kδ,m_Kδ) = linear_fit(δN_cut, Kₜ_cut)
(b_Wδ,m_Wδ) = linear_fit(δN_cut, Wₜ_cut)

@info "Plotting!"

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryU = 1:9
VaryN = 10:11

#       _____
#     K|_____ []Wₜ
#        δ²N

kp = plot(δN_cut, m_Kδ .* δN_cut .+ b_Kδ, label =  "", color=:gray50, lw = 4)
        scatter!(δN_cut[VaryN], Kₜ_cut[VaryN], marker_z=Wₜ_cut[VaryN].*1e6, label ="", markersize = 8, 
            xlabel = "δ²N [m²s⁻¹]", ylabel = "K̄ [m²s⁻³]", color=:amp,
            tickfont = 15, guidefontsize = 20, titlefont = 22, legendfont = 15, 
            colorbar_titlefontsize=15,
            bottom_margin=10mm, left_margin=10mm, right_margin = 10mm, top_margin=10mm,
            legend = :topleft, marker=:c, 
            xlims = (-1, 90), ylims = (-2e-7, 2e-5), size = (1000, 800), yformatter = :scientific,
            title =  "Tracer Weighted Buoyancy Diffusivity Scales with δ²N", #clims = (0,10),
            colorbar_title="W̄×10⁻⁶ [m²s⁻¹]")
            # simulations where U is v
        scatter!(δN_cut[VaryU], Kₜ_cut[VaryU], marker_z=Wₜ_cut[VaryU].*1e6, color=:amp, 
        label ="", markersize = 8, marker=:utriangle)
        plot!(δN, -1 .* m_Kδ .* δN .- b_Kδ, label =  @sprintf("%0.3g δ²N", m_Kδ),
                color=:gray50, lw = 4)
        scatter!(-1 .* δN_cut[VaryU] .- 10, Kₜ_cut[VaryU], color=:firebrick3, label ="Vary N", markersize = 8, marker=:c)
        scatter!(-1 .* δN_cut[VaryN] .- 10, Kₜ_cut[VaryN], color=:firebrick3, label ="Vary V₀", markersize = 8, marker=:utriangle)
        

savefig(apath * "Paper_CWtdB_v_Delta.png")


(b_Wδ,m_Wδ) = linear_fit(δ_cut, Wₜ_cut)

wp = plot(δ_cut, m_Wδ .* δ_cut .+ b_Wδ, label =  "", color=:gray50, lw = 4)
        scatter!(δ_cut[VaryN], Wₜ_cut[VaryN], color =:firebrick2, label ="Vary N", markersize = 8, 
            xlabel = "δ²N [m²s⁻¹]", ylabel = "W̄ [m²s⁻¹]", 
            tickfont = 15, guidefontsize = 20, titlefont = 22, legendfont = 15, 
            bottom_margin=10mm, left_margin=10mm, right_margin = 10mm, top_margin=10mm,
            legend = :topleft, marker=:c, 
            xlims = (0, 160), ylims = (-2.5e-7, 1.2e-5), size = (1000, 800), yformatter = :scientific,
            title =  "Tracer Weighted Buoyancy Velocity Scales with δ", )
            # simulations where U is v
        scatter!(δ_cut[VaryU], Wₜ_cut[VaryU], color =:firebrick2, 
        label ="Vary V₀", markersize = 8, marker=:utriangle)
        plot!(δ, -1 .* m_Wδ .* δ .- b_Wδ .- 10, label =  @sprintf("%0.3g δ²N", m_Wδ),
                color=:gray50, lw = 4)

savefig(apath * "Paper_CWtdBW_v_Delta.png")

VaryU = 1:8
VaryN = 9:10

deleteat!(Kₜ_cut,9)
deleteat!(Wₜ_cut,9)
deleteat!(δN_cut,9)

(b_Kδ,m_Kδ) = linear_fit(δN_cut, Kₜ_cut)


kp = plot(δN_cut, m_Kδ .* δN_cut .+ b_Kδ, label =  "", color=:gray50, lw = 4)
        scatter!(δN_cut[VaryN], Kₜ_cut[VaryN], marker_z=Wₜ_cut[VaryN].*1e6, label ="", markersize = 8, 
            xlabel = "δ²N [m²s⁻¹]", ylabel = "K̄ [m²s⁻³]", color=:amp,
            tickfont = 15, guidefontsize = 20, titlefont = 22, legendfont = 15, 
            colorbar_titlefontsize=15,
            bottom_margin=10mm, left_margin=10mm, right_margin = 10mm, top_margin=10mm,
            legend = :topleft, marker=:c, 
            xlims = (-1, 90), ylims = (-2e-7, 2e-5), size = (1000, 800), yformatter = :scientific,
            title =  "Tracer Weighted Buoyancy Diffusivity Scales with δ²N", #clims = (0,10),
            colorbar_title="W̄×10⁻⁶ [m²s⁻¹]")
            # simulations where U is v
        scatter!(δN_cut[VaryU], Kₜ_cut[VaryU], marker_z=Wₜ_cut[VaryU].*1e6, color=:amp, 
        label ="", markersize = 8, marker=:utriangle)
        plot!(δN, -1 .* m_Kδ .* δN .- b_Kδ, label =  @sprintf("%0.3g δ²N", m_Kδ),
                color=:gray50, lw = 4)
        scatter!(-1 .* δN_cut[VaryU] .- 10, Kₜ_cut[VaryU], color=:firebrick3, label ="Vary N", markersize = 8, marker=:c)
        scatter!(-1 .* δN_cut[VaryN] .- 10, Kₜ_cut[VaryN], color=:firebrick3, label ="Vary V₀", markersize = 8, marker=:utriangle)
        

savefig(apath * "Paper_CWtdB_v_Delta.png")

