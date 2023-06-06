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
filescalename1 = dpath * "DeltavCwtdB_p1.jld2"
filescalename2 = dpath * "DeltavCwtdB_p2.jld2"
filescalename3 = dpath * "DeltavCwtdB_p3.jld2"

scale_file1 = jldopen(filescalename1, "r+")
scale_file2 = jldopen(filescalename2, "r+")
scale_file3 = jldopen(filescalename3, "r+")

setnames = vcat(scale_file1["setnames"], scale_file2["setnames"], scale_file3["setnames"])
δ = vcat(scale_file1["δ"], scale_file2["δ"], scale_file3["δ"])
Ñ = vcat(scale_file1["Ñn"], scale_file2["Ñn"], scale_file3["Ñn"])
Kₜ = vcat(scale_file1["Kₜ_endAvg"], scale_file2["Kₜ_endAvg"], scale_file3["Kₜ_endAvg"])
Wₜ = vcat(scale_file1["Wₜ_endAvg"], scale_file2["Wₜ_endAvg"], scale_file3["Wₜ_endAvg"])

δN = δ.^2 .* Ñ

@info "Find Bestfit Lines"
(b_Kδ,m_Kδ) = linear_fit(δN, Kₜ)
(b_Wδ,m_Wδ) = linear_fit(δ, Wₜ)

@info "Plotting!"

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryU = 1:11
VaryN = 12:25

#       _____
#     K|_____ []Wₜ
#        δ²N

kp = plot(δN, m_Kδ .* δN .+ b_Kδ, label =  "", color=:gray50, lw = 4)
        scatter!(δN[VaryN], Kₜ[VaryN], marker_z=Wₜ[VaryN].*1e6, label ="", markersize = 8, 
            xlabel = "δ²N [m²s⁻¹]", ylabel = "K̄ [m²s⁻³]", color=:amp,
            tickfont = 15, guidefontsize = 20, titlefont = 22, legendfont = 15, 
            colorbar_titlefontsize=15,
            bottom_margin=10mm, left_margin=10mm, right_margin = 10mm, top_margin=10mm,
            legend = :topleft, marker=:c, 
            xlims = (-1, 90), ylims = (-2e-7, 3e-6), size = (1000, 800), yformatter = :scientific,
            title =  "Tracer Weighted Buoyancy Diffusivity Scales with δ²N", #clims = (0,10),
            colorbar_title="W̄×10⁻⁶ [m²s⁻¹]")
            # simulations where U is v
        scatter!(δN[VaryU], Kₜ[VaryU], marker_z=Wₜ[VaryU].*1e6, color=:amp, 
        label ="", markersize = 8, marker=:utriangle)
        plot!(δN, -1 .* m_Kδ .* δN .- b_Kδ, label =  @sprintf("%0.3g δ²N", m_Kδ),
                color=:gray50, lw = 4)
        scatter!(-1 .* δN[VaryU] .- 10, Kₜ[VaryU], color=:firebrick3, label ="Vary N", markersize = 8, marker=:c)
        scatter!(-1 .* δN[VaryN] .- 10, Kₜ[VaryN], color=:firebrick3, label ="Vary V₀", markersize = 8, marker=:utriangle)
        
savefig(apath * "Paper_CWtdB_v_Delta.png")


kp = plot(δN, m_Kδ .* δN .+ b_Kδ, label =  "", color=:gray50, lw = 4)
scatter!(δN[VaryN], Kₜ[VaryN], label ="Vary N", markersize = 8, 
    xlabel = "δ²N [m²s⁻¹]", ylabel = "K̄ [m²s⁻³]", color=:firebrick2,
    tickfont = 15, guidefontsize = 20, titlefont = 22, legendfont = 15, 
    bottom_margin=10mm, left_margin=10mm, right_margin = 10mm, top_margin=10mm,
    legend = :topright, marker=:c, 
    xlims = (-1, 90), ylims = (-2e-7, 3e-6), size = (1000, 800), yformatter = :scientific,
    title =  "Tracer Weighted Buoyancy Diffusivity vs δ²N")    # simulations where U is v
scatter!(δN[VaryU], Kₜ[VaryU], color=:dodgerblue2, 
label ="Vary V₀", markersize = 8, marker=:utriangle)
plot!(δN, -1 .* m_Kδ .* δN .- 100, label =  @sprintf("%0.3g δ²N", m_Kδ),
        color=:gray50, lw = 4)

Wₜcut = vcat(Wₜ[1:11], Wₜ[14:22])
δcut = vcat(δ[1:11], δ[14:22])
VaryU = 1:11
VaryN = 12:20

(b_Wδ,m_Wδ) = linear_fit(δcut, Wₜcut)

wp = plot(δcut, m_Wδ .* δcut .+ b_Wδ, label =  "", color=:gray50, lw = 4)
        scatter!(δcut[VaryN], Wₜcut[VaryN], color =:firebrick2, label ="Vary N", markersize = 8, 
            xlabel = "δ [m]", ylabel = "W̄ [m²s⁻¹]", 
            tickfont = 15, guidefontsize = 20, titlefont = 22, legendfont = 15, 
            bottom_margin=10mm, left_margin=10mm, right_margin = 10mm, top_margin=10mm,
            legend = :topleft, marker=:c, 
            xlims = (0, 160), ylims = (-2.5e-7, 6e-7), size = (1000, 800), yformatter = :scientific,
            title =  "Tracer Weighted Buoyancy Velocity vs δ", )
            # simulations where U is v
        scatter!(δcut[VaryU], Wₜcut[VaryU], color =:dodgerblue2, 
        label ="Vary V₀", markersize = 8, marker=:utriangle)
        plot!(δcut, -1 .* m_Wδ .* δcut .- b_Wδ .- 10, label =  @sprintf("%0.3g δ", m_Wδ),
                color=:gray50, lw = 4)

savefig(apath * "Paper_CWtdBW_v_Delta.png")

VaryU = 1:8
VaryN = 9:10

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

