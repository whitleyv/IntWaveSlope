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
filescalename_oth = dpath * "DeltavAllScale_oth.jld2"
dnames = "../SetnamesList.jld2"
scale_file = jldopen(filescalename, "r+")
scale_file_oth = jldopen(filescalename_oth, "r+")
dfile = jldopen(dnames, "r+")
setnames = vcat(scale_file["setnames"],scale_file_oth["setnames_oth"])
δ = vcat(scale_file["δ"],scale_file_oth["δ_oth"])
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
Varysg_nosub = 23:24


# simulations where N is varied to change delta
Cp = plot(δ, m_Chaδ .* δ .+ b_Chaδ, label =  "", color=:gray50, lw = 4, legend=false)
scatter!(δ[Varysg], Cheight_havg_tavg[Varysg], label ="", markersize = 8, marker=:star4, color= :darkgreen)
annotate!(δ[Varysg_sub] .+ 5, Cheight_havg_tavg[Varysg_sub] .- 5, text("Subcritical γ", 10), aspect_ratio=:true)
# varying sigma without changing criticality, changing slope instead
scatter!(δ[Varys], Cheight_havg_tavg[Varys], label ="", markersize = 8, marker=:star4, color= :darkgreen)
scatter!(δ[VaryN], Cheight_havg_tavg[VaryN], label ="", markersize = 8, color =:firebrick2, xlabel = "δ [m]", ylabel = "Intrusion Thickness [m]",
tickfont = 15, guidefontsize = 20, titlefont = 22, title = "Tracer Thickness, Lₜᵣₐ", 
bottom_margin=10mm, left_margin=10mm, legend = :topleft, marker=:c,
xlims = (0, 160), ylims = (0, 160), size = (1400, 700), yticks = 0:40:160, xticks = 0:40:160)
# simulations where U is varied to change delta
scatter!(δ[VaryU], Cheight_havg_tavg[VaryU], label ="", markersize = 8, marker=:utriangle, color = :dodgerblue2)


Np = plot(δ, m_Nhaδ .* δ .+ b_Nhaδ, label =  "", color=:gray30, lw = 4)
scatter!(δ[Varysg], Nheight_havg_tavg[Varysg], label ="", markersize = 8, marker=:star4, color =:darkgreen,)
annotate!(δ[Varysg_sub].+ 18, Nheight_havg_tavg[Varysg_sub], text("Subcritical γ", 10), aspect_ratio=:true)
# varying sigma without changing criticality, changing slope instead
scatter!(δ[Varys], Nheight_havg_tavg[Varys], label ="", markersize = 8, marker=:star4, color =:darkgreen,)
scatter!(δ[VaryN], Nheight_havg_tavg[VaryN], label ="Vary N", markersize = 8, color =:firebrick2, xlabel = "δ [m]", ylabel = "",
tickfont = 15, guidefontsize = 20, titlefont = 22, title = "Stratification Anomaly Thickness, Lₙ₂", 
bottom_margin=10mm, right_margin = 10mm, legendfont = 15, legend = :topleft, marker=:c,
xlims = (0, 160), ylims = (0, 160), size = (1400, 700), yticks = false, xticks = 0:40:160)
# simulations where U is varied to change delta
scatter!(δ[VaryU], Nheight_havg_tavg[VaryU], label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
scatter!(-1 .* δ[Varysg], Nheight_havg_tavg[Varysg], label ="Vary σ", markersize = 8, marker=:star4, color =:darkgreen)
# varying sigma to change the criticality 
plot!(δ, -1 .* m_Chaδ .* δ .- b_Chaδ, label =  @sprintf("Lₜᵣₐ = %0.1f δ + %0.1f", m_Chaδ, b_Chaδ),
        color=:gray50, lw = 4)
plot!(δ, -1 .* m_Nhaδ .* δ .- b_Nhaδ, label =  @sprintf("Lₙ₂ = %0.1f δ + %0.1f", m_Nhaδ, b_Nhaδ),
        color=:gray30, lw = 4)

bothplots = plot(Cp, Np,layout = (1,2))
y = 1:5
big_title = "Intrusion Thickness Scales Linearly with δ"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_Intrusion_v_Delta1.png")


# simulations where N is varied to change delta
Cp = plot(m_Chaδ .* δ .+ b_Chaδ, δ, label =  "", color=:gray50, lw = 4, legend=:bottomright)
scatter!(Cheight_havg_tavg[Varysg_nosub], δ[Varysg_nosub],label ="", markersize = 8, marker=:star4, color= :darkgreen)
scatter!( Cheight_havg_tavg[Varys], δ[Varys], label ="", markersize = 8, marker=:star4, color= :darkgreen)
scatter!(Cheight_havg_tavg[Varysg_sub], δ[Varysg_sub],label ="", markersize = 8, marker=:star4, color= :darkgreen, markerstrokecolor=:darkgoldenrod, markerstrokewidth=1.5)
annotate!(Cheight_havg_tavg[Varysg_sub] .- 5, δ[Varysg_sub] .+ 5, text("Subcritical γ", 10, color=:darkgoldenrod), aspect_ratio=:true)
# varying sigma without changing criticality, changing slope instead
scatter!( Cheight_havg_tavg[VaryN], δ[VaryN],label ="Vary N", markersize = 8, color =:firebrick2, ylabel = "δ [m]", xlabel = "Tracer Thickness, Lₜᵣₐ [m]",
tickfont = 15, guidefontsize = 20, titlefont = 22, legendfont = 15,title = "", 
bottom_margin=10mm, left_margin=10mm, marker=:c,
xlims = (0, 160), ylims = (0, 160), size = (1400, 700), yticks = 0:40:160, xticks = 40:40:160)
# simulations where U is varied to change delta
scatter!(Cheight_havg_tavg[VaryU], δ[VaryU], label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
scatter!(-1 .* δ[Varysg], Nheight_havg_tavg[Varysg], label ="Vary σ", markersize = 8, marker=:star4, color =:darkgreen)
# varying sigma to change the criticality 
plot!( -1 .* m_Chaδ .* δ .- b_Chaδ, δ,label =  @sprintf("Lₜᵣₐ = %0.1f δ + %0.1f", m_Chaδ, b_Chaδ),
        color=:gray50, lw = 4)
plot!( -1 .* m_Nhaδ .* δ .- b_Nhaδ, δ, label =  @sprintf("Lₙ₂ = %0.1f δ + %0.1f", m_Nhaδ, b_Nhaδ),
        color=:gray30, lw = 4)

Np = plot(m_Nhaδ .* δ .+ b_Nhaδ, δ, label =  "", color=:gray30, lw = 4)
scatter!(Nheight_havg_tavg[Varysg_nosub], δ[Varysg_nosub],label ="", markersize = 8, marker=:star4, color= :darkgreen)
scatter!(Nheight_havg_tavg[Varysg_sub], δ[Varysg_sub],label ="", markersize = 8, marker=:star4, color= :darkgreen, markerstrokecolor=:darkgoldenrod, markerstrokewidth=1.5)
annotate!(Nheight_havg_tavg[Varysg_sub[1]] .+10 ,δ[Varysg_sub[1]].+ 5, text("Subcritical γ", 10,color=:darkgoldenrod), aspect_ratio=:true)
annotate!(Nheight_havg_tavg[Varysg_sub[2]] .- 5 ,δ[Varysg_sub[2]].+ 5, text("Subcritical γ", 10,color=:darkgoldenrod), aspect_ratio=:true)
# varying sigma without changing criticality, changing slope instead
scatter!(Nheight_havg_tavg[Varys], δ[Varys], label ="", markersize = 8, marker=:star4, color =:darkgreen,)
scatter!( Nheight_havg_tavg[VaryN],δ[VaryN], label ="Vary N", markersize = 8, color =:firebrick2, xlabel = "Stratification Anomaly Thickness, Lₙ₂",
tickfont = 15, guidefontsize = 20, titlefont = 22, title = "", left_margin = -5mm,
bottom_margin=10mm, right_margin = 10mm, legendfont = 15, legend = false, marker=:c,
xlims = (0, 160), ylims = (0, 160), size = (1400, 700), yticks = false, xticks = 40:40:160)
# simulations where U is varied to change delta
scatter!( Nheight_havg_tavg[VaryU], δ[VaryU],label ="Vary V₀", markersize = 8, marker=:utriangle, color = :dodgerblue2)
scatter!(-1 .* δ[Varysg], Nheight_havg_tavg[Varysg], label ="Vary σ", markersize = 8, marker=:star4, color =:darkgreen)


bothplots = plot(Cp, Np,layout = (1,2))
y = 1:5
big_title = "Intrusion Thickness Scales Linearly with δ"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))

savefig(apath * "Paper_Intrusion_v_Delta2.png")


# simulations where N is varied to change delta
Cp = plot(m_Chaδ .* δ .+ b_Chaδ, δ, label =  "", color=:gray50, lw = 4)
scatter!( Cheight_havg_tavg[VaryN], δ[VaryN], label ="Vary N, V₀=0.25 ms⁻¹", markersize = 8, color =:firebrick2, ylabel = "δ [m]", xlabel = "Lₜᵣₐ [m]",
tickfont = 15, guidefontsize = 17, titlefont = 17, legendfont = 13, title = "Tracer Thickness, Lₜᵣₐ", 
bottom_margin=5mm, left_margin=15mm, legend = :bottomright, marker=:c,
ylims = (0, 160), xlims = (0, 160), size = (1400, 700), xticks = 0:40:160, yticks = 0:40:160)
# simulations where U is varied to change delta
scatter!( Cheight_havg_tavg[VaryU], δ[VaryU], label ="Vary V₀, N=3.5×10⁻³s⁻¹", markersize = 8, marker=:utriangle, color = :dodgerblue2)
# varying sigma to change the criticality 
plot!( -1 .* m_Chaδ .* δ .- b_Chaδ, δ, label =  @sprintf("Lₜᵣₐ = %0.1f δ + %0.1f", m_Chaδ, b_Chaδ),
        color=:gray50, lw = 4)
plot!( -1 .* m_Nhaδ .* δ .- b_Nhaδ,  δ, label =  @sprintf("Lₙ₂ = %0.1f δ + %0.1f", m_Nhaδ, b_Nhaδ),
        color=:gray50, lw = 4)

Np = plot(m_Nhaδ .* δ .+ b_Nhaδ, δ, label =  "", color=:gray50, lw = 4)
scatter!( Nheight_havg_tavg[VaryN],δ[VaryN], label ="", markersize = 8, color =:firebrick2, yticks = false, xlabel = "Lₙ₂ [m]",
tickfont = 15, guidefontsize = 17, titlefont = 17, title = "Stratification Anomaly Thickness, Lₙ₂", 
bottom_margin=10mm, right_margin = 10mm, legend = false, marker=:c,
ylims = (0, 160), xlims = (0, 160), size = (1400, 700), xticks = 40:40:160)
# simulations where U is varied to change delta
scatter!(Nheight_havg_tavg[VaryU], δ[VaryU], label ="", markersize = 8, marker=:utriangle, color = :dodgerblue2)

bothplots = plot(Cp, Np,layout = (1,2))
y = 1:5
big_title = "Mean Thickness of Interior Intrusions Scales Linearly with δ=V₀/N"
BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
Plots.text(big_title, 22)), axis=nothing, grid=false, leg=false,
foreground_color_subplot=colorant"white")

fin = plot(BigT, bothplots, layout=grid(2, 1, heights=[0.05,0.95]))
