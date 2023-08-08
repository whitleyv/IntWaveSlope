using Statistics
using Printf
using CurveFit
using JLD2
using CairoMakie
using GeometryBasics

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename = dpath * "DeltavAllScale.jld2"
filescalename_oth = dpath * "DeltavAllScale_oth.jld2"
dnames = "../SetnamesList.jld2"

scale_file_sn = jldopen(dnames, "r+")
δs = scale_file_sn["δs"][1:22]

scale_file = jldopen(filescalename, "r+")
scale_file_oth = jldopen(filescalename_oth, "r+")

setnames = vcat(scale_file["setnames"][1:22],scale_file_oth["setnames_oth"])
δ = vcat(δs,scale_file_oth["δ_oth"])

Cheight_havg_tavg =vcat(scale_file["Cheight_havg_tavg"][1:22],scale_file_oth["Cheight_havg_tavg_oth"])
Nheight_havg_tavg =vcat(scale_file["Nheight_havg_tavg"][1:22],scale_file_oth["Nheight_havg_tavg_oth"])

@info "Find Bestfit Lines"

############################### INTRUSIONS v Delta
(b_Chaδ,m_Chaδ) = linear_fit(δ, Cheight_havg_tavg)
################################ Isopycnals v Delta
(b_Nhaδ,m_Nhaδ) = linear_fit(δ, Nheight_havg_tavg)

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryN = 1:11
VaryU = 12:22
# ones that are subcritical out of sigma varying set
Varysg_sub = 25:26
# all the ones where sigma is varied in any way
Varyoth = 23:30     


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





f = Figure(resolution = (1000, 1400), fontsize=26)
ga = f[1, 1] = GridLayout()

ax1 = Axis(ga[1, 1],  xlabel = "δ [m]", ylabel = "Tracer Thickness, Lₜᵣ [m]")
ax1.xticks = 0:40:160
ax1.yticks = 0:40:160
limits!(ax1, 0, 165, 0, 165)

ax2 = Axis(ga[2, 1],  xlabel = "δ [m]", ylabel = "Stratification Anomaly Thickness, Lₙ₂ [m]")
ax2.xticks = 0:40:160
ax2.yticks = 0:40:160
limits!(ax2, 0, 165, 0, 165)

hidexdecorations!(ax1, grid = false)

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

vsp = scatter!(ax1, δ[Varyoth], Cheight_havg_tavg[Varyoth], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax1, δ[Varysg_sub], Cheight_havg_tavg[Varysg_sub], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax1, δ[VaryN], Cheight_havg_tavg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax1, δ[VaryU], Cheight_havg_tavg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

scatter!(ax2, δ[Varyoth], Nheight_havg_tavg[Varyoth], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δ[Varysg_sub], Nheight_havg_tavg[Varysg_sub], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
scatter!(ax2, δ[VaryN], Nheight_havg_tavg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δ[VaryU], Nheight_havg_tavg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend( ga[1, 1, Top()], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 15), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :bottom, orientation = :horizontal)

Label(ga[1, 1, TopLeft()], "a",
                fontsize = 30,
                font = :bold,
                padding = (-10, 5, 5, 10),
                halign = :right)
Label(ga[2, 1, TopLeft()], "b",
                fontsize = 30,
                font = :bold,
                padding = (0, 5, 5, 10),
                halign = :right)
savename = apath * "Paper_Intrusion_v_Delta_All"
save(savename * ".png", f, px_per_unit = 2)
