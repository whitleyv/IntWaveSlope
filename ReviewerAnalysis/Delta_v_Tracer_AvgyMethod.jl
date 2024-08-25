using JLD2
using Statistics
using CairoMakie
using Printf
using GeometryBasics

dpath = "ReviewerAnalysis/ReviewerData/"
apath = "ReviewerAnalysis/ReviewerPlots/"

filescalename = dpath * "DeltavAll_ExtraStats.jld2"
filesetnames =  "SetList_mp.jld2"

file_sn = jldopen(filesetnames, "r+")
file_scale = jldopen(filescalename, "r+")

Cheight_stats = file_scale["Cheight_stats"]

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
γ_σ = file_sn["γ_varyσ"]

idx_subcritical_σ = findall(γ_σ .< 1)
idx_gammachange_σ = findall(γ_σ .!= 1.9)
idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10

Cheight_Avgy = [39.1, 154.62, 190, 246.88]
idx_Avgy = [2, 5, 7, 9]
U_Avgy = U_UN[idx_Avgy]
N_Avgy = N_UN[1] .* ones(length(idx_Avgy))
δ_Avgy = δ_UN[idx_Avgy]

f = Figure(resolution = (700, 800), fontsize=26)
ga = f[1, 1] = GridLayout()

ax1 = Axis(ga[2, 1], ylabel = "Lₜᵣ",  xlabel = rich("h", subscript("w")),
xlabelsize=30, ylabelsize=30)
ax1.xticks = 0:40:300
ax1.yticks = 0:40:300
limits!(ax1, 0,160,0,280) 

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

# mean + 80/20 percentile + extrema
vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
        color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
        marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
v_Avgy = scatter!(ax1, δ_Avgy, Cheight_Avgy, markersize = 25, marker=:dtriangle, 
        color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
Legend( ga[1, 1],  [vnp, vump, vsp1, vsbp, v_Avgy], ["Vary N₀", "Vary V₀", "Vary γ", "Subcritical", "Avg y"],
        tellheight = false, tellwidth = false, rowgap = 20, colgap = 30,
        margin = (10, 10, 20, 5), framevisible = false, patchlabelgap = 7,
        labelsize = 25, nbanks = 2,
        halign = :center, valign = :top, orientation = :horizontal)

rowsize!(ga,1, Auto(0.05))       

save(apath * "Delta_v_Tracer_Avgy.png", f) 
