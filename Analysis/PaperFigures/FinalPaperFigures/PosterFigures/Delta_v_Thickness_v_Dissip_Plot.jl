using Statistics
using Printf
using CurveFit
using JLD2
using CairoMakie
using GeometryBasics

apath = "Analysis/PaperFigures/FinalPaperFigures/PosterFigures/"
dpath = "Data/"

dfilescalename_UN = dpath * "DeltavAllScale_mp.jld2"
dfilescalename_σ = dpath * "DeltavAllScale_mp_VarySigma.jld2"
efileεname_U = dpath * "DeltavAll_mtest.jld2"
efileεname_UN = dpath * "DeltavDissip_mp.jld2"
efileεname_σ = dpath * "DeltavDissip_mp_varyσ.jld2"

filesetnames =  "SetList_mp.jld2"

# getting the surrounding simulation params out
file_sn = jldopen(filesetnames, "r+")

δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
N_σ = file_sn["N_varyσ"]
γ_σ = file_sn["γ_varyσ"]

δN_UN = δ_UN.^2 .* N_UN.^3
δN_σ = δ_σ.^2 .* N_σ.^3

# determining which data sets are where:
idx_subcritical_σ = findall(γ_σ .< 1)
idx_Nsubcritical_gammachange_σ = findall((γ_σ .> 1) .& (γ_σ .< 1.9))
idx_topochange_σ = findall(γ_σ .== 1.9)
idx_gammachange_σ = findall(γ_σ .!= 1.9)
idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10
idx_VaryU_n05 = (idx_VaryU_start+1):idx_VaryU_start+10
idx_VaryN_n05 = (idx_VaryN_start+1):idx_VaryN_start+10

# pulling data 
file_scale_UN = jldopen(dfilescalename_UN, "r+")
file_scale_σ = jldopen(dfilescalename_σ, "r+")
file_ε_UN = jldopen(efileεname_UN, "r+")
file_ε_U = jldopen(efileεname_U, "r+")
file_ε_σ = jldopen(efileεname_σ, "r+")

Cheight_havg_tavg_UN =file_scale_UN["Cheight_havg_tavg"]
Cheight_havg_tavg_σ =file_scale_σ["Cheight_havg_tavg"]
eps_endAvg_UN = file_ε_UN["eps_endAvg"]
eps_endAvg_U2 = file_ε_U["eps_endAvg"]
eps_endAvg_σ = file_ε_σ["eps_endAvg"]

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

f = Figure(resolution = (1200, 1550), fontsize=40)
ga = f[1, 1] = GridLayout()

# first plot should be dissipation
ax2 = Axis(ga[2, 1],  ylabel = rich("h", subscript("w"), superscript("2"),"N", subscript("0"), superscript("3"), " [m²s", superscript("-3"),"]"),  yscale = log10,xscale = log10,
xlabel = "ϵ̄ [m²s⁻³]", xlabelsize=45, ylabelsize=45)
limits!(ax2, 1e-7, 10^-(4.25), 10^(-5.5), 10^(-2.75))

scatter!(ax2, eps_endAvg_σ[idx_gammachange_σ], δN_σ[idx_gammachange_σ], markersize = 30, marker=:star4, 
color =:darkgreen, strokewidth = 1, strokecolor = :black)
scatter!(ax2, eps_endAvg_σ[idx_subcritical_σ], δN_σ[idx_subcritical_σ], markersize = 30, 
marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
scatter!(ax2, eps_endAvg_UN[idx_VaryN], δN_UN[idx_VaryN], markersize = 30, marker = :circle, 
color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, eps_endAvg_UN[idx_VaryU], δN_UN[idx_VaryU], markersize = 30, marker=:dtriangle, 
color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

# second plot should be tracer thickness
ax1 = Axis(ga[3, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("tr")," [m]"), 
xlabelsize=45, ylabelsize=45)
ax1.xticks = 0:40:160
ax1.yticks = 0:40:160
limits!(ax1, 0, 180, 0, 180) 

vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_havg_tavg_σ[idx_gammachange_σ], markersize = 30, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_havg_tavg_σ[idx_subcritical_σ], markersize = 30, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_havg_tavg_UN[idx_VaryN], markersize = 30, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_havg_tavg_UN[idx_VaryU], markersize = 30, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend( ga[1, 1],  [vnp, vump, vsp1, vsbp], ["Vary N₀", "Vary V₀", "Vary γ", "Subcritical γ"],
            tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
            margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
            labelsize = 40,
            halign = :center, valign = :top, orientation = :horizontal)
rowsize!(ga, 1, Auto(0.05))       
colgap!(ga, 40)

Label(ga[2, 1, TopLeft()], "a",
                fontsize = 30,
                font = :bold,
                padding = (-10, 5, 5, 10),
                halign = :right)
Label(ga[3, 1, TopLeft()], "b",
                fontsize = 30,
                font = :bold,
                padding = (0, 5, 5, 10),
                halign = :right)

display(f)
savename = apath * "Paper_Dissip_v_Thickness_v_Delta_rms_log"
save(savename * ".png", f, px_per_unit = 2)
                