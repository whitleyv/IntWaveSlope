using Statistics
using Printf
using CurveFit
using JLD2
using CairoMakie
using GeometryBasics

apath = "Analysis/PaperFigures/"
dpath = "Data/"

filescalename_UN = dpath * "DeltavAllScale_mp.jld2"
filescalename_σ = dpath * "DeltavAllScale_mp_VarySigma.jld2"
filesetnames =  "SetList_mp.jld2"

# getting the surrounding simulation params out
file_sn = jldopen(filesetnames, "r+")

δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
N_σ = file_sn["N_varyσ"]
γ_σ = file_sn["γ_varyσ"]

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
file_scale_UN = jldopen(filescalename_UN, "r+")
file_scale_σ = jldopen(filescalename_σ, "r+")

Cheight_havg_tavg_UN =file_scale_UN["Cheight_havg_tavg"]
Nheight_havg_tavg_UN =file_scale_UN["Nheight_havg_tavg"]
Cheight_havg_tavg_σ =file_scale_σ["Cheight_havg_tavg"]
Nheight_havg_tavg_σ =file_scale_σ["Nheight_havg_tavg"]

#=
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  xlabel = L"$δ$ [m]", ylabel = L"Tracer Thickness, L$_{tr}$ [m]", xlabelsize=35, ylabelsize=35)
    ax1.xticks = (0:40:160, [L"0", L"40", L"80", L"120", L"160"]) 
    ax1.yticks = (0:40:160, [L"0", L"40", L"80", L"120", L"160"]) 
    limits!(ax1, 0, 180, 0, 180) 

    ax2 = Axis(ga[2, 2],  xlabel = L"$δ$ [m]", ylabel = L"Stratification Anomaly Thickness, $L_{N^2}$ [m]", xlabelsize=35, ylabelsize=35)
    ax2.xticks = (0:40:160, [L"0", L"40", L"80", L"120", L"160"]) 
    ax2.yticks = (0:40:160, [L"0", L"40", L"80", L"120", L"160"]) 
    limits!(ax2, 0, 180, 0, 180)

    #hidexdecorations!(ax1, grid = false)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    #vsp2 = scatter!(ax1, δ_σ[idx_topochange_σ], Cheight_havg_tavg_σ[idx_topochange_σ], markersize = 25, marker=:star5, 
    #            color =:chartreuse3, strokewidth = 1, strokecolor = :black)
    vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_havg_tavg_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_havg_tavg_σ[idx_subcritical_σ], markersize = 25, 
        marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_havg_tavg_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_havg_tavg_UN[idx_VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    #scatter!(ax2, δ_σ[idx_topochange_σ], Nheight_havg_tavg_σ[idx_topochange_σ], markersize = 25, marker=:star5, 
    #            color =:chartreuse3, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_σ[idx_gammachange_σ], Nheight_havg_tavg_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_σ[idx_subcritical_σ], Nheight_havg_tavg_σ[idx_subcritical_σ], markersize = 25, 
        marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, δ_UN[idx_VaryN], Nheight_havg_tavg_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_UN[idx_VaryU], Nheight_havg_tavg_UN[idx_VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend( ga[1, 1:2],  [vnp, vump, vsp1, vsbp], [L"Vary N$_0$", L"Vary V$_0$", L"Vary $γ$", L"Subcritical $γ$"],
                tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
                margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
                labelsize = 35,
                halign = :center, valign = :top, orientation = :horizontal)
    rowsize!(ga,1, Auto(0.05))       
    colgap!(ga, 40)
                
    Label(ga[2, 1, TopLeft()], "a",
                    fontsize = 30,
                    font = :bold,
                    padding = (-10, 5, 5, 10),
                    halign = :right)
    Label(ga[2, 2, TopLeft()], "b",
                    fontsize = 30,
                    font = :bold,
                    padding = (0, 5, 5, 10),
                    halign = :right)
    display(f)
savename = apath * "Paper_Intrusion_v_Delta_ALL"
save(savename * ".png", f, px_per_unit = 2)
=#

f = Figure(resolution = (1500, 800), fontsize=26)
ga = f[1, 1] = GridLayout()

ax1 = Axis(ga[2, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("tr")," [m]"), 
xlabelsize=30, ylabelsize=30)
ax1.xticks = 0:40:160
ax1.yticks = 0:40:160
limits!(ax1, 0, 180, 0, 180) 

ax2 = Axis(ga[2, 2],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("N", superscript("2")), " [m]"), 
xlabelsize=30, ylabelsize=30)
ax2.xticks = 0:40:160
ax2.yticks = 0:40:160
limits!(ax2, 0, 180, 0, 180)

#hidexdecorations!(ax1, grid = false)

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

#vsp2 = scatter!(ax1, δ_σ[idx_topochange_σ], Cheight_havg_tavg_σ[idx_topochange_σ], markersize = 25, marker=:star5, 
#            color =:chartreuse3, strokewidth = 1, strokecolor = :black)
vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_havg_tavg_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_havg_tavg_σ[idx_subcritical_σ], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_havg_tavg_UN[idx_VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_havg_tavg_UN[idx_VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

#scatter!(ax2, δ_σ[idx_topochange_σ], Nheight_havg_tavg_σ[idx_topochange_σ], markersize = 25, marker=:star5, 
#            color =:chartreuse3, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δ_σ[idx_Nsubcritical_gammachange_σ], Nheight_havg_tavg_σ[idx_Nsubcritical_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δ_UN[idx_VaryN_n05], Nheight_havg_tavg_UN[idx_VaryN_n05], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δ_UN[idx_VaryU_n05], Nheight_havg_tavg_UN[idx_VaryU_n05], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black) 

Legend( ga[1, 1:2],  [vnp, vump, vsp1, vsbp], ["Vary N₀", "Vary V₀", "Vary γ", "Subcritical γ*"],
            tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
            margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
            labelsize = 35,
            halign = :center, valign = :top, orientation = :horizontal)
rowsize!(ga,1, Auto(0.05))       
colgap!(ga, 40)
            
Label(ga[2, 1, TopLeft()], "a",
                fontsize = 30,
                font = :bold,
                padding = (-10, 5, 5, 10),
                halign = :right)
Label(ga[2, 2, TopLeft()], "b",
                fontsize = 30,
                font = :bold,
                padding = (0, 5, 5, 10),
                halign = :right)
display(f)
savename = apath * "Paper_Intrusion_v_Delta_All"
save(savename * ".png", f, px_per_unit = 2)
