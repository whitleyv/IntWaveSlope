using Statistics
using Printf
using CurveFit
using JLD2

#ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename = dpath * "DeltavAllScale.jld2"
filescalename_oth = dpath * "DeltavAllScale_oth.jld2"
filesetnames =  "SetnamesList.jld2"

scale_file_sn = jldopen(filesetnames, "r+")
δs = scale_file_sn["δs"][1:22]

scale_file = jldopen(filescalename, "r+")
scale_file_oth = jldopen(filescalename_oth, "r+")

setnames = vcat(scale_file["setnames"][1:22],scale_file_oth["setnames_oth"])
δ = vcat(δs,scale_file_oth["δ_oth"])
thorpe_hmax_tavg = vcat(scale_file["thorpe_hmax_tavg"][1:22],scale_file_oth["thorpe_hmax_tavg_oth"])

@info "Find Bestfit Lines"

(b_Ltmδ,m_Ltmδ) = linear_fit(δ, thorpe_hmax_tavg)

@info "Plotting!"
#VaryN = 12:25
# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryN = 1:11
VaryU = 12:22
# ones that are subcritical out of sigma varying set
Varysg_sub = 25:26
# all the ones where sigma is varied in any way
Varyoth = 23:30     

using CairoMakie
using GeometryBasics

f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout()

ax = Axis(ga[2, 1],  xlabel = "δ [m]", ylabel = "Lₜ [m]")
# we want a log y axis with ticks at -8 through 0
ax.xticks = 0:40:160
ax.yticks = 20:20:60
limits!(ax, 0, 165, 0, 70)

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

vsp = scatter!(ax, δ[Varyoth], thorpe_hmax_tavg[Varyoth], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax, δ[Varysg_sub], thorpe_hmax_tavg[Varysg_sub], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax, δ[VaryN], thorpe_hmax_tavg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax, δ[VaryU], thorpe_hmax_tavg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
Legend(ga[1, 1], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
                tellheight = false, tellwidth = false, 
                margin = (10, 10, 20, 5), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :top, orientation = :horizontal)
rowsize!(ga,1, Auto(0.05))                                
savename = apath * "Paper_Thorpe_v_Delta_All"
save(savename * ".png", f, px_per_unit = 2)

