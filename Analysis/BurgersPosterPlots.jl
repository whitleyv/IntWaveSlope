using Measures
using Statistics
using Printf
using CurveFit
using JLD2

#ENV["GKSwstype"] = "nul" # if on remote HPC

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename_ε = dpath * "DeltavDissip_mp.jld2"
filescalename = dpath * "DeltavAllScale_mp.jld2"
filescalename_ε2 = dpath * "DeltavDissip.jld2"
filescalename2 = dpath * "DeltavAllScale.jld2"
#filescalename_oth = dpath * "DeltavAllScale_oth.jld2"
filesetnames =  "SetList_mp.jld2"
filesetnames2 =  "SetnamesList.jld2"
filescalename3 = dpath * "DeltavAll_mtest.jld2"

scale_file_sn = jldopen(filesetnames, "r+")
scale_file_sn2 = jldopen(filesetnames2, "r+")

δ = scale_file_sn2["δs"][1:25]

δs = scale_file_sn["δs"][1:22]
Ns = scale_file_sn["Ns"][1:22]

δ4 = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.45, 0.3, 0.4, 0.5, 0.55]./ Ns[1]
δN4 = δ4.^2 .* (Ns[1])^3  

scale_file = jldopen(filescalename, "r+")
scale_file2 = jldopen(filescalename2, "r+")
scale_file3 = jldopen(filescalename3, "r+")
setname3 = vcat(scale_file3["setnames"])
sns3 = vcat(scale_file3["sns"])

scale_file_oth = jldopen(filescalename_oth, "r+")

setnames = vcat(scale_file["setnames"][1:22],scale_file_oth["setnames_oth"])
δ2 = vcat(δs,scale_file_oth["δ_oth"])
thorpe_hmax_tavg = vcat(scale_file["thorpe_hmax_tavg"][1:22],scale_file_oth["thorpe_hmax_tavg_oth"])

Cheight_havg_tavg =vcat(scale_file["Cheight_havg_tavg"][1:22],scale_file_oth["Cheight_havg_tavg_oth"])
Nheight_havg_tavg =vcat(scale_file["Nheight_havg_tavg"][1:22],scale_file_oth["Nheight_havg_tavg"])

thorpe_hmax_tavg2 = scale_file2["thorpe_hmax_tavg"]
Cheight_havg_tavg2 =scale_file2["Cheight_havg_tavg"]
Nheight_havg_tavg2 =scale_file2["Nheight_havg_tavg"]

thorpe_hmax_tavg3 = scale_file3["thorpe_hmax_tavg"]
Cheight_havg_tavg3 =scale_file3["Cheight_havg_tavg"]
eps_endAvg3 = scale_file3["eps_endAvg"]

scale_file_ε = jldopen(filescalename_ε, "r+")
scale_file_ε2 = jldopen(filescalename_ε2, "r+")

setnames = scale_file_ε2["setnames"]

include("../parameters.jl")

Ñ2 = zeros(3)
δ3 = zeros(3)
for (m,setname) in enumerate(setnames[end-2:end])
    pm = getproperty(SimParams(), Symbol(setname))
    Ñ2[m] = pm.Ñ
    δ3[m] = pm. U₀/pm.Ñ
end
Ñ2 = vcat(Ns, Ñ2)
δ3 = vcat(δs, δ3)

Ñ = Ns
eps_endAvg = scale_file_ε["eps_endAvg"]
eps_endAvg2 = scale_file_ε2["eps_endAvg"]

δN = δs.^2 .* Ñ.^3
δN2 = δ3.^2 .* Ñ2.^3

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryU = 1:11
VaryN = 12:22 #25

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryN2 = 1:11
VaryU2 = 12:22
VaryN3 = 23:25
# ones that are subcritical out of sigma varying set
Varysg_sub = 25:26
# all the ones where sigma is varied in any way
Varyoth = 23:30    

@info "Plotting..."

""" Plotting Full Plots for Avg Dissipation"""

xlimd = maximum(δN)
xmind = minimum(δN)
ylime = maximum(eps_endAvg)
ymine = minimum(eps_endAvg)

using CairoMakie
using GeometryBasics

f1 = Figure(resolution = (900, 1000), fontsize=26)
    ga = f1[1, 1] = GridLayout()

    ax = Axis(ga[2, 1],  xlabel = "δ [m]", ylabel = "Lₜ [m]")

    ax.xticks = 0:40:160
    ax.yticks = 20:20:60
    limits!(ax, 0, 165, 0, 70)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp = scatter!(ax, δ2[Varyoth], thorpe_hmax_tavg[Varyoth], markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax, δ2[Varysg_sub], thorpe_hmax_tavg[Varysg_sub], markersize = 25, 
        marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vnp = scatter!(ax, δ2[VaryN], thorpe_hmax_tavg[VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vup = scatter!(ax, δ2[VaryU], thorpe_hmax_tavg[VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
                                
    ax2 = Axis(ga[3, 1],  xlabel = "δ²N³ [m²s⁻³]", ylabel = "ε [m²s⁻³]")
    ax2.xticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
    ax2.yticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
    limits!(ax2, -1e-5, 1.15*1e-3, -3e-7, 3.4e-5)

    scatter!(ax2, δN4[7], eps_endAvg3[7], markersize = 25, marker=:utriangle, 
    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δN4[10:end], eps_endAvg3[10:end], markersize = 25, marker=:utriangle, 
    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δN[1:8], eps_endAvg[1:8], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend(ga[1, 1, Top()], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
                tellheight = false, tellwidth = false, 
                margin = (10, 10, 10, 0), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :top, orientation = :horizontal)
    rowsize!(ga,1, Auto(0.05))    
    #=
    scatter!(inset_ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(inset_ax2, δN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
    =# 
    yspace = maximum(tight_yticklabel_spacing!, [ax, ax2])
                ax.yticklabelspace = yspace
                ax2.yticklabelspace = yspace
            
savename = apath * "NOPPPoster_Dissip_v_Lt_v_Delta_All"
save(savename * ".png", f1, px_per_unit = 2)

display(f1)

f1 = Figure(resolution = (1500, 800), fontsize=26)
    ga = f1[1, 1] = GridLayout()

    ax = Axis(ga[2, 1],  xlabel = "δ [m]", ylabel = "Lₜ [m]")

    ax.xticks = 0:40:160
    ax.yticks = 20:20:60
    limits!(ax, 0, 165, 0, 70)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp = scatter!(ax, δ2[Varyoth], thorpe_hmax_tavg[Varyoth], markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax, δ2[Varysg_sub], thorpe_hmax_tavg[Varysg_sub], markersize = 25, 
        marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vnp = scatter!(ax, δ2[VaryN], thorpe_hmax_tavg[VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vup = scatter!(ax, δ2[VaryU], thorpe_hmax_tavg[VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
                                
    ax2 = Axis(ga[2, 2],  xlabel = "δ²N³ [m²s⁻³]", ylabel = "ε [m²s⁻³]")
    ax2.xticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
    ax2.yticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
    limits!(ax2, -1e-5, 1.15*1e-3, -3e-7, 3.3e-5)

scatter!(ax2, δN4[7], eps_endAvg3[7], markersize = 25, marker=:utriangle, 
    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δN4[10:end], eps_endAvg3[10:end], markersize = 25, marker=:utriangle, 
    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

scatter!(ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δN[1:8], eps_endAvg[1:8], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend(ga[1, 1:2, Top()], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
            tellheight = false, tellwidth = false, 
            margin = (10, 10, 10, 0), framevisible = false, patchlabelgap = 7,
            halign = :center, valign = :top, orientation = :horizontal)
rowsize!(ga,1, Auto(0.05))    
#=
scatter!(inset_ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(inset_ax2, δN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
=# 
yspace = maximum(tight_yticklabel_spacing!, [ax, ax2])
            ax.yticklabelspace = yspace
            ax2.yticklabelspace = yspace
            
savename = apath * "NOPPPoster_Dissip_v_Lt_v_Delta_All2"
save(savename * ".png", f1, px_per_unit = 2)

#display(f1)

""" Plotting Full Plot for Avg Intrusion Thickness"""

f2 = Figure(resolution = (1000, 700), fontsize=26)
ga = f2[1, 1] = GridLayout()

ax1 = Axis(ga[2, 1],  xlabel = "δ [m]", ylabel = "Tracer Thickness, Lₜᵣ [m]")
ax1.xticks = 0:40:160
ax1.yticks = 0:40:160
limits!(ax1, 0, 180, 0, 180)

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

vsp = scatter!(ax1, δ2[Varyoth], Cheight_havg_tavg[Varyoth], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax1, δ2[Varysg_sub], Cheight_havg_tavg[Varysg_sub], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax1, δ2[VaryN], Cheight_havg_tavg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax1, δ2[VaryU], Cheight_havg_tavg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend(ga[1, 1, Top()], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 0), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :top, orientation = :horizontal)

rowsize!(ga,1, Auto(0.05))    

#display(f2)

savename = apath * "NOPPPoster_Intrusion_v_Delta"
save(savename * ".png", f2, px_per_unit = 2)
 

f2 = Figure(resolution = (1000, 1000), fontsize=26)
ga = f2[1, 1] = GridLayout()

ax1 = Axis(ga[2, 1], ylabel = "Tracer Thickness, Lₜᵣ [m]")
ax2 = Axis(ga[3, 1],  xlabel = "δ [m]", ylabel = "Tracer Thickness, Lₜᵣ [m]")

hidexdecorations!(ax1)
ax2.xticks = 0:40:160
ax1.yticks = 0:40:160
ax2.yticks = 0:40:160

limits!(ax1, 0, 180, 0, 180)
limits!(ax2, 0, 180, 0, 120)

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

vsp = scatter!(ax1, δ2[Varyoth], Cheight_havg_tavg[Varyoth], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax1, δ2[Varysg_sub], Cheight_havg_tavg[Varysg_sub], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax1, δ2[VaryN], Cheight_havg_tavg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax1, δ2[VaryU], Cheight_havg_tavg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

vsp = scatter!(ax2, δ2[Varyoth], Nheight_havg_tavg[Varyoth], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax2, δ2[Varysg_sub], Nheight_havg_tavg[Varysg_sub], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax2, δ2[VaryN], Nheight_havg_tavg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax2, δ2[VaryU], Nheight_havg_tavg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend(ga[1, 1, Top()], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 0), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :top, orientation = :horizontal)

rowsize!(ga,1, Auto(0.05))    
rowsize!(ga, 3, Auto(0.65))
#display(f2)

savename = apath * "NOPPPoster_Intrusion_v_Strat_v_Delta"
save(savename * ".png", f2, px_per_unit = 2)
 