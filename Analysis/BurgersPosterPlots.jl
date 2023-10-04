using Measures
using Statistics
using Printf
using CurveFit
using JLD2

#ENV["GKSwstype"] = "nul" # if on remote HPC

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename_ε = dpath * "DeltavDissip.jld2"
filescalename = dpath * "DeltavAllScale.jld2"
filescalename_oth = dpath * "DeltavAllScale_oth.jld2"
filesetnames =  "SetnamesList.jld2"

scale_file_sn = jldopen(filesetnames, "r+")
δ = scale_file_sn["δs"][1:25]

δs = scale_file_sn["δs"][1:22]

scale_file = jldopen(filescalename, "r+")
scale_file_oth = jldopen(filescalename_oth, "r+")

setnames = vcat(scale_file["setnames"][1:22],scale_file_oth["setnames_oth"])
δ2 = vcat(δs,scale_file_oth["δ_oth"])
thorpe_hmax_tavg = vcat(scale_file["thorpe_hmax_tavg"][1:22],scale_file_oth["thorpe_hmax_tavg_oth"])

Cheight_havg_tavg =vcat(scale_file["Cheight_havg_tavg"][1:22],scale_file_oth["Cheight_havg_tavg_oth"])

scale_file_ε = jldopen(filescalename_ε, "r+")

setnames = scale_file_ε["setnames"]

include("../parameters.jl")

Ñ = zeros(length(setnames))
for (m,setname) in enumerate(setnames)
    pm = getproperty(SimParams(), Symbol(setname))
    Ñ[m] = pm.Ñ
end

eps_beginAvg = scale_file_ε["eps_beginAvg"]
eps_endAvg = scale_file_ε["eps_endAvg"]
eps_beginMax = scale_file_ε["eps_beginMaxAvg"]
eps_endMax = scale_file_ε["eps_endMaxAvg"]

δN = δ.^2 .* Ñ.^3

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryU = 1:11
VaryN = 12:25

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryN2 = 1:11
VaryU2 = 12:22
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

f = Figure(resolution = (900, 1000), fontsize=26)
ga = f[1, 1] = GridLayout()

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
vnp = scatter!(ax, δ2[VaryN2], thorpe_hmax_tavg[VaryN2], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax, δ2[VaryU2], thorpe_hmax_tavg[VaryU2], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
                            
ax2 = Axis(ga[3, 1],  xlabel = "δ²N³ [m²s⁻³]", ylabel = "ε [m²s⁻³]")
ax2.xticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
ax2.yticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
limits!(ax2, -3e-5, 1.15*1e-3, 0, 3.8e-5)

# inset 2
xL2 = 7e-5
xR2 = 2.5e-4
yT2 = 6e-6
yB2 = 4e-6

# where the plot is plotted
poly!(ax2, Point2f[(1e-5, 2.3e-5), (4.7e-4, 2.3e-5), (4.7e-4, 3.51e-5), (1e-5, 3.51e-5)], color = :white)
# where teh plot material comes from in the main plot
poly!(ax2, Point2f[(xL2, yB2), (xR2, yB2), (xR2, yT2), (xL2, yT2)], 
        color = :transparent, strokewidth = 3, strokecolor = (:darkgreen, 0.6))

# add box inside first plot
dx2 = round(((xR2 - xL2)/6) * 1e5)*1e-5
dy2 = round(((yT2 - yB2)/6) * 1e7)*1e-7

xtL2 = xL2 + dx2
xtR2 = xR2 - dx2
ytT2 = yB2 + dy2
ytB2 = yT2 - dy2

#860, 1101, 465, 685
#          844, 1085, 475, 695  *2 = 1688,     2170,       950,    1390      
#                                  L from L, R from L, B from B, T from B                        
inset_ax2 = Axis(f, bbox = BBox(170, 430, 313, 434), backgroundcolor = :white,
    spinewidth = 4, leftspinecolor = :darkgreen, rightspinecolor = :darkgreen, 
    topspinecolor = :darkgreen, bottomspinecolor = :darkgreen,
    yaxisposition = :right)

inset_ax2.xticks =  ([xtL2, xtR2], ["1×10⁻⁴", "2.2×10⁻⁴"])
inset_ax2.yticks =  [ytB2, ytT2] 
limits!(inset_ax2, xL2, xR2, yB2, yT2)

translate!(inset_ax2.scene, 0, 0, 50) # bring box to front

scatter!(ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend(ga[1, 1, Top()], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
            tellheight = false, tellwidth = false, 
            margin = (10, 10, 10, 0), framevisible = false, patchlabelgap = 7,
            halign = :center, valign = :top, orientation = :horizontal)
rowsize!(ga,1, Auto(0.05))    

scatter!(inset_ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(inset_ax2, δN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

yspace = maximum(tight_yticklabel_spacing!, [ax, ax2])
            ax.yticklabelspace = yspace
            ax2.yticklabelspace = yspace
            
savename = apath * "Poster_Dissip_v_Lt_v_Delta_All2"
save(savename * ".png", f, px_per_unit = 2)
display(f)

""" Plotting Full Plot for Avg Intrusion Thickness"""

f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout()

ax1 = Axis(ga[2, 1],  xlabel = "δ [m]", ylabel = "Tracer Thickness, Lₜᵣ [m]")
ax1.xticks = 0:40:160
ax1.yticks = 0:40:160
limits!(ax1, 0, 165, 0, 165)

p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

vsp = scatter!(ax1, δ2[Varyoth], Cheight_havg_tavg[Varyoth], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
vsbp = scatter!(ax1, δ2[Varysg_sub], Cheight_havg_tavg[Varysg_sub], markersize = 25, 
    marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
vnp = scatter!(ax1, δ2[VaryN2], Cheight_havg_tavg[VaryN2], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax1, δ2[VaryU2], Cheight_havg_tavg[VaryU2], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend(ga[1, 1, Top()], [vnp, vup, vsp, vsbp], ["Vary N₀", "Vary V₀", "Vary σ", "Subcritical γ"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 0), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :top, orientation = :horizontal)

rowsize!(ga,1, Auto(0.05))    

display(f)

savename = apath * "Poster_Intrusion_v_Delta"
save(savename * ".png", f, px_per_unit = 2)
 