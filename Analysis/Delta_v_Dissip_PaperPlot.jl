using Measures
using Statistics
using Printf
using CurveFit
using JLD2

#ENV["GKSwstype"] = "nul" # if on remote HPC

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename_ε = dpath * "DeltavDissip.jld2"
filescalename_Lt = dpath * "DeltavAllScale.jld2"
filesetnames =  "SetnamesList.jld2"

scale_file_sn = jldopen(filesetnames, "r+")
δ = scale_file_sn["δs"][1:25]

scale_file_ε = jldopen(filescalename_ε, "r+")
scale_file_Lt = jldopen(filescalename_Lt, "r+")

setnames = scale_file_ε["setnames"]

include("parameters.jl")

Ñ = zeros(length(setnames))
for (m,setname) in enumerate(setnames)
    pm = getproperty(SimParams(), Symbol(setname))
    Ñ[m] = pm.Ñ
end

eps_beginAvg = scale_file_ε["eps_beginAvg"]
eps_endAvg = scale_file_ε["eps_endAvg"]
eps_beginMax = scale_file_ε["eps_beginMaxAvg"]
eps_endMax = scale_file_ε["eps_endMaxAvg"]
thorpe_hmax_tavg = scale_file_Lt["thorpe_hmax_tavg"]

δN = δ.^2 .* Ñ.^3
LtN = thorpe_hmax_tavg.^2 .* Ñ.^3

@info "Find Bestfit Lines for full data"
(b_εδ,m_εδ) = linear_fit(δN, eps_endAvg)
(b_εLt,m_εLt) = linear_fit(LtN, eps_endAvg)
(b_Mεδ,m_Mεδ) = linear_fit(δN, eps_endMax)
(b_MεLt,m_MεLt) = linear_fit(LtN, eps_endMax)

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryU = 1:11
VaryN = 12:20

VaryU = 1:11
VaryN = 12:25

# removing last N and first U val from data
VaryU2 = 1:10
VaryN2 = 13:22

# indices now avai;able for two data sets after cutting out above points
VaryU_cut = 1:10
VaryN_cut = 11:20

@info "Find Bestfit Lines for data cut"

eps_endAvg_cut = vcat(eps_endAvg[VaryU2], eps_endAvg[VaryN2])
eps_endMax_cut = vcat(eps_endMax[VaryU2], eps_endMax[VaryN2])
δN_cut = vcat(δN[VaryU2], δN[VaryN2])
LtN_cut = vcat(LtN[VaryU2], LtN[VaryN2])

(b_εδ_cut,m_εδ_cut) = linear_fit(δN_cut, eps_endAvg_cut)
(b_εLt_cut,m_εLt_cut) = linear_fit(LtN_cut, eps_endAvg_cut)
(b_Mεδ_cut,m_Mεδ_cut) = linear_fit(δN_cut, eps_endMax_cut)
(b_MεLt_cut,m_MεLt_cut) = linear_fit(LtN_cut, eps_endMax_cut)

@info "Plotting..."

""" Plotting Full Plots for Avg Dissipation"""

xlimLt =maximum(LtN)
xlimd = maximum(δN)
xmind = minimum(δN)
xminLt = minimum(LtN)
ylime = maximum(eps_endAvg)
ymine = minimum(eps_endAvg)

using CairoMakie
using GeometryBasics

f = Figure(resolution = (1400, 800), fontsize=26)
ga = f[1, 1] = GridLayout()

ax1 = Axis(ga[1, 1],  xlabel = "Lₜ²N³ [m²s⁻³]", ylabel = "⟨ε⟩ [m²s⁻³]")
ax1.xticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
ax1.yticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
limits!(ax1, -3e-5, 1.2*1e-3, 0, 3.8e-5)

ax2 = Axis(ga[1, 2],  xlabel = "δ²N³ [m²s⁻³]")
ax2.xticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
ax2.yticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
limits!(ax2, -3e-5, 1.2*1e-3, 0, 3.8e-5)

# inset 1
xL1 = 0
xR1 = 6e-5
yT1 = 6e-6
yB1 = 4e-6

poly!(ax1, Point2f[(6e-4, 2.9e-6), (1.1e-3, 2.9e-6), (1.1e-3, 1.55e-5), (6e-4, 1.55e-5)], color = :white)

poly!(ax1, Point2f[(xL1, yB1), (xR1, yB1), (xR1, yT1), (xL1, yT1)], 
        color = :transparent, strokewidth = 3, strokecolor = (:darkgreen, 0.6))

# inset 2
xL2 = 7e-5
xR2 = 2.5e-4
yT2 = 6e-6
yB2 = 4e-6

poly!(ax2, Point2f[(7e-5, 2.3e-5), (5.7e-4, 2.3e-5), (5.7e-4, 3.51e-5), (7e-5, 3.51e-5)], color = :white)

poly!(ax2, Point2f[(xL2, yB2), (xR2, yB2), (xR2, yT2), (xL2, yT2)], 
        color = :transparent, strokewidth = 3, strokecolor = (:darkgreen, 0.6))

dx1 = round(((xR1 - xL1)/6) * 1e6)*1e-6
dy1 = round(((yT1 - yB1)/6) * 1e7)*1e-7

xtL1 = xL1 + dx1
xtR1 = xR1 - dx1
ytT1 = yB1 + dy1
ytB1 = yT1 - dy1

# add box inside first plot
inset_ax1 = Axis(f, bbox = BBox(448, 689, 132, 352), backgroundcolor = :white,
    spinewidth = 4, leftspinecolor = :darkgreen, rightspinecolor = :darkgreen, 
    topspinecolor = :darkgreen, bottomspinecolor = :darkgreen)

inset_ax1.xticks = [xtL1, xtR1] #(0:2e-5:4e-5, ["0", "2×10⁻⁵", "4×10⁻⁵"])
inset_ax1.yticks = [ytB1, ytT1] #(4e-6:2e-6:6e-6, ["4×10⁻⁶", "6×10⁻⁶"])
limits!(inset_ax1, xL1, xR1, yB1, yT1)

dx2 = round(((xR2 - xL2)/6) * 1e5)*1e-5
dy2 = round(((yT2 - yB2)/6) * 1e7)*1e-7

xtL2 = xL2 + dx2
xtR2 = xR2 - dx2
ytT2 = yB2 + dy2
ytB2 = yT2 - dy2

inset_ax2 = Axis(f, bbox = BBox(844, 1085, 475, 695), backgroundcolor = :white,
    spinewidth = 4, leftspinecolor = :darkgreen, rightspinecolor = :darkgreen, 
    topspinecolor = :darkgreen, bottomspinecolor = :darkgreen,
    yaxisposition = :right)

inset_ax2.xticks =  ([xtL2, xtR2], ["1×10⁻⁴", "2.2×10⁻⁴"])#(1e-4:5e-5:2e-4, ["1×10⁻⁴", "1.5×10⁻⁴", "2×10⁻⁴"])
inset_ax2.yticks =  [ytB2, ytT2] #(4e-6:2e-6:6e-6, ["4×10⁻⁶", "6×10⁻⁶"])
limits!(inset_ax2, xL2, xR2, yB2, yT2)

hideydecorations!(ax2, grid = false)

translate!(inset_ax1.scene, 0, 0, 50) # bring box to front
translate!(inset_ax2.scene, 0, 0, 50) # bring box to front

vnp = scatter!(ax1, LtN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
vup = scatter!(ax1, LtN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
#rangebars!(ax1, [5e-5], [0], [6e-6], color = :darkgreen, whiskerwidth = 10, linewidth = 4, direction = :y)

scatter!(ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

Legend( ga[1, 1, Top()], [vnp, vup], ["Vary N₀", "Vary V₀"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 15), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :bottom, orientation = :horizontal)

Label(ga[1, 1, TopLeft()], "a",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)
Label(ga[1, 2, TopLeft()], "b",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)

scatter!(inset_ax1, LtN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(inset_ax1, LtN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

scatter!(inset_ax2, δN[VaryN], eps_endAvg[VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(inset_ax2, δN[VaryU], eps_endAvg[VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

savename = apath * "Paper_Dissip_v_LtDelta_All"
save(savename * ".png", f, px_per_unit = 2)
