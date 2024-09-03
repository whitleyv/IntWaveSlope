using Statistics
using Printf
using CurveFit
using JLD2
using CairoMakie
using GeometryBasics
using CurveFit

apath = "Analysis/PaperFigures/FinalPaperFigures/FinalFiguresUsed/"
dpath = "Data/"

# setnames
filesetnames =  "SetList_mp.jld2"

# Thorpe Scale Values
filescalename_UN = dpath * "DeltavThorpeDissip.jld2"
filescalename_σ = dpath * "DeltavThorpeDissip_VarySigma.jld2"

# Dissipation Values
fileεname_U = dpath * "DeltavAll_mtest.jld2"
fileεname_UN = dpath * "DeltavDissip_mp.jld2"
fileεname_σ = dpath * "DeltavDissip_mp_varyσ.jld2"

# getting the surrounding simulation params out
file_sn = jldopen(filesetnames, "r+")

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
N_σ = file_sn["N_varyσ"]

δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
δ_U2 = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.45, 0.3, 0.4, 0.5, 0.55]./ N_UN[1]

γ_σ = file_sn["γ_varyσ"]

δN_U2 = δ_U2.^2 .* (N_UN[1]).^ 3
δN_UN = δ_UN.^2 .* N_UN.^3
δN_σ = δ_σ.^2 .* N_σ.^3

file_Lt_UN = jldopen(filescalename_UN, "r+")
file_Lt_σ = jldopen(filescalename_σ, "r+")

file_ε_UN = jldopen(fileεname_UN, "r+")
file_ε_U = jldopen(fileεname_U, "r+")
file_ε_σ = jldopen(fileεname_σ, "r+")

skeys = keys(file_Lt_UN)

thorpe_rms_UN = file_Lt_UN["thorpe_rms"]
thorpe_rms_σ = file_Lt_σ["thorpe_rms"]
Lt_hmax_tavg_U2 = file_ε_U["thorpe_hmax_tavg"]

Lt_eps_full_usingrms_UN = file_Lt_UN["Lt_eps_full_usingrms"]
Lt_eps_full_usingrms_σ = file_Lt_σ["Lt_eps_full_usingrms"]
LtN_U2 = Lt_hmax_tavg_U2.^2 .* (N_UN[1]).^ 3 

eps_endAvg_UN = file_ε_UN["eps_endAvg"]
eps_endAvg_U2 = file_ε_U["eps_endAvg"]
eps_endAvg_σ = file_ε_σ["eps_endAvg"]

idx_subcritical_σ = findall(γ_σ .< 1)
idx_topochange_σ = findall(γ_σ .== 1.9)
idx_gammachange_σ = findall(γ_σ .!= 1.9)

idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10
idx_VaryU2 = findall(δ_U2 .> 120)

# thorpe scale plot rms vs max
f = Figure(resolution = (1050, 700), fontsize=26)
    set_theme!(fonts = (; regular = "TeX Gyre Heros Makie Regular", bold = "TeX Gyre Heros Makie Bold", bold_italic = "TeX Gyre Heros Makie Bold Italic", italic = "TeX Gyre Heros Makie Italic"))

    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("T"), " [m]"), xlabelsize=35, ylabelsize=35)
    # we want a log y axis with ticks at -8 through 0
    ax1.xticks = 0:40:160
    ax1.yticks =  0:20:60
    limits!(ax1, 0, 165, 0, 65)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], thorpe_rms_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], thorpe_rms_σ[idx_subcritical_σ], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_max = scatter!(ax1, δ_U2[idx_VaryU2], Lt_hmax_tavg_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
        color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, δ_UN[idx_VaryN], thorpe_rms_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, δ_UN[idx_VaryU], thorpe_rms_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend(ga[1, 1], [vnp_rms, vump_rms, vump_max, vsp1, vsbp], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ", "Vary γ", "Subcritical γ"],
                    tellheight = false, tellwidth = false, labelsize = 35,rowgap = 30, colgap = 30,
                    margin = (10, 10, 20, 5), framevisible = false, patchlabelgap = 7,
                    halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.1))       
display(f)                         
savename = apath * "Paper_Thorpe_v_Delta_rms"
save(savename * ".png", f, px_per_unit = 2) 
#=
f = Figure(resolution = (1050, 700), fontsize=26)
    #set_theme!(fonts = (; regular = "TeX Gyre Heros Makie Regular", bold = "TeX Gyre Heros Makie Bold", bold_italic = "TeX Gyre Heros Makie Bold Italic", italic = "TeX Gyre Heros Makie Italic"))

    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1],  ylabel = rich("h", subscript("w"), " [m]"), 
        xlabel = "Lₜ [m]", xlabelsize=35, ylabelsize=35, aspect = 1)
    # we want a log y axis with ticks at -8 through 0
    ax1.yticks = 0:40:160
    ax1.xticks =  0:20:60
    limits!(ax1, 0, 65, 0, 165)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp1 = scatter!(ax1, thorpe_rms_σ[idx_gammachange_σ], δ_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, thorpe_rms_σ[idx_subcritical_σ], δ_σ[idx_subcritical_σ], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_max = scatter!(ax1, Lt_hmax_tavg_U2[idx_VaryU2], δ_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
        color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, thorpe_rms_UN[idx_VaryN], δ_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, thorpe_rms_UN[idx_VaryU], δ_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend(ga[1, 2], [vnp_rms, vump_rms, vump_max, vsp1, vsbp], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ", "Vary γ", "Subcritical γ"],
                    tellheight = false, tellwidth = false, labelsize = 35, rowgap = 30, colgap = 30,
                    framevisible = false, patchlabelgap = 7, )
                   # halign = :center, valign = :top, )#margin = (10, 10, 60, 5),

    colsize!(ga,2, Auto(0.35))     
    colgap!(ga, 5)  
display(f)                         
savename = apath * "Paper_Thorpe_v_Delta_rms_flipped_"
save(savename * ".png", f, px_per_unit = 2) 


f = Figure(resolution = (700, 1050), fontsize=26)
    #set_theme!(fonts = (; regular = "TeX Gyre Heros Makie Regular", bold = "TeX Gyre Heros Makie Bold", bold_italic = "TeX Gyre Heros Makie Bold Italic", italic = "TeX Gyre Heros Makie Italic"))

    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[3, 1],  ylabel = rich("h", subscript("w"), " [m]"), 
        xlabel = "Lₜ [m]", xlabelsize=35, ylabelsize=35, )
    # we want a log y axis with ticks at -8 through 0
    ax1.yticks = 0:40:160
    ax1.xticks =  0:20:60
    limits!(ax1, 0, 65, 0, 165)

display(f)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp1 = scatter!(ax1, thorpe_rms_σ[idx_gammachange_σ], δ_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, thorpe_rms_σ[idx_subcritical_σ], δ_σ[idx_subcritical_σ], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_max = scatter!(ax1, Lt_hmax_tavg_U2[idx_VaryU2], δ_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
        color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, thorpe_rms_UN[idx_VaryN], δ_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, thorpe_rms_UN[idx_VaryU], δ_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend(ga[1, 1], [vnp_rms, vump_rms, vump_max], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ"],
                    labelsize = 35, rowgap = 30, colgap = 30,
                    framevisible = false, patchlabelgap = 7, 
                    halign = :center, valign = :top, orientation = :horizontal)#, margin = (0, 100, 5, 5))
    #rowsize!(ga,1, Auto(0.1))       
    Legend(ga[2, 1], [vsp1, vsbp], ["Vary γ", "Subcritical γ"],
                    labelsize = 35, rowgap = 30, colgap = 30,
                    framevisible = false, patchlabelgap = 7, 
                    halign = :center, valign = :top, orientation = :horizontal)#, margin = (0, 100, 5, 5))
   rowsize!(ga,3, Auto(0.9))       

    #colgap!(ga, 5)  
display(f)                         
savename = apath * "Paper_Thorpe_v_Delta_rms_flipped2_"
save(savename * ".png", f, px_per_unit = 2) 
=#

# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = rich("L", subscript("t"), superscript("2"),"N", subscript("0"), superscript("3"), " [m²s", superscript("-3"), "]"), yscale = log10, xscale = log10,
    xlabel = "ϵ̄ [m²s⁻³]", xlabelsize=35, ylabelsize=35)
    #ax1.yticks = ( 10 .^(-5.5:.5:-3), ["0", "5×10⁻⁴", "1×10⁻³"])
    #ax1.xticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
    limits!(ax1, 1e-7, 10^-(4.25), 1e-6, 10^(-2.75))
    
    ax2 = Axis(ga[2, 2],  ylabel = rich("h", subscript("w"), superscript("2"),"N", subscript("0"), superscript("3"), " [m²s", superscript("-3"), "]"),  yscale = log10,xscale = log10,
    xlabel = "ϵ̄ [m²s⁻³]", xlabelsize=35, ylabelsize=35)
    #ax2.yticks = (1e:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
    #ax2.xticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
    limits!(ax2, 1e-7, 10^-(4.25), 1e-6, 10^(-2.75))

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp1 = scatter!(ax1, eps_endAvg_σ[idx_gammachange_σ], Lt_eps_full_usingrms_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, eps_endAvg_σ[idx_subcritical_σ], Lt_eps_full_usingrms_σ[idx_subcritical_σ], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, eps_endAvg_U2[idx_VaryU2], LtN_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryN], Lt_eps_full_usingrms_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryU], Lt_eps_full_usingrms_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
    
    scatter!(ax2, eps_endAvg_σ[idx_gammachange_σ], δN_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_σ[idx_subcritical_σ], δN_σ[idx_subcritical_σ], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, eps_endAvg_U2[idx_VaryU2], δN_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, eps_endAvg_UN[idx_VaryN], δN_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_UN[idx_VaryU], δN_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend(ga[1, 1:2], [vnp_rms, vump_rms, vump_mean, vsp1, vsbp], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ", "Vary γ", "Subcritical γ"],
                    tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
                    margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
                    labelsize = 35,
                    halign = :center, valign = :top, orientation = :horizontal)
    rowsize!(ga,1, Auto(0.05))       
    colgap!(ga, 30)
    Label(ga[2, 1, TopLeft()], "a",
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :left)
    Label(ga[2, 2, TopLeft()], "b",
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :left)
display(f)
savename = apath * "Paper_Dissip_v_Thorpe_v_Delta_rms_log_fixedbracket"
save(savename * ".png", f, px_per_unit = 2)

#=
# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = "Lₜ² N₀³ [m²s⁻³]",
    xlabel = "ϵ̄ [m²s⁻³]", xlabelsize=35, ylabelsize=35)
    ax1.yticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
    ax1.xticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
    limits!(ax1, -3e-7, 3.3e-5, -2e-5, 1.2*1e-3)
    
    ax2 = Axis(ga[2, 2],  ylabel = "δ² N₀³ [m²s⁻³]",  
    xlabel = "ϵ̄ [m²s⁻³]", xlabelsize=35, ylabelsize=35)
    ax2.yticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
    ax2.xticks = (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"])
    limits!(ax2, -3e-7, 3.3e-5, -2e-5, 1.2*1e-3)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp1 = scatter!(ax1, eps_endAvg_σ[idx_gammachange_σ], Lt_eps_full_usingrms_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, eps_endAvg_σ[idx_subcritical_σ], Lt_eps_full_usingrms_σ[idx_subcritical_σ], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, eps_endAvg_U2[idx_VaryU2], LtN_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryN], Lt_eps_full_usingrms_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryU], Lt_eps_full_usingrms_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
    
    scatter!(ax2, eps_endAvg_σ[idx_gammachange_σ], δN_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_σ[idx_subcritical_σ], δN_σ[idx_subcritical_σ], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, eps_endAvg_U2[idx_VaryU2], δN_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, eps_endAvg_UN[idx_VaryN], δN_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_UN[idx_VaryU], δN_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend(ga[1, 1:2], [vnp_rms, vump_rms, vump_mean, vsp1, vsbp], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ", "Vary γ", "Subcritical γ"],
                    tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
                    margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
                    labelsize = 35,
                    halign = :center, valign = :top, orientation = :horizontal)
    rowsize!(ga,1, Auto(0.05))       
    colgap!(ga, 30)
    Label(ga[2, 1, TopLeft()], "a",
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :left)
    Label(ga[2, 2, TopLeft()], "b",
                    fontsize = 30,
                    font = :bold,
                    padding = (5, 5, 5, 5),
                    halign = :left)
display(f)
savename = apath * "Paper_Dissip_v_Thorpe_v_Delta_rms"
save(savename * ".png", f, px_per_unit = 2)

eps_endAvg_all = vcat(vcat(eps_endAvg_σ[idx_gammachange_σ],eps_endAvg_UN ),  eps_endAvg_U2[idx_VaryU2])
δN_all = vcat(vcat(δN_σ[idx_gammachange_σ], δN_UN),  δN_U2[idx_VaryU2])

LtN_all = vcat(vcat(Lt_eps_full_usingrms_σ[idx_gammachange_σ], Lt_eps_full_usingrms_UN),  LtN_U2[idx_VaryU2])

(b, m_δ) = linear_fit(δN_all, eps_endAvg_all)
(b, m_Lt) = linear_fit(LtN_all, eps_endAvg_all)

=#