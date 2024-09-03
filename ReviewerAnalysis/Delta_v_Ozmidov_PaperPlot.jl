using Statistics
using Printf
using CurveFit
using JLD2
using CairoMakie
using GeometryBasics
using CurveFit

apath = "ReviewerAnalysis/ReviewerPlots/"
dpath = "Data/"
dpath_reviewer = "ReviewerAnalysis/ReviewerData/"

filesetnames =  "SetList_mp.jld2"

# Thorpe Scale Values
filesLtname_UN = dpath * "DeltavThorpeDissip.jld2"
fileLtname_σ = dpath * "DeltavThorpeDissip_VarySigma.jld2"

# Dissipation Values
fileεname_U = dpath * "DeltavAll_mtest.jld2"
fileεname_UN = dpath * "DeltavDissip_mp.jld2"
fileεname_σ = dpath * "DeltavDissip_mp_varyσ.jld2"

file_sn = jldopen(filesetnames, "r+")
file_Lt_UN = jldopen(filesLtname_UN, "r+")
file_Lt_σ = jldopen(fileLtname_σ, "r+")

file_ε_UN = jldopen(fileεname_UN, "r+")
file_ε_U = jldopen(fileεname_U, "r+")
file_ε_σ = jldopen(fileεname_σ, "r+")

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
N_σ = file_sn["N_varyσ"]
δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
δ_U2 = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.45, 0.3, 0.4, 0.5, 0.55]./ N_UN[1]
γ_σ = file_sn["γ_varyσ"]

thorpe_rms_UN = file_Lt_UN["thorpe_rms"]
thorpe_rms_σ = file_Lt_σ["thorpe_rms"]
thorpe_rms_U2 = file_ε_U["thorpe_hmax_tavg"]

eps_endAvg_UN = file_ε_UN["eps_endAvg"]
eps_endAvg_U2 = file_ε_U["eps_endAvg"]
eps_endAvg_σ = file_ε_σ["eps_endAvg"]

idx_subcritical_σ = findall(γ_σ .< 1)
idx_gammachange_σ = findall(γ_σ .!= 1.9)
idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10
idx_VaryU2 = findall(δ_U2 .> 120)
N_all = vcat(N_UN, N_σ[idx_gammachange_σ])

ϵ_over_N3_σ = eps_endAvg_σ ./ N_σ.^3
ϵ_over_N3_UN = eps_endAvg_UN ./ (N_UN .^ 3)
ϵ_over_N3_U2 = eps_endAvg_U2 ./ (N_UN[1]).^ 3

thorpe2_σ = thorpe_rms_σ .^ 2
thorpe2_UN = thorpe_rms_UN .^ 2
thorpe2_U2 = thorpe_rms_U2 .^ 2

# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = rich("L", subscript("T"), superscript("2")," [m²]"), yscale = log10, xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax1, 1, 1100, 10, 35000)

    ax2 = Axis(ga[2, 2],  ylabel = rich("h", subscript("w"), superscript("2")," [m²]"), yscale = log10, xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax2, 1, 1100, 10, 35000)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))
    vsp1 = scatter!(ax1, ϵ_over_N3_σ[idx_gammachange_σ], thorpe2_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, ϵ_over_N3_σ[idx_subcritical_σ], thorpe2_σ[idx_subcritical_σ], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, ϵ_over_N3_U2[idx_VaryU2], thorpe2_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, ϵ_over_N3_UN[idx_VaryN], thorpe2_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, ϵ_over_N3_UN[idx_VaryU], thorpe2_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax2, ϵ_over_N3_σ[idx_gammachange_σ], δ_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, ϵ_over_N3_σ[idx_subcritical_σ], δ_σ[idx_subcritical_σ].^2, markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, ϵ_over_N3_U2[idx_VaryU2], δ_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, ϵ_over_N3_UN[idx_VaryN], δ_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, ϵ_over_N3_UN[idx_VaryU], δ_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
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
savename = apath * "Paper_Ozmidov2_v_Thorpe2_v_Delta2_log"
save(savename * ".png", f, px_per_unit = 2)


# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = rich("L", subscript("T")," [m]"),
    xlabel = rich("L", subscript("O")," [m]"), xlabelsize=35, ylabelsize=35)
    limits!(ax1, 0, 30, 0, 170)

    ax2 = Axis(ga[2, 2],  ylabel = rich("h", subscript("w")," [m]"),
    xlabel = rich("L", subscript("O")," [m]"), xlabelsize=35, ylabelsize=35)
    limits!(ax2, 0, 30, 0, 170)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))
    vsp1 = scatter!(ax1, sqrt.(ϵ_over_N3_σ[idx_gammachange_σ]), sqrt.(thorpe2_σ[idx_gammachange_σ]), markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, sqrt.(ϵ_over_N3_σ[idx_subcritical_σ]), sqrt.(thorpe2_σ[idx_subcritical_σ]), markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, sqrt.(ϵ_over_N3_U2[idx_VaryU2]), sqrt.(thorpe2_U2[idx_VaryU2]), markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, sqrt.(ϵ_over_N3_UN[idx_VaryN]), sqrt.(thorpe2_UN[idx_VaryN]), markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, sqrt.(ϵ_over_N3_UN[idx_VaryU]), sqrt.(thorpe2_UN[idx_VaryU]), markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax2, sqrt.(ϵ_over_N3_σ[idx_gammachange_σ]), δ_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, sqrt.(ϵ_over_N3_σ[idx_subcritical_σ]), δ_σ[idx_subcritical_σ], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, sqrt.(ϵ_over_N3_U2[idx_VaryU2]), δ_U2[idx_VaryU2], markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, sqrt.(ϵ_over_N3_UN[idx_VaryN]), δ_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, sqrt.(ϵ_over_N3_UN[idx_VaryU]), δ_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
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
savename = apath * "Paper_Ozmidov_v_Thorpe_v_Delta"
save(savename * ".png", f, px_per_unit = 2)

# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  xlabel = rich("L", subscript("T")," [m]"),
    ylabel = rich("L", subscript("O")," [m]"), xlabelsize=35, ylabelsize=35)
    limits!(ax1, 0, 50, 0, 30)

    ax2 = Axis(ga[2, 2],  xlabel = rich("h", subscript("w")," [m]"),
    ylabel = rich("L", subscript("O")," [m]"), xlabelsize=35, ylabelsize=35)
    limits!(ax2, 0, 170, 0, 30)
    hideydecorations!(ax2, grid = false)
    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    vsp1 = scatter!(ax1,  sqrt.(thorpe2_σ[idx_gammachange_σ]), sqrt.(ϵ_over_N3_σ[idx_gammachange_σ]), markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, sqrt.(thorpe2_σ[idx_subcritical_σ]), sqrt.(ϵ_over_N3_σ[idx_subcritical_σ]), markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, sqrt.(thorpe2_U2[idx_VaryU2]), sqrt.(ϵ_over_N3_U2[idx_VaryU2]), markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, sqrt.(thorpe2_UN[idx_VaryN]), sqrt.(ϵ_over_N3_UN[idx_VaryN]), markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, sqrt.(thorpe2_UN[idx_VaryU]), sqrt.(ϵ_over_N3_UN[idx_VaryU]), markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax2, δ_σ[idx_gammachange_σ], sqrt.(ϵ_over_N3_σ[idx_gammachange_σ]), markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_σ[idx_subcritical_σ], sqrt.(ϵ_over_N3_σ[idx_subcritical_σ]), markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, δ_U2[idx_VaryU2], sqrt.(ϵ_over_N3_U2[idx_VaryU2]), markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, δ_UN[idx_VaryN],  sqrt.(ϵ_over_N3_UN[idx_VaryN]), markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_UN[idx_VaryU], sqrt.(ϵ_over_N3_UN[idx_VaryU]),  markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend(ga[1, 1:2], [vnp_rms, vump_rms, vump_mean, vsp1, vsbp], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ", "Vary γ", "Subcritical γ"],
                    tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
                    margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
                    labelsize = 35,
                    halign = :center, valign = :top, orientation = :horizontal)
    rowsize!(ga,1, Auto(0.05))       
    colgap!(ga, 10)
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
savename = apath * "Paper_Ozmidov_v_Thorpe_v_Delta_switchedaxes"
save(savename * ".png", f, px_per_unit = 2)

# best fit lines 
LO_all = vcat(vcat(sqrt.(ϵ_over_N3_σ[idx_gammachange_σ]), sqrt.(ϵ_over_N3_UN)), sqrt.(ϵ_over_N3_U2[idx_VaryU2]))
Lt_all = vcat(vcat(thorpe_rms_σ[idx_gammachange_σ], thorpe_rms_UN), thorpe_rms_U2[idx_VaryU2])
δ_all = vcat(vcat(δ_σ[idx_gammachange_σ],δ_UN), δ_U2[idx_VaryU2])

# linear_fit(x, y)
# Lₜ = b + m_Lt * Lₒ
(b, m_Lt) = linear_fit(LO_all, Lt_all)
# δ = b + m_δ * Lₒ
(b, m_δ) = linear_fit(LO_all, δ_all)
# Lₜ = b + m_Ltδ * δ
(b, m_Ltδ) = linear_fit(δ_all, Lt_all)

@info "Lₜ = b + $m_Lt * Lₒ"
@info "δ = b + $m_δ * Lₒ"
@info "Lₜ = b + $m_Ltδ * δ"

@info "Lₒ = b + $(1/m_Lt) * Lₜ"
@info "Lₒ = b + $(1/m_δ) * δ"
(b, m_δ2) = linear_fit(δ_all .^ 2, LO_all .^ 2)
@info "Lₒ² = b + $m_δ2 * δ²"
