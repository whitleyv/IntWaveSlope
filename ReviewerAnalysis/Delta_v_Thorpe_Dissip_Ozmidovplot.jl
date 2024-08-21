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

# setnames
filesetnames =  "SetList_mp.jld2"

# Thorpe Scale Values
filescalename_UN = dpath * "DeltavThorpeDissip.jld2"
filescalename_σ = dpath * "DeltavThorpeDissip_VarySigma.jld2"

# Dissipation Values
fileεname_U = dpath * "DeltavAll_mtest.jld2"
fileεname_UN = dpath * "DeltavDissip_mp.jld2"
fileεname_σ = dpath * "DeltavDissip_mp_varyσ.jld2"

# error bars
filescalename_extraStats = dpath_reviewer * "DeltavAll_ExtraStats.jld2"

file_sn = jldopen(filesetnames, "r+")
file_scale_UN = jldopen(filescalename_UN, "r+")
file_scale_extraStats = jldopen(filescalename_extraStats, "r+")

# getting the surrounding simulation params out
file_sn = jldopen(filesetnames, "r+")

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
N_σ = file_sn["N_varyσ"]

δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
δ_U2 = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.45, 0.3, 0.4, 0.5, 0.55]./ N_UN[1]

γ_σ = file_sn["γ_varyσ"]

file_Lt_UN = jldopen(filescalename_UN, "r+")
file_Lt_σ = jldopen(filescalename_σ, "r+")

file_ε_UN = jldopen(fileεname_UN, "r+")
file_ε_U = jldopen(fileεname_U, "r+")
file_ε_σ = jldopen(fileεname_σ, "r+")

skeys = keys(file_Lt_UN)

thorpe_rms_UN = file_Lt_UN["thorpe_rms"]
thorpe_rms_σ = file_Lt_σ["thorpe_rms"]
Lt_hmax_tavg_U2 = file_ε_U["thorpe_hmax_tavg"]

#Lt_eps_full_usingrms_UN = file_Lt_UN["Lt_eps_full_usingrms"]
#Lt_eps_full_usingrms_σ = file_Lt_σ["Lt_eps_full_usingrms"]
#LtN_U2 = Lt_hmax_tavg_U2.^2 .* (N_UN[1]).^ 3 

eps_endAvg_UN = file_ε_UN["eps_endAvg"]
eps_endAvg_U2 = file_ε_U["eps_endAvg"]
eps_endAvg_σ = file_ε_σ["eps_endAvg"]

Dissip_stats = file_scale_extraStats["Dissip_stats"]
Thorpe_stats = file_scale_extraStats["Thorpe_stats"]

# row 7 = 20th percentile, row 8 = 80th percentile
Thorpe_20thpercentile = Thorpe_stats[:, 7]
Dissip_20thpercentile = Dissip_stats[:, 7]
Thorpe_80thpercentile = Thorpe_stats[:, 8]
Dissip_80thpercentile = Dissip_stats[:, 8]

idx_subcritical_σ = findall(γ_σ .< 1)
idx_gammachange_σ = findall(γ_σ .!= 1.9)

idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10
idx_VaryU2 = findall(δ_U2 .> 120)

N_all = vcat(N_UN, N_σ[idx_gammachange_σ])

# this is error on each value of ϵ and Lₜ, but we are using ϵ/N³ and Lₜ²
Lt2_20th = Thorpe_20thpercentile .^ 2
ϵN_20th = Dissip_20thpercentile ./ (N_all .^ 3)

Lt2_80th = Thorpe_80thpercentile .^ 2 
ϵN_80th = Dissip_80thpercentile ./ (N_all .^ 3)

# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = rich("L", subscript("t"), superscript("2")," [m²]"), yscale = log10, xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax1, 3, 1000, 10, 30000)

    ax2 = Axis(ga[2, 2],  ylabel = rich("h", subscript("w"), superscript("2")," [m²]"), yscale = log10,xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax2, 3, 1000, 10, 30000)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))
    vsp1 = scatter!(ax1, eps_endAvg_σ[idx_gammachange_σ] ./ N_σ[idx_gammachange_σ].^3, thorpe_rms_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, eps_endAvg_σ[idx_subcritical_σ] ./ N_σ[idx_subcritical_σ].^3, thorpe_rms_σ[idx_subcritical_σ].^2, markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, eps_endAvg_U2[idx_VaryU2] ./ (N_UN[1]).^ 3, Lt_hmax_tavg_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryN] ./ (N_UN[idx_VaryN]).^ 3, thorpe_rms_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, thorpe_rms_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax2, eps_endAvg_σ[idx_gammachange_σ]./ (N_σ[idx_gammachange_σ].^3), δ_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_σ[idx_subcritical_σ]./ (N_σ[idx_subcritical_σ].^3), δ_σ[idx_subcritical_σ].^2, markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, eps_endAvg_U2[idx_VaryU2]./ (N_UN[1]).^ 3, δ_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, eps_endAvg_UN[idx_VaryN]  ./ (N_UN[idx_VaryN]).^ 3, δ_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, δ_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
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
savename = apath * "Paper_Ozmidov_v_Thorpe_v_Delta_rms_log"
save(savename * ".png", f, px_per_unit = 2)

# best fit lines between ϵ/N³ vs Lₜ² or h_w², ie. h_w² = m_δ ϵ/N³
epsN_endAvg_all = vcat(eps_endAvg_σ[idx_gammachange_σ] ./ N_σ[idx_gammachange_σ].^3, eps_endAvg_UN  ./ (N_UN.^ 3 ))
Lt2_all = vcat(thorpe_rms_σ[idx_gammachange_σ].^2, thorpe_rms_UN.^2)
δ2_all = vcat(δ_σ[idx_gammachange_σ].^2,δ_UN .^2)
(b, m_Lt) = linear_fit(epsN_endAvg_all, Lt2_all)
(b, m_δ) = linear_fit(epsN_endAvg_all, δ2_all)
(b, m_Ltδ) = linear_fit(δ2_all, Lt2_all)

# best fit lines between ϵ vs Lₜ²N³ or h_w²N³, ie. m_hN h_w²N³ = ϵ
eps_endAvg_all = vcat(eps_endAvg_σ[idx_gammachange_σ], eps_endAvg_UN)
Lt2N3_all = vcat(thorpe_rms_σ[idx_gammachange_σ].^2 .* N_σ[idx_gammachange_σ].^3, thorpe_rms_UN.^2 .* N_UN.^ 3)
δ2N3_all = vcat(δ_σ[idx_gammachange_σ].^2 .* N_σ[idx_gammachange_σ].^3, δ_UN .^2 .* N_UN.^ 3 )
(b, m_LtN) = linear_fit(Lt2N3_all, eps_endAvg_all)
(b, m_hN) = linear_fit(δ2N3_all, eps_endAvg_all)
# now the difference between the U and N data

Lt2N3_U =  (thorpe_rms_UN[idx_VaryU]).^2 .* (N_UN[idx_VaryU]).^ 3
Lt2N3_N =  (thorpe_rms_UN[idx_VaryN]).^2 .* (N_UN[idx_VaryN]).^ 3
δ2N3_U =  (δ_UN[idx_VaryU]).^2 .* (N_UN[idx_VaryU]).^ 3
δ2N3_N =  (δ_UN[idx_VaryN]).^2 .* (N_UN[idx_VaryN]).^ 3

(b, m_LtN_U) = linear_fit(Lt2N3_U, eps_endAvg_UN[idx_VaryU])
(b, m_LtN_N) = linear_fit(Lt2N3_N, eps_endAvg_UN[idx_VaryN])

(b, m_hN_U) = linear_fit(δ2N3_U, eps_endAvg_UN[idx_VaryU])
(b, m_hN_N) = linear_fit(δ2N3_N, eps_endAvg_UN[idx_VaryN])


# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = rich("L", subscript("t"), superscript("2")," [m²]"), yscale = log10, xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax1, 3, 1000, 10, 30000)

    ax2 = Axis(ga[2, 2],  ylabel = rich("h", subscript("w"), superscript("2")," [m²]"), yscale = log10,xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax2, 3, 1000, 10, 30000)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    # error on thorpe scale measurements
    band!(ax1, eps_endAvg_UN[idx_VaryU]./ N_UN[idx_VaryU].^3, Lt2_20th[idx_VaryU], Lt2_80th[idx_VaryU], color = (:gray70, 1.0))
    band!(ax1, eps_endAvg_UN[idx_VaryN]./ N_UN[idx_VaryN].^3, Lt2_20th[idx_VaryN], Lt2_80th[idx_VaryN], color = (:gray70, 1.0))
    # error on dissipation calculation
    pU = Polygon(Point2.(vcat(ϵN_80th[idx_VaryU], reverse!(ϵN_20th[idx_VaryU])), vcat(thorpe_rms_UN[idx_VaryU].^2, reverse!(thorpe_rms_UN[idx_VaryU].^2))))
    poly!(ax1, pU, color = (:gray70, 1.0))
    pN = Polygon(Point2.(vcat(ϵN_80th[idx_VaryN], reverse!(ϵN_20th[idx_VaryN])), vcat(thorpe_rms_UN[idx_VaryN].^2, reverse!(thorpe_rms_UN[idx_VaryN].^2))))
    poly!(ax1, pN, color = (:gray70, 1.0))

    vsp1 = scatter!(ax1, eps_endAvg_σ[idx_gammachange_σ] ./ N_σ[idx_gammachange_σ].^3, thorpe_rms_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, eps_endAvg_σ[idx_subcritical_σ] ./ N_σ[idx_subcritical_σ].^3, thorpe_rms_σ[idx_subcritical_σ].^2, markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, eps_endAvg_U2[idx_VaryU2] ./ (N_UN[1]).^ 3, Lt_hmax_tavg_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryN] ./ (N_UN[idx_VaryN]).^ 3, thorpe_rms_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, thorpe_rms_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

     # error on dissipation calculation
    pU = Polygon(Point2.(vcat(ϵN_80th[idx_VaryU], reverse!(ϵN_20th[idx_VaryU])), vcat(δ_UN[idx_VaryU].^2, reverse!(δ_UN[idx_VaryU].^2))))
        poly!(ax2, pU, color = (:gray70, 1.0))
    pN = Polygon(Point2.(vcat(ϵN_80th[idx_VaryN], reverse!(ϵN_20th[idx_VaryN])), vcat(δ_UN[idx_VaryN].^2, reverse!(δ_UN[idx_VaryN].^2))))
        poly!(ax2, pN, color = (:gray70, 1.0))

    scatter!(ax2, eps_endAvg_σ[idx_gammachange_σ]./ (N_σ[idx_gammachange_σ].^3), δ_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_σ[idx_subcritical_σ]./ (N_σ[idx_subcritical_σ].^3), δ_σ[idx_subcritical_σ].^2, markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, eps_endAvg_U2[idx_VaryU2]./ (N_UN[1]).^ 3, δ_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, eps_endAvg_UN[idx_VaryN]  ./ (N_UN[idx_VaryN]).^ 3, δ_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, δ_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
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
savename = apath * "Paper_Ozmidov_v_Thorpe_v_Delta_rms_log_errorclouds"
save(savename * ".png", f, px_per_unit = 2)



# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = rich("L", subscript("t"), superscript("2")," [m²]"), yscale = log10, xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax1, 3, 1000, 10, 30000)

    ax2 = Axis(ga[2, 2],  ylabel = rich("h", subscript("w"), superscript("2")," [m²]"), yscale = log10,xscale = log10,
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    limits!(ax2, 3, 1000, 10, 30000)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    # error on thorpe scale measurements goes in y direction
    rangebars!(ax1, eps_endAvg_σ[idx_gammachange_σ]./ N_σ[idx_gammachange_σ].^3, Lt2_20th[22 .+ idx_gammachange_σ], Lt2_80th[22 .+ idx_gammachange_σ], linewidth = 3, color = :darkgreen)
    rangebars!(ax1, eps_endAvg_UN[idx_VaryU]./ N_UN[idx_VaryU].^3, Lt2_20th[idx_VaryU], Lt2_80th[idx_VaryU], linewidth = 3, color = :dodgerblue2)
    rangebars!(ax1, eps_endAvg_UN[idx_VaryN]./ N_UN[idx_VaryN].^3, Lt2_20th[idx_VaryN], Lt2_80th[idx_VaryN], linewidth = 3, color = :firebrick2)

    # error on dissipation scale measurements goes in x direction
    rangebars!(ax1, thorpe_rms_σ[idx_gammachange_σ].^2, ϵN_20th[22 .+ idx_gammachange_σ], ϵN_80th[22 .+ idx_gammachange_σ], linewidth = 3, direction = :x, color = :darkgreen)
    rangebars!(ax1, thorpe_rms_UN[idx_VaryU].^2, ϵN_20th[idx_VaryU], ϵN_80th[idx_VaryU], linewidth = 3, direction = :x, color = :dodgerblue2)
    rangebars!(ax1, thorpe_rms_UN[idx_VaryN].^2, ϵN_20th[idx_VaryN], ϵN_80th[idx_VaryN], linewidth = 3, direction = :x, color = :firebrick2)

    vsp1 = scatter!(ax1, eps_endAvg_σ[idx_gammachange_σ] ./ N_σ[idx_gammachange_σ].^3, thorpe_rms_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, eps_endAvg_σ[idx_subcritical_σ] ./ N_σ[idx_subcritical_σ].^3, thorpe_rms_σ[idx_subcritical_σ].^2, markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, eps_endAvg_U2[idx_VaryU2] ./ (N_UN[1]).^ 3, Lt_hmax_tavg_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryN] ./ (N_UN[idx_VaryN]).^ 3, thorpe_rms_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, thorpe_rms_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    # error on dissipation scale measurements goes in x direction
    rangebars!(ax2, δ_σ[idx_gammachange_σ].^2, ϵN_20th[22 .+ idx_gammachange_σ], ϵN_80th[22 .+ idx_gammachange_σ], linewidth = 3, direction = :x, color = :darkgreen)
    rangebars!(ax2, δ_UN[idx_VaryU].^2, ϵN_20th[idx_VaryU], ϵN_80th[idx_VaryU], linewidth = 3, direction = :x, color = :dodgerblue2)
    rangebars!(ax2, δ_UN[idx_VaryN].^2, ϵN_20th[idx_VaryN], ϵN_80th[idx_VaryN], linewidth = 3, direction = :x, color = :firebrick2)
    
    scatter!(ax2, eps_endAvg_σ[idx_gammachange_σ]./ (N_σ[idx_gammachange_σ].^3), δ_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_σ[idx_subcritical_σ]./ (N_σ[idx_subcritical_σ].^3), δ_σ[idx_subcritical_σ].^2, markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, eps_endAvg_U2[idx_VaryU2]./ (N_UN[1]).^ 3, δ_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, eps_endAvg_UN[idx_VaryN]  ./ (N_UN[idx_VaryN]).^ 3, δ_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, δ_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
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
savename = apath * "Paper_Ozmidov_v_Thorpe_v_Delta_rms_log_errorbars"
save(savename * ".png", f, px_per_unit = 2)


# thorpe scale vs dissipation
f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  ylabel = rich("L", subscript("t"), superscript("2")," [m²]"),
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    #limits!(ax1, 3, 1000, 10, 30000)

    ax2 = Axis(ga[2, 2],  ylabel = rich("h", subscript("w"), superscript("2")," [m²]"),
    xlabel = rich("ϵ̄/N", subscript("0"), superscript("3"), " [m²]"), xlabelsize=35, ylabelsize=35)
    #limits!(ax2, 3, 1000, 10, 30000)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    # error on thorpe scale measurements goes in y direction
    rangebars!(ax1, eps_endAvg_σ[idx_gammachange_σ]./ N_σ[idx_gammachange_σ].^3, Lt2_20th[22 .+ idx_gammachange_σ], Lt2_80th[22 .+ idx_gammachange_σ], linewidth = 3, color = :darkgreen)
    rangebars!(ax1, eps_endAvg_UN[idx_VaryU]./ N_UN[idx_VaryU].^3, Lt2_20th[idx_VaryU], Lt2_80th[idx_VaryU], linewidth = 3, color = :dodgerblue2)
    rangebars!(ax1, eps_endAvg_UN[idx_VaryN]./ N_UN[idx_VaryN].^3, Lt2_20th[idx_VaryN], Lt2_80th[idx_VaryN], linewidth = 3, color = :firebrick2)

    # error on dissipation scale measurements goes in x direction
    rangebars!(ax1, thorpe_rms_σ[idx_gammachange_σ].^2, ϵN_20th[22 .+ idx_gammachange_σ], ϵN_80th[22 .+ idx_gammachange_σ], linewidth = 3, direction = :x, color = :darkgreen)
    rangebars!(ax1, thorpe_rms_UN[idx_VaryU].^2, ϵN_20th[idx_VaryU], ϵN_80th[idx_VaryU], linewidth = 3, direction = :x, color = :dodgerblue2)
    rangebars!(ax1, thorpe_rms_UN[idx_VaryN].^2, ϵN_20th[idx_VaryN], ϵN_80th[idx_VaryN], linewidth = 3, direction = :x, color = :firebrick2)

    vsp1 = scatter!(ax1, eps_endAvg_σ[idx_gammachange_σ] ./ N_σ[idx_gammachange_σ].^3, thorpe_rms_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, eps_endAvg_σ[idx_subcritical_σ] ./ N_σ[idx_subcritical_σ].^3, thorpe_rms_σ[idx_subcritical_σ].^2, markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_mean = scatter!(ax1, eps_endAvg_U2[idx_VaryU2] ./ (N_UN[1]).^ 3, Lt_hmax_tavg_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryN] ./ (N_UN[idx_VaryN]).^ 3, thorpe_rms_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, thorpe_rms_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    # error on dissipation scale measurements goes in x direction
    rangebars!(ax2, δ_σ[idx_gammachange_σ].^2, ϵN_20th[22 .+ idx_gammachange_σ], ϵN_80th[22 .+ idx_gammachange_σ], linewidth = 3, direction = :x, color = :darkgreen)
    rangebars!(ax2, δ_UN[idx_VaryU].^2, ϵN_20th[idx_VaryU], ϵN_80th[idx_VaryU], linewidth = 3, direction = :x, color = :dodgerblue2)
    rangebars!(ax2, δ_UN[idx_VaryN].^2, ϵN_20th[idx_VaryN], ϵN_80th[idx_VaryN], linewidth = 3, direction = :x, color = :firebrick2)
    
    scatter!(ax2, eps_endAvg_σ[idx_gammachange_σ]./ (N_σ[idx_gammachange_σ].^3), δ_σ[idx_gammachange_σ].^2, markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_σ[idx_subcritical_σ]./ (N_σ[idx_subcritical_σ].^3), δ_σ[idx_subcritical_σ].^2, markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, eps_endAvg_U2[idx_VaryU2]./ (N_UN[1]).^ 3, δ_U2[idx_VaryU2].^2, markersize = 25, marker=:utriangle, 
                color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, eps_endAvg_UN[idx_VaryN]  ./ (N_UN[idx_VaryN]).^ 3, δ_UN[idx_VaryN].^2, markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, eps_endAvg_UN[idx_VaryU] ./ (N_UN[idx_VaryU]).^ 3, δ_UN[idx_VaryU].^2, markersize = 25, marker=:dtriangle, 
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
savename = apath * "Paper_Ozmidov_v_Thorpe_v_Delta_rms_errorbars"
save(savename * ".png", f, px_per_unit = 2)
