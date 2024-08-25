using JLD2
using Statistics
using CairoMakie
using Printf
using GeometryBasics

dpath = "ReviewerAnalysis/ReviewerData/"
apath = "ReviewerAnalysis/ReviewerPlots/"

filescalename = dpath * "DeltavAll_ExtraStats.jld2"
fileerrorname = dpath * "Delta_v_all_Confint.jld2"
filesetnames =  "SetList_mp.jld2"

file_sn = jldopen(filesetnames, "r+")
file_scale = jldopen(filescalename, "r+")
file_error = jldopen(fileerrorname, "r+")

Thorpe_rms = file_scale["Thorpe_rms"]
Thorpe_stats = file_scale["Thorpe_stats"]
Nheight_stats = file_scale["Nheight_stats"]
Cheight_stats = file_scale["Cheight_stats"]
Dissip_stats = file_scale["Dissip_stats"]

lower_confint_tracer = file_error["lower_confint_tracer"]
upper_confint_tracer = file_error["upper_confint_tracer"]
lower_confint_stratanom = file_error["lower_confint_stratanom"]
upper_confint_stratanom = file_error["upper_confint_stratanom"]
lower_confint_thorpe = file_error["lower_confint_thorpe"]
upper_confint_thorpe = file_error["upper_confint_thorpe"]
lower_confint_dissiptime = file_error["lower_confint_dissiptime"]
upper_confint_dissiptime = file_error["upper_confint_dissiptime"]
lower_confint_dissipall = file_error["lower_confint_dissipall"]
upper_confint_dissipall = file_error["upper_confint_dissipall"]

# thorpe is actually exponential scaling:
λ̄ = (Thorpe_stats[:, 5] .- 2) ./ (Thorpe_stats[:, 1] .* Thorpe_stats[:, 5])
degs_freedom = 120
λ_upper =  λ̄ .* (1 - 1.96/√degs_freedom)
λ_lower =  λ̄ .* (1 + 1.96/√degs_freedom)
Thorpe_mean_upper = 1 ./ λ_upper
Thorpe_mean_lower = 1 ./ λ_lower

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
γ_σ = file_sn["γ_varyσ"]

idx_subcritical_σ = findall(γ_σ .< 1)
idx_gammachange_σ = findall(γ_σ .!= 1.9)
idx_supcritical_σ = findall((γ_σ .> 1) .& ( γ_σ .!= 1.9))
idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10

# extra_stats = [tmean, tmedian, tmax, tmin, tnum, tstd, twentieth_percentile, eightieth_percentile]

# rich("L", subscript("tr")," [m]")
#  0:40:160
function plot_stats(stats, savename, stat_string, xtix, ytix, lims)

    (xmin, xmax, ymin, ymax) = lims
    low_err20 = stats[:,1] .- stats[:, 7]
    high_err80 = stats[:, 8] .- stats[:,1] 

    stddev_error = stats[:, 6]  #./ sqrt.(stats[:, 5])

    f = Figure(resolution = (1200, 1200), fontsize=26)
        ga = f[1, 1] = GridLayout()

        ax1 = Axis(ga[2, 1], ylabel = stat_string, 
        xlabelsize=30, ylabelsize=30)
        ax1.xticks = xtix
        ax1.yticks = ytix
        limits!(ax1, xmin, xmax, ymin, ymax) 

        ax2 = Axis(ga[2, 2], 
        xlabelsize=30, ylabelsize=30)
        ax2.xticks = xtix
        ax2.yticks = ytix
        limits!(ax2, xmin, xmax, ymin, ymax) 

        ax3 = Axis(ga[3, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = stat_string, 
        xlabelsize=30, ylabelsize=30)
        ax3.xticks = xtix
        ax3.yticks = ytix
        limits!(ax3, xmin, xmax, ymin, ymax) 

        ax4 = Axis(ga[3, 2],  xlabel = rich("h", subscript("w"), " [m]"), 
        xlabelsize=30, ylabelsize=30)
        ax4.xticks = xtix
        ax4.yticks = ytix
        limits!(ax4, xmin, xmax, ymin, ymax)

        hidexdecorations!(ax1, grid = false)
        hidedecorations!(ax2, grid = false)
        hideydecorations!(ax4, grid = false)

        p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
        p_small = decompose(Point2f, Circle(Point2f(0), 0.5))
    
        # mean + 80/20 percentile + extrema
        vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
        vnp = scatter!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        vump = scatter!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        errsp = errorbars!(ax1, δ_σ[idx_gammachange_σ], stats[22 .+ idx_gammachange_σ, 1], 
                    low_err20[22 .+ idx_gammachange_σ], high_err80[22 .+ idx_gammachange_σ]; linewidth = 3, color = :darkgreen)
        errN = errorbars!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 1], 
                    low_err20[idx_VaryN], high_err80[idx_VaryN]; linewidth = 3, color = :firebrick2)
        errV = errorbars!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 1], 
                    low_err20[idx_VaryU], high_err80[idx_VaryU]; linewidth = 2, color = :dodgerblue2)
        scatter!(ax1, δ_σ[idx_gammachange_σ], stats[22 .+ idx_gammachange_σ, 3], markersize = 25, marker=:hline, 
                    color = :darkgreen, strokewidth = 1, strokecolor = :black)
        extrem = scatter!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 4], markersize = 25, marker = :hline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 4], markersize = 25, marker=:hline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 3], markersize = 25, marker = :hline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 3], markersize = 25, marker=:hline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # median
        scatter!(ax2, δ_σ[idx_gammachange_σ], stats[22 .+ idx_gammachange_σ, 2],  markersize = 25, marker=:star4, 
                 color =:darkgreen, strokewidth = 1, strokecolor = :black)
        scatter!(ax2, δ_σ[idx_subcritical_σ], stats[22 .+ idx_subcritical_σ, 2], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
        scatter!(ax2, δ_UN[idx_VaryN], stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax2, δ_UN[idx_VaryU], stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        # std dev 
        scatter!(ax3, δ_σ[idx_gammachange_σ], stats[22 .+ idx_gammachange_σ, 6],  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        scatter!(ax3, δ_σ[idx_subcritical_σ], stats[22 .+ idx_subcritical_σ, 6], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
        scatter!(ax3, δ_UN[idx_VaryN], stats[idx_VaryN, 6], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax3, δ_UN[idx_VaryU], stats[idx_VaryU, 6], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # mean + regular error bars
        vsp1 = scatter!(ax4, δ_σ[idx_gammachange_σ], stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        vsbp = scatter!(ax4, δ_σ[idx_subcritical_σ], stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
        scatter!(ax4, δ_UN[idx_VaryN], stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax4, δ_UN[idx_VaryU], stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        errsp2 = errorbars!(ax1, δ_σ[idx_gammachange_σ], stats[22 .+ idx_gammachange_σ, 1], 
                stddev_error[22 .+ idx_gammachange_σ], stddev_error[22 .+ idx_gammachange_σ]; linewidth = 3, linestyle = :dot, color = :darkgreen)
        errN2 = errorbars!(ax4, δ_UN[idx_VaryN], stats[idx_VaryN, 1], 
                stddev_error[idx_VaryN], stddev_error[idx_VaryN]; linewidth = 3, linestyle = :dot, color = :firebrick2)
        errV2 = errorbars!(ax4, δ_UN[idx_VaryU], stats[idx_VaryU, 1], 
                stddev_error[idx_VaryU], stddev_error[idx_VaryU]; linewidth = 2, linestyle = :dot, color = :dodgerblue2)

        Legend( ga[1, 1:2],  [vnp, vump, vsp1, vsbp, errN, extrem], ["Vary N₀", "Vary V₀", "Vary γ", "Subcritical", "Error", "Extrema"],
            tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
            margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
            labelsize = 35,
            halign = :center, valign = :top, orientation = :horizontal)

        rowsize!(ga,1, Auto(0.05))       

        Label(ga[2, 1, Top()], "Mean + 20/80th Percentiles",
                        fontsize = 30,
                        font = :bold,
                        padding = (-10, 5, 5, 10),
                        )
        Label(ga[2, 2, Top()], "Median",
                        fontsize = 30,
                        font = :bold,
                        padding = (0, 5, 5, 10),
                        )
        Label(ga[3, 1, Top()], "Standard Dev.",
                        fontsize = 30,
                        font = :bold,
                        padding = (-10, 5, 5, 10),
                        )
        Label(ga[3, 2, Top()], "Mean ± σ",
                        fontsize = 30,
                        font = :bold,
                        padding = (0, 5, 5, 10),
                        )

    save(savename * ".png", f) 
end

function logplot_stats(stats, savename, stat_string, indep_string, xtix, ytix, lims)

    (xmin, xmax, ymin, ymax) = lims
    #low_err20 = stats[:,1] .- stats[:, 7]
    #high_err80 = stats[:, 8] .- stats[:,1]

    stddev_error = 2*stats[:, 6]  ./ sqrt.(60)

    f = Figure(resolution = (1200, 1200), fontsize=26)
        ga = f[1, 1] = GridLayout()

        ax1 = Axis(ga[2, 1], ylabel = stat_string, 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
        #ax1.xticks = xtix
        #ax1.yticks = ytix
        limits!(ax1, ymin, ymax, xmin, xmax) 

        ax2 = Axis(ga[2, 2], 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
       #ax2.xticks = xtix
       #ax2.yticks = ytix
        limits!(ax2, ymin, ymax, xmin, xmax) 

        ax3 = Axis(ga[3, 1],  xlabel = indep_string, ylabel = stat_string, 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
        #ax3.xticks = xtix
        #ax3.yticks = ytix
        limits!(ax3, ymin, ymax, xmin, xmax) 

        ax4 = Axis(ga[3, 2],  xlabel = indep_string, 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
        #ax4.xticks = xtix
        #ax4.yticks = ytix
        limits!(ax4, ymin, ymax, xmin, xmax)

        hidexdecorations!(ax1, grid = false)
        hidedecorations!(ax2, grid = false)
        hideydecorations!(ax4, grid = false)

        p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
        p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

        # mean + 80/20 percentile + extrema
        errsp = rangebars!(ax1, δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3),
        stats[22 .+ idx_gammachange_σ, 7], stats[22 .+ idx_gammachange_σ, 8]; direction = :x, linewidth = 3, color = :darkgreen)
        errN = rangebars!(ax1, (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN],
        stats[idx_VaryN, 7], stats[idx_VaryN, 8]; linewidth = 3, direction = :x, color = :firebrick2)
        errV = rangebars!(ax1, (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU],
        stats[idx_VaryU, 7], stats[idx_VaryU, 8]; linewidth = 2, direction = :x, color = :dodgerblue2)
        vsp1 = scatter!(ax1, stats[22 .+ idx_gammachange_σ, 1], δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3),  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        vsbp = scatter!(ax1, stats[22 .+ idx_subcritical_σ, 1], δ_σ[idx_subcritical_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod)
        vnp = scatter!(ax1, stats[idx_VaryN, 1], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        vump = scatter!(ax1, stats[idx_VaryU, 1], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)        
        scatter!(ax1, stats[22 .+ idx_gammachange_σ, 3],  δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, marker=:vline, 
                    color = :darkgreen, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, stats[22 .+ idx_gammachange_σ, 4], δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, marker=:vline, 
                    color = :darkgreen, strokewidth = 1, strokecolor = :black)
        extrem = scatter!(ax1, stats[idx_VaryN, 4],  (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN],  markersize = 25, marker = :vline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1,  stats[idx_VaryU, 4], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], markersize = 25, marker=:vline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1,  stats[idx_VaryN, 3], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN], markersize = 25, marker = :vline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1,  stats[idx_VaryU, 3], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], markersize = 25, marker=:vline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # median
        scatter!(ax2, stats[22 .+ idx_gammachange_σ, 2], δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3),  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        scatter!(ax2, stats[22 .+ idx_subcritical_σ, 2], δ_σ[idx_subcritical_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
        scatter!(ax2, stats[idx_VaryN, 2], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax2, stats[idx_VaryU, 2], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        # std dev 
        scatter!(ax3, stats[22 .+ idx_gammachange_σ, 6], δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        scatter!(ax3, stats[22 .+ idx_subcritical_σ, 6], δ_σ[idx_subcritical_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
        scatter!(ax3, stats[idx_VaryN, 6],  (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax3,  stats[idx_VaryU, 6], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # mean + regular error bars
        errsp2 = errorbars!(ax4, stats[22 .+ idx_gammachange_σ, 1], δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3), 
                    stddev_error[22 .+ idx_gammachange_σ]; direction = :x, linewidth = 3, linestyle = :dot, color = :darkgreen)
        errN2 = errorbars!(ax4, stats[idx_VaryN, 1], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN],
                stddev_error[idx_VaryN]; linewidth = 3, linestyle = :dot, direction = :x, color = :firebrick2)
        errV2 = errorbars!(ax4, stats[idx_VaryU, 1], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU],
                stddev_error[idx_VaryU]; linewidth = 2, linestyle = :dot, direction = :x, color = :dodgerblue2)
        scatter!(ax4, stats[22 .+ idx_gammachange_σ, 1], δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        scatter!(ax4, stats[22 .+ idx_subcritical_σ, 1], δ_σ[idx_subcritical_σ] .^2 .* (N_UN[1] .^ 3), markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
        scatter!(ax4, stats[idx_VaryN, 1], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax4, stats[idx_VaryU, 1], (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        Legend( ga[1, 1:2],  [vnp, vump, vsp1, vsbp, extrem], ["Vary N₀", "Vary V₀", "Vary γ", "Subcritical", 
             "Extrema"],
            tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
            margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
            labelsize = 35,
            halign = :center, valign = :top, orientation = :horizontal)

        rowsize!(ga,1, Auto(0.05))       

        Label(ga[2, 1, Top()], "Mean + 20/80th Percentiles",
                        fontsize = 30,
                        font = :bold,
                        padding = (-10, 5, 5, 10),
                        )
        Label(ga[2, 2, Top()], "Median",
                        fontsize = 30,
                        font = :bold,
                        padding = (0, 5, 5, 10),
                        )
        Label(ga[3, 1, Top()], "Standard Dev.",
                        fontsize = 30,
                        font = :bold,
                        padding = (-10, 5, 5, 10),
                        )
        Label(ga[3, 2, Top()], "Mean + 2σ/√numT Error",
                        fontsize = 30,
                        font = :bold,
                        padding = (0, 5, 5, 10),
                        )

    save(savename * ".png", f) 
end

plot_stats(Thorpe_stats, apath * "Delta_v_Thorpe_ExtraStats", "Lₜ [m]", 0:40:200, 0:20:200, (0,160,0,80))

plot_stats(Nheight_stats, apath * "Delta_v_Nheight_ExtraStats", "Lₙ [m]", 0:40:200, 0:40:200, (0,160,0,300))

plot_stats(Cheight_stats, apath * "Delta_v_Cheight_ExtraStats", rich("L", subscript("tr")," [m]"), 0:40:200, 0:40:200, (0,160,0,300))


# usually it's ϵ vs hw^2N^3, we can instead do √(ϵ/ N^3) = hw
logplot_stats(Dissip_stats, apath * "Delta_v_Dissip_ExtraStats", rich("h", subscript("w"), superscript("2"),"N", subscript("0"), superscript("3"), " [m²s", superscript("-3"), "]"), "ϵ", (0:5*1e-4:1e-3),  (0:1e-5:3e-5), (1e-6, 10^(-2.75), 1e-8, 10^(-3.2)))

# plotting all medians and original means:

f = Figure(resolution = (1200, 1200), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1], ylabel =  rich("L", subscript("tr")," [m]"), 
        xlabelsize=30, ylabelsize=30)
        ax1.xticks = 0:40:200
        ax1.yticks = 0:40:200
        limits!(ax1, 0,165,0,210) 

    ax2 = Axis(ga[2, 2], ylabel = "Lₙ [m]", 
        xlabelsize=30, ylabelsize=30)
        ax2.xticks = 0:40:200
        ax2.yticks = 0:40:200
        limits!(ax2, 0,165,0,165) 

    ax3 = Axis(ga[3, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel =  "Lₜ [m]", 
        xlabelsize=30, ylabelsize=30)
        ax3.xticks = 0:40:200
        ax3.yticks = 0:10:200
        limits!(ax3, 0,165,0,45) 

    ax4 = Axis(ga[3, 2],  xlabel = rich("h", subscript("w"), superscript("2"),"N", subscript("0"), superscript("3"), " [m²s⁻³]"), 
        ylabel = "ϵ [m²s⁻³]", xscale = log10, yscale = log10,
        xlabelsize=30, ylabelsize=30)
        limits!(ax4, 1e-6, 10^(-2.75), 1e-8, 10^(-4),)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    #### TRACER
    # error bars, 80/20 percentile
    errsp = rangebars!(ax1, δ_σ[idx_gammachange_σ],  
            Cheight_stats[22 .+ idx_gammachange_σ, 7], Cheight_stats[22 .+ idx_gammachange_σ, 8]; linewidth = 3, color = :gray70)
    errN = rangebars!(ax1, δ_UN[idx_VaryN], 
            Cheight_stats[idx_VaryN, 7], Cheight_stats[idx_VaryN, 8]; linewidth = 3, color = :gray70)
    errV = rangebars!(ax1, δ_UN[idx_VaryU], 
            Cheight_stats[idx_VaryU, 7], Cheight_stats[idx_VaryU, 8]; linewidth = 2, color = :gray70)
    # mean
    scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
            color =:gray50)
    scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :gray50)
    scatter!(ax1, δ_UN[idx_VaryN], Cheight_stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:gray50)
    scatter!(ax1, δ_UN[idx_VaryU], Cheight_stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :gray50)
    # median
    vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_stats[22 .+ idx_gammachange_σ, 2],  markersize = 25, marker=:star4, 
             color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_stats[22 .+ idx_subcritical_σ, 2], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    #### STARTIFICATION

 #   rangebars!(ax2, δ_σ[idx_gammachange_σ],  
 #           Nheight_stats[22 .+ idx_gammachange_σ, 7], Nheight_stats[22 .+ idx_gammachange_σ, 8]; linewidth = 3, color = :gray70)
    rangebars!(ax2, δ_UN[idx_VaryN], 
            Nheight_stats[idx_VaryN, 7], Nheight_stats[idx_VaryN, 8]; linewidth = 3, color = :gray70)
    rangebars!(ax2, δ_UN[idx_VaryU], 
            Nheight_stats[idx_VaryU, 7], Nheight_stats[idx_VaryU, 8]; linewidth = 2, color = :gray70)
    # mean
        #    scatter!(ax2, δ_σ[idx_gammachange_σ], Nheight_stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
        #            color =:gray50)
  #  scatter!(ax2, δ_σ[idx_subcritical_σ], Nheight_stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
  #          marker=Polygon(p_big, [p_small]), color= :gray50)
    scatter!(ax2, δ_UN[idx_VaryN], Nheight_stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:gray50)
    scatter!(ax2, δ_UN[idx_VaryU], Nheight_stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :gray50)
    # median
 #   scatter!(ax2, δ_σ[idx_gammachange_σ], Nheight_stats[22 .+ idx_gammachange_σ, 2],  markersize = 25, marker=:star4, 
  #           color =:darkgreen, strokewidth = 1, strokecolor = :black)
   # scatter!(ax2, δ_σ[idx_subcritical_σ], Nheight_stats[22 .+ idx_subcritical_σ, 2], markersize = 25, 
    #        marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax2, δ_UN[idx_VaryN], Nheight_stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_UN[idx_VaryU], Nheight_stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    #### THORPE
    rangebars!(ax3, δ_σ[idx_gammachange_σ],  
            Thorpe_stats[22 .+ idx_gammachange_σ, 7], Thorpe_stats[22 .+ idx_gammachange_σ, 8]; linewidth = 3, color = :gray70)
    rangebars!(ax3, δ_UN[idx_VaryN], 
            Thorpe_stats[idx_VaryN, 7], Thorpe_stats[idx_VaryN, 8]; linewidth = 3, color = :gray70)
    rangebars!(ax3, δ_UN[idx_VaryU],  
            Thorpe_stats[idx_VaryU, 7], Thorpe_stats[idx_VaryU, 8]; linewidth = 2, color = :gray70)
    # mean
    scatter!(ax3, δ_σ[idx_gammachange_σ], Thorpe_stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
            color =:gray50)
    scatter!(ax3, δ_σ[idx_subcritical_σ], Thorpe_stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :gray50)
    scatter!(ax3, δ_UN[idx_VaryN], Thorpe_stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:gray50)
    scatter!(ax3, δ_UN[idx_VaryU], Thorpe_stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :gray50)
    # median
    scatter!(ax3, δ_σ[idx_gammachange_σ], Thorpe_stats[22 .+ idx_gammachange_σ, 2],  markersize = 25, marker=:star4, 
             color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_σ[idx_subcritical_σ], Thorpe_stats[22 .+ idx_subcritical_σ, 2], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax3, δ_UN[idx_VaryN], Thorpe_stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryU], Thorpe_stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    #### DISSIPATION

    rangebars!(ax4, δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3),  
         Dissip_stats[22 .+ idx_gammachange_σ, 7], Dissip_stats[22 .+ idx_gammachange_σ, 8]; linewidth = 3, color = :gray70)
    rangebars!(ax4, (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN],  
            Dissip_stats[idx_VaryN, 7], Dissip_stats[idx_VaryN, 8]; linewidth = 3, color = :gray70)
    rangebars!(ax4, (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], 
            Dissip_stats[idx_VaryU, 7], Dissip_stats[idx_VaryU, 8]; linewidth = 2, color = :gray70)
    # mean
    scatter!(ax4, δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3), Dissip_stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
            color =:gray50)
    scatter!(ax4, δ_σ[idx_subcritical_σ] .^2 .* (N_UN[1] .^ 3), Dissip_stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :gray50)
    scatter!(ax4, (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN], Dissip_stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:gray50)
    scatter!(ax4, (δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], Dissip_stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :gray50)
    # median
    scatter!(ax4, δ_σ[idx_gammachange_σ] .^2 .* (N_UN[1] .^ 3), Dissip_stats[22 .+ idx_gammachange_σ, 2],  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax4, δ_σ[idx_subcritical_σ] .^2 .* (N_UN[1] .^ 3), Dissip_stats[22 .+ idx_subcritical_σ, 2], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    scatter!(ax4, (δ_UN .^ 2 .* N_UN.^3)[idx_VaryN], Dissip_stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax4,(δ_UN .^ 2 .* N_UN.^3)[idx_VaryU], Dissip_stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)


    Legend( ga[1, 1:2],  [vnp, vump, vsp1, vsbp, errN], ["Vary N₀", "Vary V₀", "Vary γ", "Subcritical", 
    "20/80th Percentile"],
        tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
        margin = (10, 10, 50, -5), framevisible = false, patchlabelgap = 7,
        labelsize = 35, 
        halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.05))       
   #=
    Label(ga[2, 1, Top()], "Mean (Gray), Median (color),  20/80th Percentiles ",
    fontsize = 30,
    font = :bold,
    padding = (-10, 5, 5, 10),
    )
        =#
save(apath * "Delta_v_All_Medians.png", f) 

# plotting all means and medians and true error bars
function plotting_everything_for_stat(ax, x_UN, x_σ, yval_stats, yval_err_up, yval_err_down)


        # median
        scatter!(ax, x_σ[idx_gammachange_σ], yval_stats[22 .+ idx_gammachange_σ, 2],  markersize = 25, marker=:star4, 
                color =:gray50)
        scatter!(ax, x_σ[idx_subcritical_σ], yval_stats[22 .+ idx_subcritical_σ, 2], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :gray50)
        medn = scatter!(ax, x_UN[idx_VaryN], yval_stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                color =:gray50)
        scatter!(ax, x_UN[idx_VaryU], yval_stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                color = :gray50)

        # error bars
        errsp = rangebars!(ax, x_σ[idx_gammachange_σ],  
                yval_err_down[idx_gammachange_σ], yval_err_up[idx_gammachange_σ]; linewidth = 5, color = :darkgreen)
        errN = rangebars!(ax, x_UN[idx_VaryN], 
                yval_err_down[idx_VaryN], yval_err_up[idx_VaryN]; linewidth = 4, color = :firebrick2)
        errV = rangebars!(ax, x_UN[idx_VaryU], 
                yval_err_down[idx_VaryU], yval_err_up[idx_VaryU]; linewidth = 3, color = :dodgerblue2)

        # mean
        vsp1 = scatter!(ax, x_σ[idx_gammachange_σ], yval_stats[22 .+ idx_gammachange_σ, 1],  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
        vsbp = scatter!(ax, x_σ[idx_subcritical_σ], yval_stats[22 .+ idx_subcritical_σ, 1], markersize = 25, 
                marker=Polygon(p_big, [p_small]), color= :darkgoldenrod)
        vnp = scatter!(ax, x_UN[idx_VaryN], yval_stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
        vump = scatter!(ax, x_UN[idx_VaryU], yval_stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        
        return errsp, vsp1, vsbp, vnp, vump, medn
end

function plotting_everything_for_stat(ax, x_UN, x_σ, yval_stats, yval_err_up, yval_err_down, idx_supcritical_σ)


        # median
        scatter!(ax, x_σ[idx_supcritical_σ], yval_stats[22 .+ idx_supcritical_σ, 2],  markersize = 25, marker=:star4, 
                        color =:gray50)
        medn = scatter!(ax, x_UN[idx_VaryN], yval_stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                color =:gray50)
        scatter!(ax, x_UN[idx_VaryU], yval_stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                color = :gray50)

        # error bars
        errsp = rangebars!(ax, x_σ[idx_supcritical_σ],  
                yval_err_down[idx_supcritical_σ], yval_err_up[idx_supcritical_σ]; linewidth = 5, color = :darkgreen)
        errN = rangebars!(ax, x_UN[idx_VaryN], 
                yval_err_down[idx_VaryN], yval_err_up[idx_VaryN]; linewidth = 4, color = :firebrick2)
        errV = rangebars!(ax, x_UN[idx_VaryU], 
                yval_err_down[idx_VaryU], yval_err_up[idx_VaryU]; linewidth = 3, color = :dodgerblue2)

        # mean
        vsp1 = scatter!(ax, x_σ[idx_supcritical_σ], yval_stats[22 .+ idx_supcritical_σ, 1],  markersize = 25, marker=:star4, 
                        color =:darkgreen, strokewidth = 1, strokecolor = :black)
        vnp = scatter!(ax, x_UN[idx_VaryN], yval_stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
        vump = scatter!(ax, x_UN[idx_VaryU], yval_stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        
        return errsp, vsp1, vsbp, vnp, vump, medn
end

f = Figure(resolution = (1300, 1200), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1], ylabel =  rich("L", subscript("tr")," [m]"), 
        xlabelsize=30, ylabelsize=30)
        ax1.xticks = 0:40:200
        ax1.yticks = 0:40:200
        limits!(ax1, 0,165,0,210) 

    ax2 = Axis(ga[2, 2], ylabel = "Lₙ [m]", 
        xlabelsize=30, ylabelsize=30)
        ax2.xticks = 0:40:200
        ax2.yticks = 0:40:200
        limits!(ax2, 0,165,0,120) 

    ax3 = Axis(ga[3, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel =  "Lₜ [m]", 
        xlabelsize=30, ylabelsize=30)
        ax3.xticks = 0:40:200
        ax3.yticks = 0:10:200
        limits!(ax3, 0,165,0,40) 

    ax4 = Axis(ga[3, 2],  ylabel = "ϵ [m²s⁻³]", xlabel = rich("h", subscript("w"), superscript("2"),"N", subscript("0"), superscript("3"), " [m²s⁻³]"), 
        xscale = log10, yscale = log10,
        xlabelsize=30, ylabelsize=30)
        limits!(ax4, 10^(-5.7), 10^(-2.75), 1e-8, 10^(-4), )

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    #### TRACER
    (errsp, vsp1, vsbp, vnp, vump, medn) = plotting_everything_for_stat(ax1, δ_UN, δ_σ, Cheight_stats, upper_confint_tracer, lower_confint_tracer)
    #### STARTIFICATION

    (_, _, _, _, _, _) = plotting_everything_for_stat(ax2, δ_UN, δ_σ, Nheight_stats, upper_confint_stratanom, lower_confint_stratanom, idx_supcritical_σ)

    (_, _, _, _, _, _) = plotting_everything_for_stat(ax3, δ_UN, δ_σ, Thorpe_stats, Thorpe_mean_upper, Thorpe_mean_lower)

    (_, _, _, _, _, _) = plotting_everything_for_stat(ax4, δ_UN .^ 2 .* N_UN.^3, δ_σ.^ 2 .* (N_UN[1])^3, Dissip_stats, upper_confint_dissiptime, lower_confint_dissiptime)


    Legend( ga[1, 1:2],  [vnp, vump, vsp1, vsbp, medn, errsp], ["Vary N₀", "Vary V₀", "Vary γ", "Subcritical", 
        "Median", "95% Confidence"],
        tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
        margin = (10, 10, 50, -5), framevisible = false, patchlabelgap = 7,
        labelsize = 30, 
        halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.05))      

save(apath * "Delta_v_All_MeansMediansError.png", f) 