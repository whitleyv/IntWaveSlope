using JLD2
using Statistics
using CairoMakie
using Printf

dpath = "ReviewerAnalysis/ReviewerData/"
apath = "ReviewerAnalysis/"

filescalename_UN = dpath * "DeltavDissipThorpe.jld2"
filesetnames =  "SetList_mp.jld2"

file_sn = jldopen(filesetnames, "r+")
file_scale_UN = jldopen(filescalename_UN, "r+")

Thorpe_rms = file_scale_UN["Thorpe_rms"]
Thorpe_stats = file_scale_UN["Thorpe_stats"]
Nheight_stats = file_scale_UN["Nheight_stats"]
Dissip_stats = file_scale_UN["Dissip_stats"]

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
δ_UN = file_sn["δs"]

idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10

# extra_stats = [tmean, tmedian, tmax, tmin, tnum, tstd, twentieth_percentile, eightieth_percentile]

# rich("L", subscript("tr")," [m]")
#  0:40:160
function plot_stats(stats, δ_UN, savename, idx_VaryN, idx_VaryU, stat_string, xtix, ytix, lims)

    (xmin, xmax, ymin, ymax) = lims
    low_err20 = stats[:,1] .- stats[:, 7]
    high_err80 = stats[:, 8] .- stats[:,1] 

    stddev_error = 2 .* stats[:, 6]  #./ sqrt.(stats[:, 5])

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

        # mean + 80/20 percentile + extrema
        vnp = scatter!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        vump = scatter!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        errN = errorbars!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 1], 
                    low_err20[idx_VaryN], high_err80[idx_VaryN]; linewidth = 3, color = :firebrick2)
        errV = errorbars!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 1], 
                    low_err20[idx_VaryU], high_err80[idx_VaryU]; linewidth = 2, color = :dodgerblue2)
        extrem = scatter!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 4], markersize = 25, marker = :hline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 4], markersize = 25, marker=:hline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, δ_UN[idx_VaryN], stats[idx_VaryN, 3], markersize = 25, marker = :hline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, δ_UN[idx_VaryU], stats[idx_VaryU, 3], markersize = 25, marker=:hline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # median
        scatter!(ax2, δ_UN[idx_VaryN], stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax2, δ_UN[idx_VaryU], stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        # std dev 
        scatter!(ax3, δ_UN[idx_VaryN], stats[idx_VaryN, 6], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax3, δ_UN[idx_VaryU], stats[idx_VaryU, 6], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # mean + regular error bars
        scatter!(ax4, δ_UN[idx_VaryN], stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax4, δ_UN[idx_VaryU], stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        errN2 = errorbars!(ax4, δ_UN[idx_VaryN], stats[idx_VaryN, 1], 
                stddev_error[idx_VaryN], stddev_error[idx_VaryN]; linewidth = 3, linestyle = :dash, color = :firebrick2)
        errV2 = errorbars!(ax4, δ_UN[idx_VaryU], stats[idx_VaryU, 1], 
                stddev_error[idx_VaryU], stddev_error[idx_VaryU]; linewidth = 2, linestyle = :dash, color = :dodgerblue2)

        Legend( ga[1, 1:2],  [vnp, vump, errN, errN2, extrem], ["Vary N₀", "Vary V₀", "20/80th percentile", 
            "2σ", "Extrema"],
            tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
            margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
            labelsize = 35,
            halign = :center, valign = :top, orientation = :horizontal)

        rowsize!(ga,1, Auto(0.05))       

        Label(ga[2, 1, Top()], "Mean + Percentiles",
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
        Label(ga[3, 2, Top()], "Mean + Std Dev Error",
                        fontsize = 30,
                        font = :bold,
                        padding = (0, 5, 5, 10),
                        )

    save(savename * ".png", f) 
end

function logplot_stats(stats, indep, savename, idx_VaryN, idx_VaryU, stat_string, indep_string, xtix, ytix, lims)

    (xmin, xmax, ymin, ymax) = lims
    low_err20 = abs.(stats[:,1] .- stats[:, 7])
    high_err80 = abs.(stats[:, 8] .- stats[:,1] )

    stddev_error = 2 .* stats[:, 6]  #./ sqrt.(stats[:, 5])

    f = Figure(resolution = (1200, 1200), fontsize=26)
        ga = f[1, 1] = GridLayout()

        ax1 = Axis(ga[2, 1], ylabel = stat_string, 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
        #ax1.xticks = xtix
       # ax1.yticks = ytix
        limits!(ax1, xmin, xmax, ymin, ymax) 

        ax2 = Axis(ga[2, 2], 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
       # ax2.xticks = xtix
       # ax2.yticks = ytix
        limits!(ax2, xmin, xmax, ymin, ymax) 

        ax3 = Axis(ga[3, 1],  xlabel = indep_string, ylabel = stat_string, 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
        #ax3.xticks = xtix
        #ax3.yticks = ytix
        limits!(ax3, xmin, xmax, ymin, ymax) 

        ax4 = Axis(ga[3, 2],  xlabel = indep_string, 
        xlabelsize=30, ylabelsize=30, yscale = log10, xscale = log10)
        #ax4.xticks = xtix
        #ax4.yticks = ytix
        limits!(ax4, xmin, xmax, ymin, ymax)

        hidexdecorations!(ax1, grid = false)
        hidedecorations!(ax2, grid = false)
        hideydecorations!(ax4, grid = false)

        # mean + 80/20 percentile + extrema
        vnp = scatter!(ax1, indep[idx_VaryN], stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        vump = scatter!(ax1, indep[idx_VaryU], stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        errN = errorbars!(ax1, indep[idx_VaryN], stats[idx_VaryN, 1], 
                    low_err20[idx_VaryN], high_err80[idx_VaryN]; linewidth = 3, color = :firebrick2)
        errV = errorbars!(ax1, indep[idx_VaryU], stats[idx_VaryU, 1], 
                    low_err20[idx_VaryU], high_err80[idx_VaryU]; linewidth = 2, color = :dodgerblue2)
        extrem = scatter!(ax1, indep[idx_VaryN], stats[idx_VaryN, 4], markersize = 25, marker = :hline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, indep[idx_VaryU], stats[idx_VaryU, 4], markersize = 25, marker=:hline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, indep[idx_VaryN], stats[idx_VaryN, 3], markersize = 25, marker = :hline, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax1, indep[idx_VaryU], stats[idx_VaryU, 3], markersize = 25, marker=:hline, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # median
        scatter!(ax2, indep[idx_VaryN], stats[idx_VaryN, 2], markersize = 25, marker = :circle, 
                    color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax2, indep[idx_VaryU], stats[idx_VaryU, 2], markersize = 25, marker=:utriangle, 
                    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        # std dev 
        scatter!(ax3, indep[idx_VaryN], stats[idx_VaryN, 6], markersize = 25, marker = :circle, 
            color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax3, indep[idx_VaryU], stats[idx_VaryU, 6], markersize = 25, marker=:utriangle, 
            color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

        # mean + regular error bars
        scatter!(ax4, indep[idx_VaryN], stats[idx_VaryN, 1], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
        scatter!(ax4, indep[idx_VaryU], stats[idx_VaryU, 1], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        errN2 = errorbars!(ax4, indep[idx_VaryN], stats[idx_VaryN, 1], 
            0, stddev_error[idx_VaryN]; linewidth = 3, linestyle = :dash, color = :firebrick2)
        errV2 = errorbars!(ax4, indep[idx_VaryU], stats[idx_VaryU, 1], 
            0, stddev_error[idx_VaryU]; linewidth = 2, linestyle = :dash, color = :dodgerblue2)

        Legend( ga[1, 1:2],  [vnp, vump, errN, errN2, extrem], ["Vary N₀", "Vary V₀", "20/80th percentile", 
             "2σ", "Extrema"],
            tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
            margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
            labelsize = 35,
            halign = :center, valign = :top, orientation = :horizontal)

        rowsize!(ga,1, Auto(0.05))       

        Label(ga[2, 1, Top()], "Mean + Percentiles",
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
        Label(ga[3, 2, Top()], "Mean + Std Dev Error",
                        fontsize = 30,
                        font = :bold,
                        padding = (0, 5, 5, 10),
                        )

    save(savename * ".png", f) 
end

plot_stats(Thorpe_stats, δ_UN, apath * "Delta_v_Thorpe_ExtraStats", idx_VaryN, idx_VaryU, "Lₜ", 0:40:200, 0:20:200, (0,160,0,80))

plot_stats(Nheight_stats, δ_UN, apath * "Delta_v_Nheight_ExtraStats", idx_VaryN, idx_VaryU, "Lₙ", 0:40:200, 0:40:200, (0,160,0,300))

# usually it's ϵ vs hw^2N^3, we can instead do √(ϵ/ N^3) = hw
logplot_stats(Dissip_stats, δ_UN .^ 2 .* N_UN.^3 , apath * "Delta_v_Dissip_ExtraStats", idx_VaryN, idx_VaryU, "ϵ", rich("h", subscript("w"), superscript("2"),"N", subscript("0"), superscript("3"), " [m²s", superscript("-3"), "]"), (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"]),  (0:1e-5:3e-5, ["0", "1×10⁻⁵", "2×10⁻⁵", "3×10⁻⁵"]), (1e-6, 10^(-2.75), 1e-6, 10^-(3.25)))

