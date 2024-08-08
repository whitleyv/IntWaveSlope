using JLD2
using Statistics
using CairoMakie
using Printf

dpath = "ReviewerAnalysis/ReviewerData/"
apath = "ReviewerAnalysis/"

filescalename_UN = dpath * "DeltavDye_ExtraStats.jld2"
filesetnames =  "SetList_mp.jld2"

file_sn = jldopen(filesetnames, "r+")
file_scale_UN = jldopen(filescalename_UN, "r+")

Cheight_htmean = file_scale_UN["Cheight_htmean"]
Cheight_htmedian = file_scale_UN["Cheight_htmedian"]
Cheight_htmiddle = file_scale_UN["Cheight_htmiddle"]
Cheight_htmin = file_scale_UN["Cheight_htmin"]
Cheight_htmax = file_scale_UN["Cheight_htmax"]
Cheight_htnum = file_scale_UN["Cheight_htnum"]
Cheight_htstd = file_scale_UN["Cheight_htstd"]
Cheight_ht90p = file_scale_UN["Cheight_ht90p"]
Cheight_ht10p = file_scale_UN["Cheight_ht10p"]
Cheight_ht80p = file_scale_UN["Cheight_ht80p"]
Cheight_ht20p = file_scale_UN["Cheight_ht20p"]

low_err10 = Cheight_htmean .- Cheight_ht10p
low_err20 = Cheight_htmean .- Cheight_ht20p
high_err80 = Cheight_ht80p .- Cheight_htmean
high_err90 = Cheight_ht90p .- Cheight_htmean

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
δ_UN = file_sn["δs"]

idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10

f = Figure(resolution = (1200, 1200), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1], ylabel = rich("L", subscript("tr")," [m]"), 
    xlabelsize=30, ylabelsize=30)
    ax1.xticks = 0:40:200
    ax1.yticks = 0:40:200
    limits!(ax1, 0, 200, 0, 200) 

    ax2 = Axis(ga[2, 2], 
    xlabelsize=30, ylabelsize=30)
    ax2.xticks = 0:40:200
    ax2.yticks = 0:40:200
    limits!(ax2, 0, 200, 0, 200) 

    ax3 = Axis(ga[3, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("tr")," [m]"), 
    xlabelsize=30, ylabelsize=30)
    ax3.xticks = 0:40:300
    ax3.yticks = 0:40:300
    #limits!(ax3, 0, 300, 0, 200) 

    ax4 = Axis(ga[3, 2],  xlabel = rich("h", subscript("w"), " [m]"), 
    xlabelsize=30, ylabelsize=30)
    #ax4.xticks = 0:40:200
    #ax4.yticks = 0:40:200
    #limits!(ax4, 0, 200, 0, 200)

    hidexdecorations!(ax1, grid = false)
    hidedecorations!(ax2, grid = false)
    #hideydecorations!(ax4, grid = false)

    vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_htmean[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_htmean[idx_VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax2, δ_UN[idx_VaryN], Cheight_htmedian[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_UN[idx_VaryU], Cheight_htmedian[idx_VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax3, δ_UN[idx_VaryN], Cheight_htmiddle[idx_VaryN], markersize = 25, marker = :circle, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryU], Cheight_htmiddle[idx_VaryU], markersize = 25, marker=:utriangle, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryN], Cheight_htmin[idx_VaryN], markersize = 25, marker = :hline, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryU], Cheight_htmin[idx_VaryU], markersize = 25, marker=:hline, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        scatter!(ax3, δ_UN[idx_VaryN], Cheight_htmax[idx_VaryN], markersize = 25, marker = :hline, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryU], Cheight_htmax[idx_VaryU], markersize = 25, marker=:hline, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)


    scatter!(ax4, δ_UN[idx_VaryN], Cheight_htnum[idx_VaryN], markersize = 25, marker = :circle, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax4, δ_UN[idx_VaryU], Cheight_htnum[idx_VaryU], markersize = 25, marker=:utriangle, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend( ga[1, 1:2],  [vnp, vump], ["Vary N₀", "Vary V₀"],
        tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
        margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
        labelsize = 35,
        halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.05))       
    #colgap!(ga, 40)

    Label(ga[2, 1, Top()], "Mean",
                    fontsize = 30,
                    font = :bold,
                    padding = (-10, 5, 5, 10),
                    )
    Label(ga[2, 2, Top()], "Median",
                    fontsize = 30,
                    font = :bold,
                    padding = (0, 5, 5, 10),
                    )
    Label(ga[3, 1, Top()], "Extrema/ Middle",
                    fontsize = 30,
                    font = :bold,
                    padding = (-10, 5, 5, 10),
                    )
    Label(ga[3, 2, Top()], "Number of Obs",
                    fontsize = 30,
                    font = :bold,
                    padding = (0, 5, 5, 10),
                    )

savename = apath * "Delta_v_TracerThickness_extraStats"
save(savename * ".png", f) 


f = Figure(resolution = (1200, 1200), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1], ylabel = rich("L", subscript("tr")," [m]"), 
    xlabelsize=30, ylabelsize=30)
    ax1.xticks = 0:40:200
    ax1.yticks = 0:40:200
    limits!(ax1, 0, 200, 0, 200) 

    ax2 = Axis(ga[2, 2], 
    xlabelsize=30, ylabelsize=30)
    ax2.xticks = 0:40:200
    ax2.yticks = 0:40:200
    limits!(ax2, 0, 200, 0, 200) 

    ax3 = Axis(ga[3, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("tr")," [m]"), 
    xlabelsize=30, ylabelsize=30)
    ax3.xticks = 0:40:300
    ax3.yticks = 0:40:300
    limits!(ax3, 0, 200, 0, 200) 

    ax4 = Axis(ga[3, 2],  xlabel = rich("h", subscript("w"), " [m]"), 
    xlabelsize=30, ylabelsize=30)
    ax4.xticks = 0:40:200
    ax4.yticks = 0:40:200
    limits!(ax4, 0, 200, 0, 200)

    hidexdecorations!(ax1, grid = false)
    hidedecorations!(ax2, grid = false)
    hideydecorations!(ax4, grid = false)

    vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_htmean[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_htmean[idx_VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax2, δ_UN[idx_VaryN], Cheight_htmedian[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_UN[idx_VaryU], Cheight_htmedian[idx_VaryU], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    scatter!(ax3, δ_UN[idx_VaryN], Cheight_htmiddle[idx_VaryN], markersize = 25, marker = :circle, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryU], Cheight_htmiddle[idx_VaryU], markersize = 25, marker=:utriangle, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryN], Cheight_htmin[idx_VaryN], markersize = 25, marker = :hline, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryU], Cheight_htmin[idx_VaryU], markersize = 25, marker=:hline, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)
        scatter!(ax3, δ_UN[idx_VaryN], Cheight_htmax[idx_VaryN], markersize = 25, marker = :hline, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax3, δ_UN[idx_VaryU], Cheight_htmax[idx_VaryU], markersize = 25, marker=:hline, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)


    scatter!(ax4, δ_UN[idx_VaryN], Cheight_htstd[idx_VaryN], markersize = 25, marker = :circle, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax4, δ_UN[idx_VaryU], Cheight_htstd[idx_VaryU], markersize = 25, marker=:utriangle, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    Legend( ga[1, 1:2],  [vnp, vump], ["Vary N₀", "Vary V₀"],
        tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
        margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
        labelsize = 35,
        halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.05))       
    #colgap!(ga, 40)

    Label(ga[2, 1, Top()], "Mean",
                    fontsize = 30,
                    font = :bold,
                    padding = (-10, 5, 5, 10),
                    )
    Label(ga[2, 2, Top()], "Median",
                    fontsize = 30,
                    font = :bold,
                    padding = (0, 5, 5, 10),
                    )
    Label(ga[3, 1, Top()], "Extrema/ Middle",
                    fontsize = 30,
                    font = :bold,
                    padding = (-10, 5, 5, 10),
                    )
    Label(ga[3, 2, Top()], "Std Dev",
                    fontsize = 30,
                    font = :bold,
                    padding = (0, 5, 5, 10),
                    )

savename = apath * "Delta_v_TracerThickness_extraStats2"
save(savename * ".png", f) 

# mean + percentiles

f = Figure(resolution= (800, 900), fontsize = 26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1], ylabel = rich("L", subscript("tr")," [m]"), 
    xlabel = rich("h", subscript("w"), " [m]"),
    xlabelsize=30, ylabelsize=30)
    #ax1.xticks = 0:40:200
    #ax1.yticks = 0:40:200
    #limits!(ax1, 0, 200, 0, 200) 

    vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_htmean[idx_VaryN], markersize = 25, marker = :circle, 
    color =:firebrick2, strokewidth = 1, strokecolor = :black)

    vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_htmean[idx_VaryU], markersize = 25, marker=:utriangle, 
    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    errN = errorbars!(ax1, δ_UN[idx_VaryN], Cheight_htmean[idx_VaryN], 
        low_err10[idx_VaryN], high_err90[idx_VaryN]; linewidth = 3, color = :firebrick2)

    errorbars!(ax1, δ_UN[idx_VaryU], Cheight_htmean[idx_VaryU], 
        low_err10[idx_VaryU], high_err90[idx_VaryU]; linewidth = 2, color = :dodgerblue2, linestyel= :dash)

     Legend( ga[1, 1],  [vnp, vump, errN], ["Vary N₀", "Vary V₀", "90/10th percentile"],
        tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
        margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
        labelsize = 35,
        halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.05))    

    Label(ga[2, 1, Top()], "Mean",
                    fontsize = 30,
                    font = :bold,
                    padding = (-10, 5, 5, 10),
                    )
savename = apath * "Delta_v_TracerThickness_mean_percentiles.png"

save(savename * ".png", f) 
                    

f = Figure(resolution= (800, 900), fontsize = 26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1], ylabel = rich("L", subscript("tr")," [m]"), 
    xlabel = rich("h", subscript("w"), " [m]"),
    xlabelsize=30, ylabelsize=30)
    #ax1.xticks = 0:40:200
    #ax1.yticks = 0:40:200
    #limits!(ax1, 0, 200, 0, 200) 

    vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_htmean[idx_VaryN], markersize = 25, marker = :circle, 
    color =:firebrick2, strokewidth = 1, strokecolor = :black)

    vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_htmean[idx_VaryU], markersize = 25, marker=:utriangle, 
    color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    err80 = errorbars!(ax1, δ_UN[idx_VaryN], Cheight_htmean[idx_VaryN], 
        low_err20[idx_VaryN], high_err80[idx_VaryN]; linewidth = 3, color = :firebrick2)
    errorbars!(ax1, δ_UN[idx_VaryU], Cheight_htmean[idx_VaryU], 
        low_err20[idx_VaryU], high_err80[idx_VaryU]; linewidth = 2, color = :dodgerblue2)
        
     Legend( ga[1, 1],  [vnp, vump, err80], ["Vary N₀", "Vary V₀", "80/20th percentile"],
        tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
        margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7,
        labelsize = 35,
        halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.05))    

    Label(ga[2, 1, Top()], "Mean",
                    fontsize = 30,
                    font = :bold,
                    padding = (-10, 5, 5, 10),
                    )
savename = apath * "Delta_v_TracerThickness_mean_percentiles2.png"

save(savename * ".png", f) 
 
