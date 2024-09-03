using Statistics
using Printf
using JLD2
using CairoMakie
using GeometryBasics

dpath = "ReviewerAnalysis/ReviewerData/"
apath = "ReviewerAnalysis/ReviewerPlots/"
dpath2 = "Data/"

#filestatsname = dpath * "DeltavAll_ExtraStats.jld2"
#fileerrorname = dpath * "Delta_v_all_Confint.jld2"
fileerrorname = dpath * "Delta_v_tracstrat_Confint.jld2" # has medians, mean, and error bars
filesetnames =  "SetList_mp.jld2"

#filescalename_UN = dpath2 * "DeltavAllScale_mp.jld2"
#filescalename_σ = dpath2 * "DeltavAllScale_mp_VarySigma.jld2"
#filescalename_U = dpath2 * "DeltavAll_mtest.jld2"

# getting the surrounding simulation params out
file_sn = jldopen(filesetnames, "r+")

δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
δ_V0B = file_sn["δ_varyV0B"]

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
N_σ = file_sn["N_varyσ"]
γ_σ = file_sn["γ_varyσ"]
N_V0B = file_sn["N_varyV0B"]

# determining which data sets are where:
idx_subcritical_σ = findall(γ_σ .< 1)
idx_supcritical_σ = findall((γ_σ .> 1) .& ( γ_σ .!= 1.9))
idx_gammachange_σ = findall(γ_σ .!= 1.9)
idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10
idx_VaryU_n05 = (idx_VaryU_start+1):idx_VaryU_start+10
idx_VaryN_n05 = (idx_VaryN_start+1):idx_VaryN_start+10
idx_VaryV0B = 26 .+ (1:3)

# pulling data 
#file_scale_UN = jldopen(filescalename_UN, "r+")
#file_scale_σ = jldopen(filescalename_σ, "r+")
#file_scale_U2 = jldopen(filescalename_U, "r+")
#file_stat = jldopen(filestatsname, "r+")
file_error = jldopen(fileerrorname, "r+")

#Cheight_havg_tavg_UN =file_scale_UN["Cheight_havg_tavg"]
#Nheight_havg_tavg_UN =file_scale_UN["Nheight_havg_tavg"]
#Cheight_havg_tavg_σ =file_scale_σ["Cheight_havg_tavg"]
#Nheight_havg_tavg_σ =file_scale_σ["Nheight_havg_tavg"]

Nheight_mean = file_error["mean_stratanom"]
Nheight_median = file_error["median_stratanom"]

Cheight_mean = file_error["mean_tracer"]
Cheight_median = file_error["median_tracer"]

lower_confint_tracer = file_error["lower_confint_tracer"]
upper_confint_tracer = file_error["upper_confint_tracer"]
lower_confint_stratanom = file_error["lower_confint_stratanom"]
upper_confint_stratanom = file_error["upper_confint_stratanom"]

# fixing error bars to be around the correct mean:

f = Figure(resolution = (1500, 800), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1], xlabel = rich("h", subscript("w"), " [m]"), ylabel =  rich("L", subscript("tr")," [m]"), 
        xlabelsize=30, ylabelsize=30)
        ax1.xticks = 0:40:200
        ax1.yticks = 0:40:200
        limits!(ax1, 0,165,0,210) 

    ax2 = Axis(ga[2, 2], xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("N", superscript("2")), " [m]"), 
        xlabelsize=30, ylabelsize=30)
        ax2.xticks = 0:40:200
        ax2.yticks = 0:40:200
        limits!(ax2, 0,165,0,210) 

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    elem_3 = [LineElement(linewidth = 6, color =  (:firebrick2, 0.5),
        points = Point2f[(0, 0), (0, 1)]), LineElement(linewidth = 6, color =  (:dodgerblue2, 0.5),
        points = Point2f[(0.5,0), (0.5, 1)]), LineElement(linewidth = 6, color =  (:darkgreen, 0.5),
        points = Point2f[(1, 0), (1, 1)])]

    # error bars
    errsp = rangebars!(ax1, δ_σ[idx_gammachange_σ],  
        lower_confint_tracer[22 .+ idx_gammachange_σ], upper_confint_tracer[22 .+ idx_gammachange_σ]; linewidth = 6, color =  (:darkgreen, 0.5))
    rangebars!(ax1, δ_V0B, 
        lower_confint_tracer[idx_VaryV0B], upper_confint_tracer[idx_VaryV0B]; linewidth = 6, color =  (:dodgerblue2, 0.5))
    rangebars!(ax1, δ_UN[idx_VaryN], 
        lower_confint_tracer[idx_VaryN], upper_confint_tracer[idx_VaryN]; linewidth = 6, color =  (:firebrick2, 0.5))
    rangebars!(ax1, δ_UN[idx_VaryU], 
        lower_confint_tracer[idx_VaryU], upper_confint_tracer[idx_VaryU]; linewidth = 6, color =  (:dodgerblue2, 0.5))

    # median
    scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_median[22 .+ idx_gammachange_σ],  markersize = 25, marker=:star4, 
        color =(:gray30, 0.7),  strokewidth = 1, strokecolor = (:black, 0.7))
    scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_median[22 .+ idx_subcritical_σ], markersize = 25, 
        marker=Polygon(p_big, [p_small]), color= (:gray30, 0.7))
    scatter!(ax1, δ_V0B, Cheight_median[idx_VaryV0B], markersize = 25, marker=:dtriangle, 
        color = (:white, 0.7),  strokewidth = 1, strokecolor = color = (:gray30, 0.7))
    medn = scatter!(ax1, δ_UN[idx_VaryN], Cheight_median[idx_VaryN], markersize = 25, marker = :circle, 
        color =(:gray30, 0.7),  strokewidth = 1, strokecolor = (:black, 0.7))
    scatter!(ax1, δ_UN[idx_VaryU], Cheight_median[idx_VaryU], markersize = 25, marker=:utriangle, 
        color = (:gray30, 0.7),  strokewidth = 1, strokecolor = (:black, 0.7))

        # mean
    vuexd = scatter!(ax1, δ_V0B, Cheight_mean[idx_VaryV0B], markersize = 25, marker=:dtriangle, 
        color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], Cheight_mean[22 .+ idx_gammachange_σ],  markersize = 25, marker=:star4, 
        color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], Cheight_mean[22 .+ idx_subcritical_σ], markersize = 25, 
        marker=Polygon(p_big, [p_small]), color= :darkgoldenrod)
    vnp = scatter!(ax1, δ_UN[idx_VaryN], Cheight_mean[idx_VaryN], markersize = 25, marker = :circle, 
        color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump = scatter!(ax1, δ_UN[idx_VaryU], Cheight_mean[idx_VaryU], markersize = 25, marker=:utriangle, 
        color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    #### STARTIFICATION

        # error bars
    rangebars!(ax2, δ_σ[idx_supcritical_σ],  
        lower_confint_stratanom[22 .+ idx_supcritical_σ], upper_confint_stratanom[22 .+ idx_supcritical_σ]; linewidth = 6, color = (:darkgreen, 0.5))
    rangebars!(ax2, δ_V0B, 
        lower_confint_stratanom[idx_VaryV0B], upper_confint_stratanom[idx_VaryV0B]; linewidth = 6, color = (:dodgerblue2, 0.5))
    rangebars!(ax2, δ_UN[idx_VaryN_n05], 
        lower_confint_stratanom[idx_VaryN_n05], upper_confint_stratanom[idx_VaryN_n05]; linewidth = 6, color = (:firebrick2, 0.5))
    rangebars!(ax2, δ_UN[idx_VaryU_n05], 
        lower_confint_stratanom[idx_VaryU_n05], upper_confint_stratanom[idx_VaryU_n05]; linewidth = 6, color = (:dodgerblue2, 0.5))

        # median
    scatter!(ax2, δ_σ[idx_supcritical_σ], Nheight_median[22 .+ idx_supcritical_σ],  markersize = 25, marker=:star4, 
        color =(:gray30, 0.7),  strokewidth = 1, strokecolor = (:black, 0.7))
    scatter!(ax2, δ_V0B, Nheight_median[idx_VaryV0B], markersize = 25, marker=:dtriangle, 
        color = (:white, 0.7),  strokewidth = 3, strokecolor = (:gray30, 0.7))
    scatter!(ax2, δ_UN[idx_VaryN_n05], Nheight_median[idx_VaryN_n05], markersize = 25, marker = :circle, 
        color =(:gray30, 0.7),  strokewidth = 1, strokecolor = (:black, 0.7))
    scatter!(ax2, δ_UN[idx_VaryU_n05], Nheight_median[idx_VaryU_n05], markersize = 25, marker=:utriangle, 
        color = (:gray30, 0.7),  strokewidth = 1, strokecolor = (:black, 0.7))
        
        # mean
    scatter!(ax2, δ_V0B, Nheight_mean[idx_VaryV0B], markersize = 25, marker=:dtriangle, 
               color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    scatter!(ax2, δ_σ[idx_supcritical_σ], Nheight_mean[22 .+ idx_supcritical_σ],  markersize = 25, marker=:star4, 
                color =:darkgreen, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_UN[idx_VaryN_n05], Nheight_mean[idx_VaryN_n05], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    scatter!(ax2, δ_UN[idx_VaryU_n05], Nheight_mean[idx_VaryU_n05], markersize = 25, marker=:utriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)


    Legend( ga[1, 1:2],  [vnp, vump, vuexd, vsp1, vsbp, medn, elem_3], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ", "Vary γ", "Subcritical*", 
        "Median", "95% Confidence"],
        tellheight = false, tellwidth = false, rowgap = 30, colgap = 30,
        margin = (10, 10, 50, -5), framevisible = false, patchlabelgap = 7,
        labelsize = 30, 
        halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.05))   
    colgap!(ga, 40)

    Label(ga[2, 1, TopLeft()], "a",
        fontsize = 30,
        font = :bold,
        padding = (-10, 5, 5, 10),
        halign = :right)
    Label(ga[2, 2, TopLeft()], "b",
        fontsize = 30,
        font = :bold,
        padding = (0, 5, 5, 10),
        halign = :right)

save(apath * "Paper_Intrusion_v_Delta_v_All_wMediansError2.png", f, px_per_unit = 2) 