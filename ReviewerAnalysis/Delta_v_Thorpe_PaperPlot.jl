using Statistics
using Printf
using CurveFit
using JLD2
using CairoMakie
using GeometryBasics
using CurveFit

apath = "ReviewerAnalysis/ReviewerPlots/"
dpath = "Data/"
dpath_review = "ReviewerAnalysis/ReviewerData/"
# setnames
filesetnames =  "SetList_mp.jld2"

# Thorpe Scale Values
filescalename_UN = dpath * "DeltavThorpeDissip.jld2"
filescalename_σ = dpath * "DeltavThorpeDissip_VarySigma.jld2"
filescalename_V0B = dpath_review * "Delta_v_tracstrat_Confint.jld2"

fileerrorname = dpath_review * "DeltavAll_ExtraStats.jld2"

fileerrorname_V0B = filescalename_V0B

file_sn = jldopen(filesetnames, "r+")

U_UN = file_sn["Us"]
N_UN = file_sn["Ns"]
N_σ = file_sn["N_varyσ"]

δ_UN = file_sn["δs"]
δ_σ = file_sn["δ_varyσ"]
δ_V0B = file_sn["δ_varyV0B"]

γ_σ = file_sn["γ_varyσ"]

idx_subcritical_σ = findall(γ_σ .< 1)
idx_topochange_σ = findall(γ_σ .== 1.9)
idx_gammachange_σ = findall(γ_σ .!= 1.9)

idx_VaryN_start = findfirst(N_UN .> 3.51e-3)
idx_VaryN = idx_VaryN_start:idx_VaryN_start+10
idx_VaryU_start = findfirst(U_UN .< 0.25)
idx_VaryU = idx_VaryU_start:idx_VaryU_start+10

file_Lt_UN = jldopen(filescalename_UN, "r+")
file_Lt_σ = jldopen(filescalename_σ, "r+")
file_Lt_V0B = jldopen(filescalename_V0B, "r+")

thorpe_rms_UN = file_Lt_UN["thorpe_rms"]
thorpe_rms_σ = file_Lt_σ["thorpe_rms"]
thorpe_rms_V0B = file_Lt_V0B["mean_thorpe"]

# error bars:
file_Lt_error = jldopen(fileerrorname,"r+")
Thorpe_stats = file_Lt_error["Thorpe_stats"]

# thorpe is exponential scaling:

# λ̄ = (Num - 2) / mean * Num = (Num - 2)/sum 
thorpe_mean_UNσ = vcat(thorpe_rms_UN, thorpe_rms_σ[1:end-4])
λ̄ = (Thorpe_stats[:, 5] .- 2) ./ (thorpe_mean_UNσ .* Thorpe_stats[:, 5])
degs_freedom = 120
λ_upper =  λ̄ .* (1 - 1.96/√degs_freedom)
λ_lower =  λ̄ .* (1 + 1.96/√degs_freedom)

thorpe_upper_UNσ = 1 ./ λ_upper
thorpe_lower_UNσ = 1 ./ λ_lower
thorpe_upper_V0B = file_Lt_V0B["upper_confint_thorpe"]
thorpe_lower_V0B = file_Lt_V0B["lower_confint_thorpe"]


# thorpe scale plot rms vs max
f = Figure(resolution = (1050, 800), fontsize=26)
    #set_theme!(fonts = (; regular = "TeX Gyre Heros Makie Regular", bold = "TeX Gyre Heros Makie Bold", bold_italic = "TeX Gyre Heros Makie Bold Italic", italic = "TeX Gyre Heros Makie Italic"))

    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[2, 1],  xlabel = rich("h", subscript("w"), " [m]"), ylabel = rich("L", subscript("T"), " [m]"), xlabelsize=35, ylabelsize=35)
    # we want a log y axis with ticks at -8 through 0
    ax1.xticks = 0:40:160
    ax1.yticks =  0:20:60
    limits!(ax1, 0, 165, 0, 65)

    p_big = decompose(Point2f, Circle(Point2f(0), 0.6))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.5))

    elem_3 = [LineElement(linewidth = 6, color =  (:firebrick2, 0.5),
    points = Point2f[(0, 0), (0, 1)]), LineElement(linewidth = 6, color =  (:dodgerblue2, 0.5),
    points = Point2f[(0.5,0), (0.5, 1)]), LineElement(linewidth = 6, color =  (:darkgreen, 0.5),
    points = Point2f[(1, 0), (1, 1)])]

    # error bars
    errsp = rangebars!(ax1, δ_σ[idx_gammachange_σ],  
        thorpe_lower_UNσ[22 .+ idx_gammachange_σ], thorpe_upper_UNσ[22 .+ idx_gammachange_σ]; linewidth = 6, color =  (:darkgreen, 0.5))
    rangebars!(ax1, δ_V0B, 
        thorpe_lower_V0B, thorpe_upper_V0B; linewidth = 6, color =  (:dodgerblue2, 0.5))
    rangebars!(ax1, δ_UN[idx_VaryN], 
        thorpe_lower_UNσ[idx_VaryN], thorpe_upper_UNσ[idx_VaryN]; linewidth = 6, color =  (:firebrick2, 0.5))
    rangebars!(ax1, δ_UN[idx_VaryU], 
        thorpe_lower_UNσ[idx_VaryU], thorpe_upper_UNσ[idx_VaryU]; linewidth = 6, color =  (:dodgerblue2, 0.5))

    vsp1 = scatter!(ax1, δ_σ[idx_gammachange_σ], thorpe_rms_σ[idx_gammachange_σ], markersize = 25, marker=:star4, 
            color =:darkgreen, strokewidth = 1, strokecolor = :black)
    vsbp = scatter!(ax1, δ_σ[idx_subcritical_σ], thorpe_rms_σ[idx_subcritical_σ], markersize = 25, 
            marker=Polygon(p_big, [p_small]), color= :darkgoldenrod, )
    vump_max = scatter!(ax1, δ_V0B, thorpe_rms_V0B, markersize = 25, marker=:utriangle, 
        color = :white, strokewidth = 3, strokecolor = :dodgerblue2)
    vnp_rms = scatter!(ax1, δ_UN[idx_VaryN], thorpe_rms_UN[idx_VaryN], markersize = 25, marker = :circle, 
                color =:firebrick2, strokewidth = 1, strokecolor = :black)
    vump_rms = scatter!(ax1, δ_UN[idx_VaryU], thorpe_rms_UN[idx_VaryU], markersize = 25, marker=:dtriangle, 
                color = :dodgerblue2, strokewidth = 1, strokecolor = :black)

    #Legend(ga[1, 1], [vnp_rms, vump_rms, vump_max, vsp1, vsbp, elem_3], ["Vary N₀", "Vary V₀", "Vary V₀ᴮ", "Vary γ", "Subcritical γ", "95% Confidence"],
    Legend(ga[1, 1], [vnp_rms, vump_max, vump_rms, vsp1, elem_3, vsbp], ["Vary N₀", "Vary V₀ᴮ",  "Vary V₀", "Vary γ", "95% Confidence", "Subcritical γ"],
                    tellheight = false, tellwidth = false, labelsize = 35, rowgap = 15, colgap = 30,
                    margin = (10, 10, 10, 5), framevisible = false, patchlabelgap = 7, nbanks = 2,
                    halign = :center, valign = :top, orientation = :horizontal)

    rowsize!(ga,1, Auto(0.17))       
display(f)                  

savename = apath * "Paper_Thorpe_v_Delta_rms_witherror"
save(savename * ".png", f, px_per_unit = 2) 
