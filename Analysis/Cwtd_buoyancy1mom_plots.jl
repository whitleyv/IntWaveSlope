using Statistics
using Printf
using JLD2
using CairoMakie

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename = dpath * "cwtdb1mom.jld2"
filescalename2 = dpath * "cwtdb1mom_small_U250N100Lz100g100.jld2"
scale_file = jldopen(filescalename, "r+")
scale_file2 = jldopen(filescalename2, "r+")
c_weighted_bavg_mp = vcat(scale_file["c_weighted_bavg_mp"][1:4,:], scale_file["c_weighted_bavg_mp"][7:9,:])
c_weighted_bavg_mn = scale_file["c_weighted_bavg_mn"]
setnames = scale_file["setnames"]
setnames2 = scale_file["setnames2"]

wave_times = 0:1/15:10+10/15
######################
#                  b CENTROID PLOT FOR ALL GAUSSIAN
######################

f = Figure(resolution = (700, 1000), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "Δb̄ [ms⁻²]")
ax3 = Axis(ga[2, 1], ylabel = "Δb̄ [ms⁻²]")
ax2 = Axis(ga[3, 1], ylabel = "Δb̄ [ms⁻²]", xlabel = "Wave Periods [Tσ]")

limits!(ax1, 0, 11, -1.2e-4, 1.2e-4)
limits!(ax3, 0, 11, -1.2e-4, 1.2e-4)
limits!(ax2, 0, 11, -2.7e-4, 2.7e-4)

hidexdecorations!(ax1, grid = false)
hidexdecorations!(ax3, grid = false)

dbmp = lines!(ax1, wave_times, (c_weighted_bavg_mp[1,:] .- c_weighted_bavg_mp[1,1]), color = :dodgerblue1, label = "28 m", linewidth = 5)
    lines!(ax1, wave_times, (c_weighted_bavg_mp[2,:] .- c_weighted_bavg_mp[2,1]), color = :dodgerblue2, label = "42 m", linewidth = 5)
    lines!(ax1, wave_times, (c_weighted_bavg_mp[3,:] .- c_weighted_bavg_mp[3,1]), color = :dodgerblue3, label = "57 m", linewidth = 5)
    lines!(ax1, wave_times, (c_weighted_bavg_mp[4,:] .- c_weighted_bavg_mp[4,1]), color = :dodgerblue4, label ="71 m", linewidth = 5)

lines!(ax2, wave_times, (c_weighted_bavg_mp[5,:] .- c_weighted_bavg_mp[5,1]), color = :dodgerblue1, label ="85 m", linewidth = 5)
    lines!(ax2, wave_times, (c_weighted_bavg_mp[6,:] .- c_weighted_bavg_mp[6,1]), color = :dodgerblue3, label ="114 m", linewidth = 5)
    lines!(ax2, wave_times, (c_weighted_bavg_mp[7,:] .- c_weighted_bavg_mp[7,1]), color = :dodgerblue4, label ="142 m", linewidth = 5)

dbmn = lines!(ax3, wave_times, (c_weighted_bavg_mn[1,:] .- c_weighted_bavg_mn[1,1]), color = :firebrick1, label = "28 m", linewidth = 5)
    lines!(ax3, wave_times, (c_weighted_bavg_mn[2,:] .- c_weighted_bavg_mn[2,1]), color = :firebrick2, label = "42 m", linewidth = 5)
    lines!(ax3, wave_times, (c_weighted_bavg_mn[3,:] .- c_weighted_bavg_mn[3,1]), color = :firebrick3,  label = "57 m", linewidth = 5)
    lines!(ax3, wave_times, (c_weighted_bavg_mn[4,:] .- c_weighted_bavg_mn[4,1]), color = :firebrick4, label ="71 m", linewidth = 5)

axislegend(ax2, [dbmp, dbmn], ["m > 0", "m < 0"], position = :rt)
axislegend(ax1, position = :lb, orientation = :horizontal)
axislegend(ax2, position = :lt, orientation = :horizontal)
axislegend(ax3, position = :lt, orientation = :horizontal)

savename = "cwtd_bbar_all_split" 

save(apath * savename * ".png", f)

###################

f = Figure(resolution = (1000, 800), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "Δb̄ [ms⁻²]", xlabel = "Wave Periods [Tσ]")

limits!(ax1, 0, 11, -2e-4, 2e-4)
ax1.xticks = 0:2:10

dbmp = lines!(ax1, wave_times, (c_weighted_bavg_mp[2,:] .- c_weighted_bavg_mp[2,1]), color = :dodgerblue2, label = "δ = 42 m, m > 0", linewidth = 5)
dbmn = lines!(ax1, wave_times, (c_weighted_bavg_mn[2,:] .- c_weighted_bavg_mn[2,1]), color = :dodgerblue4, label = "δ = 42 m, m < 0", linewidth = 5)
dbmp2 = lines!(ax1, wave_times, (c_weighted_bavg_mp[5,:] .- c_weighted_bavg_mp[5,1]), color = :firebrick2, label = "δ = 85 m, m > 0", linewidth = 5)

axislegend(ax1, position = :lt)

savename = "cwtd_bbar_representative" 

save(apath * savename * ".png", f)

####################
b̄_Cg20 = scale_file["c_weighted_bavg_mp"][4,:]
b̄_Cg15 = scale_file2["b̄_Cg15"]
b̄_Cg05 = scale_file2["b̄_Cg05"]

f = Figure(resolution = (1000, 800), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "Δb̄ [ms⁻²]", xlabel = "Wave Periods [Tσ]")

limits!(ax1, 0, 11, -1.8e-4, 1.8e-4)
ax1.xticks = 0:2:10

dbmp = lines!(ax1, wave_times, (b̄_Cg20 .- b̄_Cg20[1]), color = :dodgerblue2, label = "Concentration of 10⁻⁴ at  1.10 δ", linewidth = 5)
dbmn = lines!(ax1, wave_times, (b̄_Cg15 .- b̄_Cg15[1]), color = :firebrick2, label = "Concentration of 10⁻⁴ at  0.84 δ", linewidth = 5)
dbmp2 = lines!(ax1, wave_times, (b̄_Cg05 .- b̄_Cg05[1]), color = :gray30, label = "Concentration of 10⁻⁴ at  0.28 δ", linewidth = 5)

axislegend(ax1, position = :lb)

savename = "cwtd_bbar_U250N100Lz100g100" 

save(apath * savename * ".png", f)

