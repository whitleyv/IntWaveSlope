using Statistics
using Printf
using Measures
using JLD2
using CairoMakie

path_name = "Data/"

sn1 = "U150N100Lz100g100" 
sn2 = "U250N100Lz100g100" 
sn3 = "U300N100Lz100g100" 

sn = sn1

include("../parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

        
filescalename1 = path_name * "BFluxTrip_rAvg_prof_mp_" * sn1  * ".jld2"
filescalename2 = path_name * "BFluxTrip_rAvg_prof_mp_" * sn2  * ".jld2"
filescalename3 = path_name * "BFluxTrip_rAvg_prof_mp_" * sn3  * ".jld2"

scale_file1 = jldopen(filescalename1, "r+")
scale_file2 = jldopen(filescalename2, "r+")
scale_file3 = jldopen(filescalename3, "r+")

skeys = keys(scale_file1)

∇_rWavg_prof150 = -1 .* (scale_file1[skeys[4]] .+ scale_file1[skeys[5]]) ; 
∇_rWavg_prof250 = -1 .* (scale_file2[skeys[4]] .+ scale_file2[skeys[5]]) ; 
∇_rWavg_prof300 = -1 .* (scale_file3[skeys[4]] .+ scale_file3[skeys[5]]) ; 

c_rWavg_prof150 = scale_file1[skeys[2]];
c_rWavg_prof250 = scale_file2[skeys[2]];
c_rWavg_prof300 = scale_file3[skeys[2]];

zL = 75
HAB = 2:2:zL*2
dHAB = HAB[2:end]

# center of mass of profile:
c_rWavg_prof_CoM150 = sum(c_rWavg_prof150 .* HAB, dims = 1)[ 1, :] ./ sum(c_rWavg_prof150, dims = 1)[ 1, :]
c_rWavg_prof_CoM250 = sum(c_rWavg_prof250 .* HAB, dims = 1)[ 1, :] ./ sum(c_rWavg_prof250, dims = 1)[ 1, :]
c_rWavg_prof_CoM300 = sum(c_rWavg_prof300 .* HAB, dims = 1)[ 1, :] ./ sum(c_rWavg_prof300, dims = 1)[ 1, :]

rolling_phase_times1 = scale_file1[skeys[12]];
rolling_phase_times2 = scale_file1[skeys[13]];

fsc3 = path_name * "BFluxTripb_rAvg_prof_mp.jld2"
b_file = jldopen(fsc3, "r+")

b_rWavg_prof150 = b_file["b_rWavg_prof150"];
b_rWavg_prof250 = b_file["b_rWavg_prof250"];
b_rWavg_prof300 = b_file["b_rWavg_prof300"];

Δb150 = b_rWavg_prof150 .- b_rWavg_prof150[:,1]
Δb250 = b_rWavg_prof250 .- b_rWavg_prof250[:,1]
Δb300 = b_rWavg_prof300 .- b_rWavg_prof300[:,1]


# v            c               
# wbdz        vbdy          -nablaph
# wbdz        vbdy          -nablat
# SGS   -nabturb-nabph       db/dt

nabmax = round(maximum(abs.(∇_rWavg_prof250))*1e7*0.65)*1e-7
dbmax = round(maximum(abs.(Δb250))*1e4)*1e-4

nablims = (-nabmax, nabmax)
dblims = (-dbmax, dbmax)
blims = (-3.2e-3, -1.2e-3)
delta150 = 0.15/pm.Ñ
delta250 = 0.25/pm.Ñ
delta300 = 0.3/pm.Ñ

f1 = Figure(resolution = (2000, 1500), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[1, 1], ylabel = "z' [m]") #vh
    ax2 = Axis(ga[1, 2], ) 
    ax3 = Axis(ga[1, 3], )

    ax1c = Axis(ga[2, 1], ylabel = "z' [m]") #vh
    ax2c = Axis(ga[2, 2], ) 
    ax3c = Axis(ga[2, 3], )

    ax1b = Axis(ga[3, 1], ylabel = "z' [m]", xlabel = "Tσ") #vh
    ax2b = Axis(ga[3, 2], xlabel = "Tσ" ) 
    ax3b = Axis(ga[3, 3], xlabel = "Tσ") 

    gcb1 = ga[1, 4] = GridLayout()
    gcb2 = ga[2, 4] = GridLayout()
    gcb3 = ga[3, 4] = GridLayout()       

    ax1.yticks = [50, 100]
    ax1b.yticks = [50, 100]

    ax1b.xticks = 2:2:10
    ax2b.xticks = 2:2:10
    ax3b.xticks = 2:2:10

    limits!(ax1, 1,10, 0, 150)
    limits!(ax2, 1,10, 0, 150)
    limits!(ax3, 1,10, 0, 150)
    limits!(ax1b, 1,10, 0, 150)
    limits!(ax2b, 1,10, 0, 150)
    limits!(ax3b, 1,10, 0, 150)

    hidedecorations!(ax2)
    hidedecorations!(ax3)
    hidexdecorations!(ax1)
    hideydecorations!(ax2b)
    hideydecorations!(ax3b)

    hb = heatmap!(ax1, rolling_phase_times1, HAB, b_rWavg_prof150', colormap = :thermal, colorrange = blims)
    contour!(ax1, rolling_phase_times1, HAB, b_rWavg_prof150', color = :white, linewidth = 2, levels = -0.003:0.0005:-0.001, alpha = 0.5)
    text!(ax1, Point.(1.5, 130), text = "b̄", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax1, 1:0.5:4, delta150.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax1, rolling_phase_times1, c_rWavg_prof_CoM150, color = :black, linewidth = 4)

    heatmap!(ax2, rolling_phase_times1, HAB, b_rWavg_prof250', colormap = :thermal, colorrange = blims)
    contour!(ax2, rolling_phase_times1, HAB, b_rWavg_prof250', color = :white, linewidth = 2, levels = -0.003:0.0005:-0.001, alpha = 0.5)
    text!(ax2, Point.(1.5, 130), text = "b̄", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax2, 1:0.5:4, delta250.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax2, rolling_phase_times1, c_rWavg_prof_CoM250, color = :black, linewidth = 4)

    heatmap!(ax3, rolling_phase_times1, HAB, b_rWavg_prof300', colormap = :thermal, colorrange = blims)
    contour!(ax3, rolling_phase_times1, HAB, b_rWavg_prof300', color = :white, linewidth = 2, levels = -0.003:0.0005:-0.001, alpha = 0.5)
    text!(ax3, Point.(1.5, 130), text = "b̄", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax3, 1:0.5:4, delta300.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax3, rolling_phase_times1, c_rWavg_prof_CoM300, color = :black, linewidth = 4)

    hdb = heatmap!(ax1c, rolling_phase_times1, HAB, Δb150', colormap = :balance, colorrange = dblims)
    text!(ax1c, Point.(1.5, 130), text = "b̄ - b̄₀", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax1c, 1:0.5:4, delta150.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax1c, rolling_phase_times1, c_rWavg_prof_CoM150, color = :black, linewidth = 4)

    heatmap!(ax2c, rolling_phase_times1, HAB, Δb250', colormap = :balance, colorrange = dblims)
    text!(ax2c, Point.(1.5, 130), text = "b̄ - b̄₀", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax2c, 1:0.5:4, delta250.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax2c, rolling_phase_times1, c_rWavg_prof_CoM250, color = :black, linewidth = 4)

    heatmap!(ax3c, rolling_phase_times1, HAB, Δb300', colormap = :balance, colorrange = dblims)
    text!(ax3c, Point.(1.5, 130), text = "b̄ - b̄₀", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax3c, 1:0.5:4, delta300.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax3c, rolling_phase_times1, c_rWavg_prof_CoM300, color = :black, linewidth = 4)

    hnabp = heatmap!(ax1b, rolling_phase_times2, dHAB, ∇_rWavg_prof150', colormap = :balance, colorrange = nablims)
    text!(ax1b, Point.(1.5, 130), text = "-∇⋅(u⃗̃b̃)-∇⋅(u⃗'b')", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax1b, 1:0.5:4, delta150.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax1b, rolling_phase_times1, c_rWavg_prof_CoM150, color = :black, linewidth = 4)

    heatmap!(ax2b, rolling_phase_times2, dHAB, ∇_rWavg_prof250', colormap = :balance, colorrange = nablims)
    text!(ax2b, Point.(1.5, 130), text = "-∇⋅(u⃗̃b̃)-∇⋅(u⃗'b')", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax2b, 1:0.5:4, delta250.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax2b, rolling_phase_times1, c_rWavg_prof_CoM250, color = :black, linewidth = 4)

    heatmap!(ax3b, rolling_phase_times2, dHAB, ∇_rWavg_prof300', colormap = :balance, colorrange = nablims)
    text!(ax3b, Point.(1.5, 130), text = "-∇⋅(u⃗̃b̃)-∇⋅(u⃗'b')", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(ax3b, 1:0.5:4, delta300.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax3b, rolling_phase_times1, c_rWavg_prof_CoM300, color = :black, linewidth = 4)

    cb2 = Colorbar(gcb1[1,1], hb, ticks = (-3e-3:5e-4:1e-3), size =35, label = "Buoyancy")
    cb1 = Colorbar(gcb2[1,1], hdb, ticks = (-dbmax/2:dbmax/2:dbmax/2), size =35, label = "Change in Buoyancy")
    cb6 = Colorbar(gcb3[1,2], hnabp, ticks = (-nabmax/2:nabmax/2:nabmax/2), size =35, label = "∂b/∂t")

    rd1 = round(delta150)
    rd2 = round(10*delta250)* 0.1
    rd3 = round(10*delta300)* 0.1
    
    Label(ga[1, 1, Top()], "δ = $rd1 m",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)
    Label(ga[1, 2, Top()], "δ = $rd2 m",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)    
    Label(ga[1, 3, Top()], "δ = $rd3 m",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)
display(f1)
savename = "Analysis/Plots/NewerAnalysisFluxesDye/SlopeNormalHovmoller_deltacomp"

save(savename * ".png", f1)


f1 = Figure(resolution = (2000, 800), fontsize=35)
    ga = f1[1, 1] = GridLayout()    

    ax1b = Axis(ga[1, 1], ylabel = "z' [m]", xlabel = "Tσ") #vh
    ax2b = Axis(ga[1, 2], xlabel = "Tσ" ) 
    ax3b = Axis(ga[1, 3], xlabel = "Tσ") 

    gcb1 = ga[1, 4] = GridLayout()

    ax1b.yticks = [50, 100]

    ax1b.xticks = 2:2:10
    ax2b.xticks = 2:2:10
    ax3b.xticks = 2:2:10

    limits!(ax1b, 2,9, 2, 150)
    limits!(ax2b, 2,9, 2, 150)
    limits!(ax3b, 2,9, 2, 150)

    hideydecorations!(ax2b)
    hideydecorations!(ax3b)

    hnabp = heatmap!(ax1b, rolling_phase_times2, dHAB, ∇_rWavg_prof150', colormap = :balance, colorrange = nablims)
    lines!(ax1b, 1:0.5:10, delta150.*ones(length( 1:0.5:10)), color = :black, linewidth = 2, linestyle = :dash)
    #lines!(ax1b, rolling_phase_times1, c_rWavg_prof_CoM150, color = :black, linewidth = 4)
    contour!(ax1b, rolling_phase_times1, HAB, log10.(c_rWavg_prof150'), color = :black, linewidth = 2, levels = -4:1:-1)
    #contour!(ax1b, rolling_phase_times1, HAB, c_rWavg_prof150', color = :green, linewidth = 2, levels = 0:.03:.1)
    #contour!(ax1b, rolling_phase_times1, HAB, c_rWavg_prof150', color = :green, linewidth = 2, levels = 0:.003:.01)
    #contour!(ax1b, rolling_phase_times1, HAB, c_rWavg_prof150', color = :green, linewidth = 2, levels = 0:.0003:.001)
    text!(ax1b, Point.(2.9, 61), text = "10⁻²", align = (:right, :center), color = :black, font = :bold, fontsize = 26)
    text!(ax1b, Point.(2.9, 86), text = "10⁻³", align = (:right, :center), color = :black, font = :bold, fontsize = 26)
    text!(ax1b, Point.(2.9, 108), text = "10⁻⁴", align = (:right, :center), color = :black, font = :bold, fontsize = 26)


    heatmap!(ax2b, rolling_phase_times2, dHAB, ∇_rWavg_prof250', colormap = :balance, colorrange = nablims)
    lines!(ax2b, 1:0.5:10, delta250.*ones(length( 1:0.5:10)), color = :black, linewidth = 2, linestyle = :dash)
    #lines!(ax2b, rolling_phase_times1, c_rWavg_prof_CoM250, color = :black, linewidth = 4)
    contour!(ax2b, rolling_phase_times1, HAB, log10.(c_rWavg_prof250'), color = :black, linewidth = 2, levels = -4:1:-1)
    #contour!(ax2b, rolling_phase_times1, HAB, c_rWavg_prof250', color = :green, linewidth = 2, levels = 0:.03:.1)
    #contour!(ax2b, rolling_phase_times1, HAB, c_rWavg_prof250', color = :green, linewidth = 2, levels = 0:.003:.01)
    #contour!(ax2b, rolling_phase_times1, HAB, c_rWavg_prof250', color = :green, linewidth = 2, levels = 0:.0003:.001)
    text!(ax2b, Point.(2.9, 54), text = "10⁻²", align = (:right, :center), color = :black, font = :bold, fontsize = 26)
    text!(ax2b, Point.(3.93, 135), text = "10⁻³", align = (:right, :center), color = :black, font = :bold, fontsize = 26)

    heatmap!(ax3b, rolling_phase_times2, dHAB, ∇_rWavg_prof300', colormap = :balance, colorrange = nablims)
    lines!(ax3b, 1:0.5:10, delta300.*ones(length( 1:0.5:10)), color = :black, linewidth = 2, linestyle = :dash)
    #lines!(ax3b, rolling_phase_times1, c_rWavg_prof_CoM300, color = :black, linewidth = 4)
    contour!(ax3b, rolling_phase_times1, HAB, log10.(c_rWavg_prof300'), color = :black, linewidth = 2, levels = -4:1:-1)
    #contour!(ax3b, rolling_phase_times1, HAB, c_rWavg_prof300', color = :green, linewidth = 2, levels = 0:.03:.1)
    #contour!(ax3b, rolling_phase_times1, HAB, c_rWavg_prof300', color = :green, linewidth = 2, levels = 0:.003:.01)
    #contour!(ax3b, rolling_phase_times1, HAB, c_rWavg_prof300', color = :green, linewidth = 2, levels = 0:.0003:.001)
    text!(ax3b, Point.(3.95, 35), text = "10⁻²", align = (:right, :center), color = :black, font = :bold, fontsize = 26)
    text!(ax3b, Point.(3.15, 130), text = "10⁻³", align = (:right, :center), color = :black, font = :bold, fontsize = 26)

    cb6 = Colorbar(gcb1[1,1], hnabp, ticks = (-nabmax:nabmax/2:nabmax), size =35, label = "-∇⋅(u⃗̃b̃)-∇⋅(u⃗'b')")

    rd1 = round(delta150)
    rd2 = round(10*delta250)* 0.1
    rd3 = round(10*delta300)* 0.1
    
    Label(ga[1, 1, Top()], "δ = $rd1 m",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)
    Label(ga[1, 2, Top()], "δ = $rd2 m",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)    
    Label(ga[1, 3, Top()], "δ = $rd3 m",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)
display(f1)
savename = "Analysis/Plots/NewerAnalysisFluxesDye/SlopeNormalHovmoller_deltacomp_dbdtonly"

save(savename * ".png", f1)

