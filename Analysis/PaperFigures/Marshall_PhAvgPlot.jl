using Statistics
using Printf
using JLD2
using CairoMakie

apath =  "Analysis/PaperFigures/"
dpath = "Data/"
dpath = "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/Analysis/"

include("../../parameters.jl")

setname = "U350N100Lz100g100"

pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
m = π/pm.Lz,
l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
Tf = 2*π/pm.f, 
Tσ = 2*π/pm.σ))

Phavg_filename = "VolumeMVelocityDissipN2_nb_25_110_Phavg_" * setname

Phf = jldopen(dpath * Phavg_filename * ".jld2")
skeys = keys(Phf)

phase_times  = Phf[skeys[1]];
b_bins_bot  = Phf[skeys[2]];
b_bins_bot_velocity  = Phf[skeys[3]];
∂tV_inb_Phavg  = Phf[skeys[4]];
∫M_inb_Phavg  = Phf[skeys[5]];
Adif_Phavg = Phf[skeys[6]];
vavg_inb_Phavg = Phf[skeys[7]];
log_eavg_inb_Phavg = Phf[skeys[8]];
N2avg_inb_Phavg = Phf[skeys[9]];
M2avg_inb_Phavg = Phf[skeys[11]];
∇bavg_inb_Phavg= Phf[skeys[10]];

ini_Nisos = 25
ini_Nisos_velocity = 50 # double resolution for velocity contours
volstart =  2
volend = ini_Nisos - 2 
ini_Nisos_range = volstart:volend
ini_Nisos_range_velocity = volstart*2:ini_Nisos_velocity-4

Δb = -500*pm.Ñ^2/ini_Nisos
Δb_velocity = -500*pm.Ñ^2/ini_Nisos_velocity

b_bins_bot = (Δb .* (1:ini_Nisos))[ini_Nisos_range]
b_bins_bot_velocity = (Δb_velocity .* (1:ini_Nisos_velocity))[ini_Nisos_range_velocity]
#(1500, 700)

f1 = Figure(resolution = (2000, 1400), fontsize=26)
f1 = Figure(size = (2000, 1400), fontsize=26)

    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[2, 2]) 
    ax3 = Axis(ga[2, 3]) 
    ax4 = Axis(ga[4, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax5 = Axis(ga[4, 2], xlabel = "Tσ" ) 
    ax6 = Axis(ga[4, 3], xlabel = "Tσ" ) 

    gcb1 = ga[1, 1:3] = GridLayout()
    gcb2 = ga[3, 1:3] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3
    ax4.yticks = -5e-3:1e-3:-1e-3

    ax4.xticks = .2:.2:1
    ax5.xticks = .2:.2:1
    ax6.xticks = .2:.2:1

    limits!(ax1, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax2, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax3, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax4, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax5, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax6, -0.03, 0.99, -5.2e-3, -0.8e-3)

    hidedecorations!(ax2)
    hidedecorations!(ax3)
    hidexdecorations!(ax1)
    hideydecorations!(ax5)
    hideydecorations!(ax6)

    scaling_levelsn = -.1:.1:.1
    scaling_levelsp = 0:.1:.1

    N2scaling_levelsn = -5e-6:5e-6:0
    N2scaling_levelsp = 0:1e-5:1e-5

    hV = heatmap!(ax1, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 25, labelfont = :bold)    
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 25, labelfont = :bold)    
    text!(ax1, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax1, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax1, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax1, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax1, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax1, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    heatmap!(ax2, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 25, labelfont = :bold)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 25, labelfont = :bold)    
    text!(ax2, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    hA = heatmap!(ax3, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-400,400))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 25, labelfont = :bold)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 25, labelfont = :bold)    
    text!(ax3, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    hv = heatmap!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg', colormap = :balance, colorrange = (-0.2,0.2))
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 25, labelfont = :bold)    
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 25, labelfont = :bold)    
    text!(ax4, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    hE = heatmap!(ax5, phase_times, b_bins_bot_velocity, log_eavg_inb_Phavg', colormap = :thermal, colorrange = (-7.5,-4.5))
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 25, labelfont = :bold)    
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 25, labelfont = :bold)    
    #contour!(ax5, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg'; levels = N2scaling_levelsp, color = :white , linewidth = 3)    
    #contour!(ax5, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg'; levels = N2scaling_levelsn, color = :white , linewidth = 3, linestyle = :dot)    
  text!(ax5, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    hN = heatmap!(ax6, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg', colormap = :balance, colorrange = (-2e-5,2e-5))
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 25, labelfont = :bold)    
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 25, labelfont = :bold)    
    #contour!(ax6, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg'; levels = N2scaling_levelsp, color = :white , linewidth = 3)    
    #contour!(ax6, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg'; levels = N2scaling_levelsn, color = :white , linewidth = 3, linestyle = :dot)    
  text!(ax6, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-500:250:500), size =35, label = "∂V(b,t)/∂t [m³s⁻¹]", vertical = false, labelsize = 30)
    cb2 = Colorbar(gcb1[1,2], hV, ticks = (-500:250:500), size =35, label = "- M(b,t) [m³s⁻¹]", vertical = false, labelsize = 30)
    cb3 = Colorbar(gcb1[1,3], hA, ticks = (-300:100:300), size =35, label = "A(b+Δb,t) - A(b,t) [m³s⁻¹]", vertical = false, labelsize = 30)
    cb4 = Colorbar(gcb2[1,1], hv, ticks = (-0.2:0.1:0.2), size =35, label = "v̅ [ms⁻¹]", vertical = false, labelsize = 30)
    cb4 = Colorbar(gcb2[1,2], hE, ticks = (-7:1:-5, ["10⁻⁷", "10⁻⁸", "10⁻⁹"]), size =35, label = "ε̅ [m³s⁻²]", vertical = false, labelsize = 30)
    cb4 = Colorbar(gcb2[1,3], hN, ticks = (-1e-5:1e-5:1e-5), size =35, label = "N²' [s⁻²]", vertical = false, labelsize = 30)

    rowsize!(ga, 1, Relative(0.05))
    rowsize!(ga, 3, Relative(0.05))

savename = apath * Phavg_filename * "_Hovmoller" 
display(f1)
save(savename * ".png", f1)


f2 = Figure(resolution = (2100, 1300), fontsize=26)
    gb = f2[1, 1] = GridLayout()    

    ax1 = Axis(gb[1, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(gb[1, 3]) 
    ax3 = Axis(gb[1, 5]) 
    ax4 = Axis(gb[2, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax5 = Axis(gb[2, 3], xlabel = "Tσ" ) 
    ax6 = Axis(gb[2, 5], xlabel = "Tσ" ) 

    gcb1 = gb[1:2, 2] = GridLayout()
    gcb2 = gb[1:2, 4] = GridLayout()
    gcb3 = gb[1:2, 6] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3
    ax4.yticks = -5e-3:1e-3:-1e-3

    ax4.xticks = .2:.2:1
    ax5.xticks = .2:.2:1
    ax6.xticks = .2:.2:1

    limits!(ax1, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax2, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax3, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax4, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax5, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax6, -0.03, 0.99, -5.2e-3, -0.8e-3)

    hidedecorations!(ax2)
    hidedecorations!(ax3)
    hidexdecorations!(ax1)
    hideydecorations!(ax5)
    hideydecorations!(ax6)

    scaling_levelsn = -.1:.1:.1
    scaling_levelsp = 0:.1:.1

    N2scaling_levelsn = -5e-6:5e-6:0
    N2scaling_levelsp = 0:1e-5:1e-5

    hV = heatmap!(ax1, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-650,650))
        contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
        contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
        text!(ax1, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
        text!(ax1, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
        text!(ax1, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
        text!(ax1, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
        text!(ax1, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
        text!(ax1, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    heatmap!(ax2, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    text!(ax2, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax2, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    heatmap!(ax3, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    text!(ax3, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax3, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    hv = heatmap!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg', colormap = :balance, colorrange = (-0.2,0.2))
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    text!(ax4, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax4, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    hE = heatmap!(ax5, phase_times, b_bins_bot_velocity, log_eavg_inb_Phavg', colormap = :thermal, colorrange = (-7.5,-4.5))
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    text!(ax5, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax5, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    hN = heatmap!(ax6, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg', colormap = :balance, colorrange = (-2e-5,2e-5))
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3)    
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dot)    
    text!(ax6, Point.(.34, -9.85e-4), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.6, -1e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.865, -1.7e-3), text = "0.00", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.19, -2.4e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.53, -3.05e-3), text = "0.10", align = (:left, :center), color = :black, fontsize = 25)
    text!(ax6, Point.(.89, -3.9e-3), text = "-0.10", align = (:left, :center), color = :black, fontsize = 25)

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-500:250:500), size =35)
    cb2 = Colorbar(gcb2[1,1], hV, ticks = (-500:250:500), size =35)
    cb3 = Colorbar(gcb3[1,1], hV, ticks = (-500:250:500), size =35)
    cb4 = Colorbar(gcb1[2,1], hv, ticks = (-0.2:0.1:0.2), size =35)
    cb4 = Colorbar(gcb2[2,1], hE, ticks = (-7:1:-5, ["10⁻⁷", "10⁻⁸", "10⁻⁹"]), size =35)
    cb4 = Colorbar(gcb3[2,1], hN, ticks = (-1e-5:1e-5:1e-5), size =35)

    Label(gb[1, 1, Top()], "∂V(b,t)/∂t [m³s⁻¹]", fontsize = 30, font = :bold, padding = (5, 5, 5, 5), halign = :center)
    Label(gb[1, 3, Top()], "- M(b,t) [m³s⁻¹]", fontsize = 30, font = :bold, padding = (5, 5, 5, 5), halign = :center)
    Label(gb[1, 5, Top()], "A(b+Δb,t) - A(b,t) [m³s⁻¹]", fontsize = 30, font = :bold, padding = (5, 5, 5, 5), halign = :center)
    Label(gb[2, 1, Top()], "v [ms⁻¹]", fontsize = 30, font = :bold, padding = (5, 5, 5, 5), halign = :center)
    Label(gb[2, 3, Top()], "ε [m³s⁻²]", fontsize = 30, font = :bold, padding = (5, 5, 5, 5), halign = :center)
    Label(gb[2, 5, Top()], "ΔN² [s⁻²]", fontsize = 30, font = :bold, padding = (5, 5, 5, 5), halign = :center)

    Label(ga[2, 1, TopLeft()], "a", fontsize = 30, font = :bold, padding = (5, 5, 5, 5), halign = :left)

    colsize!(gb, 2, Relative(0.05))
    colsize!(gb, 4, Relative(0.05))
    colsize!(gb, 6, Relative(0.05))
    colgap!(gb, 5)
    rowgap!(gb, 10)

savename = apath * Phavg_filename * "_Hovmollerflipped" 
display(f2)
save(savename * ".png", f2)


### using updated julia you can do contour labels for real!
f1 = Figure(size = (2000, 1400), fontsize=26)

    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[2, 2]) 
    ax3 = Axis(ga[2, 3]) 
    ax4 = Axis(ga[4, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax5 = Axis(ga[4, 2], xlabel = "Tσ" ) 
    ax6 = Axis(ga[4, 3], xlabel = "Tσ" ) 

    gcb1 = ga[1, 1:3] = GridLayout()
    gcb2 = ga[3, 1:3] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3
    ax4.yticks = -5e-3:1e-3:-1e-3

    ax4.xticks = .2:.2:1
    ax5.xticks = .2:.2:1
    ax6.xticks = .2:.2:1

    limits!(ax1, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax2, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax3, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax4, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax5, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax6, -0.03, 0.99, -5.2e-3, -0.8e-3)

    hidedecorations!(ax2)
    hidedecorations!(ax3)
    hidexdecorations!(ax1)
    hideydecorations!(ax5)
    hideydecorations!(ax6)

    scaling_levelsn = -.1:.5:0
    scaling_levelsp = 0:.1:.1

    N2scaling_levelsn = -5e-6:5e-6:-.1e-6
    N2scaling_levelsp = 0:1e-5:1e-5

    hV = heatmap!(ax1, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    heatmap!(ax2, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hA = heatmap!(ax3, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hv = heatmap!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg', colormap = :balance, colorrange = (-0.2,0.2))
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hE = heatmap!(ax5, phase_times, b_bins_bot_velocity, log_eavg_inb_Phavg', colormap = :thermal, colorrange = (-7.25,-4.25))
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hN = heatmap!(ax6, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg', colormap = :balance, colorrange = (-2e-5,2e-5))
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    cb1 = Colorbar(gcb1[1,1], hV, ticks = (-500:250:500), size =35, label = "∂V(b,t)/∂t [m³s⁻¹]", vertical = false, labelsize = 30)
    cb2 = Colorbar(gcb1[1,2], hV, ticks = (-500:250:500), size =35, label = "- M(b,t) [m³s⁻¹]", vertical = false, labelsize = 30)
    cb3 = Colorbar(gcb1[1,3], hA, ticks = (-500:250:500), size =35, label = "A(b+Δb,t) - A(b,t) [m³s⁻¹]", vertical = false, labelsize = 30)
    cb4 = Colorbar(gcb2[1,1], hv, ticks = (-0.2:0.1:0.2), size =35, label = "v [ms⁻¹]", vertical = false, labelsize = 30)
    cb4 = Colorbar(gcb2[1,2], hE, ticks = (-7:1:-5, ["10⁻⁷", "10⁻⁸", "10⁻⁹"]), size =35, label = "ε [m³s⁻²]", vertical = false, labelsize = 30)
    cb4 = Colorbar(gcb2[1,3], hN, ticks = (-1e-5:1e-5:1e-5), size =35, label = "N²' [s⁻²]", vertical = false, labelsize = 30)

    rowsize!(ga, 1, Relative(0.05))
    rowsize!(ga, 3, Relative(0.05))

savename = apath * Phavg_filename * "_Hovmoller" 
save(savename * ".png", f1)

# inverting the order
f1 = Figure(size = (2000, 1400), fontsize=26)

    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[2, 1], ylabel = "b [ms⁻²]") #vh
    ax2 = Axis(ga[2, 2]) 
    ax3 = Axis(ga[2, 3]) 
    ax4 = Axis(ga[4, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax5 = Axis(ga[4, 2], xlabel = "Tσ" ) 
    ax6 = Axis(ga[4, 3], xlabel = "Tσ" ) 

    gcb1 = ga[1, 1:3] = GridLayout()
    gcb2 = ga[3, 1:3] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3
    ax4.yticks = -5e-3:1e-3:-1e-3

    ax4.xticks = .2:.2:1
    ax5.xticks = .2:.2:1
    ax6.xticks = .2:.2:1

    limits!(ax1, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax2, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax3, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax4, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax5, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax6, -0.03, 0.99, -5.2e-3, -0.8e-3)

    hidedecorations!(ax2)
    hidedecorations!(ax3)
    hidexdecorations!(ax1)
    hideydecorations!(ax5)
    hideydecorations!(ax6)

    scaling_levelsn = -.1:.5:0
    scaling_levelsp = 0:.1:.1

    N2scaling_levelsn = -5e-6:5e-6:-.1e-6
    N2scaling_levelsp = 0:1e-5:1e-5

    hV = heatmap!(ax4, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    heatmap!(ax5, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax5, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hA = heatmap!(ax6, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax6, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hv = heatmap!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg', colormap = :balance, colorrange = (-0.22,0.22))
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hE = heatmap!(ax2, phase_times, b_bins_bot_velocity, log_eavg_inb_Phavg', colormap = :thermal, colorrange = (-7.25,-4.25))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    hN = heatmap!(ax3, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg', colormap = :balance, colorrange = (-2e-5,2e-5))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 3, labels = true, labelsize = 30, labelfont = :bold)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 3, linestyle = :dash, labels = true, labelsize = 30, labelfont = :bold)    

    cb1 = Colorbar(gcb2[1,1], hV, ticks = (-500:250:500), size =35, label = "∂V(b,t)/∂t [m³s⁻¹]", vertical = false, labelsize = 30)
    cb2 = Colorbar(gcb2[1,2], hV, ticks = (-500:250:500), size =35, label = "- M(b,t) [m³s⁻¹]", vertical = false, labelsize = 30)
    cb3 = Colorbar(gcb2[1,3], hA, ticks = (-500:250:500), size =35, label = "A(b+Δb,t) - A(b,t) [m³s⁻¹]", vertical = false, labelsize = 30)
    cb4 = Colorbar(gcb1[1,1], hv, ticks = (-0.2:0.1:0.2), size =35, label = "v [ms⁻¹]", vertical = false, labelsize = 30)
    cb5 = Colorbar(gcb1[1,2], hE, ticks = (-7:1:-5, ["10⁻⁷", "10⁻⁸", "10⁻⁹"]), size =35, label = "ε [m³s⁻²]", vertical = false, labelsize = 30)
    cb6 = Colorbar(gcb1[1,3], hN, ticks = (-1e-5:1e-5:1e-5), size =35, label = "N²' [s⁻²]", vertical = false, labelsize = 30)

    cb1.alignmode = Mixed(top = 0)
    cb2.alignmode = Mixed(top = 0)
    cb3.alignmode = Mixed(top = 0)

    cb4.alignmode = Mixed(top = 0)
    cb5.alignmode = Mixed(top = 0)
    cb6.alignmode = Mixed(top = 0)

    rowsize!(ga, 1, Relative(0.1))
    rowsize!(ga, 3, Relative(0.1))

    colgap!(ga, 30)
    rowgap!(ga, 30)
    colgap!(gcb2, 30)
    rowgap!(gcb2, 30)
    colgap!(gcb1, 30)
    rowgap!(gcb1, 30)
    
    Label(gcb1[1, 1, TopLeft()], "a",fontsize = 30, font = :bold, padding = (0, 5, 5, 5), halign = :left)
    Label(gcb1[1, 2, TopLeft()], "b",fontsize = 30, font = :bold, padding = (0, 5, 5, 5), halign = :left)
    Label(gcb1[1, 3, TopLeft()], "c",fontsize = 30, font = :bold, padding = (0, 5, 5, 5), halign = :left)
    Label(gcb2[1, 1, TopLeft()], "d",fontsize = 30, font = :bold, padding = (0, 5, 5, 5), halign = :left)
    Label(gcb2[1, 2, TopLeft()], "e",fontsize = 30, font = :bold, padding = (0, 5, 5, 5), halign = :left)
    Label(gcb2[1, 3, TopLeft()], "f",fontsize = 30, font = :bold, padding = (0, 5, 5, 5), halign = :left)


savename = apath * Phavg_filename * "_Hovmoller_reverse" 
save(savename * ".png", f1, px_per_unit = 2)

### PLOT FOR INITIAL STORY:
# v, e grda b?

f1 = Figure(size = (2000, 800), fontsize=35)

    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[2, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax2 = Axis(ga[2, 2], xlabel = "Tσ") 
    ax3 = Axis(ga[2, 3], xlabel = "Tσ") 

    gcb1 = ga[1, 1:3] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3

    ax1.xticks = .2:.2:1
    ax2.xticks = .2:.2:1
    ax3.xticks = .2:.2:1

    limits!(ax1, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax2, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax3, -0.03, 0.99, -5.2e-3, -0.8e-3)

    hideydecorations!(ax2)
    hideydecorations!(ax3)

    scaling_levelsn = -.1:.5:0
    scaling_levelsp = 0:.1:.1

    hv = heatmap!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg', colormap = :balance, colorrange = (-0.2,0.2))
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    hE = heatmap!(ax2, phase_times, b_bins_bot_velocity, log_eavg_inb_Phavg', colormap = :thermal, colorrange = (-7.5,-4.5))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    hN = heatmap!(ax3, phase_times, b_bins_bot_velocity, N2avg_inb_Phavg', colormap = :balance, colorrange = (-2e-5,2e-5))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    cb4 = Colorbar(gcb1[1,1], hv, ticks = (-0.2:0.1:0.2), size =45, label = "v [ms⁻¹]", vertical = false, labelsize = 40)
    cb4 = Colorbar(gcb1[1,2], hE, ticks = (-7:1:-5, ["10⁻⁷", "10⁻⁶", "10⁻⁵"]), size =45, label = "ε [m³s⁻²]", vertical = false, labelsize = 40)
    cb4 = Colorbar(gcb1[1,3], hN, ticks = (-1e-5:1e-5:1e-5), size =45, label = "ΔN² [s⁻²]", vertical = false, labelsize = 40)

    rowsize!(ga, 1, Relative(0.07))
    colgap!(ga, 15)

Ph_savename = "ven2_nb_25_110_Phavg_U350N100Lz100g100"
savename = apath * Ph_savename * "_Hovmoller" 
save(savename * ".png", f1, px_per_unit = 2)


f1 = Figure(size = (2000, 800), fontsize=35)

    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[2, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax2 = Axis(ga[2, 2], xlabel = "Tσ") 
    ax3 = Axis(ga[2, 3], xlabel = "Tσ") 

    gcb1 = ga[1, 1:3] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3

    ax1.xticks = .2:.2:1
    ax2.xticks = .2:.2:1
    ax3.xticks = .2:.2:1

    limits!(ax1, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax2, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax3, -0.03, 0.99, -5.2e-3, -0.8e-3)

    hideydecorations!(ax2)
    hideydecorations!(ax3)

    scaling_levelsn = -.1:.5:0
    scaling_levelsp = 0:.1:.1

    hv = heatmap!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg', colormap = :balance, colorrange = (-0.2,0.2))
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    hE = heatmap!(ax2, phase_times, b_bins_bot_velocity, log_eavg_inb_Phavg', colormap = :thermal, colorrange = (-7.5,-4.5))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    hN = heatmap!(ax3, phase_times, b_bins_bot_velocity, M2avg_inb_Phavg', colormap = :balance, colorrange = (-4e-6,4e-6))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    cb4 = Colorbar(gcb1[1,1], hv, ticks = (-0.2:0.1:0.2), size =45, label = "v [ms⁻¹]", vertical = false, labelsize = 40)
    cb4 = Colorbar(gcb1[1,2], hE, ticks = (-7:1:-5, ["10⁻⁷", "10⁻⁶", "10⁻⁵"]), size =45, label = "ε [m³s⁻²]", vertical = false, labelsize = 40)
    cb4 = Colorbar(gcb1[1,3], hN, ticks = (-2e-6:2e-6:2e-6), size =45, label = "M² [s⁻²]", vertical = false, labelsize = 40)

    rowsize!(ga, 1, Relative(0.07))
    colgap!(ga, 15)

Ph_savename = "vem2_nb_25_110_Phavg_U350N100Lz100g100"
savename = apath * Ph_savename * "_Hovmoller" 
save(savename * ".png", f1, px_per_unit = 2)


###
#  v dV/dt, M, A
###

f1 = Figure(size = (2800, 800), fontsize=35)

    ga = f1[1, 1] = GridLayout()    

    ax1 = Axis(ga[2, 1], ylabel = "b [ms⁻²]", xlabel = "Tσ") #vh
    ax2 = Axis(ga[2, 3], xlabel = "Tσ") 
    ax3 = Axis(ga[2, 4], xlabel = "Tσ") 
    ax4 = Axis(ga[2, 5], xlabel = "Tσ") 

    gcb1 = ga[1, 1:5] = GridLayout()

    ax1.yticks = -5e-3:1e-3:-1e-3

    ax1.xticks = .2:.2:1
    ax2.xticks = .2:.2:1
    ax3.xticks = .2:.2:1
    ax4.xticks = .2:.2:1

    limits!(ax1, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax2, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax3, -0.03, 0.99, -5.2e-3, -0.8e-3)
    limits!(ax4, -0.03, 0.99, -5.2e-3, -0.8e-3)

    hideydecorations!(ax2)
    hideydecorations!(ax3)
    hideydecorations!(ax4)

    scaling_levelsn = -.1:.5:0
    scaling_levelsp = 0:.1:.1

    hv = heatmap!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg', colormap = :balance, colorrange = (-0.2,0.2))
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax1, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    hdV = heatmap!(ax2, phase_times, b_bins_bot, ∂tV_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax2, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    heatmap!(ax3, phase_times, b_bins_bot, -1 .* ∫M_inb_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax3, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    hA = heatmap!(ax4, phase_times, b_bins_bot, Adif_Phavg', colormap = :balance, colorrange = (-650,650))
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsp, color = :black , linewidth = 4, labels = true, labelsize = 35, labelfont = :bold)    
    contour!(ax4, phase_times, b_bins_bot_velocity, vavg_inb_Phavg'; levels = scaling_levelsn, color = :black , linewidth = 4, linestyle = :dash, labels = true, labelsize = 35, labelfont = :bold)    

    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-0.2:0.1:0.2), size =45, label = "v [ms⁻¹]", vertical = false, labelsize = 40)
    cb2 = Colorbar(gcb1[1,3], hdV, ticks = (-500:250:500), size =45, label = "∂V(b,t)/∂t [m³s⁻¹]", vertical = false, labelsize = 40)
    cb3 = Colorbar(gcb1[1,4], hdV, ticks = (-500:250:500), size =45, label = "- M(b,t) [m³s⁻¹]", vertical = false, labelsize = 40)
    cb4 = Colorbar(gcb1[1,5], hA, ticks = (-500:250:500), size =45, label = "A(b+Δb,t) - A(b,t) [m³s⁻¹]", vertical = false, labelsize = 40)
    
    rowsize!(ga, 1, Relative(0.05))
    colsize!(ga, 2, Relative(0.07))
    colsize!(gcb1, 2, Relative(0.07))

    colgap!(ga, 15)

Ph_savename = "vdVMA_nb_25_110_Phavg_U350N100Lz100g100"
savename = apath * Ph_savename * "_Hovmoller" 
save(savename * ".png", f1, px_per_unit = 2)

