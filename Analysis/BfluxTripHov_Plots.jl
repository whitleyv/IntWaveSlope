using Statistics
using Printf
using Measures
using JLD2
using CairoMakie

path_name = "Data/"

sn1 = "U300N100Lz100g100" #"U150N100Lz100g100"

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

        
filescalename = path_name * "BFluxTrip_rAvg_prof_mp_" * sn1  * ".jld2"
fsc2 = path_name * "BFluxTripdbdt_rAvg_prof_mp.jld2"
scale_file = jldopen(filescalename, "r+")
dbdt_file = jldopen(fsc2, "r+")

skeys = keys(scale_file)

SGS_rWavg_prof = scale_file[skeys[1]];
c_rWavg_prof = scale_file[skeys[2]];
v_rWavg_prof = scale_file[skeys[3]];
∇_phasedep_rWavg_prof = -1 .* scale_file[skeys[4]];
∇_turb_rWavg_prof = -1 .* scale_file[skeys[5]];
wb_phasedep_rWavg_dz_prof = -1 .* scale_file[skeys[6]];
vb_phasedep_rWavg_dy_prof = -1 .* scale_file[skeys[7]];
wb_turb_rWavg_dz_prof = -1 .* scale_file[skeys[8]];
vb_turb_rWavg_dy_prof = -1 .* scale_file[skeys[9]];
rolling_phase_times1 = scale_file[skeys[12]];
rolling_phase_times2 = scale_file[skeys[13]];

dt_b_rWavg_prof150 = dbdt_file["dt_b_rWavg_prof150"];
dt_b_rWavg_prof250 = dbdt_file["dt_b_rWavg_prof250"];
dt_b_rWavg_prof300 = dbdt_file["dt_b_rWavg_prof300"];

∇_rWavg_prof = ∇_phasedep_rWavg_prof .+ ∇_turb_rWavg_prof;
dbdt_rWavg_prof = ∇_rWavg_prof .+ SGS_rWavg_prof[2:end, 16:end-15];

zL = 75
HAB = 2:2:zL*2
dHAB = HAB[2:end]
tlength = 135
# center of mass of profile:
c_rWavg_prof_CoM = sum(c_rWavg_prof .* HAB, dims = 1)[ 1, :] ./ sum(c_rWavg_prof, dims = 1)[ 1, :]

# v            c               
# wbdz        vbdy          -nablaph
# wbdz        vbdy          -nablat
# SGS   -nabturb-nabph       db/dt

vmax = round(maximum(v_rWavg_prof)*1e2)*1e-2
cmax = round(maximum(c_rWavg_prof)*1e2)*0.5*1e-2
cmax = 0.06
sgsmax = round(maximum(abs.(SGS_rWavg_prof))*1e8)*1e-8
dmaxp = round(maximum(abs.(vb_phasedep_rWavg_dy_prof))*1e7)*1e-7
nabmaxp = round(maximum(abs.(∇_phasedep_rWavg_prof))*1e7)*1e-7
dmaxt = round(maximum(abs.(vb_turb_rWavg_dy_prof))*1e8)*1e-8
nabmaxt = round(maximum(abs.(∇_turb_rWavg_prof))*1e8)*1e-8
nabmax = round(maximum(abs.(∇_rWavg_prof))*1e7)*1e-7

vlims = (-vmax, vmax)
clims = (0, cmax)
dplims = (-dmaxp,dmaxp)
dtlims = (-dmaxt,dmaxt)
nabplims = (-nabmaxp, nabmaxp)
nabtlims = (-nabmaxt, nabmaxt)
nablims = (-nabmax, nabmax)
sgslims = (-sgsmax, sgsmax)

delta = pm.U₀./pm.Ñ

f1 = Figure(resolution = (2000, 2000), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    axv = Axis(ga[1, 1], ylabel = "z' [m]") #vh
    axc = Axis(ga[1, 2], ) 
    axrealdbdt = Axis(ga[1, 3], )

    axdzp = Axis(ga[2, 1], ylabel = "z' [m]") #vh
    axdyp = Axis(ga[2, 2], ) 
    axnabp = Axis(ga[2, 3], ) 

    axdzt = Axis(ga[3, 1], ylabel = "z' [m]") #vh
    axdyt = Axis(ga[3, 2], ) 
    axnabt = Axis(ga[3, 3], ) 

    axsgs = Axis(ga[4, 1], ylabel = "z' [m]", xlabel = "Tσ") #vh
    axnab = Axis(ga[4, 2], xlabel = "Tσ" ) 
    axdbdt = Axis(ga[4, 3], xlabel = "Tσ") 

    gcb1 = ga[1, 4] = GridLayout()
    gcb2 = ga[2, 4] = GridLayout()       
    gcb3 = ga[3, 4] = GridLayout()       
    gcb4 = ga[4, 4] = GridLayout()       

    axv.yticks = [50, 100]
    axdzp.yticks = [50, 100]
    axdzt.yticks = [50, 100]
    axsgs.yticks = [50, 100]

    axdbdt.xticks = 2:2:10
    axnab.xticks = 2:2:10
    axsgs.xticks = 2:2:10

    limits!(axv, 1,10, 0, 150)
    limits!(axc, 1,10, 0, 150)
    limits!(axdzp, 1,10, 0, 150)
    limits!(axdyp, 1,10, 0, 150)
    limits!(axnabp, 1,10, 0, 150)
    limits!(axdzt, 1,10, 0, 150)
    limits!(axdyt, 1,10, 0, 150)
    limits!(axnabt, 1,10, 0, 150)
    limits!(axsgs, 1,10, 0, 150)
    limits!(axnab, 1,10, 0, 150)
    limits!(axdbdt,1,10, 0, 150)
    limits!(axrealdbdt,1,10, 0, 150)

    hidedecorations!(axc)
    hidedecorations!(axdyp)
    hidedecorations!(axnabp)
    hidedecorations!(axdyt)
    hidedecorations!(axnabt)
    hidexdecorations!(axv)
    hidexdecorations!(axdzp)
    hidexdecorations!(axdzt)
    hideydecorations!(axnab)
    hideydecorations!(axdbdt)
    hidedecorations!(axrealdbdt)

    hv = heatmap!(axv, rolling_phase_times1, HAB, v_rWavg_prof', colormap = :balance, colorrange = vlims)
    text!(axv, Point.(1.5, 130), text = "⟨v⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axv, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(axv, rolling_phase_times1, c_rWavg_prof_CoM, color = :dodgerblue2, linewidth = 4)

    hc = heatmap!(axc, rolling_phase_times1, HAB, c_rWavg_prof', colormap = :thermal, colorrange = clims)
    text!(axc, Point.(1.5, 130), text = "⟨c⟩", align = (:left, :center), color = :white, font = :bold, fontsize = 26)
    lines!(axc, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(axc, rolling_phase_times1, c_rWavg_prof_CoM, color = :dodgerblue2, linewidth = 4)

    hrdbdt = heatmap!(axrealdbdt, rolling_phase_times1, HAB, dt_b_rWavg_prof300', colormap = :balance, colorrange = nablims)
    text!(axrealdbdt, Point.(1.5, 130), text = "∂⟨b⟩/∂t", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axrealdbdt, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(axrealdbdt, rolling_phase_times1, c_rWavg_prof_CoM, color = :dodgerblue2, linewidth = 4)

    hwbp = heatmap!(axdzp, rolling_phase_times2, dHAB, wb_phasedep_rWavg_dz_prof', colormap = :balance, colorrange = dplims)
    text!(axdzp, Point.(1.5, 130), text = "-∂z⟨w̃b̃⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axdzp, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)

    hvbp = heatmap!(axdyp, rolling_phase_times2, dHAB, vb_phasedep_rWavg_dy_prof', colormap = :balance, colorrange = dplims)
    text!(axdyp, Point.(1.5, 130), text = "-∂y⟨ṽb̃⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axdyp, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)

    hnabp = heatmap!(axnabp, rolling_phase_times2, dHAB, ∇_phasedep_rWavg_prof', colormap = :balance, colorrange = nabplims)
    text!(axnabp, Point.(1.5, 130), text = "-∇⋅⟨u⃗̃b̃⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axnabp, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(axnabp, rolling_phase_times1, c_rWavg_prof_CoM, color = :dodgerblue2, linewidth = 4)

    hwbt = heatmap!(axdzt, rolling_phase_times2, dHAB, wb_turb_rWavg_dz_prof', colormap = :balance, colorrange = dtlims)
    text!(axdzt, Point.(1.5, 130), text = "-∂z⟨w'b'⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axdzt, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)

    hvbt = heatmap!(axdyt, rolling_phase_times2, dHAB, vb_turb_rWavg_dy_prof', colormap = :balance, colorrange = dtlims)
    text!(axdyt, Point.(1.5, 130), text = "-∂y⟨v'b'⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axdyt, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)

    hnabt = heatmap!(axnabt, rolling_phase_times2, dHAB, ∇_turb_rWavg_prof', colormap = :balance, colorrange = nabtlims)
    text!(axnabt, Point.(1.5, 130), text = "-∇⋅⟨u⃗'b'⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axnabt, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(axnabt, rolling_phase_times1, c_rWavg_prof_CoM, color = :dodgerblue2, linewidth = 4)

    hsgs = heatmap!(axsgs, rolling_phase_times1, HAB, SGS_rWavg_prof', colormap = :balance, colorrange = sgslims)
    text!(axsgs, Point.(1.5, 130), text = "⟨∇⋅κ∇b⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)    
    lines!(axsgs, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)

    hnab = heatmap!(axnab, rolling_phase_times2, dHAB, ∇_rWavg_prof', colormap = :balance, colorrange = nablims)
    text!(axnab, Point.(1.5, 130), text = "-∇⋅⟨u⃗̃b̃⟩-∇⋅⟨u⃗'b'⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axnab, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)

    hdbdt = heatmap!(axdbdt, rolling_phase_times2, dHAB, dbdt_rWavg_prof', colormap = :balance, colorrange = nablims)
    text!(axdbdt, Point.(1.5, 130), text = "-∇⋅⟨u⃗̃b̃⟩-∇⋅⟨u⃗'b'⟩+⟨∇⋅κ∇b⟩", align = (:left, :center), color = :black, font = :bold, fontsize = 26)
    lines!(axdbdt, 1:0.5:4, delta.*ones(length( 1:0.5:4)), color = :black, linewidth = 2, linestyle = :dash)
    lines!(axdbdt, rolling_phase_times1, c_rWavg_prof_CoM, color = :dodgerblue2, linewidth = 4)

    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax:vmax:vmax), size =35, flipaxis=false, label = "velocity")
    cb2 = Colorbar(gcb1[1,2], hc, ticks = (0:cmax/2:cmax), size =35, label = "tracer")

    cb3 = Colorbar(gcb2[1,1], hwbp, ticks = (-dmaxp/2:dmaxp/2:dmaxp/2), size =35, flipaxis=false, label = "Derivatives")
    cb6 = Colorbar(gcb2[1,2], hnabp, ticks = (-nabmaxp/2:nabmaxp/2:nabmaxp/2), size =35, label = "Sum")

    cb4 = Colorbar(gcb3[1,1], hwbt, ticks = (-dmaxt/2:dmaxt/2:dmaxt/2), size =35, flipaxis=false, label = "Derivatives")
    cb6 = Colorbar(gcb3[1,2], hnabt, ticks = (-nabmaxt/2:nabmaxt/2:nabmaxt/2), size =35, label = "Sum")

    cb5 = Colorbar(gcb4[1,2], hdbdt, ticks = (-nabmax/2:nabmax/2:nabmax/2), size =35, label = "db/dt")
    cb6 = Colorbar(gcb4[1,1], hsgs, ticks = (-sgsmax/2:sgsmax/2:sgsmax/2), size =35, flipaxis=false, label = "SGS")

    rd = round(delta)
    Label(ga[1, 1:3, Top()], "Buoyancy Evolution Terms for δ = $rd m",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)

#display(f1)
savename = "Analysis/Plots/NewerAnalysisFluxesDye/SlopeNormalHovmoller_wcCOm" * sn

save(savename * ".png", f1)
