using Statistics
using Printf
using Measures
using JLD2
using CairoMakie

dpath = "Data/"
apath = "Analysis/Plots/"

sn = "U350N100Lz100g100" #"U250N100Lz100g100"

filescalename = dpath * "BFluxTrip_full_mp_noS_xavg_" * sn * ".jld2" #"BFluxTrip_fullxavg_mp_" * sn * ".jld2"

include("../../parameters.jl") 
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       nz = round(Int,pm.Lz/2),
                       m = -π/pm.Lz,
                       l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                       Tf = 2*π/pm.f, 
                       Tσ = 2*π/pm.σ))

scale_file = jldopen(filescalename, "r+")

skeys = keys(scale_file)

∇_phasedepWavg_xavg  = scale_file[skeys[1]];
∇_turbWavg_xavg  = scale_file[skeys[2]];
vb_phasedepWavg_dy_xavg  = scale_file[skeys[3]];
wb_phasedepWavg_dz_xavg  = scale_file[skeys[4]];
vb_turbWavg_dy_xavg  = scale_file[skeys[5]];
wb_turbWavg_dz_xavg = scale_file[skeys[6]];
SGS_Wavg_xavg = scale_file[skeys[7]];
u∇b_Wavg_xavg = scale_file[skeys[8]];

yb = scale_file[skeys[9]];
zb = scale_file[skeys[10]];
#phase_times = scale_file[skeys[10]];

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ

Sp_extra = ifelse(z_start>0, 250.0, 0.0)
Ly = pm.Lyˢ+Sp_extra
ny = round(Int,Ly/4)
slope_end = pm.Lzˢ/pm.Tanα

pm = merge(pm, (;Ly=Ly,ny=ny, slope_end=slope_end, Sp_extra=Sp_extra))

# if slope is in different spot than usual, need to move the curved part too!
const zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

land = curvedslope.(yb)
land_pdel = (curvedslope.(yb) .+ pm.U₀/pm.Ñ)[1:382]

zlength = length(zb)
ylength = 880

nabmaxp = 4e-7 #round(maximum(abs.(wb_phasedepWavg_dz_xavg))*1e7/4)*1e-7

nabmaxt = 4e-7 #round(maximum(abs.(wb_turbWavg_dz_xavg))*1e8/5)*1e-8
umaxt = 3.2e-7 #round(maximum(abs.(u∇b_Wavg_xavg))*1e8/4)*1e-8
sgsmax = 1e-8 
#dmax =  2.4e-7#1e-7


    ###[ v'b'p  ][ w'b'p  ][ -∇p        ]
    ###[ v'b't  ][ w'b't  ][ -∇t        ]

    phlims = (-nabmaxp, nabmaxp)
    tlims = (-nabmaxt, nabmaxt)
    ulims = (-umaxt, umaxt)
    sgslims = (-sgsmax, sgsmax)
        
    yb_cut = yb[1:ylength]
        
f1 = Figure(resolution = (1800, 700), fontsize=26)
        ga = f1[1:2, 1] = GridLayout()    
        gb = f1[1:2, 2] = GridLayout()
        gc = f1[1:2, 3] = GridLayout()
    
        gcb1 = f1[1, 4] = GridLayout()       
        gcb2 = f1[2, 4] = GridLayout()       
    
        axvbp = Axis(ga[1, 1], ylabel = "z [m]") #v
        axvbt = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
        
        axvbt.xticks = 500:500:1000

        axvbp.yticks = [-250, 0]
        axvbt.yticks = [-250, 0]
    
        axwbp = Axis(gb[1, 1]) #vh'b'
        axwbt = Axis(gb[2, 1], xlabel = "y [m]") #dy v'b'
        
        axgradp = Axis(gc[1, 1]) #wh'b'
        axgradt = Axis(gc[2, 1], xlabel = "y [m]") #dz w'b'
        axgradt.xticks = 500:500:1000
    
        limits!(axvbp, 0, 1500, -500, 0)
        limits!(axvbt, 0, 1500, -500, 0)
        limits!(axwbp, 0, 1500, -500, 0)
        limits!(axwbt, 0, 1500, -500, 0)
        limits!(axgradp, 0, 1500, -500, 0)
        limits!(axgradt, 0, 1500, -500, 0)
    
        hidedecorations!(axwbp)
        hidedecorations!(axgradp)
        hidexdecorations!(axvbp)
        hideydecorations!(axwbt)
        hideydecorations!(axgradt)
    
        colsize!(f1.layout, 4, Relative(0.05))
    
        hv = heatmap!(axvbp, yb_cut[2:end], zb[2:end], -1 .* vb_phasedepWavg_dy_xavg, colormap = :balance, colorrange = phlims)
            band!(axvbp, yb, land,-500, color=:black)
            lines!(axvbp, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
            #text!(axvbp, Point.(50, -450), text = "- ∂y(ṽb̃)", align = (:left, :center), color = :black, 
            #font = :bold, fontsize = 26)
    
        hw = heatmap!(axwbp, yb_cut, zb[2:end], -1 .* wb_phasedepWavg_dz_xavg, colormap = :balance, colorrange = phlims)
            band!(axwbp, yb, land,-500, color=:black)
            lines!(axwbp, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
         #   text!(axwbp, Point.(50, -450), text = "- ∂z(w̃b̃)", align = (:left, :center), color = :black, 
         #   font = :bold, fontsize = 26)

        hnab = heatmap!(axgradp, yb_cut[2:end], zb[2:end], -1 .* ∇_phasedepWavg_xavg, colormap = :balance, colorrange = phlims)
            band!(axgradp, yb, land,-500, color=:black)
            lines!(axgradp, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
         #   text!(axgradp, Point.(50, -450), text = "- ∇⋅(ũ⃗b̃)", align = (:left, :center), color = :black, 
         #   font = :bold, fontsize = 26)

        heatmap!(axvbt, yb_cut[2:end], zb[2:end],  -1 .* vb_turbWavg_dy_xavg, colormap = :balance, colorrange = tlims)
        band!(axvbt, yb, land,-500, color=:black)
        lines!(axvbt, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
         #   text!(axvbt, Point.(50, -450), text = "- ∂y(v'b')", align = (:left, :center), color = :black, 
         #   font = :bold, fontsize = 26)
            
        heatmap!(axwbt, yb_cut[2:end], zb[2:end],  -1 .* wb_turbWavg_dz_xavg, colormap = :balance, colorrange = tlims)
        band!(axwbt, yb, land,-500, color=:black)
        lines!(axwbt, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
         #   text!(axwbt, Point.(50, -450), text = "- ∂z(w'b')", align = (:left, :center), color = :black, 
         #   font = :bold, fontsize = 26)

        hn = heatmap!(axgradt, yb_cut[2:end], zb[2:end], -1 .* ∇_turbWavg_xavg, colormap = :balance, colorrange = tlims)
        band!(axgradt, yb, land,-500, color=:black)
        lines!(axgradt, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
          #  text!(axgradt, Point.(50, -450), text = "- ∇⋅(u⃗'b')", align = (:left, :center), color = :black, 
         #   font = :bold, fontsize = 26)
    
        # create colorbars the size of the whole data set
        cb1 = Colorbar(gcb1[1,1], hv, ticks = (-nabmaxp/2:nabmaxp/2:nabmaxp/2), size =35, label = "Wave Terms [ms⁻³]")
        cb2 = Colorbar(gcb2[1,1], hn, ticks = (-nabmaxt/2:nabmaxt/2:nabmaxt/2), size =35, label = "Turbulent Terms [ms⁻³]" )

        colgap!(ga, 15)
        rowgap!(ga, 5)
        colgap!(gb, 15)
        rowgap!(gb, 5)
        colgap!(gc, 15)
        rowgap!(gc, 5)
    
        #display(f1)
save(apath * "BFluxTrip_full_mp_noS_2_" * sn * ".png", f1, px_per_unit = 2)

 
f1 = Figure(resolution = (1800, 1000), fontsize=26)
        ga = f1[1:3, 1] = GridLayout()    
        gb = f1[1:3, 2] = GridLayout()
        gc = f1[1:3, 3] = GridLayout()
    
        gcb1 = f1[1, 4] = GridLayout()       
        gcb2 = f1[2, 4] = GridLayout()       
        gcb3 = f1[3, 4] = GridLayout()       
    
        axvbp = Axis(ga[1, 1], ylabel = "z [m]") #v
        axvbt = Axis(ga[2, 1], ylabel = "z [m]") #vh
        axu = Axis(ga[3, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
        
        axu.xticks = 500:500:1000

        axvbp.yticks = [-250, 0]
        axvbt.yticks = [-250, 0]
        axu.yticks = [-250, 0]
    
        axwbp = Axis(gb[1, 1]) #vh'b'
        axwbt = Axis(gb[2, 1]) #dy v'b'
        axnab = Axis(gb[3, 1], xlabel = "y [m]") #dy v'b'

        axgradp = Axis(gc[1, 1]) #wh'b'
        axgradt = Axis(gc[2, 1]) #dz w'b'
        axdif = Axis(gc[3, 1], xlabel = "y [m]") #dz w'b'

        axdif.xticks = 500:500:1000
        axnab.xticks = 500:500:1000
    
        limits!(axvbp, 0, 1500, -500, 0)
        limits!(axvbt, 0, 1500, -500, 0)
        limits!(axwbp, 0, 1500, -500, 0)
        limits!(axwbt, 0, 1500, -500, 0)
        limits!(axgradp, 0, 1500, -500, 0)
        limits!(axgradt, 0, 1500, -500, 0)
        limits!(axu, 0, 1500, -500, 0)
        limits!(axnab, 0, 1500, -500, 0)
        limits!(axdif, 0, 1500, -500, 0)
    
        hidedecorations!(axwbp)
        hidedecorations!(axgradp)
        hidexdecorations!(axvbp)
        hidedecorations!(axwbt)
        hidedecorations!(axgradt)
        hidexdecorations!(axvbt)
        hideydecorations!(axnab)
        hideydecorations!(axdif)
    
        colsize!(f1.layout, 4, Relative(0.05))
    
        hv = heatmap!(axvbp, yb_cut[2:end], zb[2:end], -1 .* vb_phasedepWavg_dy_xavg, colormap = :balance, colorrange = ulims)
            lines!(axvbp, yb, land, color=:black, lw = 4)
            lines!(axvbp, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axvbp, Point.(50, -450), text = "- ∂y(ṽb̃)", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)
    
        hw = heatmap!(axwbp, yb_cut, zb[2:end], -1 .* wb_phasedepWavg_dz_xavg, colormap = :balance, colorrange = ulims)
            lines!(axwbp, yb, land, color=:black, lw = 4)
            lines!(axwbp, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axwbp, Point.(50, -450), text = "- ∂z(w̃b̃)", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        hnab = heatmap!(axgradp, yb_cut[2:end], zb[2:end], -1 .* ∇_phasedepWavg_xavg, colormap = :balance, colorrange = ulims)
            lines!(axgradp, yb, land, color=:black, lw = 4)
            lines!(axgradp, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axgradp, Point.(50, -450), text = "- ∇⋅(u⃗̃b̃)", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        heatmap!(axvbt, yb_cut[2:end], zb[2:end],  -1 .* vb_turbWavg_dy_xavg, colormap = :balance, colorrange = tlims)
            lines!(axvbt, yb, land, color=:black, lw = 4)
            lines!(axvbt, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axvbt, Point.(50, -450), text = "- ∂y(v'b')", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)
    
        heatmap!(axwbt, yb_cut[2:end], zb[2:end],  -1 .* wb_turbWavg_dz_xavg, colormap = :balance, colorrange = tlims)
            lines!(axwbt, yb, land, color=:black, lw = 4)
            lines!(axwbt, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axwbt, Point.(50, -450), text = "- ∂z(w'b')", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        hn = heatmap!(axgradt, yb_cut[2:end], zb[2:end], -1 .* ∇_turbWavg_xavg, colormap = :balance, colorrange = tlims)
            lines!(axgradt, yb, land, color=:black, lw = 4)
            lines!(axgradt, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axgradt, Point.(50, -450), text = "- ∇⋅(u⃗'b')", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        hu = heatmap!(axu, yb_cut[2:end], zb[2:end],  u∇b_Wavg_xavg, colormap = :balance, colorrange = ulims)
            lines!(axu, yb, land, color=:black, lw = 4)
            lines!(axu, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axu, Point.(50, -450), text = "u⃗̄⋅∇b̄", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)
    
        heatmap!(axnab, yb_cut[2:end], zb[2:end], SGS_Wavg_xavg[2:end,2:end] .- ∇_phasedepWavg_xavg .- ∇_turbWavg_xavg, colormap = :balance, colorrange = ulims)
            lines!(axnab, yb, land, color=:black, lw = 4)
            lines!(axnab, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axnab, Point.(30, -450), text = "- ∇⋅(u⃗̃b̃) - ∇⋅(u⃗'b') + ∇(κ⋅∇b)", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        hd = heatmap!(axdif, yb_cut[2:end], zb[2:end], SGS_Wavg_xavg[2:end,2:end], colormap = :balance, colorrange = (-5e-9,5e-9))
            lines!(axdif, yb, land, color=:black, lw = 4)
            lines!(axdif, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axdif, Point.(30, -450), text = "∇(κ⋅∇b)", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)
        # create colorbars the size of the whole data set
        cb1 = Colorbar(gcb1[1,1], hv, ticks = (-umaxt/2:umaxt/2:umaxt/2), size =35, label = "Wave Terms [ms⁻³]")
        cb2 = Colorbar(gcb2[1,1], hn, ticks = (-nabmaxt/2:nabmaxt/2:nabmaxt/2), size =35, label = "Turbulent Terms [ms⁻³]" )
        cb3 = Colorbar(gcb3[1,1], hu, ticks = (-umaxt/2:umaxt/2:umaxt/2), size =35, label = "Mean Advection [ms⁻³]" )#,flipaxis=false,)
        #cb4 = Colorbar(gcb3[1,1], hd, ticks = (-dmax/2:dmax/2:dmax/2), size =35, label = "Difference [ms⁻³]" )

        colgap!(ga, 15)
        rowgap!(ga, 5)
        colgap!(gb, 15)
        rowgap!(gb, 5)
        colgap!(gc, 15)
        rowgap!(gc, 5)
    
        #display(f1)
save(apath * "BFluxTrip_full_mp_noS_wu_" * sn * ".png", f1, px_per_unit = 2)
 

f1 = Figure(resolution = (1800, 600), fontsize=26)
        ga = f1[1, 1] = GridLayout()    
        
        axu = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
        
        axu.xticks = 500:500:1000

        axu.yticks = [-500, -250, 0]
    
        axnab = Axis(ga[2,2], xlabel = "y [m]") #dy v'b'

        axdif = Axis(ga[2, 3], xlabel = "y [m]") #dz w'b'

        axdif.xticks = 500:500:1000
        axnab.xticks = 500:500:1000
    
        limits!(axu, 0, 1500, -500, 0)
        limits!(axnab, 0, 1500, -500, 0)
        limits!(axdif, 0, 1500, -500, 0)
    
        hideydecorations!(axnab)
        hideydecorations!(axdif)
    
        hu = heatmap!(axu, yb_cut[2:end], zb[2:end],  u∇b_Wavg_xavg, colormap = :balance, colorrange = ulims)
            lines!(axu, yb, land, color=:black, lw = 4)
            lines!(axu, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axu, Point.(50, -450), text = "u⃗̄⋅∇b̄", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)
    
        hf = heatmap!(axnab, yb_cut[2:end], zb[2:end], -1 .*  ∇_phasedepWavg_xavg .- ∇_turbWavg_xavg, colormap = :balance, colorrange = ulims)
            lines!(axnab, yb, land, color=:black, lw = 4)
            lines!(axnab, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axnab, Point.(30, -450), text = "- ∇⋅(u⃗̃b̃) - ∇⋅(u⃗'b')", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        hd = heatmap!(axdif, yb_cut[2:end], zb[2:end],  SGS_Wavg_xavg[2:end,2:end] , colormap = :balance, colorrange = sgslims)
            lines!(axdif, yb, land, color=:black, lw = 4)
            lines!(axdif, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axdif, Point.(30, -450), text = "∇(κ⋅∇b)", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        # create colorbars the size of the whole data set
        cb1 = Colorbar(ga[1,1], hu, ticks = (-umaxt/2:umaxt/2:umaxt/2), size =35, label = "Advective Term [ms⁻³]", vertical = false,)
        cb2 = Colorbar(ga[1,2], hf, ticks = (-umaxt/2:umaxt/2:umaxt/2), size =35, label = "Flux Terms [ms⁻³]" , vertical = false,)
        cb3 = Colorbar(ga[1,3], hd, ticks = (-sgsmax/2:sgsmax/2:sgsmax/2), size =35, label = "Subgrid-Scale Term [ms⁻³]" , vertical = false,) 
 
        rowsize!(ga, 1, Relative(0.05))

        colgap!(ga, 15)
        rowgap!(ga, 15)
    
        #display(f1)
#save(apath * "BFluxTrip_full_mp_noS_onlytotal_" * sn * ".png", f1, px_per_unit = 2)
save(apath * "BFluxTrip_full_mp_noS_onlytotal_Res_" * sn * ".png", f1, px_per_unit = 2)
 
    


f1 = Figure(resolution = (1800, 600), fontsize=26)
        ga = f1[1, 1] = GridLayout()    
        
        axu = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
        
        axu.xticks = 500:500:1000

        axu.yticks = [-500, -250, 0]
    
        axnab = Axis(ga[2,2], xlabel = "y [m]") #dy v'b'

        axdif = Axis(ga[2, 3], xlabel = "y [m]") #dz w'b'

        axdif.xticks = 500:500:1000
        axnab.xticks = 500:500:1000
    
        limits!(axu, 0, 1500, -500, 0)
        limits!(axnab, 0, 1500, -500, 0)
        limits!(axdif, 0, 1500, -500, 0)
    
        hideydecorations!(axnab)
        hideydecorations!(axdif)
    
        hu = heatmap!(axu, yb_cut, zb,  dbdt_first_Wavg, colormap = :balance, colorrange = (-2e-8,2e-8))
            lines!(axu, yb, land, color=:black, lw = 4)
            lines!(axu, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axu, Point.(50, -450), text = "db/dt", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        hf = heatmap!(axnab, yb_cut[2:end], zb[2:end], dbdt, colormap = :balance, colorrange = (-1e-8,1e-8))
            lines!(axnab, yb, land, color=:black, lw = 4)
            lines!(axnab, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axnab, Point.(30, -450), text = "∂b̄/∂t", align = (:left, :center), color = :black, 
            font = :bold, fontsize = 26)

        hd = heatmap!(axdif, yb_cut[2:end], zb[2:end],  -1 .*  ∇_phasedepWavg_xavg .- ∇_turbWavg_xavg .+ SGS_Wavg_xavg[2:end,2:end] .- u∇b_Wavg_xavg, colormap = :balance, colorrange = ulims)#colorrange = sgslims)
            lines!(axdif, yb, land, color=:black, lw = 4)
            lines!(axdif, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
            text!(axdif, Point.(30, -450), text = "Residual", align = (:left, :center), color = :black, #text = "∇(κ⋅∇b)"
            font = :bold, fontsize = 26)

        # create colorbars the size of the whole data set
        cb1 = Colorbar(ga[1,1], hu, ticks =  (-1e-8:1e-8:1e-8), size =35, label = "db/dt then wave avg [ms⁻³]", vertical = false,)
        cb2 = Colorbar(ga[1,2], hf, ticks = (-5e-9:5e-9:5e-9), size =35, label = "db̄/dt [ms⁻³]" , vertical = false,)
        #cb3 = Colorbar(ga[1,3], hd, ticks = (-sgsmax/2:sgsmax/2:sgsmax/2), size =35, label = "Subgrid-Scale Term [ms⁻³]" , vertical = false,) 
        cb3 = Colorbar(ga[1,3], hd, ticks = (-umaxt/2:umaxt/2:umaxt/2), size =35, label = "Residual Term [ms⁻³]" , vertical = false,) # label = "Subgrid-Scale Term [ms⁻³]" 
 
        rowsize!(ga, 1, Relative(0.05))

        colgap!(ga, 15)
        rowgap!(ga, 15)
    
        #display(f1)
#save(apath * "BFluxTrip_full_mp_noS_onlytotal_" * sn * ".png", f1, px_per_unit = 2)
save(apath * "BFluxTrip_full_mp_noS_onlytotal_Res_" * sn * ".png", f1, px_per_unit = 2)
 
    

