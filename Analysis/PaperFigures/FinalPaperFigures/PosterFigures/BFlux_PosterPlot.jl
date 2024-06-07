using Statistics
using Printf
using Measures
using JLD2
using CairoMakie

dpath = "Data/"
apath = "Analysis/PaperFigures/FinalPaperFigures/PosterFigures/"

sn = "U350N100Lz100g100" #"U250N100Lz100g100"

filescalename = dpath * "BFluxTrip_full_mp_noS_xavg_" * sn * ".jld2" #"BFluxTrip_fullxavg_mp_" * sn * ".jld2"

include("../../../../parameters.jl") 
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

### [-v'b'p]
### [-w'b'p]
### [-∇p]

phlims = (-nabmaxp, nabmaxp)
        
yb_cut = yb[1:ylength]
        
f1 = Figure(resolution = (1000, 1300), fontsize=36)
        ga = f1[1, 1] = GridLayout()    
        gb = f1[2, 1] = GridLayout()
        gc = f1[3, 1] = GridLayout()
    
        gcb1 = f1[1:3, 2] = GridLayout()       
    
        axvbp = Axis(ga[1, 1], ylabel = "z [m]") #v
        axwbp = Axis(gb[1, 1], ylabel = "z [m]") #v
        axgradp = Axis(gc[1, 1], ylabel = "z [m]", xlabel = "y [m]") #v
        
        axgradp.xticks = 500:500:1000

        axvbp.yticks = [-250, 0]
        axwbp.yticks = [-250, 0]
        axgradp.yticks = [-250, 0]
        
        limits!(axvbp, 0, 1500, -500, 0)
        limits!(axwbp, 0, 1500, -500, 0)
        limits!(axgradp, 0, 1500, -500, 0)
    
        hidexdecorations!(axvbp)
        hidexdecorations!(axwbp)
    
        colsize!(f1.layout, 2, Relative(0.05))
    
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

         # create colorbars the size of the whole data set
        cb1 = Colorbar(gcb1[1,1], hv, ticks = (-nabmaxp/2:nabmaxp/2:nabmaxp/2), size =45, label = "Wave Flux Terms [ms⁻³]")

        colgap!(ga, 15)
        rowgap!(ga, 5)
        colgap!(gb, 15)
        rowgap!(gb, 5)
        colgap!(gc, 15)
        rowgap!(gc, 5)
    
save(apath * "BFluxTrip_waveterms_" * sn * ".png", f1, px_per_unit = 2)

  
f1 = Figure(resolution = (2000, 700), fontsize=36)
        ga = f1[2, 1] = GridLayout()    
        gb = f1[2, 2] = GridLayout()
        gc = f1[2, 3] = GridLayout()
    
        gcb1 = f1[1, 1:3] = GridLayout()       
    
        axvbp = Axis(ga[1, 1], ylabel = "z [m]", xlabel = "y [m]") #v
        axwbp = Axis(gb[1, 1], xlabel = "y [m]") #v
        axgradp = Axis(gc[1, 1], xlabel = "y [m]") #v
        
        axvbp.xticks = 500:500:1000
        axwbp.xticks = 500:500:1000
        axgradp.xticks = 500:500:1000

        axvbp.yticks = [-250, 0]
        axwbp.yticks = [-250, 0]
        axgradp.yticks = [-250, 0]
        
        limits!(axvbp, 0, 1500, -500, 0)
        limits!(axwbp, 0, 1500, -500, 0)
        limits!(axgradp, 0, 1500, -500, 0)
    
        hideydecorations!(axwbp)
        hideydecorations!(axgradp)
    
        rowsize!(f1.layout, 1, Relative(0.05))
    
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

         # create colorbars the size of the whole data set
        cb1 = Colorbar(gcb1[1,1], hv, ticks = (-nabmaxp/2:nabmaxp/2:nabmaxp/2), size =45, label = "Wave Buoyancy Flux Divergence [ms⁻³]", vertical = false)

        colgap!(ga, 15)
        rowgap!(ga, 5)
        colgap!(gb, 15)
        rowgap!(gb, 5)
        colgap!(gc, 15)
        rowgap!(gc, 5)
    
save(apath * "BFluxTrip_waveterms_horiz_" * sn * ".png", f1, px_per_unit = 2)


f1 = Figure(resolution = (2500, 650), fontsize=36)
        ga = f1[1, 1] = GridLayout()    
        gb = f1[1, 2] = GridLayout()
        gc = f1[1, 3] = GridLayout()
    
        gcb1 = f1[1, 4] = GridLayout()       
    
        axvbp = Axis(ga[1, 1], ylabel = "z [m]", xlabel = "y [m]") #v
        axwbp = Axis(gb[1, 1], xlabel = "y [m]") #v
        axgradp = Axis(gc[1, 1], xlabel = "y [m]") #v
        
        axvbp.xticks = 500:500:1000
        axwbp.xticks = 500:500:1000
        axgradp.xticks = 500:500:1000

        axvbp.yticks = [-250, 0]
        axwbp.yticks = [-250, 0]
        axgradp.yticks = [-250, 0]
        
        limits!(axvbp, 0, 1500, -500, 0)
        limits!(axwbp, 0, 1500, -500, 0)
        limits!(axgradp, 0, 1500, -500, 0)
    
        hideydecorations!(axwbp)
        hideydecorations!(axgradp)
    
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

         # create colorbars the size of the whole data set
        cb1 = Colorbar(gcb1[1,1], hv, ticks = (-nabmaxp/2:nabmaxp/2:nabmaxp/2), size =45, label = "Buoyancy Flux Divergence [ms⁻³]")

        colgap!(ga, 15)
        rowgap!(ga, 5)
        colgap!(gb, 15)
        rowgap!(gb, 5)
        colgap!(gc, 15)
        rowgap!(gc, 5)
    
save(apath * "BFluxTrip_waveterms_horiz2_" * sn * ".png", f1, px_per_unit = 2)
