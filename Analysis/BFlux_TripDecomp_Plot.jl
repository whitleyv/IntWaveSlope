using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using ArgParse
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/Analysis/"
sn = "U300N100Lz100g100"

filescalename1 = path_name * "BFluxTrip_beg_" * sn * ".jld2"
filescalename2 = path_name * "BFluxTrip_mid_" * sn * ".jld2"
filescalename3 = path_name * "BFluxTrip_end_" * sn * ".jld2"

include("parameters.jl") 
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       nz = round(Int,pm.Lz/2),
                       m = -π/pm.Lz,
                       l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                       Tf = 2*π/pm.f, 
                       Tσ = 2*π/pm.σ))

scale_file1 = jldopen(filescalename1, "r+")
scale_file2 = jldopen(filescalename2, "r+")
scale_file3 = jldopen(filescalename3, "r+")
                       
skeys = keys(scale_file1)

FluxDiv_turbWavg_beg  = scale_file1[skeys[1]];
FluxDiv_phasedepWavg_beg  = scale_file1[skeys[2]];
PhaseDependentVals_beg  = scale_file1[skeys[3]];
TurbVals_beg  = scale_file1[skeys[4]];
WaveAveragedVals_beg  = scale_file1[skeys[5]];

yb = scale_file1[skeys[6]];
zb = scale_file1[skeys[7]];
phase_times = scale_file1[skeys[8]];

FluxDiv_turbWavg_mid = scale_file2[skeys[1]];
FluxDiv_phasedepWavg_mid = scale_file2[skeys[2]];
PhaseDependentVals_mid = scale_file2[skeys[3]];
TurbVals_mid = scale_file2[skeys[4]];
WaveAveragedVals_mid = scale_file2[skeys[5]];

FluxDiv_turbWavg_end  = scale_file3[skeys[1]];
FluxDiv_phasedepWavg_end  = scale_file3[skeys[2]];
PhaseDependentVals_end  = scale_file3[skeys[3]];
TurbVals_end  = scale_file3[skeys[4]];
WaveAveragedVals_end  = scale_file3[skeys[5]];

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

vbmax = round(maximum(WaveAveragedVals_beg.vb_xpertxavgWavg)*1e6)*1e-6*0.5
wbmax = round(maximum(abs.(WaveAveragedVals_beg.wb_xpertxavgWavg))*1e6)*1e-6*0.5
nabmax = round(maximum(abs.(FluxDiv_Wavg_beg.vb_xpertxavgWavg_dy))*1e7)*1e-7*0.25
SGSmax = round(maximum(FluxDiv_Wavg_beg.wb_xpertxavgWavg_dz)*1e7)*1e-7*0.25

function WaveAverageFluxPlot(savename, apath , ylength, vb_phasedepWavg, wb_phasedepWavg, ∇_phasedepWavg, 
    vb_turbWavg, wb_turbWavg, ∇_turbWavg, SGS_Wavg, c_Wavg, vbmax, wbmax, nabmax, SGSmax)

    ###[ v'b'p  ][ w'b'p  ][ ∇p        ]
    ###[ v'b't  ][ w'b't  ][ ∇t        ]
    ###[ ∇p+∇t  ][ SGS    ][-∇p-∇t+SGS ]
    
    vblims = (-vbmax,vbmax)
    wblims = (-wbmax,wbmax)
    gradblims = (-nabmax, nabmax)
    SGSlims = (-SGSmax, SGSmax)

    Full∇_Wavg = ∇_turbWavg .+ ∇_phasedepWavg
    Fulldbdt_Wavg = Full∇_Wavg .+ SGS_Wavg

    yb_cut = yb[1:ylength]

    clog = log10.(clamp.(c_Wavg, 1e-8, 1))

    f1 = Figure(resolution = (1300, 1100), fontsize=26)
    ga = f1[2:4, 1] = GridLayout()
    gb = f1[2:4, 2] = GridLayout()
    gc = f1[2:4, 3] = GridLayout()

    gcb1 = f1[1, 1] = GridLayout()
    gcb2 = f1[1, 2] = GridLayout()       
    gcb3 = f1[1, 3] = GridLayout()       

    axv = Axis(ga[1, 1], ylabel = "z [m]") #v
    axvh = Axis(ga[2, 1], ylabel = "z [m]") #v
    axgrad = Axis(ga[3, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
    axgrad.xticks = 500:500:1000
    axv.yticks = [-250, 0]
    axvh.yticks = [-250, 0]
    axgrad.yticks = [-250, 0]

    axvb = Axis(gb[1, 1]) #v'b'
    axvbh = Axis(gb[2, 1]) #vh'b'
    axdyvb = Axis(gb[3, 1], xlabel = "y [m]") #dy v'b'

    axdyvb.xticks = 500:500:1000

    axwb = Axis(gc[1, 1]) #w'b'
    axwbh = Axis(gc[2, 1]) #wh'b'
    axdzwb = Axis(gc[3, 1], xlabel = "y [m]") #dz w'b'
    axdzwb.xticks = 500:500:1000

    limits!(axv, 0, 1500, -500, 0)
    limits!(axvh, 0, 1500, -500, 0)
    limits!(axvb, 0, 1500, -500, 0)
    limits!(axvbh, 0, 1500, -500, 0)
    limits!(axwb, 0, 1500, -500, 0)
    limits!(axwbh, 0, 1500, -500, 0)
    limits!(axgrad, 0, 1500, -500, 0)
    limits!(axdyvb, 0, 1500, -500, 0)
    limits!(axdzwb, 0, 1500, -500, 0)

    hidedecorations!(axvb)
    hidedecorations!(axwb)
    hidedecorations!(axvbh)
    hidedecorations!(axwbh)
    hidexdecorations!(axv)
    hidexdecorations!(axvh)
    hideydecorations!(axdyvb)
    hideydecorations!(axdzwb)

    rowsize!(f1.layout, 1, Relative(0.05))

    hv = heatmap!(axv, yb_cut, zb, vb_phasedepWavg, colormap = :balance, colorrange = vblims)
        lines!(axv, yb, land, color=:black, lw = 4)
        lines!(axv, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axv, Point.(50, -450), text = "⟨ṽb̃⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hvhat = heatmap!(axvh, yb_cut, zb, vb_turbWavg, colormap = :balance, colorrange = vblims)
        lines!(axvh, yb, land, color=:black, lw = 4)
        lines!(axvh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvh, Point.(50, -450), text = "⟨v'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hnab = heatmap!(axgrad, yb_cut[2:end], zb[2:end], Full∇_Wavg, colormap = :balance, colorrange = gradblims)
        lines!(axgrad, yb, land, color=:black, lw = 4)
        lines!(axgrad, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axgrad, Point.(50, -450), text = "∇⋅⟨u⃗'b'⟩ + ∇⋅⟨u⃗̃b̃⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hvb = heatmap!(axvb, yb_cut, zb, wb_phasedepWavg, colormap = :balance, colorrange = wblims)
        lines!(axvb, yb, land, color=:black, lw = 4)
        lines!(axvb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvb, Point.(50, -450), text = "⟨w̃b̃⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        contour!(axvb, yb_cut, zb, clog, color = :gray30, lw = 6, levels = -5:1:0, labels=true)

    heatmap!(axvbh, yb_cut, zb, wb_turbWavg, colormap = :balance, colorrange = wblims)
        lines!(axvbh, yb, land, color=:black, lw = 4)
        lines!(axvbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvbh, Point.(50, -450), text = "⟨w'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        
    hdyvb = heatmap!(axdyvb, yb_cut, zb, SGS_Wavg, colormap = :balance, colorrange = SGSlims)
        lines!(axdyvb, yb, land, color=:black, lw = 4)
        lines!(axdyvb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdyvb, Point.(50, -450), text = "SGS", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hwb = heatmap!(axwb, yb_cut, zb, ∇_phasedepWavg, colormap = :balance, colorrange = gradblims)
        lines!(axwb, yb, land, color=:black, lw = 4)
        lines!(axwb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axwb, Point.(50, -450), text = "∇⋅⟨u⃗̃b̃⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        contour!(axwb, yb_cut, zb, clog, color = :gray30, lw = 6, levels = -5:1:0)

    heatmap!(axwbh, yb_cut, zb, ∇_turbWavg, colormap = :balance, colorrange = gradblims)
        lines!(axwbh, yb, land, color=:black, lw = 4)
        lines!(axwbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axwbh, Point.(50, -450), text = "∇⋅⟨u⃗'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hdzwb = heatmap!(axdzwb, yb_cut, zb, Fulldbdt_Wavg, colormap = :balance, colorrange = gradblims)
        lines!(axdzwb, yb, land, color=:black, lw = 4)
        lines!(axdzwb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdzwb, Point.(50, -450), text = "∂ₜ⟨b⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vbmax:vbmax:vbmax), size =35,vertical = false)
    cb2 = Colorbar(gcb2[1,1], hvb, ticks = (-wbmax/2:wbmax/2:wbmax/2), size =35,vertical = false, )
    cb3 = Colorbar(gcb3[1,1], hwb, ticks = (-nabmax/2:nabmax/2:nabmax/2), size =35,vertical = false, )

    cb3 = Colorbar(gcb1[1,1], hnab, ticks = (-nabmax/2:nabmax/2:nabmax/2),
            size =35,vertical = false, flipaxis=false)
    cb4 = Colorbar(gcb2[1,1], hdyvb, ticks = (-SGSmax/2:SGSmax/2:SGSmax/2),
            size =35,vertical = false, flipaxis=false)
    cb5 = Colorbar(gcb3[1,1], hdzwb, ticks = (-nabmax/2:nabmax/2:nabmax/2),
            size =35,vertical = false, flipaxis=false)

    colgap!(ga, 15)
    rowgap!(ga, 5)
    colgap!(gb, 15)
    rowgap!(gb, 5)
    colgap!(gc, 15)
    rowgap!(gc, 5)

    save(apath * savename * ".png", f1)
end

vmax = round(maximum(WaveAveragedVals_beg.v_Wavg)*1e2)*1e-2
nabmax = round(maximum(abs.(FluxDiv_phasedepWavg_beg.vb_phasedepWavg_dy))*1e7)*1e-7*0.25
SGSmax = round(maximum(FluxDiv_Wavg_beg.wb_xpertxavgWavg_dz)*1e7)*1e-7*0.25

function WaveAverageFluxPlot(savename, apath , ylength, v_Wavg, ∇_phasedepWavg, 
    ∇_turbWavg, SGS_Wavg, c_Wavg, vmax, nabmax, SGSmax)

    ###[ v      ][ ∇p     ][ ∇t        ]
    ###[ ∇p+∇t  ][ SGS    ][-∇p-∇t+SGS ]
    
    vlims = (-vmax,vmax)
    gradblims = (-nabmax, nabmax)
    SGSlims = (-SGSmax, SGSmax)

    Full∇_Wavg = ∇_turbWavg + ∇_phasedepWavg
    Fulldbdt_Wavg = (-1 .* Full∇_Wavg) + SGS_Wavg[2:end,2:end]

    yb_cut = yb[1:ylength]

    clog = log10.(clamp.(c_Wavg, 1e-8, 1))

    f1 = Figure(resolution = (1300, 1100), fontsize=26)
    ga = f1[2:3, 1] = GridLayout()    
    gb = f1[2:3, 2] = GridLayout()
    gc = f1[2:3, 3] = GridLayout()

    gcb1 = f1[1, 1] = GridLayout()       
    gcb2 = f1[1, 2] = GridLayout()       
    gcb3 = f1[1, 3] = GridLayout()       

    axvh = Axis(ga[1, 1], ylabel = "z [m]") #v
    axgrad = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
    axgrad.xticks = 500:500:1000
    axvh.yticks = [-250, 0]
    axgrad.yticks = [-250, 0]

    axvbh = Axis(gb[1, 1]) #vh'b'
    axdyvb = Axis(gb[2, 1], xlabel = "y [m]") #dy v'b'

    axdyvb.xticks = 500:500:1000

    axwbh = Axis(gc[1, 1]) #wh'b'
    axdzwb = Axis(gc[2, 1], xlabel = "y [m]") #dz w'b'
    axdzwb.xticks = 500:500:1000

    limits!(axvh, 0, 1500, -500, 0)
    limits!(axvbh, 0, 1500, -500, 0)
    limits!(axwbh, 0, 1500, -500, 0)
    limits!(axgrad, 0, 1500, -500, 0)
    limits!(axdyvb, 0, 1500, -500, 0)
    limits!(axdzwb, 0, 1500, -500, 0)

    hidedecorations!(axvbh)
    hidedecorations!(axwbh)
    hidexdecorations!(axvh)
    hideydecorations!(axdyvb)
    hideydecorations!(axdzwb)

    rowsize!(f1.layout, 1, Relative(0.05))

    hv = heatmap!(axvh, yb_cut, zb, v_Wavg, colormap = :balance, colorrange = vlims)
        lines!(axvh, yb, land, color=:black, lw = 4)
        lines!(axvh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvh, Point.(50, -450), text = "⟨v⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        contour!(axvh, yb_cut, zb, clog, color = :gray30, linewidth = 2, levels = -5:1:0)

    hnab = heatmap!(axgrad, yb_cut[2:end], zb[2:end], Full∇_Wavg, colormap = :balance, colorrange = gradblims)
        lines!(axgrad, yb, land, color=:black, lw = 4)
        lines!(axgrad, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axgrad, Point.(50, -450), text = "∇⋅⟨u⃗'b'⟩ + ∇⋅⟨u⃗̃b̃⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hnab1 = heatmap!(axvbh, yb_cut, zb, ∇_phasedepWavg, colormap = :balance, colorrange = gradblims)
        lines!(axvbh, yb, land, color=:black, lw = 4)
        lines!(axvbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvbh, Point.(50, -450), text = "∇⋅⟨u⃗̃b̃⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        contour!(axvbh, yb_cut, zb, clog, color = :gray30, linewidth = 2, levels = -5:1:0)

    hdyvb = heatmap!(axdyvb, yb_cut, zb, SGS_Wavg, colormap = :balance, colorrange = SGSlims)
        lines!(axdyvb, yb, land, color=:black, lw = 4)
        lines!(axdyvb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdyvb, Point.(50, -450), text = "SGS", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)        

    hnab2 = heatmap!(axwbh, yb_cut, zb, ∇_turbWavg, colormap = :balance, colorrange = gradblims)
        lines!(axwbh, yb, land, color=:black, lw = 4)
        lines!(axwbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axwbh, Point.(50, -450), text = "∇⋅⟨u⃗'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
        contour!(axwbh, yb_cut, zb, clog, color = :gray30, linewidth = 2, levels = -5:1:0)


    hdzwb = heatmap!(axdzwb, yb_cut, zb, Fulldbdt_Wavg, colormap = :balance, colorrange = gradblims)
        lines!(axdzwb, yb, land, color=:black, lw = 4)
        lines!(axdzwb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdzwb, Point.(50, -450), text = "∂ₜ⟨b⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax:vmax:vmax), size =35,vertical = false)
    cb2 = Colorbar(gcb2[1,1], hnab1, ticks = (-nabmax/2:nabmax/2:nabmax/2), size =35,vertical = false, )
    cb3 = Colorbar(gcb3[1,1], hnab2, ticks = (-nabmax/2:nabmax/2:nabmax/2), size =35,vertical = false, )

    cb3 = Colorbar(gcb1[1,1], hnab, ticks = (-nabmax/2:nabmax/2:nabmax/2),
            size =35,vertical = false, flipaxis=false)
    cb4 = Colorbar(gcb2[1,1], hdyvb, ticks = (-SGSmax/2:SGSmax/2:SGSmax/2),
            size =35,vertical = false, flipaxis=false)
    cb5 = Colorbar(gcb3[1,1], hdzwb, ticks = (-nabmax/2:nabmax/2:nabmax/2),
            size =35,vertical = false, flipaxis=false)

    colgap!(ga, 15)
    rowgap!(ga, 5)
    colgap!(gb, 15)
    rowgap!(gb, 5)
    colgap!(gc, 15)
    rowgap!(gc, 5)

    save(apath * savename * ".png", f1)
end

v_Wavg_begxavg = mean(WaveAveragedVals_beg.v_Wavg, dims = 1)[1,:,:];
∇_phasedepWavg_begxavg = mean(FluxDiv_phasedepWavg_beg.∇_phasedepWavg, dims = 1)[1,:,:];
∇_turbWavg_begxavg = mean(FluxDiv_turbWavg_beg.∇_turbWavg, dims = 1)[1,:,:];
c_Wavg_begxavg = mean(WaveAveragedVals_beg.c_Wavg, dims = 1)[1,:,:];
SGS_Wavg_begxavg = mean(WaveAveragedVals_beg.SGS_Wavg, dims = 1)[1,:,:];

v_Wavg_endxavg = mean(WaveAveragedVals_end.v_Wavg, dims = 1)[1,:,:];
∇_phasedepWavg_endxavg = mean(FluxDiv_phasedepWavg_end.∇_phasedepWavg, dims = 1)[1,:,:];
∇_turbWavg_endxavg = mean(FluxDiv_turbWavg_end.∇_turbWavg, dims = 1)[1,:,:];
c_Wavg_endxavg = mean(WaveAveragedVals_end.CGl_Wavg, dims = 1)[1,1:ylength,:];
SGS_Wavg_endxavg = mean(WaveAveragedVals_end.SGS_Wavg, dims = 1)[1,:,:];

v_Wavg_midxavg = mean(WaveAveragedVals_mid.v_Wavg, dims = 1)[1,:,:];
∇_phasedepWavg_midxavg = mean(FluxDiv_phasedepWavg_mid.∇_phasedepWavg, dims = 1)[1,:,:];
∇_turbWavg_midxavg = mean(FluxDiv_turbWavg_mid.∇_turbWavg, dims = 1)[1,:,:];
c_Wavg_midxavg = mean(WaveAveragedVals_mid.CGl_Wavg, dims = 1)[1,1:ylength,:];
SGS_Wavg_midxavg = mean(WaveAveragedVals_mid.SGS_Wavg, dims = 1)[1,:,:];

vmax = round(maximum(v_Wavg_begxavg)*1e2)*1e-2
nabmax1 = maximum((maximum(abs.(∇_phasedepWavg_begxavg)), maximum(abs.(∇_turbWavg_begxavg))))
nabmax = round(maximum(abs.(nabmax1))*1e7)*1e-7*0.25
SGSmax = round(maximum(SGS_Wavg_begxavg)*1e8)*1e-8*0.25

WaveAverageFluxPlot("blfuxtrip_Wavg_beg" * sn,  path_name, ylength, 
v_Wavg_begxavg, ∇_phasedepWavg_begxavg, ∇_turbWavg_begxavg,
SGS_Wavg_begxavg, c_Wavg_begxavg, vmax, nabmax, SGSmax)

WaveAverageFluxPlot("blfuxtrip_Wavg_end" * sn,  path_name, ylength, 
    v_Wavg_endxavg, ∇_phasedepWavg_endxavg, ∇_turbWavg_endxavg,
    SGS_Wavg_endxavg, c_Wavg_endxavg,
    vmax, nabmax, SGSmax)

WaveAverageFluxPlot("blfuxtrip_Wavg_mid" * sn, path_name, ylength, 
v_Wavg_midxavg, ∇_phasedepWavg_midxavg, ∇_turbWavg_midxavg,
SGS_Wavg_midxavg, c_Wavg_midxavg,
    vmax, nabmax, SGSmax)
    

#####################
# Taking a Volume Integral
#####################

zlength = length(zb)
yslopelength = 334

Ygrid = reshape(repeat(yb[2:ylength], zlength-1), ylength-1, zlength-1);
Zgrid = permutedims(reshape(repeat(zb[2:zlength], ylength-1), zlength-1, ylength-1), [2,1]);
SlopeGridY = curvedslope.(Ygrid);
LinSlopeGridY = linslope.(Ygrid);

Ygrid2 = reshape(repeat(yb[1:ylength], zlength), ylength, zlength);
Zgrid2 = permutedims(reshape(repeat(zb[1:zlength], ylength), zlength, ylength), [2,1]);
SlopeGridY2 = curvedslope.(Ygrid2);
LinSlopeGridY2 = linslope.(Ygrid2);
################
# \   \
#  \   \
#   \   \
#    \   \
#     \   \
#      \   \
#       \   \
#________\ __\_
#################

endslopedel = (500 +  pm.U₀/pm.Ñ)/pm.Tanα
endslope2del = (500 +  2*pm.U₀/pm.Ñ)/pm.Tanα

# all the values greater than slope
boolZY_aboveslope = ((SlopeGridY .+ 4) .<= zb[2:end]') .& (Ygrid .<= endslopedel) .& (Zgrid .<= 0) ;
boolZY_belowdelta = boolZY_aboveslope .& ((LinSlopeGridY .+ pm.U₀/pm.Ñ) .>= zb[2:end]');
boolZY_abovedelta = (((LinSlopeGridY .+ pm.U₀/pm.Ñ) .<= zb[2:end]') .& ((LinSlopeGridY .+ 2*pm.U₀/pm.Ñ) .>= zb[2:end]')) .& (Ygrid .<= endslope2del) .& (Zgrid .<= 0);

boolZY_aboveslope2 = ((SlopeGridY2 .+ 4) .<= zb[1:end]') .& (Ygrid2 .<= endslopedel) .& (Zgrid2 .<= 0) ;
boolZY_belowdelta2 = boolZY_aboveslope2 .& ((LinSlopeGridY2 .+ pm.U₀/pm.Ñ) .>= zb[1:end]');
boolZY_abovedelta2 = (((LinSlopeGridY2 .+ pm.U₀/pm.Ñ) .<= zb[1:end]') .& ((LinSlopeGridY2 .+ 2*pm.U₀/pm.Ñ) .>= zb[1:end]')) .& (Ygrid2 .<= endslope2del) .& (Zgrid2 .<= 0);

################
# \  \
#  \  \
#   \__\
#    \  \
#-250 \  \** Gaussian Released
#      \__\
#       \  \
#________\ _\_
#################

# on the slope
boolZY_bd_bot = boolZY_belowdelta .& (Zgrid .< -334) 
boolZY_bd_mid = boolZY_belowdelta .& ((Zgrid .>= -334) .& (Zgrid .< -168))
boolZY_bd_top = boolZY_belowdelta .& ((Zgrid .>= -168) .& (Zgrid .< -0))

# off teh slope
boolZY_ad_bot = boolZY_abovedelta .& (Zgrid .< -334) 
boolZY_ad_mid = boolZY_abovedelta .& ((Zgrid .>= -334) .& (Zgrid .< -168))
boolZY_ad_top = boolZY_abovedelta .& ((Zgrid .>= -168) .& (Zgrid .< -0))

# on the slope
boolZY_bd_bot2 = boolZY_belowdelta2 .& (Zgrid2 .< -334) 
boolZY_bd_mid2 = boolZY_belowdelta2 .& ((Zgrid2 .>= -334) .& (Zgrid2 .< -168))
boolZY_bd_top2 = boolZY_belowdelta2 .& ((Zgrid2 .>= -168) .& (Zgrid2 .< -0))

# off teh slope
boolZY_ad_bot2 = boolZY_abovedelta2 .& (Zgrid2 .< -334) 
boolZY_ad_mid2 = boolZY_abovedelta2 .& ((Zgrid2 .>= -334) .& (Zgrid2 .< -168))
boolZY_ad_top2 = boolZY_abovedelta2 .& ((Zgrid2 .>= -168) .& (Zgrid2 .< -0))

function volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_xpertxavgWavg)
    int_∇_xpertxavgWavg_bot = sum(∇_xpertxavgWavg[boolZY_bd_bot])
    int_∇_xpertxavgWavg_mid = sum(∇_xpertxavgWavg[boolZY_bd_mid])
    int_∇_xpertxavgWavg_top = sum(∇_xpertxavgWavg[boolZY_bd_top])

    boolZY_all = (boolZY_bd_top .* int_∇_xpertxavgWavg_top) .+ (boolZY_bd_mid .* int_∇_xpertxavgWavg_mid) .+ (boolZY_bd_bot .* int_∇_xpertxavgWavg_bot)

    vert_split_int_∇_xpertxavgWavg = (; int_∇_xpertxavgWavg_bot, int_∇_xpertxavgWavg_mid, int_∇_xpertxavgWavg_top)
    return vert_split_int_∇_xpertxavgWavg, boolZY_all
end

# on the slope
int_∇_phasedepWavgxavg_belowdelta_beg = sum(∇_phasedepWavg_begxavg[boolZY_belowdelta])
int_∇_phasedepWavgxavg_belowdelta_mid = sum(∇_phasedepWavg_midxavg[boolZY_belowdelta])
int_∇_phasedepWavgxavg_belowdelta_end = sum(∇_phasedepWavg_endxavg[boolZY_belowdelta])

int_∇_turbWavgxavg_belowdelta_beg = sum(∇_turbWavg_begxavg[boolZY_belowdelta])
int_∇_turbWavgxavg_belowdelta_mid = sum(∇_turbWavg_midxavg[boolZY_belowdelta])
int_∇_turbWavgxavg_belowdelta_end = sum(∇_turbWavg_endxavg[boolZY_belowdelta])

int_SGS_Wavg_belowdelta_beg = sum(SGS_Wavg_begxavg[boolZY_belowdelta2])
int_SGS_Wavg_belowdelta_mid = sum(SGS_Wavg_midxavg[boolZY_belowdelta2])
int_SGS_Wavg_belowdelta_end = sum(SGS_Wavg_endxavg[boolZY_belowdelta2])


vert_split_int_∇_phasedepWavgxavg_begT_bd, phasedep_boolZY_all_begT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_phasedepWavg_begxavg);
vert_split_int_∇_phasedepWavgxavg_midT_bd, phasedep_boolZY_all_midT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_phasedepWavg_midxavg);
vert_split_int_∇_phasedepWavgxavg_endT_bd, phasedep_boolZY_all_endT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_phasedepWavg_endxavg);

vert_split_int_∇_turbWavgxavg_begT_bd, turb_boolZY_all_begT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_turbWavg_begxavg);
vert_split_int_∇_turbWavgxavg_midT_bd, turb_boolZY_all_midT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_turbWavg_midxavg);
vert_split_int_∇_turbWavgxavg_endT_bd, turb_boolZY_all_endT_bd = volume_vertsplit(boolZY_bd_bot, boolZY_bd_mid, boolZY_bd_top, ∇_turbWavg_endxavg);

vert_split_int_SGS_Wavg_begT_bd, SGS_boolZY_all_begT_bd = volume_vertsplit(boolZY_bd_bot2, boolZY_bd_mid2, boolZY_bd_top2, SGS_Wavg_begxavg);
vert_split_int_SGS_Wavg_midT_bd, SGS_boolZY_all_midT_bd = volume_vertsplit(boolZY_bd_bot2, boolZY_bd_mid2, boolZY_bd_top2, SGS_Wavg_midxavg);
vert_split_int_SGS_Wavg_endT_bd, SGS_boolZY_all_endT_bd = volume_vertsplit(boolZY_bd_bot2, boolZY_bd_mid2, boolZY_bd_top2, SGS_Wavg_endxavg);

#_______________

vert_split_int_∇_phasedepWavgxavg_begT_ad, phasedep_boolZY_all_begT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, ∇_phasedepWavg_begxavg);
vert_split_int_∇_phasedepWavgxavg_midT_ad, phasedep_boolZY_all_midT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, ∇_phasedepWavg_midxavg);
vert_split_int_∇_phasedepWavgxavg_endT_ad, phasedep_boolZY_all_endT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, ∇_phasedepWavg_endxavg);

vert_split_int_∇_turbWavgxavg_begT_ad, turb_boolZY_all_begT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, ∇_turbWavg_begxavg);
vert_split_int_∇_turbWavgxavg_midT_ad, turb_boolZY_all_midT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, ∇_turbWavg_midxavg);
vert_split_int_∇_turbWavgxavg_endT_ad, turb_boolZY_all_endT_ad = volume_vertsplit(boolZY_ad_bot, boolZY_ad_mid, boolZY_ad_top, ∇_turbWavg_endxavg);

vert_split_int_SGS_Wavg_begT_ad, SGS_boolZY_all_begT_ad = volume_vertsplit(boolZY_ad_bot2, boolZY_ad_mid2, boolZY_ad_top2, SGS_Wavg_begxavg);
vert_split_int_SGS_Wavg_midT_ad, SGS_boolZY_all_midT_ad = volume_vertsplit(boolZY_ad_bot2, boolZY_ad_mid2, boolZY_ad_top2, SGS_Wavg_midxavg);
vert_split_int_SGS_Wavg_endT_ad, SGS_boolZY_all_endT_ad = volume_vertsplit(boolZY_ad_bot2, boolZY_ad_mid2, boolZY_ad_top2, SGS_Wavg_endxavg);

#_______________

phasedep_boolZY_all_begT = phasedep_boolZY_all_begT_bd .+ phasedep_boolZY_all_begT_ad
phasedep_boolZY_all_midT = phasedep_boolZY_all_midT_bd .+ phasedep_boolZY_all_midT_ad
phasedep_boolZY_all_endT = phasedep_boolZY_all_endT_bd .+ phasedep_boolZY_all_endT_ad

turb_boolZY_all_begT = turb_boolZY_all_begT_bd .+ turb_boolZY_all_begT_ad
turb_boolZY_all_midT = turb_boolZY_all_midT_bd .+ turb_boolZY_all_midT_ad
turb_boolZY_all_endT = turb_boolZY_all_endT_bd .+ turb_boolZY_all_endT_ad

SGS_boolZY_all_begT = SGS_boolZY_all_begT_bd .+ SGS_boolZY_all_begT_ad
SGS_boolZY_all_midT = SGS_boolZY_all_midT_bd .+ SGS_boolZY_all_midT_ad
SGS_boolZY_all_endT = SGS_boolZY_all_endT_bd .+ SGS_boolZY_all_endT_ad

divmax = maximum((maximum(phasedep_boolZY_all_begT), maximum(phasedep_boolZY_all_midT), maximum(phasedep_boolZY_all_endT)))
divmin = minimum((minimum(phasedep_boolZY_all_begT), minimum(phasedep_boolZY_all_midT), minimum(phasedep_boolZY_all_endT)))

divext = maximum((divmax, abs(divmin)))

divext = 1.2e-4

clog_beg = log10.(clamp.(c_Wavg_begxavg, 1e-8, 1))
clog_mid = log10.(clamp.(c_Wavg_midxavg, 1e-8, 1))
clog_end = log10.(clamp.(c_Wavg_endxavg, 1e-8, 1))

f1 = Figure(resolution = (2000, 1200), fontsize=26)
ga = f1[1, 1] = GridLayout() # beg
gb = f1[2, 1] = GridLayout() # mid
gc = f1[3, 1] = GridLayout() # end
gcb = f1[1:3, 2] = GridLayout()

axb_p = Axis(ga[1, 1], ylabel = "z [m]") #v
axm_p = Axis(gb[1, 1], ylabel = "z [m]") #vh
axe_p = Axis(gc[1, 1], ylabel = "z [m]", xlabel = "y [m]") #vh

axb_t = Axis(ga[1, 2]) #v
axm_t = Axis(gb[1, 2]) #vh
axe_t = Axis(gc[1, 2], xlabel = "y [m]") #vh

axb_tp = Axis(ga[1, 3]) #v
axm_tp = Axis(gb[1, 3]) #vh
axe_tp = Axis(gc[1, 3], ylabel = "z [m]", xlabel = "y [m]") #vh

axb_tps = Axis(ga[1, 4]) #v
axm_tps = Axis(gb[1, 4]) #vh
axe_tps = Axis(gc[1, 4], ylabel = "z [m]", xlabel = "y [m]") #vh


axe_p.xticks = 500:500:1500
axe_t.xticks = 500:500:1500
axe_tp.xticks = 500:500:1500
axe_tps.xticks = 500:500:1500

axb_p.yticks = [-250, 0]
axm_p.yticks = [-250, 0]
axe_p.yticks = [-250, 0]

limits!(axb_p, 0, 2000, -500, 0)
limits!(axm_p, 0, 2000, -500, 0)
limits!(axe_p, 0, 2000, -500, 0)
limits!(axb_t, 0, 2000, -500, 0)
limits!(axm_t, 0, 2000, -500, 0)
limits!(axe_t, 0, 2000, -500, 0)
limits!(axb_tp, 0, 2000, -500, 0)
limits!(axm_tp, 0, 2000, -500, 0)
limits!(axe_tp, 0, 2000, -500, 0)
limits!(axb_tps, 0, 2000, -500, 0)
limits!(axm_tps, 0, 2000, -500, 0)
limits!(axe_tps, 0, 2000, -500, 0)

hidexdecorations!(axb_p)
hidexdecorations!(axm_p)
hidedecorations!(axb_t)
hidedecorations!(axb_tp)
hidedecorations!(axb_tps)
hidedecorations!(axm_t)
hidedecorations!(axm_tp)
hidedecorations!(axm_tps)
hideydecorations!(axe_t)
hideydecorations!(axe_tp)
hideydecorations!(axe_tps)

colsize!(f1.layout, 2, Relative(0.05))

SGS_boolZY_all_begT_cut = SGS_boolZY_all_begT[2:end, 2:end]
SGS_boolZY_all_midT_cut = SGS_boolZY_all_midT[2:end, 2:end]
SGS_boolZY_all_endT_cut = SGS_boolZY_all_endT[2:end, 2:end]

hv = heatmap!(axb_p, yb[2:ylength], zb[2:zlength], -1 .* phasedep_boolZY_all_begT, colormap = :balance, colorrange = (-divext, divext))
    text!(axb_p, Point.(50, -450), text = "1-4 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)
    contour!(axb_p, yb[1:ylength], zb, clog_beg, color = :black, linewidth = 3, levels = -5:1:0)

hvhat = heatmap!(axm_p, yb[2:ylength], zb[2:zlength], -1 .* phasedep_boolZY_all_midT, colormap = :balance, colorrange = (-divext, divext))
    text!(axm_p, Point.(50, -450), text = "4-7 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)
    contour!(axm_p, yb[1:ylength], zb, clog_mid, color = :black, linewidth = 3, levels = -5:1:0)

hnab = heatmap!(axe_p, yb[2:ylength], zb[2:zlength], -1 .* phasedep_boolZY_all_endT, colormap = :balance, colorrange = (-divext, divext))
    text!(axe_p, Point.(50, -450), text = "7-11 Tσ", align = (:left, :center), color = :black, 
    font = :bold, fontsize = 26)
    contour!(axe_p, yb[1:ylength], zb, clog_end, color = :black, linewidth = 3, levels = -5:1:0)

heatmap!(axb_t, yb[2:ylength], zb[2:zlength], -1 .* turb_boolZY_all_begT, colormap = :balance, colorrange = (-divext, divext))
heatmap!(axm_t, yb[2:ylength], zb[2:zlength], -1 .* turb_boolZY_all_midT, colormap = :balance, colorrange = (-divext, divext))
heatmap!(axe_t, yb[2:ylength], zb[2:zlength], -1 .* turb_boolZY_all_endT, colormap = :balance, colorrange = (-divext, divext))

heatmap!(axb_tp, yb[2:ylength], zb[2:zlength], -1 .* (turb_boolZY_all_begT .+ phasedep_boolZY_all_begT), colormap = :balance, colorrange = (-divext, divext))
heatmap!(axm_tp, yb[2:ylength], zb[2:zlength], -1 .* (turb_boolZY_all_midT .+ phasedep_boolZY_all_midT), colormap = :balance, colorrange = (-divext, divext))
heatmap!(axe_tp, yb[2:ylength], zb[2:zlength], -1 .* (turb_boolZY_all_endT .+ phasedep_boolZY_all_endT), colormap = :balance, colorrange = (-divext, divext))

heatmap!(axb_tps, yb[2:ylength], zb[2:zlength], -1 .* (turb_boolZY_all_begT .+ phasedep_boolZY_all_begT) .+ SGS_boolZY_all_begT_cut, colormap = :balance, colorrange = (-divext, divext))
heatmap!(axm_tps, yb[2:ylength], zb[2:zlength], -1 .* (turb_boolZY_all_midT .+ phasedep_boolZY_all_midT) .+ SGS_boolZY_all_begT_cut, colormap = :balance, colorrange = (-divext, divext))
heatmap!(axe_tps, yb[2:ylength], zb[2:zlength], -1 .* (turb_boolZY_all_endT .+ phasedep_boolZY_all_endT) .+ SGS_boolZY_all_begT_cut, colormap = :balance, colorrange = (-divext, divext))

# create colorbars the size of the whole data set
cb1 = Colorbar(gcb[1,1], hnab, ticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0",  "-8×10⁻⁵", "-4×10⁻⁵"]), size =35, label = "∂ₜb terms")

colgap!(ga, 15)
colgap!(gb, 15)
colgap!(gc, 15)

Label(ga[1, 1, Top()], "-∇⋅⟨u⃗̃b̃⟩",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)
Label(ga[1, 2, Top()], "-∇⋅⟨u⃗'b'⟩",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)
Label(ga[1, 3, Top()], "-∇⋅⟨u⃗̃b̃⟩-∇⋅⟨u⃗'b'⟩",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)
Label(ga[1, 4, Top()], "-∇⋅⟨u⃗̃b̃⟩-∇⋅⟨u⃗'b'⟩+ ⟨∇⋅κ∇b⟩",
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :left)

savename = path_name * "IntegratesFlux_TripDecomp_" * sn
save(savename * ".png", f1)


#################
#          PHASE BREAKDOWN
################

FluxDiv_turbWavg_beg  = scale_file1[skeys[1]];
FluxDiv_phasedepWavg_beg  = scale_file1[skeys[2]];

v_Phdep_begxavg = mean(PhaseDependentVals_beg.v_Phdep, dims = 1)[1,:,:,:];
w_Phdep_begxavg = mean(PhaseDependentVals_beg.w_Phdep, dims = 1)[1,:,:,:];
b_Phdep_begxavg = mean(PhaseDependentVals_beg.b_Phdep, dims = 1)[1,:,:,:];
c_Phdep_begxavg = mean(PhaseDependentVals_beg.c_Phdep, dims = 1)[1,:,:,:];
∇_phasedepWavg_begxavg = mean(FluxDiv_phasedepWavg_beg.∇_phasedepWavg, dims = 1)[1,:,:];
∇_turbWavg_begxavg = mean(FluxDiv_turbWavg_beg.∇_turbWavg, dims = 1)[1,:,:];

v_turb_begxavg = mean(TurbVals_beg.v_turb, dims = 1)[1,:,:,:, :];
w_turb_begxavg = mean(TurbVals_beg.w_turb, dims = 1)[1,:,:,:, :];
b_turb_begxavg = mean(TurbVals_beg.b_turb, dims = 1)[1,:,:,:, :];

vb_Phdep = v_Phdep_begxavg .* b_Phdep_begxavg;
wb_Phdep = w_Phdep_begxavg .* b_Phdep_begxavg;

vb_turb = mean(v_turb_begxavg .* b_turb_begxavg, dims = 4);
wb_turb = mean(w_turb_begxavg .* b_turb_begxavg, dims = 4);

zlength = length(zb)
ylength = 880

Ygrid = reshape(repeat(yb[2:ylength], zlength-1), ylength-1, zlength-1)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

Δy = yb[1]-yb[2]
Δz = zb[1]-zb[2]

wb_turb_dz = ((wb_turb[2:end,1:end-1, :] .- wb_turb[2:end,2:end, :])./Δz) .* boolZY;
vb_turb_dy = ((vb_turb[1:end-1,2:end, :] .- vb_turb[2:end,2:end, :])./Δy) .* boolZY;
∇_turb_Phavg = wb_turb_dz .+ vb_turb_dy;

wb_Phdep_dz = ((wb_Phdep[2:end,1:end-1, :] .- wb_Phdep[2:end,2:end, :])./Δz) .* boolZY;
vb_Phdep_dy = ((vb_Phdep[1:end-1,2:end, :] .- vb_Phdep[2:end,2:end, :])./Δy) .* boolZY;
∇_Phdep = wb_Phdep_dz .+ vb_Phdep_dy;

∇_PhDepFluxTerm =  ∇_phasedepWavg_begxavg .- ∇_Phdep;
∇_turbFluxTerm =  ∇_turb_Phavg .- ∇_turbWavg_begxavg;
∇_RHS = ∇_PhDepFluxTerm .- ∇_turbFluxTerm;

tindxs = vcat(vcat(1:2:7, 8),9:2:13)
yb_cut = yb[1:ylength]

function PhaseAverageFluxPlot(savename, apath, tindxs, v_Phdep_begxavg, w_Phdep_begxavg, ∇_PhDepFluxTerm, 
    ∇_turbFluxTerm, ∇_RHS, b_Phdep_begxavg, c_Phdep_begxavg)

    f = Figure(resolution = (2200, 1400), fontsize=26)
    ga = f[1, 1:8] = GridLayout() # v
    gb = f[2, 1:8] = GridLayout() # w
    gc = f[3, 1:8] = GridLayout() # grad
    gbd = f[4, 1:8] = GridLayout() # grad
    gcd = f[5, 1:8] = GridLayout() # grad
    ge = f[6, 1:8] = GridLayout() # b
    gd = f[7, 1:8] = GridLayout() # c

    # v'b'
    gcb1 = f[1:2, 9] = GridLayout()
    # w'b'
    gcb2 = f[3:5, 9] = GridLayout()
    gcb3 = f[6:7, 9] = GridLayout()

    gplots = [ga, gb, gc, gbd, gcd, ge, gd]

    for g in gplots
        if g == gd
            for m = 1:8
                if m==1
                    axi = Axis(g[1, m], ylabel = "z [m]", xlabel = "y [m]") 
                    global ax = axi
                    axi.xticks = 500:500:1000
                    axi.yticks = [-250, 0]
                else
                    axi = Axis(g[1,m], xlabel = "y [m]")
                    hideydecorations!(axi)
                    global ax = hcat(ax, axi)
                    axi.xticks = 500:500:1000
                end
            end
            global axs = vcat(axs, ax)
        else
            for m = 1:8
                if m==1
                    axi = Axis(g[1, m], ylabel = "z [m]")
                    global ax = axi
                    hidexdecorations!(axi)
                    axi.yticks = [-250, 0]
                else
                    axi = Axis(g[1, m])
                    global ax = hcat(ax, axi)
                    hidedecorations!(axi)
                end
            end
            if g == ga
                global axs = ax
            else
                global axs = vcat(axs, ax)
            end
        end
        rowgap!(g, 5)
        colgap!(g, 5)
    end

    for j = 1:length(axs)
        limits!(axs[j], 0, 1500, -500, 0)
    end

    vmax = 0.3
    ∇max = round(maximum(abs.(∇_PhDepFluxTerm))*0.1*1e6)*1e-6
    bmax = 0.0011

    vlims = (-vmax,vmax)
    nabblims = (-∇max, ∇max)
    blims = (-bmax, bmax)

    time_pre = "t = "
    time_post = " Tσ"

    zb_cut = zb[2:end]
    yb_cut2 = yb_cut[2:end]

    for (m,idx) in enumerate(tindxs)
        t = phase_times[idx]
        phaselabel = time_pre * @sprintf("%0.2f", t) * time_post

        v_Phdep_begxavgi = v_Phdep_begxavg[:,:,idx]
        w_Phdep_begxavgi = w_Phdep_begxavg[:,:,idx]
        ∇_PhDepFluxTermi = ∇_PhDepFluxTerm[:,:,idx]
        ∇_turbFluxTermi = ∇_turbFluxTerm[:,:,idx]
        ∇_RHSi = ∇_RHS[:,:,idx]
        b_Phdep_begxavgi = b_Phdep_begxavg[:,:,idx]
        c_Phdep_begxavgi  = log10.(clamp.(c_Phdep_begxavg[:,:,idx], 1e-12,1))
        
        global hv = heatmap!(axs[1, m], yb_cut, zb, v_Phdep_begxavgi, colormap = :balance, colorrange = vlims)
        lines!(axs[1, m], yb, land, color=:black, lw = 4)
        lines!(axs[1, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)

        global hw = heatmap!(axs[2, m], yb_cut, zb, w_Phdep_begxavgi, colormap = :balance, colorrange = vlims)
        lines!(axs[2, m], yb, land, color=:black, lw = 4)
        lines!(axs[2, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        
        global hwb = heatmap!(axs[3, m], yb_cut2, zb_cut, ∇_PhDepFluxTermi, colormap = :balance, colorrange = nabblims)
        lines!(axs[3, m], yb, land, color=:black, lw = 4)
        lines!(axs[3, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)

        global hdyvb = heatmap!(axs[4, m], yb_cut2, zb_cut, ∇_turbFluxTermi, colormap = :balance, colorrange = nabblims)
        lines!(axs[4, m], yb, land, color=:black, lw = 4)
        lines!(axs[4, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)

        global hdzwb = heatmap!(axs[5, m], yb_cut2, zb_cut, ∇_RHSi, colormap = :balance, colorrange = nabblims)
        lines!(axs[5, m], yb, land, color=:black, lw = 4)
        lines!(axs[5, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)

        global hb = heatmap!(axs[6, m], yb_cut, zb, b_Phdep_begxavgi, colormap = :balance, colorrange = blims)
        lines!(axs[6, m], yb, land, color=:black, lw = 4)
        lines!(axs[6, m], yb[1:382], land_pdel, color=:white, lw = 4, linestyle = :dash)

        global hc = heatmap!(axs[7, m], yb_cut, zb, c_Phdep_begxavgi, colormap = :thermal, colorrange = (-6,0))
        lines!(axs[end, m], yb, land, color=:white, lw = 4)
        lines!(axs[end, m], yb[1:382], land_pdel, color=:white, lw = 4, linestyle = :dash)
        
        Label(f[1, m, Top()], phaselabel, valign = :bottom,
        font = :bold, 
        padding = (0, 0, 5, 0))
    end

    text!(axs[1, 1], Point.(50, -350), text = "ṽ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[2, 1], Point.(50, -350), text = "w̃", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[3, 1], Point.(50, -350), text = "∇⋅([̄ũb̃]̄ - ũb̃)", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[4, 1], Point.(50, -350), text = "∇⋅(⟨u'b'⟩ - [̄u'b']̄)", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[5, 1], Point.(50, -350), text = "RHS", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[6, 1], Point.(50, -350), text = "b̃", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[7, 1], Point.(50, -350), text = "log₁₀c̃", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax/2:vmax/2:vmax/2, ["-0.15", "0", "0.15"] ),
    size =25, label = "velocity")
    cb2 = Colorbar(gcb2[1,1], hwb, ticks = (-∇max/2:∇max/2:∇max/2,  ["-5×10⁻⁵", "0", "5×10⁻⁵"] ), 
    size =25, label = "Fluxes")
    cb3 = Colorbar(gcb3[1,1], hb, ticks = (-bmax/2:bmax/2:bmax/2,  ["-5.5×10⁻⁴", "0", "5.5×10⁻⁴"] ), 
    size =25, label = "Buoyancy")


    colsize!(f.layout, 9, Relative(0.02))

    Label(f[1, 1:8, Top()],"Phase Triple Decomp, Beginning", valign = :bottom,
        font = :bold, fontsize = 30,
        padding = (0, 0, 50, 0))

    save(apath * savename * ".png", f)

    save(path_name * "TripPhaseDecomp_beg_" * sn * ".png", f )
end
