using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using ArgParse
using CairoMakie

dpath = "Data/"
apath = "Analysis/Plots/"

sn = "U300N100Lz100g100"
#filescalename1 = dpath * "BFluxFull_" * sn * ".jld2"
#filescalename2 = dpath * "cPhavg_" * sn * ".jld2"
filescalename1 = dpath * "BFlux_beg_" * sn * ".jld2"
filescalename1 = dpath * "BFlux_mid_" * sn * ".jld2"
filescalename1 = dpath * "BFlux_end_" * sn * ".jld2"

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

skeys = keys(scale_file1)

Phavg_idx = findall(contains("Phavg"), skeys)

vb_WpertPhavg  = scale_file1[skeys[Phavg_idx[1]]];
wb_WpertPhavg  = scale_file1[skeys[Phavg_idx[2]]];
v̂b_WpertPhavg  = scale_file1[skeys[Phavg_idx[3]]];
ŵb_WpertPhavg  = scale_file1[skeys[Phavg_idx[4]]];

v_Phavg  = scale_file1[skeys[Phavg_idx[5]]];
v̂_Phavg  = scale_file1[skeys[Phavg_idx[6]]];
ŵ_Phavg  = scale_file1[skeys[Phavg_idx[7]]];
b_Phavg  = scale_file1[skeys[Phavg_idx[8]]];
e_Phavg  = scale_file1[skeys[Phavg_idx[9]]];

#c_Phavg_xavg = scale_file2["c_Phavg_xavg"]
#c_Phavg_xavg_b = scale_file2["c_Phavg_xavg_b"]

yb=scale_file1[skeys[end-2]]
zb=scale_file1[skeys[end-1]]
phase_times= scale_file1[skeys[end]]

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

""" Phase averaging """
v̂_Phavg_xavg = mean(v̂_Phavg, dims = 1)[1,:,:,:];
v_Phavg_xavg = mean(v_Phavg, dims = 1)[1,:,:,:];
ŵ_Phavg_xavg = mean(ŵ_Phavg, dims = 1)[1,:,:,:];

vb_WpertPhavg_xavg = mean(vb_WpertPhavg, dims = 1)[1,:,:,:];
wb_WpertPhavg_xavg = mean(wb_WpertPhavg, dims = 1)[1,:,:,:];
ŵb_WpertPhavg_xavg = mean(ŵb_WpertPhavg, dims = 1)[1,:,:,:];
v̂b_WpertPhavg_xavg = mean(v̂b_WpertPhavg, dims = 1)[1,:,:,:];
b_Phavg_xavg = mean(b_Phavg, dims =1)[1,:,:,:];
e_Phavg_xavg = mean(e_Phavg, dims =1)[1,:,:,:];

zlength = length(zb)
Ygrid = reshape(repeat(yb[2:ylength], zlength-1), ylength-1, zlength-1)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

wb_WpertPhavg_xavg_dz = ((wb_WpertPhavg_xavg[2:end,1:end-1,:] .- wb_WpertPhavg_xavg[2:end,2:end,:])./(zb[1]-zb[2])).* boolZY;
vb_WpertPhavg_xavg_dy = ((vb_WpertPhavg_xavg[1:end-1,2:end,:] .- wb_WpertPhavg_xavg[2:end,2:end,:])./(yb[1]-yb[2])).* boolZY;
∇WpertPhavg_xavg = wb_WpertPhavg_xavg_dz .+ vb_WpertPhavg_xavg_dy;

(_,_,Lw) = size(b_Phavg_xavg)

yb_cut = yb[1:ylength]

function PhaseAverageFluxPlot(savename, apath, tindxs, v_Phavg_xavg, vb_WpertPhavg_xavg, 
    wb_WpertPhavg_xavg, vb_WpertPhavg_xavg_dy, wb_WpertPhavg_xavg_dz, ∇WpertPhavg_xavg, 
    b_Phavg_xavg)

    #f = Figure(resolution = (2500, 1800), fontsize=26)
    f = Figure(resolution = (2200, 1400), fontsize=26)
    ga = f[1, 1:8] = GridLayout() # vhat
    gb = f[2, 1:8] = GridLayout() # what
    gc = f[3, 1:8] = GridLayout() #<v'b'>
    gbd = f[4, 1:8] = GridLayout() #<w'b'>
    gcd = f[5, 1:8] = GridLayout() # grad
    ge = f[6, 1:8] = GridLayout() # diss
    gd = f[end+1, 1:8] = GridLayout()

    # v'b'
    gcb1 = f[1:3, 9] = GridLayout()
    # w'b'
    #gcb2 = f[5:end, 9] = GridLayout()
    gcb2 = f[4:end, 9] = GridLayout()

    gplots = [ga, gb, gc, gbd, gcd, gd]

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

    vmax = pm.U₀
    vbmax = round(maximum(abs.(vb_WpertPhavg_xavg))*0.5*1e5)*1e-5
    wbmax = round(maximum(abs.(wb_WpertPhavg_xavg))*0.5*1e5)*1e-5
    nabmax = round(maximum(abs.(vb_WpertPhavg_xavg_dy))*0.5*1e5)*1e-5
    dzmax = maximum(wb_WpertPhavg_xavg_dz)*0.1 < 1e-6 ? round(maximum(wb_WpertPhavg_xavg_dz)*1e7*0.1)*1e-7 : round(maximum(wb_WpertPhavg_xavg_dz)*1e6*0.1)*1e-6

    vbmax2 = string(round(maximum(abs.(vb_WpertPhavg_xavg))*0.5*1e5)/2)
    wbmax2 = string(round(maximum(abs.(wb_WpertPhavg_xavg))*0.5*1e5)/2)
    nabmax2 = string(round(maximum(abs.(vb_WpertPhavg_xavg_dy))*0.5*1e5)/2)
    dzmax2 = maximum(wb_WpertPhavg_xavg_dz)*0.1 < 1e-6 ? string(round(maximum(wb_WpertPhavg_xavg_dz)*1e7*0.1)/2) : string(round(maximum(wb_WpertPhavg_xavg_dz)*1e6*0.1)/2)
    dzmaxender = maximum(wb_WpertPhavg_xavg_dz)*0.1 < 1e-6 ? "10⁻⁷" : "10⁻⁶"

    vblims = (-vbmax,vbmax)
    wblims = (-wbmax,wbmax)
    gradblims = (-nabmax, nabmax)
    dwblims = (-dzmax, dzmax)

    time_pre = "t = "
    time_post = " Tσ"

    zb_cut = zb[2:end]
    yb_cut2 = yb_cut[2:end]

    bmin = round(pm.Ñ^2*500*1e3)*1e-3
    bstep = bmin/8

    for (m,idx) in enumerate(tindxs)
        t = phase_times[idx]
        phaselabel = time_pre * @sprintf("%0.2f", t) * time_post

        v_Phavg_xavgi = v_Phavg_xavg[:,:,idx]
        vb_WpertPhavg_xavgi = vb_WpertPhavg_xavg[:,:,idx]
        wb_WpertPhavg_xavgi = wb_WpertPhavg_xavg[:,:,idx]
        vb_WpertPhavg_xavg_dyi = vb_WpertPhavg_xavg_dy[:,:,idx]
        wb_WpertPhavg_xavg_dzi = wb_WpertPhavg_xavg_dz[:,:,idx]
        b_Phavg_xavgi = b_Phavg_xavg[:,:,idx]
        ∇WpertPhavg_xavgi  = ∇WpertPhavg_xavg[:,:,idx]
        #e_Phavg_xavgi = e_Phavg_xavg[:,:,idx]
        
        hv = heatmap!(axs[1, m], yb_cut, zb, v_Phavg_xavgi, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axs[1, m], yb, land, color=:black, lw = 4)
        lines!(axs[1, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[1, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hvb = heatmap!(axs[2, m], yb_cut, zb, vb_WpertPhavg_xavgi, colormap = :balance, colorrange = vblims)
        lines!(axs[2, m], yb, land, color=:black, lw = 4)
        lines!(axs[2, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[2, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)
        
        global hwb = heatmap!(axs[3, m], yb_cut, zb, wb_WpertPhavg_xavgi, colormap = :balance, colorrange = wblims)
        lines!(axs[3, m], yb, land, color=:black, lw = 4)
        lines!(axs[3, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[3, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hdyvb = heatmap!(axs[4, m], yb_cut2, zb_cut, vb_WpertPhavg_xavg_dyi, colormap = :balance, colorrange = gradblims)
        lines!(axs[4, m], yb, land, color=:black, lw = 4)
        lines!(axs[4, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[4, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hdzwb = heatmap!(axs[5, m], yb_cut2, zb_cut, wb_WpertPhavg_xavg_dzi, colormap = :balance, colorrange = dwblims)
        lines!(axs[5, m], yb, land, color=:black, lw = 4)
        lines!(axs[5, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[5, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        #heps = heatmap!(axs[6, m], yb_cut, zb, log10.(clamp.(e_Phavg_xavgi, 1e-14, 1)), colormap = :thermal, colorrange = (-8,-4))
        #lines!(axs[6, m], yb, land, color=:black, lw = 4)
        #lines!(axs[6, m], yb[1:382], land_pdel, color=:white, lw = 4, linestyle = :dash)

        hgradub = heatmap!(axs[end, m], yb_cut2, zb_cut, ∇WpertPhavg_xavgi, colormap = :balance, colorrange = gradblims)
        lines!(axs[end, m], yb, land, color=:black, lw = 4)
        lines!(axs[end, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[end, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)
        
        Label(f[1, m, Top()], phaselabel, valign = :bottom,
        font = :bold, 
        padding = (0, 0, 5, 0))
    end

    text!(axs[1, 1], Point.(50, -350), text = "v", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[2, 1], Point.(50, -350), text = "⟨v'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[3, 1], Point.(50, -350), text = "⟨w'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[4, 1], Point.(50, -350), text = "∂y⟨v'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axs[5, 1], Point.(50, -350), text = "∂z⟨w'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    #text!(axs[6, 1], Point.(50, -350), text = "log₁₀⟨ε⟩", align = (:left, :center), color = :white, 
    #    font = :bold, fontsize = 26)
    text!(axs[end, 1], Point.(50, -350), text = "∇⋅⟨u⃗'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hvb, ticks = (-vbmax/2:vbmax/2:vbmax/2, ["-" * vbmax2, "0", vbmax2] ),
    size =25, flipaxis=false, label = "⟨v'b'⟩ × 10⁻⁵")
    cb2 = Colorbar(gcb1[1,1], hwb, ticks = (-wbmax/2:wbmax/2:wbmax/2, ["-" * wbmax2, "0", wbmax2] ), 
    size =25, label = "⟨w'b'⟩ × 10⁻⁵")
    cb3 = Colorbar(gcb2[1,1], hdyvb, ticks = (-nabmax/2:nabmax/2:nabmax/2,  ["-" * nabmax2, "0", nabmax2] ), 
    size =25, flipaxis=false, label = "∂y⟨v'b'⟩, ∇⋅⟨u⃗'b'⟩ × 10⁻⁵")
    cb3 = Colorbar(gcb2[1,1], hdzwb, ticks = (-dzmax/2:dzmax/2:dzmax/2,  ["-" * dzmax2, "0", dzmax2] ),
    size =25, label = "∂z⟨w'b'⟩ × " * dzmaxender)

    colsize!(f.layout, 9, Relative(0.02))

    #Label(f[1, 1:8, Top()],big_titlep, valign = :bottom,
    #    font = :bold, fontsize = 30,
    #    padding = (0, 0, 50, 0))

    save(apath * savename * ".png", f)

end

PhaseAverageFluxPlot("blfux_Phavg_" * sn, path_name * "Analysis/", vcat(vcat(1:2:7, 8),9:2:13), v_Phavg_xavg, vb_WpertPhavg_xavg, 
    wb_WpertPhavg_xavg, vb_WpertPhavg_xavg_dy, wb_WpertPhavg_xavg_dz, ∇WpertPhavg_xavg, b_Phavg_xavg)

PhaseAverageFluxPlot("blfux_Phavg_b_" * sn, path_name * "Analysis/", vcat(vcat(1:2:7, 8),9:2:13), v_Phavg_xavg_b, vb_WpertPhavg_xavg_b, 
    wb_WpertPhavg_xavg_b, vb_WpertPhavg_xavg_dy_b, wb_WpertPhavg_xavg_dz_b, ∇WpertPhavg_xavg_b, b_Phavg_xavg_b)

     
function PhaseAverageFluxPaperPlot(savename, apath, tindxs, v̂_Phavg_xavg, ŵ_Phavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, vb_WpertPhavg_xavg_dy, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg)

    #f = Figure(resolution = (2500, 1800), fontsize=26)
    f = Figure(resolution = (2200, 1400), fontsize=26)
    ga = f[1, 1:5] = GridLayout() # vhat
    gb = f[2, 1:5] = GridLayout() # what
    gc = f[3, 1:5] = GridLayout() #<v'b'>
    gd = f[4, 1:5] = GridLayout() #<w'b'>
    ge = f[5, 1:5] = GridLayout() # grad
    gf = f[6, 1:5] = GridLayout() # diss

    # v
    gcb1 = f[1:2, 6] = GridLayout()
    # w

    #v'b', grad
    gcb2 = f[3:4, 6] = GridLayout()
    #w'b'

    # dissip
    gcb3 = f[5:6, 6] = GridLayout()

    gplots = [ga, gb, gc, gd, ge, gf]

    for g in gplots
        if g == gf
            for m = 1:5
                if m==1
                    axi = Axis(g[1, m], ylabel = "z [m]", xlabel = "y [m]") 
                    global ax = axi
                    axi.xticks = 500:500:1000
                    axi.yticks = [-250, 0]
                else
                    axi = Axis(g[1,m])
                    hideydecorations!(axi)
                    global ax = hcat(ax, axi)
                    axi.xticks = 500:500:1000
                end
            end
            global axs = vcat(axs, ax)
        else
            for m = 1:5
                if m==1
                    axi = Axis(g[1, m])
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

    vmax = pm.U₀
    wmax = round(maximum(abs.(ŵ_Phavg))*0.8*100)*1e-2

    vbmax = round(maximum(abs.(vb_WpertPhavg_xavg))*0.5*1e5)*1e-5
    wbmax = round(maximum(abs.(wb_WpertPhavg_xavg))*0.5*1e5)*1e-5
    nabmax = round(maximum(abs.(vb_WpertPhavg_xavg_dy))*0.5*1e5)*1e-5

    vmax2 = string(round(maximum(abs.(0.2))*100)*1e-2/2)
    wmax2 = string(round(maximum(abs.(ŵ_Phavg))*0.8*100)*1e-2/2)
    vbmax2 = string(round(maximum(abs.(vb_WpertPhavg_xavg))*0.5*1e5)/2)
    wbmax2 = string(round(maximum(abs.(wb_WpertPhavg_xavg))*0.5*1e5)/2)
    nabmax2 = string(round(maximum(abs.(vb_WpertPhavg_xavg_dy))*0.5*1e5)/2)

    vblims = (-vbmax,vbmax)
    wblims = (-wbmax,wbmax)
    gradblims = (-nabmax, nabmax)

    time_pre = "t = "
    time_post = " Tσ"

    zb_cut = zb[2:end]
    yb_cut2 = yb_cut[2:end]

    bmin = round(pm.Ñ^2*500*1e3)*1e-3
    bstep = bmin/8

    for (m,idx) in enumerate(tindxs)
        t = phase_times[idx]
        phaselabel = time_pre * @sprintf("%0.2f", t) * time_post

        v̂_Phavg_xavgi = v̂_Phavg_xavg[:,:,idx]
        ŵ_Phavgi = ŵ_Phavg[:,:,idx]

        vb_WpertPhavg_xavgi = vb_WpertPhavg_xavg[:,:,idx]
        wb_WpertPhavg_xavgi = wb_WpertPhavg_xavg[:,:,idx]
        b_Phavg_xavgi = b_Phavg_xavg[:,:,idx]
        ∇WpertPhavg_xavgi  = ∇WpertPhavg_xavg[:,:,idx]
        e_Phavg_xavgi = e_Phavg_xavg[:,:,idx]
        
        global hv = heatmap!(axs[1, m], yb_cut, zb, v̂_Phavg_xavgi, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axs[1, m], yb, land, color=:black, lw = 4)
        lines!(axs[1, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[1, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hw = heatmap!(axs[2, m], yb_cut, zb, ŵ_Phavgi, colormap = :balance, colorrange = (-wmax, wmax))
        lines!(axs[2, m], yb, land, color=:black, lw = 4)
        lines!(axs[2, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[2, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hvb = heatmap!(axs[3, m], yb_cut, zb, vb_WpertPhavg_xavgi, colormap = :balance, colorrange = vblims)
        lines!(axs[3, m], yb, land, color=:black, lw = 4)
        lines!(axs[3, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[3, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)
        
        global hwb = heatmap!(axs[4, m], yb_cut, zb, wb_WpertPhavg_xavgi, colormap = :balance, colorrange = wblims)
        lines!(axs[4, m], yb, land, color=:black, lw = 4)
        lines!(axs[4, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[4, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hgradub = heatmap!(axs[5, m], yb_cut2, zb_cut, ∇WpertPhavg_xavgi, colormap = :balance, colorrange = gradblims)
        lines!(axs[5, m], yb, land, color=:black, lw = 4)
        lines!(axs[5, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[5, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global heps = heatmap!(axs[6, m], yb_cut, zb, log10.(clamp.(e_Phavg_xavgi, 1e-14, 1)), colormap = :thermal, colorrange = (-8,-5))
        lines!(axs[6, m], yb, land, color=:black, lw = 4)
        lines!(axs[6, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        
        Label(f[1, m, Top()], phaselabel, valign = :bottom,
        font = :bold, fontsize = 36,
        padding = (0, 0, 5, 0))
    end

    text!(axs[1, 1], Point.(50, -350), text = "⟨v̂⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[2, 1], Point.(50, -350), text = "⟨ŵ⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[3, 1], Point.(50, -350), text = "⟨v'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[4, 1], Point.(50, -350), text = "⟨w'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[5, 1], Point.(50, -350), text = "∇⋅⟨u⃗'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[6, 1], Point.(50, -350), text = "⟨ε⟩ₚₕ", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 36)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax/2:vmax/2:vmax/2, ["-" * vmax2, "0", vmax2] ),
    size =25, flipaxis=false, label = "⟨v̂⟩ₚₕ", labelsize = 30)
    cb2 = Colorbar(gcb1[1,2], hw, ticks = (-wmax/2:wmax/2:wmax/2, ["-" * wmax2, "0", wmax2] ), 
    size =25, label = "⟨ŵ⟩ₚₕ", labelsize = 30)

    cb3 = Colorbar(gcb2[1,1], hvb, ticks = (-vbmax/2:vbmax/2:vbmax/2, ["-" * vbmax2, "0", vbmax2] ),
    size =25, flipaxis=false, label = "⟨v'b'⟩ × 10⁻⁵", labelsize = 30)
    cb4 = Colorbar(gcb2[1,2], hwb, ticks = (-wbmax/2:wbmax/2:wbmax/2, ["-" * wbmax2, "0", wbmax2] ), 
    size =25, label = "⟨w'b'⟩ × 10⁻⁵", labelsize = 30)

    cb5 = Colorbar(gcb3[1,1], hgradub, ticks = (-nabmax/2:nabmax/2:nabmax/2,  ["-" * nabmax2, "0", nabmax2] ), 
    size =25, flipaxis=false, label = "∇⋅⟨u⃗'b'⟩ₚₕ  × 10⁻⁵", labelsize = 30)
    cb6 = Colorbar(gcb3[1,2], heps, ticks = (-8:1:-6, ["10⁻⁸", "10⁻⁷", "10⁻⁶"] ), 
    size =25, label = "⟨ε⟩ₚₕ" , labelsize = 30)

    colgap!(gcb1, 10)
    colgap!(gcb2, 10)
    colgap!(gcb3, 10)

    colsize!(f.layout, 6, Relative(0.02))

    save(apath * savename * ".png", f)
end

PhaseAverageFluxPaperPlot("blfux_Phavg_pap_" * sn, path_name * "Analysis/", 7:2:15, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, vb_WpertPhavg_xavg_dy, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg)

PhaseAverageFluxPaperPlot("blfux_Phavg_pap_b_" * sn, path_name * "Analysis/", 7:2:15, v̂_Phavg_xavg_b, ŵ_Phavg_xavg_b, vb_WpertPhavg_xavg_b, 
        wb_WpertPhavg_xavg_b, vb_WpertPhavg_xavg_dy_b, ∇WpertPhavg_xavg_b, e_Phavg_xavg_b, b_Phavg_xavg_b)

    
vmax = pm.U₀
wmax = 0.1
vbmax=21e-5
wbmax=10e-5
nabmax=7e-5

vmax2= "0.15"
wmax2= "0.05"
vbmax2= "10.5"
wbmax2= "5"
nabmax2= "3.5"

function PhaseAverageFluxClosePlot(savename, apath, tindxs, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
    wb_WpertPhavg_xavg, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg,
    vmax, wmax, vbmax, nabmax, vmax2, wmax2,vbmax2,wbmax2, nabmax2)

    f = Figure(resolution = (2200, 1400), fontsize=26)
    ga = f[1, 1:5] = GridLayout() # vhat
    gb = f[2, 1:5] = GridLayout() # what
    gc = f[3, 1:5] = GridLayout() #<v'b'>
    gd = f[4, 1:5] = GridLayout() #<w'b'>
    ge = f[5, 1:5] = GridLayout() # grad
    gf = f[6, 1:5] = GridLayout() # diss

    # v
    gcb1 = f[1:2, 6] = GridLayout()
    # w

    #v'b', grad
    gcb2 = f[3:4, 6] = GridLayout()
    #w'b'

    # dissip
    gcb3 = f[5:6, 6] = GridLayout()

    gplots = [ga, gb, gc, gd, ge, gf]

    for g in gplots
        if g == gf
            for m = 1:5
                if m==1
                    axi = Axis(g[1, m], ylabel = "z [m]", xlabel = "y [m]") 
                    global ax = axi
                    axi.xticks = 500:300:1000
                    axi.yticks = [-400, -300]
                else
                    axi = Axis(g[1,m])
                    hideydecorations!(axi)
                    global ax = hcat(ax, axi)
                    axi.xticks = 500:300:1000
                end
            end
            global axs = vcat(axs, ax)
        else
            for m = 1:5
                if m==1
                    axi = Axis(g[1, m])
                    global ax = axi
                    hidexdecorations!(axi)
                    axi.yticks = [-400, -300]
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
        limits!(axs[j], 700, 1500, -450, -250)
    end

    vblims = (-vbmax,vbmax)
    wblims = (-wbmax,wbmax)
    gradblims = (-nabmax, nabmax)

    time_pre = "t = "
    time_post = " Tσ"

    zb_cut = zb[2:end]
    yb_cut2 = yb_cut[2:end]

    bmin = round(pm.Ñ^2*500*1e3)*1e-3
    bstep = bmin/16

    for (m,idx) in enumerate(tindxs)
        t = phase_times[idx]
        phaselabel = time_pre * @sprintf("%0.2f", t) * time_post

        v̂_Phavg_xavgi = v̂_Phavg_xavg[:,:,idx]
        ŵ_Phavg_xavgi = ŵ_Phavg_xavg[:,:,idx]

        vb_WpertPhavg_xavgi = vb_WpertPhavg_xavg[:,:,idx]
        wb_WpertPhavg_xavgi = wb_WpertPhavg_xavg[:,:,idx]
        b_Phavg_xavgi = b_Phavg_xavg[:,:,idx]
        ∇WpertPhavg_xavgi  = ∇WpertPhavg_xavg[:,:,idx]
        e_Phavg_xavgi = e_Phavg_xavg[:,:,idx]
        
        global hv = heatmap!(axs[1, m], yb_cut, zb, v̂_Phavg_xavgi, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axs[1, m], yb, land, color=:black, lw = 4)
        lines!(axs[1, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[1, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hw = heatmap!(axs[2, m], yb_cut, zb, ŵ_Phavg_xavgi, colormap = :balance, colorrange = (-wmax, wmax))
        lines!(axs[2, m], yb, land, color=:black, lw = 4)
        lines!(axs[2, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[2, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hvb = heatmap!(axs[3, m], yb_cut, zb, vb_WpertPhavg_xavgi, colormap = :balance, colorrange = vblims)
        lines!(axs[3, m], yb, land, color=:black, lw = 4)
        lines!(axs[3, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[3, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)
        
        global hwb = heatmap!(axs[4, m], yb_cut, zb, wb_WpertPhavg_xavgi, colormap = :balance, colorrange = wblims)
        lines!(axs[4, m], yb, land, color=:black, lw = 4)
        lines!(axs[4, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[4, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hgradub = heatmap!(axs[5, m], yb_cut2, zb_cut, ∇WpertPhavg_xavgi, colormap = :balance, colorrange = gradblims)
        lines!(axs[5, m], yb, land, color=:black, lw = 4)
        lines!(axs[5, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[5, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global heps = heatmap!(axs[6, m], yb_cut, zb, log10.(clamp.(e_Phavg_xavgi, 1e-14, 1)), colormap = :thermal, colorrange = (-8,-5))
        lines!(axs[6, m], yb, land, color=:black, lw = 4)
        lines!(axs[6, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        
        Label(f[1, m, Top()], phaselabel, valign = :bottom,
        font = :bold, fontsize = 36,
        padding = (0, 0, 5, 0))
    end

    text!(axs[1, 1], Point.(750, -400), text = "⟨v̂⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[2, 1], Point.(750, -400), text = "⟨ŵ⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[3, 1], Point.(750, -400), text = "⟨v'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[4, 1], Point.(750, -400), text = "⟨w'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[5, 1], Point.(750, -400), text = "∇⋅⟨u⃗'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[6, 1], Point.(750, -400), text = "⟨ε⟩ₚₕ", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 36) 

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax/2:vmax/2:vmax/2, ["-" * vmax2, "0", vmax2] ),
    size =25, flipaxis=false, label = "⟨v̂⟩ₚₕ", labelsize = 30)
    cb2 = Colorbar(gcb1[1,2], hw, ticks = (-wmax/2:wmax/2:wmax/2, ["-" * wmax2, "0", wmax2] ), 
    size =25, label = "⟨ŵ⟩ₚₕ", labelsize = 30)

    cb3 = Colorbar(gcb2[1,1], hvb, ticks = (-vbmax/2:vbmax/2:vbmax/2, ["-" * vbmax2, "0", vbmax2] ),
    size =25, flipaxis=false, label = "⟨v'b'⟩ × 10⁻⁵", labelsize = 30)
    cb4 = Colorbar(gcb2[1,2], hwb, ticks = (-wbmax/2:wbmax/2:wbmax/2, ["-" * wbmax2, "0", wbmax2] ), 
    size =25, label = "⟨w'b'⟩ × 10⁻⁵", labelsize = 30)

    cb5 = Colorbar(gcb3[1,1], hgradub, ticks = (-nabmax/2:nabmax/2:nabmax/2,  ["-" * nabmax2, "0", nabmax2] ), 
    size =25, flipaxis=false, label = "∇⋅⟨u⃗'b'⟩ₚₕ  × 10⁻⁵", labelsize = 30)
    cb6 = Colorbar(gcb3[1,2], heps, ticks = (-8:1:-6, ["10⁻⁸", "10⁻⁷", "10⁻⁶"] ), 
    size =25, label = "⟨ε⟩ₚₕ" , labelsize = 30)

    colgap!(gcb1, 10)
    colgap!(gcb2, 10)
    colgap!(gcb3, 10)

    colsize!(f.layout, 6, Relative(0.02))

    save(apath * savename * ".png", f)
end

PhaseAverageFluxClosePlot("blfux_Phavg_sm_mid_" * sn, "Analysis/Plots/WMToptions/", 7:2:15, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg,
        vmax, wmax, vbmax, nabmax, vmax2, wmax2,vbmax2,wbmax2,nabmax2)

PhaseAverageFluxClosePlot("blfux_Phavg_sm_beg_" * sn, "Analysis/Plots/WMToptions/", 7:2:15, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg,
        vmax, wmax, vbmax, nabmax, vmax2, wmax2,vbmax2,wbmax2,nabmax2)

PhaseAverageFluxClosePlot("blfux_Phavg_sm_end_" * sn, "Analysis/Plots/WMToptions/", 7:2:15, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg,
        vmax, wmax, vbmax, nabmax, vmax2, wmax2,vbmax2,wbmax2,nabmax2)

 
vmax = pm.U₀
wmax = round(maximum(abs.(ŵ_Phavg_xavg))*100)*1e-2

vbmax = round(maximum(abs.(vb_WpertPhavg_xavg))*1e5)*1e-5
wbmax = round(maximum(abs.(wb_WpertPhavg_xavg))*1e5)*1e-5
nabmax = round(maximum(abs.(vb_WpertPhavg_xavg_dy))*1e5)*1e-5

vmax2 = string(round(maximum(abs.(pm.U₀/2))*1000)*1e-3)
wmax2 = string(round(maximum(abs.(ŵ_Phavg_xavg)/2)*100)*1e-2)
vbmax2 = string(round(maximum(abs.(vb_WpertPhavg_xavg)/2)*1e5))
wbmax2 = string(round(maximum(abs.(wb_WpertPhavg_xavg)/2)*1e5))
nabmax2 = string(round(maximum(abs.(vb_WpertPhavg_xavg_dy)/2)*1e5))

function PhaseAverageFluxClosePlot_wDye(savename, apath, tindxs, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg, c_Phavg_xavg,
        vmax, wmax, vbmax, nabmax, vmax2, wmax2,vbmax2,wbmax2,nabmax2)

    f = Figure(resolution = (2200, 1400), fontsize=26)
    ga = f[1, 1:5] = GridLayout() # vhat
    gb = f[2, 1:5] = GridLayout() # what
    gc = f[3, 1:5] = GridLayout() #<v'b'>
    gd = f[4, 1:5] = GridLayout() #<w'b'>
    ge = f[5, 1:5] = GridLayout() # grad
    gf = f[6, 1:5] = GridLayout() # diss

    # v
    gcb1 = f[1:2, 6] = GridLayout()
    # w

    #v'b', grad
    gcb2 = f[3:4, 6] = GridLayout()
    #w'b'

    # dissip
    gcb3 = f[5:6, 6] = GridLayout()

    gplots = [ga, gb, gc, gd, ge, gf]

    for g in gplots
        if g == gf
            for m = 1:5
                if m==1
                    axi = Axis(g[1, m], ylabel = "z [m]", xlabel = "y [m]") 
                    global ax = axi
                    axi.xticks = 500:300:1000
                    axi.yticks = [-400, -300]
                else
                    axi = Axis(g[1,m])
                    hideydecorations!(axi)
                    global ax = hcat(ax, axi)
                    axi.xticks = 500:300:1000
                end
            end
            global axs = vcat(axs, ax)
        else
            for m = 1:5
                if m==1
                    axi = Axis(g[1, m])
                    global ax = axi
                    hidexdecorations!(axi)
                    axi.yticks = [-400, -300]
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
        limits!(axs[j], 700, 1500, -450, -250)
    end

    vblims = (-vbmax,vbmax)
    wblims = (-wbmax,wbmax)
    gradblims = (-nabmax, nabmax)

    time_pre = "t = "
    time_post = " Tσ"

    zb_cut = zb[2:end]
    yb_cut2 = yb_cut[2:end]

    bmin = round(pm.Ñ^2*500*1e3)*1e-3
    bstep = bmin/16

    for (m,idx) in enumerate(tindxs)
        t = phase_times[idx]
        phaselabel = time_pre * @sprintf("%0.2f", t) * time_post

        v̂_Phavg_xavgi = v̂_Phavg_xavg[:,:,idx]
        ŵ_Phavg_xavgi = ŵ_Phavg_xavg[:,:,idx]

        vb_WpertPhavg_xavgi = vb_WpertPhavg_xavg[:,:,idx]
        wb_WpertPhavg_xavgi = wb_WpertPhavg_xavg[:,:,idx]
        b_Phavg_xavgi = b_Phavg_xavg[:,:,idx]
        #c_Phavg_xavgi = log10.(clamp.(c_Phavg_xavg[:,:, idx], 1e-8, 1))
        c_Phavg_xavgi = clamp.(c_Phavg_xavg[:,:, idx], 1e-8, 1)
        ∇WpertPhavg_xavgi  = ∇WpertPhavg_xavg[:,:,idx]
        e_Phavg_xavgi = e_Phavg_xavg[:,:,idx]
        
        global hv = heatmap!(axs[1, m], yb_cut, zb, v̂_Phavg_xavgi, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axs[1, m], yb, land, color=:black, lw = 4)
        lines!(axs[1, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[1, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hw = heatmap!(axs[2, m], yb_cut, zb, ŵ_Phavg_xavgi, colormap = :balance, colorrange = (-wmax, wmax))
        lines!(axs[2, m], yb, land, color=:black, lw = 4)
        lines!(axs[2, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[2, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hvb = heatmap!(axs[3, m], yb_cut, zb, vb_WpertPhavg_xavgi, colormap = :balance, colorrange = vblims)
        lines!(axs[3, m], yb, land, color=:black, lw = 4)
        lines!(axs[3, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[3, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)
        
        global hwb = heatmap!(axs[4, m], yb_cut, zb, wb_WpertPhavg_xavgi, colormap = :balance, colorrange = wblims)
        lines!(axs[4, m], yb, land, color=:black, lw = 4)
        lines!(axs[4, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[4, m], yb_cut, zb, b_Phavg_xavgi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        global hgradub = heatmap!(axs[5, m], yb_cut2, zb_cut, ∇WpertPhavg_xavgi, colormap = :balance, colorrange = gradblims)
        lines!(axs[5, m], yb, land, color=:black, lw = 4)
        lines!(axs[5, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        contour!(axs[5, m], yb_cut, zb, c_Phavg_xavgi, color = :black, lw = 6, levels = 1e-3:1e-2:0.1, alpha = 0.3) # levels = -3:0.5:-1, alpha = 0.3)

        global heps = heatmap!(axs[6, m], yb_cut, zb, c_Phavg_xavgi, colormap = :thermal, colorrange = (1e-3,0.1))
        #log10.(clamp.(e_Phavg_xavgi, 1e-14, 1)), colormap = :thermal, colorrange = (-8,-5))
        lines!(axs[6, m], yb, land, color=:black, lw = 4)
        lines!(axs[6, m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        
        Label(f[1, m, Top()], phaselabel, valign = :bottom,
        font = :bold, fontsize = 36,
        padding = (0, 0, 5, 0))
    end

    text!(axs[1, 1], Point.(750, -400), text = "⟨v̂⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[2, 1], Point.(750, -400), text = "⟨ŵ⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[3, 1], Point.(750, -400), text = "⟨v'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[4, 1], Point.(750, -400), text = "⟨w'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[5, 1], Point.(750, -400), text = "∇⋅⟨u⃗'b'⟩ₚₕ", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 36)
    text!(axs[6, 1], Point.(750, -400), text = "⟨c⟩ₚₕ", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 36) # "⟨ε⟩ₚₕ"

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax/2:vmax/2:vmax/2, ["-" * vmax2, "0", vmax2] ),
    size =25, flipaxis=false, label = "⟨v̂⟩ₚₕ", labelsize = 30)
    cb2 = Colorbar(gcb1[1,2], hw, ticks = (-wmax/2:wmax/2:wmax/2, ["-" * wmax2, "0", wmax2] ), 
    size =25, label = "⟨ŵ⟩ₚₕ", labelsize = 30)

    cb3 = Colorbar(gcb2[1,1], hvb, ticks = (-vbmax/2:vbmax/2:vbmax/2, ["-" * vbmax2, "0", vbmax2] ),
    size =25, flipaxis=false, label = "⟨v'b'⟩ × 10⁻⁵", labelsize = 30)
    cb4 = Colorbar(gcb2[1,2], hwb, ticks = (-wbmax/2:wbmax/2:wbmax/2, ["-" * wbmax2, "0", wbmax2] ), 
    size =25, label = "⟨w'b'⟩ × 10⁻⁵", labelsize = 30)

    cb5 = Colorbar(gcb3[1,1], hgradub, ticks = (-nabmax/2:nabmax/2:nabmax/2,  ["-" * nabmax2, "0", nabmax2] ), 
    size =25, flipaxis=false, label = "∇⋅⟨u⃗'b'⟩ₚₕ  × 10⁻⁵", labelsize = 30)
    #cb6 = Colorbar(gcb3[1,2], heps, ticks = (-8:1:-6, ["10⁻⁸", "10⁻⁷", "10⁻⁶"] ), 
    cb6 = Colorbar(gcb3[1,2], heps, ticks = (-3:1:-1, ["10⁻³", "10⁻²", "10⁻¹"] ), 
    size =25, label = "⟨c⟩ₚₕ" , labelsize = 30)

    colgap!(gcb1, 10)
    colgap!(gcb2, 10)
    colgap!(gcb3, 10)

    colsize!(f.layout, 6, Relative(0.02))

    save(apath * savename * ".png", f)
end

PhaseAverageFluxClosePlot_wDye("blfux_Phavg_sm_" * sn, "Analysis/Plots/WMToptions/", 7:2:15, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg, c_Phavg_xavg,
        vmax, wmax, vbmax, nabmax, vmax2, wmax2,vbmax2,wbmax2,nabmax2)

PhaseAverageFluxClosePlot_wDye("blfux_Phavg_sm_b_" * sn, "Analysis/Plots/WMToptions/", 7:2:15, v̂_Phavg_xavg_b, ŵ_Phavg_xavg_b, vb_WpertPhavg_xavg_b, 
        wb_WpertPhavg_xavg_b, ∇WpertPhavg_xavg_b, e_Phavg_xavg_b, b_Phavg_xavg_b, c_Phavg_xavg_b,
        vmax, wmax, vbmax, nabmax, vmax2, wmax2,vbmax2,wbmax2,nabmax2)
