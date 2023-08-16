using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using ArgParse
using CairoMakie

dpath = "Data/"

sn = "U250N100Lz100g100"
filescalename1 = dpath * "BFluxFull_" * sn * ".jld2"
#filescalenameEps = dpath * "dissipwave_" * sn * ".jld2"

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
#eps_file = jldopen(filescalenameEps, "r+")
skeys = keys(scale_file1)

vb_WpertPhavg  = scale_file1[skeys[1]]
wb_WpertPhavg  = scale_file1[skeys[2]]
v̂b_WpertPhavg  = scale_file1[skeys[3]]
ŵb_WpertPhavg  = scale_file1[skeys[4]]
vb_WpertWavg  = scale_file1[skeys[5]]
wb_WpertWavg  = scale_file1[skeys[6]]
v̂b_WpertWavg  = scale_file1[skeys[7]]
ŵb_WpertWavg  = scale_file1[skeys[8]]
v_Phavg  = scale_file1[skeys[9]]
v_Wavg  = scale_file1[skeys[10]]
v̂_Phavg  = scale_file1[skeys[11]]
v̂_Wavg  = scale_file1[skeys[12]]
b_Phavg  = scale_file1[skeys[13]]

ŵ_Phavg  = whatfile["ŵ_Phavg_xavg"]

#e_Phavg_xavg = eps_file["e_Phavg_xavg"]
yb=scale_file1[skeys[14]]
zb=scale_file1[skeys[15]]
phase_times= scale_file1[skeys[16]]

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

big_titlep = @sprintf("Phase Averaged, δ=%0.1f m, N = %0.1f × 10⁻³", pm.U₀/pm.Ñ, pm.Ñ*1e3)
big_titlew = @sprintf("Wave Averaged, δ=%0.1f m, N = %0.1f × 10⁻³", pm.U₀/pm.Ñ, pm.Ñ*1e3)

zlength = length(zb)
ylength = 880

Ygrid = reshape(repeat(yb[2:ylength], zlength-1), ylength-1, zlength-1)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

vb_WpertWavg_xavg = mean(vb_WpertWavg, dims = 1)[1,:,:];
wb_WpertWavg_xavg = mean(wb_WpertWavg, dims = 1)[1,:,:];
ŵb_WpertWavg_xavg = mean(ŵb_WpertWavg, dims = 1)[1,:,:];
v̂b_WpertWavg_xavg = mean(v̂b_WpertWavg, dims = 1)[1,:,:];
v̂_Wavg_xavg = mean(v̂_Wavg, dims =1)[1,:,:];
v_Wavg_xavg = mean(v_Wavg, dims =1)[1,:,:];

vb_WpertWavg_xavg_b = mean(vb_WpertWavg_b, dims = 1)[1,:,:];
wb_WpertWavg_xavg_b = mean(wb_WpertWavg_b, dims = 1)[1,:,:];
ŵb_WpertWavg_xavg_b = mean(ŵb_WpertWavg_b, dims = 1)[1,:,:];
v̂b_WpertWavg_xavg_b = mean(v̂b_WpertWavg_b, dims = 1)[1,:,:];
v̂_Wavg_xavg_b = mean(v̂_Wavg_b, dims =1)[1,:,:];
v_Wavg_xavg_b = mean(v_Wavg_b, dims =1)[1,:,:];

wb_WpertWavg_xavg_dz = ((wb_WpertWavg_xavg[2:end,1:end-1] .- wb_WpertWavg_xavg[2:end,2:end])./(zb[1]-zb[2])).* boolZY;
vb_WpertWavg_xavg_dy = ((vb_WpertWavg_xavg[1:end-1,2:end] .- wb_WpertWavg_xavg[2:end,2:end])./(yb[1]-yb[2])).* boolZY;
∇WpertWavg_xavg = wb_WpertWavg_xavg_dz .+ vb_WpertWavg_xavg_dy;

wb_WpertWavg_xavg_dz_b = ((wb_WpertWavg_xavg_b[2:end,1:end-1] .- wb_WpertWavg_xavg_b[2:end,2:end])./(zb[1]-zb[2])).* boolZY;
vb_WpertWavg_xavg_dy_b = ((vb_WpertWavg_xavg_b[1:end-1,2:end] .- wb_WpertWavg_xavg_b[2:end,2:end])./(yb[1]-yb[2])).* boolZY;
∇WpertWavg_xavg_b = wb_WpertWavg_xavg_dz_b .+ vb_WpertWavg_xavg_dy_b;


""" Wave Average Plotting """
function WaveAverageFluxPlot(savename, apath , ylength, v_Wavg_xavg, vb_WpertWavg_xavg, wb_WpertWavg_xavg, v̂_Wavg_xavg, v̂b_WpertWavg_xavg, ŵb_WpertWavg_xavg, vb_WpertWavg_xavg_dy, wb_WpertWavg_xavg_dz, ∇WpertWavg_xavg)
    vmax =  round(maximum(v_Wavg_xavg)*100)*1e-2
    vbmax = round(maximum(vb_WpertWavg_xavg)*1e5)*1e-5
    wbmax = round(maximum(abs.(wb_WpertWavg_xavg))*1e6)*1e-6

    nabmax = round(maximum(abs.(vb_WpertWavg_xavg_dy))*1e6)*1e-6
    dymax = round(maximum(abs.(vb_WpertWavg_xavg_dy))*1e6)*1e-6
    dzmax = round(maximum(wb_WpertWavg_xavg_dz)*1e7*0.75)*1e-7

    vblims = (-vbmax,vbmax)
    wblims = (-wbmax,wbmax)
    gradblims = (-nabmax, nabmax)
    dvblims = (-dymax, dymax)
    dwblims = (-dzmax, dzmax)

    yb_cut = yb[1:ylength]

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

    hv = heatmap!(axv, yb_cut, zb, v_Wavg_xavg, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axv, yb, land, color=:black, lw = 4)
        lines!(axv, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axv, Point.(50, -450), text = "v", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hvhat = heatmap!(axvh, yb_cut, zb, v̂_Wavg_xavg, colormap = :balance, colorrange = (-vmax, vmax))
        lines!(axvh, yb, land, color=:black, lw = 4)
        lines!(axvh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvh, Point.(50, -450), text = "v̂", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hnab = heatmap!(axgrad, yb_cut[2:end], zb[2:end], ∇WpertWavg_xavg, colormap = :balance, colorrange = gradblims)
        lines!(axgrad, yb, land, color=:black, lw = 4)
        lines!(axgrad, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axgrad, Point.(50, -450), text = "∇⋅⟨u⃗'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hvb = heatmap!(axvb, yb_cut, zb, vb_WpertWavg_xavg, colormap = :balance, colorrange = vblims)
        lines!(axvb, yb, land, color=:black, lw = 4)
        lines!(axvb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvb, Point.(50, -450), text = "⟨v'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    heatmap!(axvbh, yb_cut, zb, v̂b_WpertWavg_xavg, colormap = :balance, colorrange = vblims)
        lines!(axvbh, yb, land, color=:black, lw = 4)
        lines!(axvbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axvbh, Point.(50, -450), text = "⟨v̂'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hdyvb = heatmap!(axdyvb, yb_cut, zb, vb_WpertWavg_xavg_dy, colormap = :balance, colorrange = dvblims)
        lines!(axdyvb, yb, land, color=:black, lw = 4)
        lines!(axdyvb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdyvb, Point.(50, -450), text = "∂y⟨v'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hwb = heatmap!(axwb, yb_cut, zb, wb_WpertWavg_xavg, colormap = :balance, colorrange = wblims)
        lines!(axwb, yb, land, color=:black, lw = 4)
        lines!(axwb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axwb, Point.(50, -450), text = "⟨w'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    heatmap!(axwbh, yb_cut, zb, ŵb_WpertWavg_xavg, colormap = :balance, colorrange = vblims)
        lines!(axwbh, yb, land, color=:black, lw = 4)
        lines!(axwbh, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axwbh, Point.(50, -450), text = "⟨ŵ'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    hdzwb = heatmap!(axdzwb, yb_cut, zb, wb_WpertWavg_xavg_dz, colormap = :balance, colorrange = dwblims)
        lines!(axdzwb, yb, land, color=:black, lw = 4)
        lines!(axdzwb, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(axdzwb, Point.(50, -450), text = "∂z⟨w'b'⟩", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = (-vmax:vmax:vmax), size =35,vertical = false)
    cb2 = Colorbar(gcb2[1,1], hvb, ticks = (-vbmax/2:vbmax/2:vbmax/2), size =35,vertical = false, )
    cb3 = Colorbar(gcb3[1,1], hwb, ticks = (-wbmax/2:wbmax/2:wbmax/2), size =35,vertical = false, )

    cb3 = Colorbar(gcb1[1,1], hnab, ticks = (-nabmax/2:nabmax/2:nabmax/2),
            size =35,vertical = false, flipaxis=false)
    cb4 = Colorbar(gcb2[1,1], hdyvb, ticks = (-dymax/2:dymax/2:dymax/2),
            size =35,vertical = false, flipaxis=false)
    cb5 = Colorbar(gcb3[1,1], hdzwb, ticks = (-dzmax/2:dzmax/2:dzmax/2),
            size =35,vertical = false, flipaxis=false)

    #Label(f1[1, 1:3, Top()],big_titlew, valign = :bottom,
    #    font = :bold,
    #    padding = (0, 0, 50, 0))

    colgap!(ga, 15)
    rowgap!(ga, 5)
    colgap!(gb, 15)
    rowgap!(gb, 5)
    colgap!(gc, 15)
    rowgap!(gc, 5)

    save(apath * savename * ".png", f1)
end

WaveAverageFluxPlot("blfux_Wavg_" * sn, path_name * "Analysis/", ylength, v_Wavg_xavg, 
vb_WpertWavg_xavg, wb_WpertWavg_xavg, v̂_Wavg_xavg, v̂b_WpertWavg_xavg,
 ŵb_WpertWavg_xavg, vb_WpertWavg_xavg_dy, wb_WpertWavg_xavg_dz, ∇WpertWavg_xavg)

WaveAverageFluxPlot("blfux_Wavg_b_" * sn, path_name * "Analysis/", ylength, v_Wavg_xavg_b, 
vb_WpertWavg_xavg_b, wb_WpertWavg_xavg_b, v̂_Wavg_xavg_b, v̂b_WpertWavg_xavg_b,
 ŵb_WpertWavg_xavg_b, vb_WpertWavg_xavg_dy_b, wb_WpertWavg_xavg_dz_b, ∇WpertWavg_xavg_b)

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

v̂_Phavg_xavg_b = mean(v̂_Phavg_b, dims = 1)[1,:,:,:];
v_Phavg_xavg_b = mean(v_Phavg_b, dims = 1)[1,:,:,:];
vb_WpertPhavg_xavg_b = mean(vb_WpertPhavg_b, dims = 1)[1,:,:,:];
wb_WpertPhavg_xavg_b = mean(wb_WpertPhavg_b, dims = 1)[1,:,:,:];
ŵb_WpertPhavg_xavg_b = mean(ŵb_WpertPhavg_b, dims = 1)[1,:,:,:];
v̂b_WpertPhavg_xavg_b = mean(v̂b_WpertPhavg_b, dims = 1)[1,:,:,:];
b_Phavg_xavg_b = mean(b_Phavg_b, dims =1)[1,:,:,:];
e_Phavg_xavg_b = mean(e_Phavg_b, dims =1)[1,:,:,:];
ŵ_Phavg_xavg_b = mean(ŵ_Phavg_b, dims = 1)[1,:,:,:];

zlength = length(zb)
Ygrid = reshape(repeat(yb[2:ylength], zlength-1), ylength-1, zlength-1)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

wb_WpertPhavg_xavg_dz = ((wb_WpertPhavg_xavg[2:end,1:end-1,:] .- wb_WpertPhavg_xavg[2:end,2:end,:])./(zb[1]-zb[2])).* boolZY;
vb_WpertPhavg_xavg_dy = ((vb_WpertPhavg_xavg[1:end-1,2:end,:] .- wb_WpertPhavg_xavg[2:end,2:end,:])./(yb[1]-yb[2])).* boolZY;
∇WpertPhavg_xavg = wb_WpertPhavg_xavg_dz .+ vb_WpertPhavg_xavg_dy;

wb_WpertPhavg_xavg_dz_b = ((wb_WpertPhavg_xavg_b[2:end,1:end-1,:] .- wb_WpertPhavg_xavg_b[2:end,2:end,:])./(zb[1]-zb[2])).* boolZY;
vb_WpertPhavg_xavg_dy_b = ((vb_WpertPhavg_xavg_b[1:end-1,2:end,:] .- wb_WpertPhavg_xavg_b[2:end,2:end,:])./(yb[1]-yb[2])).* boolZY;
∇WpertPhavg_xavg_b = wb_WpertPhavg_xavg_dz_b .+ vb_WpertPhavg_xavg_dy_b;

(_,_,Lw) = size(b_Phavg_xavg_b)
tindxs = vcat(vcat(1:2:7, 8),9:2:13)
#tindxs = vcat(vcat(1:3:10, 11),12:3:Lw)

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

PhaseAverageFluxPlot("blfux_Phavg_" * sn, path_name * "Analysis/", tindxs, v_Phavg_xavg, vb_WpertPhavg_xavg, 
    wb_WpertPhavg_xavg, vb_WpertPhavg_xavg_dy, wb_WpertPhavg_xavg_dz, ∇WpertPhavg_xavg, b_Phavg_xavg)

PhaseAverageFluxPlot("blfux_Phavg_b_" * sn, path_name * "Analysis/", tindxs, v_Phavg_xavg_b, vb_WpertPhavg_xavg_b, 
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

tindxs = 7:2:15

PhaseAverageFluxPaperPlot("blfux_Phavg_pap_" * sn, path_name * "Analysis/", tindxs, v̂_Phavg_xavg, ŵ_Phavg_xavg, vb_WpertPhavg_xavg, 
        wb_WpertPhavg_xavg, vb_WpertPhavg_xavg_dy, ∇WpertPhavg_xavg, e_Phavg_xavg, b_Phavg_xavg)

PhaseAverageFluxPaperPlot("blfux_Phavg_pap_b_" * sn, path_name * "Analysis/", tindxs, v̂_Phavg_xavg_b, ŵ_Phavg_xavg_b, vb_WpertPhavg_xavg_b, 
        wb_WpertPhavg_xavg_b, vb_WpertPhavg_xavg_dy_b, ∇WpertPhavg_xavg_b, e_Phavg_xavg_b, b_Phavg_xavg_b)