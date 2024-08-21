using Statistics
using Printf
using Oceananigans
using Measures
using JLD2
using CairoMakie

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"
sn = "U350N100Lz100g100"

ENV["GKSwstype"] = "nul" # if on remote HPC

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

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

filesetnames = "SetList_mp.jld2"
scale_file = jldopen(filesetnames, "r+")
sns = scale_file["setnames"]
sfiles = scale_file["setfilenames"]


name_prefix =  sfiles[7] * sns[7]

filepath = path_name * name_prefix * ".jld2"
apath = path_name * "Analysis/"

b_timeseries = FieldTimeSeries(filepath, "b");
N2_timeseries = FieldTimeSeries(filepath, "N2");

xb, yb, zb = nodes(b_timeseries) #CCC
xn, yn, zn = nodes(N2_timeseries) #CCC

tlength = 161

zlength = length(zb)
ylength = 880
xlength = length(xb)

bi = interior(b_timeseries)[:,1:ylength,:,:];
N2i = interior(N2_timeseries)[:,1:ylength,:,:];

# derivative 
oneover_Δy = 1 / (yb[2]-yb[1])
M2i = zeros(size(bi))
M2i[:,1,:,:] = (bi[:, 2, :,:] .- bi[:, 1, :,:] ) .* oneover_Δy
M2i[:,end,:,:] = (bi[:, end, :,:] .- bi[:, end-1, :,:] ) .* oneover_Δy
M2i[:,2:end-1,:,:] = (bi[:, 3:end, :,:] .- bi[:, 1:end-2, :,:] ) .* (0.5*oneover_Δy);

# remove slope adjacent grid points
Ygrid = reshape(repeat(yb[1:ylength], xlength*(zlength)), ylength, zlength, xlength)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb';  # all the values greater than slope

boolZYX = permutedims(boolZY, [3, 1, 2])

M2i_noslope = M2i .* boolZYX
N2i_noslope = N2i .* boolZYX

# diapycnal direction
M2_over_N2 = M2i_noslope ./ N2i_noslope

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl

fullwaves = wave_info.WavePeriods[:,4:end] # t = 4T : t = 7T

# phase indexing
M2_Ph = M2i_noslope[:,:,:,fullwaves]
N2_Ph = N2i_noslope[:,:,:,fullwaves]
M2_over_N2_Ph = M2_over_N2[:,:,:,fullwaves];
b_Ph = bi[:,:,:,fullwaves];

# wave and x averaging [y,z]
M2_Wavg = mean(M2_Ph, dims = (1,4,5))[1,:,:,1,1];
N2_Wavg = mean(N2_Ph, dims = (1,4,5))[1,:,:,1,1];
M2_over_N2_Wavg = mean(M2_over_N2_Ph, dims = (1,4,5))[1,:,:,1,1];
b_Wavg = mean(b_Ph, dims = (1,4,5))[1,:,:,1,1];

M2_Wavg_over_N2_Wavg = M2_Wavg ./ N2_Wavg

N2_Phavg = mean(N2_Ph, dims = (1,5))[1,:,:,:,1];
M2_Phavg = mean(M2_Ph, dims = (1,5))[1,:,:,:,1];
b_Phavg = mean(b_Ph, dims = (1,5))[1,:,:,:,1];

M2_Phdep = M2_Phavg .- M2_Wavg;
N2_Phdep = N2_Phavg .- N2_Wavg;
b_Phdep = b_Phavg .- b_Wavg;
M2_Phdep_over_N2_Phdep = M2_Phdep ./ N2_Phdep


land = curvedslope.(yb)
land_pdel = (curvedslope.(yb) .+ pm.U₀/pm.Ñ)[1:382]
yb_cut = yb[1:ylength]

f1 = Figure(resolution = (800, 700), fontsize=26)
    ga = f1[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
    ax1.xticks = 500:500:1000
    ax1.yticks = [-250, 0]

    limits!(ax1, 0, 1500, -500, 0)

    hm = heatmap!(ax1, yb_cut, zb, N2_over_M2_Wavg, colormap = :balance, colorrange = (-5, 5))
        lines!(ax1, yb, land, color=:black, lw = 4)
        band!(ax1, yb, -500 .* ones(length(land)), land, color=:black)
        lines!(ax1, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(ax1, Point.(50, -450), text = "N²/M²", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)

    cb1 = Colorbar(ga[1,2], hm, ticks = (-5:1:5), size =35, label = "Diapycnal Direction")
    colsize!(ga, 2, Relative(0.05))

save(apath * "DiapycnalDirection_Wavg_" * sn * ".png", f1, px_per_unit = 2)

f1 = Figure(resolution = (800, 700), fontsize=26)
    ga = f1[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
    ax1.xticks = 500:500:1000
    ax1.yticks = [-250, 0]

    limits!(ax1, 0, 1500, -500, 0)

    hm = heatmap!(ax1, yb_cut, zb, M2_Wavg_over_N2_Wavg, colormap = :balance, colorrange = (-0.25, 0.25))
        lines!(ax1, yb, land, color=:black, lw = 4)
        band!(ax1, yb, -500 .* ones(length(land)), land, color=:black)
        lines!(ax1, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(ax1, Point.(50, -450), text = "M̄²/N̄²", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)

    cb1 = Colorbar(ga[1,2], hm, ticks = (-0.2:0.1:0.2), size =35, label = "Diapycnal Direction")
    colsize!(ga, 2, Relative(0.05))

save(apath * "DiapycnalDirection_Wavgfirst_Fixed_" * sn * ".png", f1, px_per_unit = 2)

f1 = Figure(resolution = (800, 700), fontsize=26)
    ga = f1[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
    ax1.xticks = 500:500:1000
    ax1.yticks = [-250, 0]

    limits!(ax1, 0, 1500, -500, 0)

    hm = heatmap!(ax1, yb_cut, zb, (180/π) .* (atan.(M2_Wavg_over_N2_Wavg)), colormap = :balance, colorrange = (-5, 5))
        lines!(ax1, yb, land, color=:black, lw = 4)
        band!(ax1, yb, -500 .* ones(length(land)), land, color=:black)
        lines!(ax1, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(ax1, Point.(50, -450), text = "M̄²/N̄²", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)

    cb1 = Colorbar(ga[1,2], hm, ticks = (-4:2:4), size =35, label = "Diapycnal Direction (∘)")
    colsize!(ga, 2, Relative(0.05))

save(apath * "DiapycnalDirection_Wavgfirst_FixedDegs_" * sn * ".png", f1, px_per_unit = 2)

f1 = Figure(resolution = (800, 700), fontsize=26)
    ga = f1[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1], ylabel = "z [m]", xlabel = "y [m]") #vh
    ax1.xticks = 500:500:1000
    ax1.yticks = [-250, 0]

    limits!(ax1, 0, 1500, -500, 0)

    hm = heatmap!(ax1, yb_cut, zb, abs.(M2_Wavg_over_N2_Wavg), colormap = :amp, colorrange = (0, 0.2))
        lines!(ax1, yb, land, color=:black, lw = 4)
        band!(ax1, yb, -500 .* ones(length(land)), land, color=:black)
        lines!(ax1, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(ax1, Point.(50, -450), text = "|M̄²/N̄²|", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)

    cb1 = Colorbar(ga[1,2], hm, ticks = (0:0.1:0.2), size =35, label = "Diapycnal Direction")
    colsize!(ga, 2, Relative(0.05))

save(apath * "DiapycnalDirection_Wavgfirst_FixedAmp_" * sn * ".png", f1, px_per_unit = 2)

tindxs = vcat(vcat(1:2:7, 8),9:2:13)



f = Figure(resolution = (2300, 1000), fontsize=26)
    ga = f[1:2, 1:4] = GridLayout() # vhat
    gcb1 = f[1:2, 5] = GridLayout()

    for m = 1:4
        if m == 1 
            ax1 = Axis(ga[1, m], ylabel = "z [m]") 
                    global ax = ax1
                    ax1.yticks = [-250, 0]
            hidexdecorations!(ax1)
            ax2 = Axis(ga[2, m], ylabel = "z [m]", xlabel = "y [m]") 
                    global ax = ax2
                    ax2.xticks = 500:500:1000
                    ax2.yticks = [-250, 0]
            global ax = [ax1; ax2]
        else
            axt = Axis(ga[1, m])
            axb = Axis(ga[2,m], xlabel = "y [m]")
            hidedecorations!(axt)
            hideydecorations!(axb)
            axb.xticks = 500:500:1000
            global ax = hcat(ax, [axt; axb])
        end
    end

    rowgap!(ga, 30)
    #colgap!(ga, 20)
    
    for j = 1:length(ax)
        limits!(ax[j], 0, 1500, -500, 0)
    end

    time_pre = "t = "
    time_post = " Tσ"

    bmin = round(pm.Ñ^2*500*1e3)*1e-3
    bstep = bmin/8
    phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

    for (m, idx) in enumerate(tindxs)
        t = phase_times[idx]
        phaselabel = time_pre * @sprintf("%0.2f", t) * time_post

        M2_Phdep_over_N2_Phdepi = M2_Phdep_over_N2_Phdep[:,:,idx]
        b_Phdepi = b_Phdep[:,:,idx]
        
        global hv = heatmap!(ax[m], yb_cut, zb, M2_Phdep_over_N2_Phdepi, colormap = :balance, colorrange = (-1.5, 1.5))
        lines!(ax[m], yb, land, color=:black, lw = 4)
        band!(ax[m], yb, -500, land, color=:black)
        lines!(ax[m], yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        #contour!(ax[m], yb_cut, zb, b_Phdepi, color = :black, lw = 6, levels = -bmin:bstep:0, alpha = 0.3)

        if m < 5
            Label(f[1, m, Top()], phaselabel, valign = :bottom, font = :bold, padding = (0, 0, 5, 0))
        else
            Label(f[2, m - 4, Top()], phaselabel, valign = :bottom, font = :bold, padding = (0, 0, 5, 0))
        end   
    end

    text!(ax[1, 1], Point.(50, -350), text = "M̃²/Ñ²", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hv, ticks = -1:0.25:1, size =35, label = "Phase Dep. Diapycnal Direction")

    colsize!(f.layout, 5, Relative(0.1))


save(apath * "DiapycnalDirection_Phdepfirst_" * sn * ".png", f)

