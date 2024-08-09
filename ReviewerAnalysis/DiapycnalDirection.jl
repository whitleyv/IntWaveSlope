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
N2_over_M2 = N2i_noslope ./ M2i_noslope


include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
Wl = wave_info.Wl

fullwaves = wave_info.WavePeriods[:,4:end] # t = 4T : t = 7T

# phase indexing
M2_Ph = M2i_noslope[:,:,:,fullwaves]
N2_Ph = N2i_noslope[:,:,:,fullwaves]
N2_over_M2_Ph = N2_over_M2[:,:,:,fullwaves];

# wave and x averaging [y,z]
M2_Wavg = mean(M2_Ph, dims = (1,4,5))[1,:,:,1,1];
N2_Wavg = mean(N2_Ph, dims = (1,4,5))[1,:,:,1,1];
N2_over_M2_Wavg = mean(N2_over_M2_Ph, dims = (1,4,5))[1,:,:,1,1];

N2_Wavg_over_M2_Wavg = N2_Wavg ./ M2_Wavg

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

    hm = heatmap!(ax1, yb_cut, zb, N2_Wavg_over_M2_Wavg, colormap = :balance, colorrange = (-100, 100))
        lines!(ax1, yb, land, color=:black, lw = 4)
        band!(ax1, yb, -500 .* ones(length(land)), land, color=:black)
        lines!(ax1, yb[1:382], land_pdel, color=:black, lw = 4, linestyle = :dash)
        text!(ax1, Point.(50, -450), text = "N̄²/M̄²", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)

    cb1 = Colorbar(ga[1,2], hm, ticks = (-75:25:75), size =35, label = "Diapycnal Direction")
    colsize!(ga, 2, Relative(0.05))

save(apath * "DiapycnalDirection_Wavgfirst_" * sn * ".png", f1, px_per_unit = 2)
#=
# phase averaging [x,y,z,ph)]
M2_Phavg = mean(PhaseOrientVals.N2_Ph, dims = (5))[:,:,:,:,1];
N2_Phavg = mean(PhaseOrientVals.M2_Ph, dims = (5))[:,:,:,:,1];

# phase dependent [x,y,z,ph)]
M2_Phdep = M2_Phavg .- M2_Wavg
N2_Phdep = N2_Phavg .- N2_Wavg

# turbulence [x,y,z,ph,w)]
M2_turb = M2_Ph .- M2_Phavg
N2_turb = N2_Ph .- N2_Phavg

=#