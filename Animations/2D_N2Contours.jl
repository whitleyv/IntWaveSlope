using Measures
using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

path_name = "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"
apath = path_name * "Analysis/"

filesetnames =  "SetList_mp.jld2"
file_sn = jldopen(filesetnames, "r+")
γ_σ = file_sn["γ_varyσ"]
idx_subcritical_σ = findall(γ_σ .< 1)
δ_σ = file_sn["δ_varyσ"]

setname_idx = idx_subcritical_σ[1]
setname = file_sn["setnames_varyσ"][setname_idx]

@info "Loading in parameters..."

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
    Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
    m = -π/pm.Lz,
    l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
    Tf = 2*π/pm.f, 
    Tσ = 2*π/pm.σ))
zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame
    
@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

@info "Pulling Data..."
name_prefix = file_sn["setfilenames_varyσ"][setname_idx] * setname
filepath = path_name * name_prefix * ".jld2"

N_timeseries = FieldTimeSeries(filepath, "N2");
xn, yn, zn = nodes(N_timeseries) #CCC

@info "Smoothing Stratification Values..."
# rolling averages width
rolWidy = 20
rolWidz = 3

# cutting off top and bottom to avoid boundaries:
lastH = 450 #start at z = -50 (used to be -250)
firstH = 50 #end at z = -450
# indices to start and stop 
z_st = round(Int, lastH/2) # start at z = -50
z_en = round(Int, firstH/2) # end at z = -450 (smaller indices = deeper)
y_en = round(Int, 2500/4) # just choosing this y value to include most of dye excursions
zlength_sm = length(z_en:z_st)
zlength = 250
z_sm = zn[z_en:z_st]

y_st = round(Int, ((pm.Lz-lastH)/pm.Tanα)/4) # find corresponding y value on slope when z = 250
y_st = (y_st > rolWidy) ?  y_st : (rolWidy+1)
ylength_sm = length(y_st:y_en)
y_sm = yn[y_st:y_en]

tlength = length(N_timeseries.times)

include("WaveValues.jl")
wave_info=get_wave_indices(N_timeseries, pm, tlength)
# starts at 7Tσ
W7length = wave_info.Wl * 7  + 1 
# thorpe starts at 3Tσ
W3length = wave_info.Wl * 3  + 1 
# total number of data points
Wtlength = wave_info.Wl * wave_info.nTσ
# 4 waves long
cutWtlength = 4*wave_info.Wl
# ends at 10.9Tσ
W11length = wave_info.Wl * 11
#for rolling wave avg, can't use the last half of wave if no more waves after that...
hWl = floor(Int64, wave_info.Wl/2)

# - 1/2 wave if doesn't fit since rolling 
if Wtlength > (W11length + hWl)
    # if room to go to end of "last" wave, then go!
    rWtlength = cutWtlength
else 
    rWtlength = cutWtlength - hWl
end

# averaging in x
ñ2i = mean(interior(N_timeseries), dims=1)[1,:,:,:]

# perturbations
ñ2i_init = ñ2i[:,:,1]
ñ2i_pert = ñ2i .- ñ2i_init
# initial N^2 should just be pm.N^2

# rolling y avg
ñ2i_pert_ry = zeros(ylength_sm, zlength, tlength)
for (yi, y_md) in enumerate(y_st:y_en)
    ñ2i_pert_ry[yi,:, :] = mean(ñ2i_pert[y_md-rolWidy:rolWidy+y_md, 1:zlength, :], dims = 1)
end

# rolling z avg
ñ2i_pert_ryrz = zeros(ylength_sm, zlength_sm, tlength)
# rolling vertical average over 7 grid points, to smooth things out
for (ki, zk) in enumerate(z_en:z_st)
    ñ2i_pert_ryrz[:,ki,:] = mean(ñ2i_pert_ry[:,zk-rolWidz:zk+rolWidz,:], dims =2)[:,1,:]
end

# rolling wave avg
ñ2i_pert_ryrzrW = zeros(ylength_sm, zlength_sm, rWtlength)
for (i,l) in enumerate(W7length:W7length+rWtlength-1)
    ñ2i_pert_ryrzrW[:,:,i] = mean(ñ2i_pert_ryrz[:,:,wave_info.WavePeriods[l-hWl:l+hWl]], dims=3)[:,:,1]
end

land = curvedslope.(yn) 

delta = pm.U₀/pm.Ñ

##############
#   CAIRO PLOT MOVIE
##############
times = N_timeseries.times[W7length:W7length+rWtlength-1]

n = Observable(1)

ñ = @lift ñ2i_pert_ryrzrW[:,:,$n];

title = @lift @sprintf("Internal Wave Breaking, t = %.2f hrs, Tσ = %.2f", times[$n]/3600, times[$n]/pm.Tσ)

f = Figure(resolution = (2000, 1000),fontsize=35) 

ga = f[1, 1] = GridLayout()
axv = Axis(ga[1, 1], ylabel = "z [m]")

axv.xticks = 500:500:2500
axv.yticks = [-400, -250, -50]

limits!(axv, 135, 2500, -451, -50)

land_pdel = (linslope.(yn) .+ delta)[1:382]

global hmv = heatmap!(axv, y_sm, z_sm, ñ, colormap = :balance, colorrange = (-1.2e-5, 1.2e-5))
band!(axv, yn, land, -500, color=:black)
lines!(axv, yn[1:382], land_pdel, color=:black, linewidth = 2, linestyle = :dash)

Label(ga[1, 1:2, Top()], title,
                fontsize = 30,
                font = :bold,
                padding = (5, 5, 5, 5),
                halign = :center)

cb1 = Colorbar(ga[1,2], hmv, ticks = (-1e-5:5e-6:1e-5), size =35, label = "Ñ²' Smoothed [s⁻²]")

colsize!(ga, 2, Relative(0.05))

savename = "presentation_n2contours_mp_" * setname

frames = 1:length(times)
record(f, apath * savename * ".mp4", frames, framerate=6) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end

