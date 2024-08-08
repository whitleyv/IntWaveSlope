using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"

apath = path_name * "Analysis/"


sn = "U250Nfd250Lz100g100"

include("parameters.jl")

pm = getproperty(SimParams(), Symbol(sn))

resS = 1.0
dzr = pm.dz * resS
dhr = pm.dh * resS

pm = merge(pm, (; dzr=dzr, dhr=dhr, Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2), nx = round(Int, pm.Lx/dhr),
                m = π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ
Sp_extra = ifelse(z_start>0, 250.0, 0.0)
Ly = pm.Lyˢ+Sp_extra
ny = round(Int,Ly/pm.dhr)
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


# how many times were saved?
zlength = pm.nz
ylength = 425 # roughly out to y = 1700 (is this enough still??)
xlength = pm.nx

# cutting off top and bottom to avoid boundaries:

lastH = 450 #start at z = -50 (used to be -250)
firstH = 50 #end at z = -450
# indices to start and stop 
z_st = round(Int, lastH/2) # start at z = -50
z_en = round(Int, firstH/2) # end at z = -450 (smaller indices = deeper)
y_st = round(Int, ((pm.Lz-lastH)/pm.Tanα)/4) # find corresponding y value on slope when z = 250
y_en = round(Int, 2500/4) # just choosing this y value to include most of dye excursions

zlength_sm = length(z_en:z_st)
ylength_sm = length(y_st:y_en)

filesetnames = "SetList_mp.jld2"

scale_file = jldopen(filesetnames, "r+")

sns = scale_file["setnames"]
sfiles = scale_file["setfilenames"]

setnames = sns[1:22]
Lvals = length(setnames)

cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns
isemp(A) = length(A) < 1
isNemp(A) = length(A) >= 1
thresh(z) = z>10^(-6)

rolWidy = 20
rolWidz = 3

# number of individual results to log at each (y,t) location
numRes = 25

m = 5
setname  = setnames[m]

#name_prefix = "vIntWave_" * setname
name_prefix = sfiles[m] * setname
filepath = path_name * name_prefix * ".jld2"
@info "getting data from: " * setname

c_timeseries = FieldTimeSeries(filepath, "Cs");
xc, yc, zc = nodes(c_timeseries) #CCC
tlength = length(c_timeseries.times)
pm2 = getproperty(SimParams(), Symbol(setname))

pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
            Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
            nz = round(Int,pm2.Lz/2),
            m = -π/pm2.Lz,
            l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
            Tf = 2*π/pm2.f, 
            Tσ = 2*π/pm2.σ))

@info "Calculating Wave Indices..."
include("WaveValues.jl")
wave_info=get_wave_indices(c_timeseries, pm2, tlength)

# total number of data points
Wtlength = wave_info.Wl * wave_info.nTσ
# 4 waves long
cutWtlength = 4*wave_info.Wl
# starts at 7Tσ
W7length = wave_info.Wl * 7  + 1 
# thorpe starts at 3Tσ
W3length = wave_info.Wl * 3  + 1 
# ends at 10.9Tσ
W11length = wave_info.Wl * 11
LtcutWtlength = 8*wave_info.Wl
#for rolling wave avg, can't use the last half of wave if no more waves after that...
hWl = floor(Int64, wave_info.Wl/2)

# - 1/2 wave if doesn't fit since rolling 
if Wtlength > (W11length + hWl)
    # if room to go to end of "last" wave, then go!
    rWtlength = cutWtlength
else 
    rWtlength = cutWtlength - hWl
end

All_Cheights = zeros(cutWtlength, numRes)

function findroots(dcdz)
    roots = []
    for k = 2:length(dcdz)-2
        if dcdz[k] <= 0 && dcdz[k+1] >0
            roots = [roots;k]
        end
    end
    return roots
end

function find_intrusions(roots, cmi)
    if isNemp(roots)
        rt_widths = []
        rt_zst = []
        rt_zen = []
        w=1
        for p = 2:length(roots)
            maxc = maximum(cmi[roots[p-1]:roots[p]])
            minr = maxc*.75
            if maxc > 10^(-4) && cmi[roots[p]] <= minr && cmi[roots[p-1]] <= minr
                rng = roots[p-1]:roots[p]
                fr = findfirst(thresh, cmi[rng])
                lr = findfirst(thresh, reverse(cmi[rng]))
                newr1 = roots[p-1] + fr -1
                newr2 = roots[p] - lr
                z1 = zc[newr1 + z_en]
                z2 = zc[newr2 + z_en]
                new_wid = z2-z1
                rt_widths=[rt_widths; new_wid]
                rt_zst = [rt_zst; z1]
                rt_zen = [rt_zen; z2]
                w +=1
            end
        end
    end
    return rt_widths, rt_zst, rt_zen
end

apath = path_name * "Analysis/"

for (j, i) in enumerate((wave_info.WavePeriods[W7length:W11length]))
    # pulling c values within range
    c = c_timeseries[i]
    #c = c_timeseries[:,:,:,i]
    ci = mean(interior(c)[:,:, z_en:z_st], dims =1)[1,:,:]
    #ci = mean(c[:,:, z_en:z_st], dims =1)[1,:,:]

    # instead of a rolling window, average each row only if the value has any dye,   
    # see what that looks like compared to just straight averaging
    cmi_thresh = zeros(zlength_sm)
    cmi_nothresh = zeros(zlength_sm)
    cmi_nodelta_nothresh = zeros(zlength_sm)
    delta  = round(Int, (pm2.U₀/ pm2.Ñ) / 4)

    for (k, depth_idx) in enumerate(z_en:z_st) # at each depth
        y_slope_val = -zc[depth_idx]./ pm2.Tanα
        y_slope_idx = round(Int,y_slope_val/4)
        y_cutstart = maximum([y_st, y_slope_idx+4])
        y_delstart = maximum([y_st, y_slope_idx + delta])

        ci_k = ci[y_cutstart:y_en,k]
        ci_k_del = ci[y_delstart:y_en,k]

        cmi_thresh[k] =  mean(ci_k[ci_k .> 1e-16])
        cmi_nothresh[k] = mean(ci_k)
        cmi_nodelta_nothresh[k] = mean(ci_k_del)

    end

    Δc_thresh = cmi_thresh[1:end-1]-cmi_thresh[2:end]
    Δc_nothresh = cmi_nothresh[1:end-1]-cmi_nothresh[2:end]
    Δc_nodelta_nothresh = cmi_nodelta_nothresh[1:end-1]-cmi_nodelta_nothresh[2:end]

    dcdz_thresh = Δc_thresh /-2
    dcdz_nothresh = Δc_nothresh /-2
    dcdz_nodelta_nothresh = Δc_nodelta_nothresh /-2

    roots_thresh = findroots(dcdz_thresh)
    roots_nothresh = findroots(dcdz_nothresh)
    roots_nodelta_nothresh = findroots(dcdz_nodelta_nothresh)

    (intru_thresh, rt_zst_thresh, rt_zen_thresh)= find_intrusions(roots_thresh, cmi_thresh)
    (intru_nothresh, rt_zst_nothresh, rt_zen_nothresh)= find_intrusions(roots_nothresh, cmi_nothresh)
    (intru_nodelta_nothresh, rt_zst_nodelta_nothresh, rt_zen_nodelta_nothresh) = find_intrusions(roots_nodelta_nothresh, cmi_nodelta_nothresh)

    All_Cheights[j, 1:length(intru_nothresh)] = intru_nothresh

    
    f = Figure(resolution= (1300, 700), fontsize = 26)
        ga = f[1, 1] = GridLayout()

        ax1 = Axis(ga[1, 1], ylabel = "z", xlabel = "y")
            ax1.xticks = 500:1000:3000
            ax1.yticks = [-450, -250, -50]

            limits!(ax1, 0, 3000, -450, -50)
            
            heatmap!(ax1, yc, zc[z_en:z_st], log10.(clamp.(ci, 1e-8, 1)), colormap = :thermal, colorrange = (-4, 0))

        ax2 = Axis(ga[1, 2], xlabel = "c")
            hideydecorations!(ax2)
            lines!(ax2, cmi_nothresh, zc[z_en:z_st], linewidth = 3, color = :gray50)
            lines!(ax2, cmi_thresh, zc[z_en:z_st], linewidth = 3, color = :black)
            lines!(ax2, cmi_nodelta_nothresh, zc[z_en:z_st], linewidth = 3, color = :gray30, linestyle = :dash)

            limits!(ax2, 0.0, 0.09, -450, -50)
            colsize!(ga, 2, Relative(0.3))
            scatter!(ax2, cmi_thresh[roots_thresh], zc[z_en .+ roots_thresh], markersize = 10, marker = :circle, 
            color =:firebrick2)
            scatter!(ax2, cmi_nothresh[roots_nothresh], zc[z_en .+ roots_nothresh], markersize = 15, marker = :circle, 
            color =:dodgerblue2)
            scatter!(ax2, cmi_nodelta_nothresh[roots_nodelta_nothresh], zc[z_en .+ roots_nodelta_nothresh], markersize = 15, marker = :utriangle, 
            color =:dodgerblue4)
            hspan!(ax2, rt_zst_nodelta_nothresh, rt_zen_nodelta_nothresh, color = (:dodgerblue2, 0.2))
            hspan!(ax2, rt_zst_nothresh, rt_zen_nothresh, color = (:firebrick2, 0.2))
            hspan!(ax2, rt_zst_thresh, rt_zen_thresh, color = (:gray, 0.2))

    save(apath * "ContourTest$i" * setname * ".png", f)

    
end
