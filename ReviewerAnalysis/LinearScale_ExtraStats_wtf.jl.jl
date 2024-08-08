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

Cheight_htmean = zeros(Lvals)
#median value if in a sorted array kind of
Cheight_htmedian = zeros(Lvals)
# finding the extrema and computing the mean
Cheight_htmiddle = zeros(Lvals)
Cheight_htmax = zeros(Lvals)
Cheight_htmin = zeros(Lvals)
Cheight_htnum = zeros(Lvals)
Cheight_htstd = zeros(Lvals)
Cheight_ht10p = zeros(Lvals)
Cheight_ht90p = zeros(Lvals)

for (m, setname) in enumerate(setnames)
    
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

    All_Cheights = zeros(cutWtlength, ylength_sm, numRes)

    for (j, i) in enumerate(wave_info.WavePeriods[W7length:W11length])

        # pulling c values within range
        c = c_timeseries[i]
        ci = mean(interior(c)[:,:, z_en:z_st], dims =1)[1,:,:]
    
        # averaging rolling window
        cmi = zeros(ylength_sm, zlength_sm)

        # yi = y index in cut range
        # y_md y index from total range
        for (yi, y_md) in enumerate(y_st:y_en)
            cmi[yi,:] = mean(ci[y_md-rolWidy:rolWidy+y_md, :], dims = 1)
        end
        Δc = cmi[:,1:end-1]-cmi[:,2:end]
        dcdz = Δc /-2
    
        # finding all "roots"
        for (yi, y_md) in enumerate(y_st:y_en)
            roots = []
            for k = 2:size(dcdz)[2]-2
                if dcdz[yi,k] <= 0 && dcdz[yi,k+1] >0
                    roots = [roots;k]
                end
            end
    
            if isNemp(roots)
    
                top_cond(z) = z>=curvedslope(yc[y_md])
                # first index above the slope, fy : end in the fluid
                fy = findfirst(top_cond, zc)
                # cutting off any roots that are within 4 indices of bottom (6m)
                cutbot = sum( roots .< fy + 3)
                roots = roots[cutbot+1:end]
    
                rt_widths = []
    
                w=1
    
                for p = 2:length(roots)
                    maxc = maximum(cmi[yi,roots[p-1]:roots[p]])
                    minr = maxc*.5 
                    if maxc > 10^(-4) && cmi[yi, roots[p]] <= minr && cmi[yi, roots[p-1]] <= minr
                        rng = roots[p-1]:roots[p]
                        fr = findfirst(thresh, cmi[yi,rng])
                        lr = findfirst(thresh, reverse(cmi[yi,rng]))
                        newr1 = roots[p-1] + fr -1
                        newr2 = roots[p] - lr
                        z1 = zc[newr1 + z_en]
                        z2 = zc[newr2 + z_en]
                        new_wid = z2-z1
                        rt_widths=[rt_widths; new_wid]
                        All_Cheights[j, yi, w] = new_wid
                        w +=1
                    end
                end
            end    
            
        end
    
    end
    
    Fin_numRes_idx = sum(cond_gt0.(maximum(All_Cheights, dims = (1,2))[1,1,:]))
    Cut_Cheights = All_Cheights[:,:, 1:Fin_numRes_idx]

    @info "Calculating Dye Height stats in time"
    
    non0Cheights = Cut_Cheights[Cut_Cheights .> 2]

    @info "Sorting and percentiles..."

    non0Cheights_sorted = sort(non0Cheights)
    non0Cheights_length = length(non0Cheights_sorted)
    ninetieth_percentile_idx = round(Int, 0.9 * non0Cheights_length)
    tenth_percentile_idx = round(Int, 0.1 * non0Cheights_length)
    ninetieth_percentile = non0Cheights_sorted[ninetieth_percentile_idx]
    tenth_percentile = non0Cheights_sorted[tenth_percentile_idx]

    eightieth_percentile_idx = round(Int, 0.8 * non0Cheights_length)
    twentieth_percentile_idx = round(Int, 0.2 * non0Cheights_length)
    eightieth_percentile = non0Cheights_sorted[eightieth_percentile_idx]
    twentieth_percentile = non0Cheights_sorted[twentieth_percentile_idx]

    if isNemp(non0Cheights)
        Cheight_htmean[m] = mean(non0Cheights)
        Cheight_htmedian[m] = median(non0Cheights)
        Cheight_htmiddle[m] = middle(non0Cheights)
        Cheight_htmax[m] = maximum(non0Cheights)
        Cheight_htmin[m] = minimum(non0Cheights)
        Cheight_htnum[m] = length(non0Cheights)
        Cheight_htstd[m] = std(non0Cheights)
        Cheight_ht10p[m] = tenth_percentile
        Cheight_ht90p[m] = ninetieth_percentile
    end

end

filescalename = apath * "DeltavDye_ExtraStats.jld2"

jldsave(filescalename; setnames, 
Cheight_htmean, Cheight_htmedian,
Cheight_htmiddle, Cheight_htmax,
Cheight_htmin, Cheight_htnum, Cheight_htstd, 
Cheight_ht90p, Cheight_ht10p)