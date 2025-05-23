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
ylength_dissip = 880
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
sns_σ = scale_file["setnames_varyσ"][1:4]
sfiles = scale_file["setfilenames"]
sfiles_σ = scale_file["setfilenames_varyσ"][1:4]

setnames = vcat(sns, sns_σ)
setfiles = vcat(sfiles, sfiles_σ)
Lvals = length(setnames)

cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns
isemp(A) = length(A) < 1
isNemp(A) = length(A) >= 1
thresh(z) = z>10^(-6)

dissipthresh = 1e-8
rolWidy = 20
rolWidz = 3

# number of individual results to log at each (y,t) location
numRes = 25
numstats = 8
Thorpe_rms = zeros(Lvals)
Thorpe_stats = zeros(Lvals, numstats)
Nheight_stats = zeros(Lvals, numstats)
Dissip_stats = zeros(Lvals, numstats)
Cheight_stats = zeros(Lvals, numstats)

function ExtraStats(no0s)
    tmean = mean(no0s)
    tmedian = median(no0s)
    tmax = maximum(no0s)
    tmin = minimum(no0s)
    tnum = length(no0s)
    tstd = std(no0s)

    tsorted = sort(no0s)
    sorted_length = length(tsorted)
    eightieth_percentile_idx = round(Int, 0.8 * sorted_length)
    twentieth_percentile_idx = round(Int, 0.2 * sorted_length)
    eightieth_percentile = tsorted[eightieth_percentile_idx]
    twentieth_percentile = tsorted[twentieth_percentile_idx]

    extra_stats = [tmean, tmedian, tmax, tmin, tnum, tstd, twentieth_percentile, eightieth_percentile]
    return extra_stats
end

for (m, setname) in enumerate(setnames)
    
    #name_prefix = "vIntWave_" * setname
    name_prefix = setfiles[m] * setname
    filepath = path_name * name_prefix * ".jld2"
    @info "getting data from: " * setname

    b_timeseries = FieldTimeSeries(filepath,"b");
    xb, yb, zb = nodes(b_timeseries) #CCC
    c_timeseries = FieldTimeSeries(filepath, "Cs");
    xc, yc, zc = nodes(c_timeseries) #CCC
    N_timeseries = FieldTimeSeries(filepath, "N2");
    xn, yn, zn = nodes(N_timeseries) #CCC
    e_timeseries = FieldTimeSeries(filepath, "ϵ");
    xe, ye, ze = nodes(e_timeseries) #CCC

    tlength = length(b_timeseries.times)

    pm2 = getproperty(SimParams(), Symbol(setname))

    pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
                Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
                nz = round(Int,pm2.Lz/2),
                m = -π/pm2.Lz,
                l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
                Tf = 2*π/pm2.f, 
                Tσ = 2*π/pm2.σ))

    bnoise = 2*2*pm2.Ñ^2

    include("WaveValues.jl")
    wave_info=get_wave_indices(b_timeseries, pm2, tlength)
            
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
    #for rolling wave avg, can't use the last half of wave if no more waves after that...
    hWl = floor(Int64, wave_info.Wl/2)
    
    # - 1/2 wave if doesn't fit since rolling 
    if Wtlength > (W11length + hWl)
        # if room to go to end of "last" wave, then go!
        rWtlength = cutWtlength
    else 
        rWtlength = cutWtlength - hWl
    end

    LtcutWtlength = length(wave_info.WavePeriods[W3length:W11length])

    @info "Calculating Tracer Thicknesses..."

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
    
    @info "Calculating Thorpe Displacements..."
    # calculating RMS value over only the over turn at every time and y position
    All_Toverturns = zeros(LtcutWtlength, ylength-1, numRes)
    Displace_Lt = zeros(ylength-1, zlength, LtcutWtlength)

    b = interior(b_timeseries)[:,1:ylength, :, :];
    b_xavg = mean(b, dims = 1)[1,:,:,:];

    for (ni, i) in enumerate(wave_info.WavePeriods[W3length:W11length])

        #@info "time $ni/$LtcutWtlength..."
        b_xavgi = b_xavg[:,:, i]

        for (j,y) in enumerate(yb[2:ylength])
            #@info "y $j/$ylength..."
            top_cond(z) = z>=curvedslope(y)
            # first index above the slope, fy : end in the fluid
            fy = findfirst(top_cond, zb)
            # concatenate vertical profile and z values
            bz = cat(b_xavgi[j+1,:], zb, dims=2)
            # cut arrays to only include relevant z values after topo removes
            bz_cut = bz[fy:zlength,:]

            # sort both rows by buoyancy values
            bz_std = bz_cut[sortperm(bz_cut[:,1]),:]

            # z in bz_cut is actually their new loc, while bz_std stored old val. at idx
            # taking the new poisition of the parcel - old
            # + displacement means a downward move to stabilize
            #              d_o[k] =  z new[k] .- z old[k]  for k = slope:sfc
            Displace_Lt[j,fy:zlength, ni] =  bz_cut[:,2] .- bz_std[:,2]

            # now take the sum of these displacement from the sfc to each depth
            Sum_Disp = zeros(zlength)
            for k = fy:zlength
                Sum_Disp[k] = sum(Displace_Lt[j,k:end,ni])
            end

            # find all the values above threshold, defining overturns
            O_idxs = findall(cond_gt0, Sum_Disp)
            
            if isNemp(O_idxs)
                # find all the jumps in values 
                sO_idxs = findall(cond_gt1, O_idxs[2:end] - O_idxs[1:end-1])
                # find idx for 1st value that changed, (beginning of overturn)
                b_over_idxs = isNemp(sO_idxs) ? [O_idxs[1];O_idxs[sO_idxs .+ 1]] : O_idxs[1]
                e_over_idxs =  isNemp(sO_idxs) ? [O_idxs[sO_idxs]; O_idxs[end]] :  O_idxs[end]

                for o in 1:length(b_over_idxs)
                    bo = b_over_idxs[o]
                    eo = e_over_idxs[o]
                    overturn = Displace_Lt[j,bo-1:eo, ni]
                    overturnL = length(overturn)
                    # in the sorted buoyancy profile, bo will be more negative (denser) 
                    # and eo lighter (more buoyant) so b[eo] - b[bo] > 0 value
                    # compare to positive noise check
                    noisecheck = bz_std[eo-fy+1,1] .- bz_std[bo-fy,1] 
                    if noisecheck > bnoise
                        All_Toverturns[ni, j, o] = sqrt(sum((overturn).^2)/overturnL)
                    end
                end
            end

        end
        
    end
    
    @info "Smoothing Stratification Values..."
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

    @info "Calculating Unstable Heights..."
    All_Nheights = zeros(rWtlength, ylength_sm, numRes)

    for (i,l) in enumerate(W7length:W7length+rWtlength-1)

        for (yi, y_md) in enumerate(y_st:y_en)
            #@info "$yi"
            # find all the negative values of N²' 
            # indices will be based on smaller size domain with 50 cut off each end
            negs = []
            negs = findall(cond_lt0, ñ2i_pert_ryrzrW[yi,:,i])

            # assuming there was at least one value that was negative
            if isNemp(negs)
                # assuming there was at least one value that was negative
                top_cond(z) = z>=curvedslope(yn[y_md])
                if top_cond(zn[z_en])
                    # first index above the slope, fy : end in the fluid
                    fy = findfirst(top_cond, zn[z_en:z_st])
                    # cutting off any negative values that are within 4 indices of bottom (6m)
                    cutbot = sum( negs .< fy + 3)
                    negs = negs[cutbot+1:end]
                    # taking the "derivativ" of N^2' in this profile
                    Nanom_dif = ñ2i_pert_ryrzrW[yi,2:end,i] .- ñ2i_pert_ryrzrW[yi,1:end-1,i] 
                    dnegs = negs[2:end] .- negs[1:end-1]
                    # find all the places there was a jump in values
                    skip_negs = findall(cond_gt1, dnegs)
                    # the starting points will be the 1st neg + all the values after the jump
                    if isNemp(negs)
                        st_negs = [negs[1] ; negs[skip_negs.+1]]
                        en_negs = [negs[skip_negs] ; negs[end]]

                        for p = 1:length(st_negs)
                           # @info "$p"
                            # if the last negative isn't the last avaialble option
                            # and the first isn't the first available option
                            if (en_negs[p] < zlength_sm) & (st_negs[p] > 1)
                                #anom_en = en_negs[p] + findfirst(cond_gt0, banom[en_negs[p]+1:end])
                                anom_en = findfirst(cond_lte0, Nanom_dif[en_negs[p]+1:end])
                                #anom_st = findlast(cond_lte0,banom[1:st_negs[p]-1])
                                anom_st = findlast(cond_gt0,Nanom_dif[1:st_negs[p]-1])
                                if (typeof(anom_st) == Int64) & (typeof(anom_en) == Int64) 
                                    All_Nheights[i, yi, p] = zn[anom_en+en_negs[p]] .- zn[anom_st]
                                end
                            end
                        end

                        #neg_hts = zn[en_negs] .- zn[st_negs]
                        #for p = 1:length(neg_hts)
                        #    All_Nheights[i, yi, p] = neg_hts[p]
                        #end
                    end

                end   
            end 
            
        end
    end

    @info "Calculating dissipation Averages..."
    ei = interior(e_timeseries)[:,1:ylength_dissip,1:zlength,W7length:tlength];
    e_xavg = mean(ei, dims=1)[1,:,:,:];
    Ygrid = reshape(repeat(ye[1:ylength_dissip], zlength), ylength_dissip, zlength)
    SlopeGrid = curvedslope.(Ygrid)
    boolZ = (SlopeGrid .+ 4) .> ze[1:zlength]'
    # above slope values at each time step
    e_fvals = e_xavg[boolZ,:]
    e_abovethresh = e_fvals[e_fvals .> dissipthresh]

    @info "Calculating Dye Height Stats..."
    non0Cheights = All_Cheights[All_Cheights .> 2]
    Cheight_stats[m, :] = ExtraStats(non0Cheights)

    @info "Calculating Thorpe Stats..."
    overturns_overthreshold = All_Toverturns[cond_gt4.(All_Toverturns)]
    Thorpe_rms[m] = sqrt(mean(overturns_overthreshold.^2))
    Thorpe_stats[m, :] = ExtraStats(overturns_overthreshold)

    @info "Calculating start anom stats..."
    strat_anoms = All_Nheights[cond_gt0.(All_Nheights)]
    Nheight_stats[m, :] = ExtraStats(strat_anoms)

    @info "Calculating dissip stats..."
    Dissip_stats[m, :] = ExtraStats(e_abovethresh)

end

filescalename = apath * "DeltavAll_ExtraStats.jld2"

jldsave(filescalename; setnames, 
Thorpe_rms, Thorpe_stats, Nheight_stats, 
Dissip_stats, Cheight_stats)
