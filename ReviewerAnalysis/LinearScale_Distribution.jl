using Statistics
using Printf
using Oceananigans
using CurveFit
using ArgParse
using JLD2
using CairoMakie
###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

#apath = "/glade/scratch/whitleyv/NewAdvection/Parameters/"
#apath = "Plots"
path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"

#path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
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

#### ------------Each Simulation and Array Info----------###
@info "Setting up Save Arrays "

#filesetnames =  "../SetnamesList.jld2"
filesetnames = "SetList_mp.jld2"

scale_file = jldopen(filesetnames, "r+")

sns = scale_file["setnames"]

setnames = sns[1:22]

cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns
isemp(A) = length(A) < 1
isNemp(A) = length(A) >= 1
thresh(z) = z>10^(-6)

start_time = time_ns()

# rolling averages width
rolWidy = 20
rolWidz = 3

# number of individual results to log at each (y,t) location
numRes = 25

sfiles = scale_file["setfilenames"]

m = 5
setname = setnames[m]
    
#name_prefix = "vIntWave_" * setname
name_prefix = sfiles[m] * setname
filepath = path_name * name_prefix * ".jld2"
@info "getting data from: " * setname

c_timeseries = FieldTimeSeries(filepath, "Cs");
xc, yc, zc = nodes(c_timeseries) #CCC

b_timeseries = FieldTimeSeries(filepath,"b");
xb, yb, zb = nodes(b_timeseries) #CCC

N_timeseries = FieldTimeSeries(filepath, "N2");
xn, yn, zn = nodes(N_timeseries) #CCC

# move tlength here since some will be longer than 161 as usual...
tlength = length(N_timeseries.times)

# need to recalculate parameters each time because N, and wave length will change!
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
wave_info=get_wave_indices(N_timeseries, pm2, tlength)

# gather the last 4 waves in a usual sim
# this is 7Tσ - 10.9Tσ or waves 7-10
# in the array this is columns 8:11

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

### Calculation ###########################################################################

All_Cheights = zeros(cutWtlength,ylength_sm, numRes)

@info "Calculating Dye Heights at every time..."
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

# remove rows with nothing in them:
Fin_numRes_idx = sum(cond_gt0.(maximum(All_Cheights, dims = (1,2))[1,1,:]))
Cut_Cheights = All_Cheights[:,:, 1:Fin_numRes_idx]

@info "Calculating Dye Height stats..."
Cht_havg = zeros(cutWtlength)
Cht_hstd = zeros(cutWtlength)
Cht_hvar = zeros(cutWtlength)
length_intime = zeros(cutWtlength)
length_intime1 = zeros(cutWtlength)

for i = 1:cutWtlength
    n0_ct_tim = count(cond_lte0, Cut_Cheights[i,:,:], dims=1)
    Δz_std = sort(Cut_Cheights[i,:,:],dims=1)
    Δz_time_i = []
    # collecting all the widths at each time step individually
    for q = 1:Fin_numRes_idx
        Δz_time_i = [Δz_time_i; Δz_std[n0_ct_tim[q]+1:end,1]]
    end
    # only including overturn in rows that are non-zero
    Cht_havg[i] = isemp(Δz_time_i) ? 0 : mean(Δz_time_i)
    Cht_hstd[i] = isemp(Δz_time_i) ? 0 : std(Δz_time_i)
    Cht_hvar[i] = isemp(Δz_time_i) ? 0 : var(Δz_time_i)
    length_intime[i] = length(Δz_time_i)
end

wave_idxs_cut = wave_info.WavePeriods[W7length:W11length]
waves_cut = (c_timeseries.times[wave_idxs_cut])./pm2.Tσ

nonzero_idxs = findall(cond_gt0, Cut_Cheights)
Cut_Cheights_nonzero = Cut_Cheights[nonzero_idxs]
Cut_Cheights_mean = mean(Cut_Cheights_nonzero)

Cht_errorbars = 2 .* Cht_hstd ./ sqrt.(length_intime)
Cht_errorbars1 = 2 .* Cht_hstd1 ./ sqrt.(length_intime1)

# plot average in time:
using CairoMakie

f = Figure(resolution = (700, 700), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1],  ylabel = rich("L", subscript("tr"), " [m]"), 
    xlabel = "Tσ", xlabelsize=35, ylabelsize=35, )
    # we want a log y axis with ticks at -8 through 0
    lines!(ax1, waves_cut, Cht_havg, label = "42 m")
    lines!(ax1, waves_cut, Cht_havg1, label = "71 m")
    errorbars!(waves_cut, Cht_havg, Cht_errorbars, color = :red) 
    errorbars!(waves_cut, Cht_havg1, Cht_errorbars1, color = :red) 
    axislegend()
    savename = apath * "AvgTracerThickness_intime_werr"
save(savename * ".png", f) 

BINs = 10:10:250
# plot an entire distribution:
f = Figure(resolution = (700, 700), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1],  ylabel = rich("L", subscript("tr"), " [m]"), 
    xlabel = "Number", xlabelsize=25, ylabelsize=25, )
    #limits!(ax1, 0, 250, 0, 2000)
    # we want a log y axis with ticks at -8 through 0
    hist!(ax1, Cut_Cheights_nonzero, bins = BINs, strokewidth = 1, strokecolor = :black, label = "42 m")
    hist!(ax1, Cut_Cheights_nonzero1, bins = BINs, strokewidth = 1, strokecolor = :black, label = "71 m")
    axislegend()
    savename = apath * "TotalTracerThicknessDistribution"
save(savename * ".png", f) 

std(Cut_Cheights_nonzero)
mean(Cut_Cheights_nonzero)
var(Cut_Cheights_nonzero)
ErrorBar = 2 * std(Cut_Cheights_nonzero) / sqrt(length(Cut_Cheights_nonzero))


@info "Calculating Thorpe Displacements..."
# calculating RMS value over only the over turn at every time and y position
All_Toverturns = zeros(LtcutWtlength, ylength-1, numRes)

for (ni, i) in enumerate(wave_info.WavePeriods[W3length:W11length])

    Displace_Lt = zeros(ylength-1, zlength)
    @info "time $i/$tlength..."
    b = b_timeseries[i]
    bi = interior(b)[:, 1:ylength, :]

    # average in x you get in terms of y and z
    b_xavg = mean(bi, dims = 1)[1,:,:]

    for (j,y) in enumerate(yb[2:ylength])
        top_cond(z) = z>=curvedslope(y)
        # first index above the slope, fy : end in the fluid
        fy = findfirst(top_cond, zb)
        # concatenate vertical profile and z values
        bz = cat(b_xavg[j+1,:], zb, dims=2)
        # cut arrays to only include relevant z values after topo removes
        bz_cut = bz[fy:zlength,:]

        # sort both rows by buoyancy values
        bz_std = bz_cut[sortperm(bz_cut[:,1]),:]

        # z in bz_cut is actually their new loc, while bz_std stored old val. at idx
        # taking the new poisition of the parcel - old
        # + displacement means a downward move to stabilize
        #              d_o[k] =  z new[k] .- z old[k]  for k = slope:sfc
        Displace_Lt[j,fy:zlength] =  bz_cut[:,2] .- bz_std[:,2]

        # now take the sum of these displacement from the sfc to each depth
        Sum_Disp = zeros(zlength)
        for k = fy:zlength
            Sum_Disp[k] = sum(Displace_Lt[j,k:end])
        end

        # find all the values above threshold, defining overturns
        O_idxs = findall(cond_gt0, Sum_Disp)
        
        if isNemp(O_idxs)
            # find all the jumps in values 
            sO_idxs = findall(cond_gt1, O_idxs[2:end] - O_idxs[1:end-1])
            # find idx for 1st value that changed, (beginning of overturn)
            b_over_idxs = isNemp(sO_idxs) ? [O_idxs[1];O_idxs[sO_idxs .+ 1]] : O_idxs[1]
            e_over_idxs =  isNemp(sO_idxs) ? [O_idxs[sO_idxs]; O_idxs[end]] :  O_idxs[end]

            # if nothing then just set it to 1 so that it takes an index at all
            #Overturns_rmsL = ones(length(b_over_idxs))
            for o in 1:length(b_over_idxs)
                bo = b_over_idxs[o]
                eo = e_over_idxs[o]
                overturn = Displace_Lt[j,bo-1:eo]
                overturnL = length(overturn)
                #Overturns_rmsL[o] = sqrt(sum((overturn).^2)/overturnL)
                All_Toverturns[ni, j, o] = sqrt(sum((overturn).^2)/overturnL)

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

### STATISTICS ###########################################################################

@info "Calculating Dye Height stats..."
Cht_havg = zeros(cutWtlength)
Cht_hrms = zeros(cutWtlength)

for i = 1:cutWtlength
    n0_ct_tim = count(cond_lte0, All_Cheights[i,:,:], dims=1)
    Δz_std = sort(All_Cheights[i,:,:],dims=1)
    Δz_time_i = []
    # collecting all the widths at each time step individually
    for q = 1:numRes
        Δz_time_i = [Δz_time_i; Δz_std[n0_ct_tim[q]+1:end,1]]
    end
    # only including overturn in rows that are non-zero
    Cht_havg[i] = isemp(Δz_time_i) ? 0 : mean(Δz_time_i)
    Cht_hrms[i] = isemp(Δz_time_i) ? 0 : sqrt(sum(Δz_time_i.^2)/length(Δz_time_i))
end

non0_idx = findall(cond_gt0, Cht_havg)
# averaging in time to get one value per delta
Cheight_havg_tavg[m] = isemp(non0_idx) ? 0 : mean(Cht_havg[non0_idx])
Cheight_hrms_trms[m] = isemp(non0_idx) ? 0 : sqrt.(sum(Cht_hrms[non0_idx].^2)/length(non0_idx))

@info "Calculating Thorpe Stats..."
Lt_havg = zeros(LtcutWtlength)
Lt_hmax = zeros(LtcutWtlength)
Lt_med = zeros(LtcutWtlength)
for i = 1:LtcutWtlength
    n0_ct_tim = count(cond_lte0, All_Toverturns[i,:,:], dims=1)
    Lt_std = sort(All_Toverturns[i,:,:],dims=1)
    Lt_time_i = []
    # collecting all the widths at each time step individually
    for q = 1:numRes
        Lt_time_i = [Lt_time_i; Lt_std[n0_ct_tim[q]+1:end,1]]
    end
    Lt_havg[i] = isemp(Lt_time_i) ? 0 : mean(Lt_time_i)
    Lt_hmax[i] = isemp(Lt_time_i) ? 0 : maximum(Lt_time_i)
end

for i = 1:LtcutWtlength
    # time step
    overturns_t = All_Toverturns[i,:,:]
    # all overturns above 4m threshold
    overturns_overthreshold = overturns_t[cond_gt4.(overturns_t)]
    Lt_havg[i] = mean(overturns_overthreshold)
    Lt_hmax[i] = maximum(overturns_overthreshold)
    Lt_med[i] = median(overturns_overthreshold)
end

non0_idx = findall(cond_gt0, Lt_havg[wave_info.T_Tσs[3]:end]) .+ wave_info.T_Tσs[3] .- 1
# averaging in time to get one value per delta
# only averaging over past the 3rd wave
thorpe_havg_tavg[m] = isemp(non0_idx) ? 0 : mean(Lt_havg[non0_idx])
thorpe_hmax_tavg[m] = isemp(non0_idx) ? 0 : mean(Lt_hmax[non0_idx])

@info "Calculating Unstable Strat Height stats..."
Nht_havg = zeros(rWtlength)
Nht_hrms = zeros(rWtlength)

# for each time
for i = 1:rWtlength
    # count up the number of indices at each height save spot (1:25)
    n0_ct_tim = count(cond_lte0, All_Nheights[i,:,:], dims=1)
    # sort the y value heights under each save spot so that these values are at the beginning
    Δz_std = sort(All_Nheights[i,:,:],dims=1)
    Δz_time_i = []
    # collecting all the heights height spot individually
    for q = 1:numRes
        # at each new height possibility collect the non-zero results pver all y places 
        Δz_time_i = [Δz_time_i; Δz_std[n0_ct_tim[q]+1:end,1]]
    end
    # only including heights in rows that are non-zero
    # averaging over all heights saved in that time step
    Nht_havg[i] = isemp(Δz_time_i) ? 0 : mean(Δz_time_i)
    Nht_hrms[i] = isemp(Δz_time_i) ? 0 : sqrt(sum(Δz_time_i.^2)/length(Δz_time_i))
end

non0_idx = findall(cond_gt0, Nht_havg)
# averaging in time to get one value per delta
Nheight_havg_tavg[m] = isemp(non0_idx) ? 0 : mean(Nht_havg[non0_idx])
Nheight_hrms_trms[m] = isemp(non0_idx) ? 0 : sqrt.(sum(Nht_hrms[non0_idx].^2)/length(non0_idx))

