using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using CurveFit
using ArgParse
using JLD2

ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSame = 1336.6                           # point where planar and curved corner math up the best    
 
apath = "/glade/scratch/whitleyv/NewAdvection/Parameters/"

sn = "U250Nfd250Lz100g100"

include("parameters.jl")
#include("../parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ
y_start = -(pm.Lz - pm.Lzˢ)/pm.Tanα
Ly = pm.Lyˢ-y_start
ny = round(Int,Ly/4)

pm = merge(pm, (;Ly=Ly,ny=ny))

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center, gausT_width)) * heaviside(y-ySlopeSame)

# how many times were saved?
zlength = pm.nz
ylength = 425 # roughly out to y = 1700 (is this enough still??)
xlength = 38

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

#### ------------Each Simulation and Array Info----------####
global setnames=[]
for u = 50:50:550
    global setnames = [setnames ; @sprintf("U250Nfd%dLz100g100", u)]
end
for u = 50:50:550
    global setnames = [setnames ; @sprintf("U%dN100Lz100g100", u)]
end
@info "Setting up Save Arrays "

Lvals = length(setnames)

################################ INTRUSIONS

# averaging over all times and space vs delta
Cheight_havg_tavg = zeros(Lvals)
Cheight_hrms_trms = zeros(Lvals)

################################  THORPE SCALE

# mean and max over all values for each delta
thorpe_havg_tavg = zeros(Lvals)
thorpe_hmax_tavg = zeros(Lvals)

################################ Isopycnals

# averaging over all times and space vs delta
Nheight_havg_tavg = zeros(Lvals)
Nheight_hrms_trms = zeros(Lvals)

cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns
isemp(A) = length(A) < 1
isNemp(A) = length(A) > 1
isnada(val) = val == nothing
thresh(z) = z>10^(-6)

start_time = time_ns()

# rolling averages width
rolWidy = 20
rolWidz = 3

# number of individual results to log at each (y,t) location
numRes = 25

for (m, setname) in enumerate(setnames)

    if m > 11
        path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU0/"
    else
        path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryNetc/"
    end
    
    name_prefix = "vIntWave_" * setname
    filepath = path_name * name_prefix * ".jld2"
    @info "getting data from: " * setname

    c_timeseries = FieldTimeSeries(filepath, "c")
    xc, yc, zc = nodes(c_timeseries) #CCC

    b_timeseries = FieldTimeSeries(filepath,"b")
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

    @info "Calculating Thorpe Displacements..."
    # calculating RMS value over only the over turn at every time and y position
    All_Toverturns = zeros(LtcutWtlength, ylength-1, numRes)

    for (ni, i) in enumerate(wave_info.WavePeriods[W3length:W11length])

        Displace_Lt = zeros(ylength-1, zlength)
        #@info "time $i/$tlength..."
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
            O_idxs = findall(cond_gt4, Sum_Disp)
            
            if isNemp(O_idxs)
                # find all the jumps in values 
                sO_idxs = findall(cond_gt1, O_idxs[2:end] - O_idxs[1:end-1])
                # find idx for 1st value that changed, (beginning of overturn)
                b_over_idxs = isNemp(sO_idxs) ? [O_idxs[1];O_idxs[sO_idxs .+ 1]] : O_idxs[1]
                e_over_idxs =  isNemp(sO_idxs) ? [O_idxs[sO_idxs]; O_idxs[end]] :  O_idxs[end]
    
                # if nothing then just set it to 1 so that it takes an index at all
                Overturns_rmsL = ones(length(b_over_idxs))
                for o in 1:length(b_over_idxs)
                    bo = b_over_idxs[o]
                    eo = e_over_idxs[o]
                    overturn = Displace_Lt[j,bo:eo]
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
            # find all the negative values of N²'
            negs = []
            negs = findall(cond_lt0, ñ2i_pert_ryrzrW[yi,:,i])

            # assuming there was at least one value that was negative
            if isNemp(negs)
                # assuming there was at least one value that was negative
                top_cond(z) = z>=curvedslope(yn[y_md])
                # first index above the slope, fy : end in the fluid
                fy = findfirst(top_cond, zn)
                # cutting off any negative values that are within 4 indices of bottom (6m)
                cutbot = sum( negs .< fy + 3)
                negs = negs[cutbot+1:end]

                dnegs = negs[2:end] .- negs[1:end-1]
                # find all the places there was a jump in values
                skip_negs = findall(cond_gt1, dnegs)
                # the starting points will be the 1st neg + all the values after the jump
                if isNemp(skip_negs)
                    st_negs = [negs[1] ; negs[skip_negs.+1]]
                    en_negs = [negs[skip_negs] ; negs[end]]
                    neg_hts = zn[en_negs] .- zn[st_negs]
                    for p = 1:length(neg_hts)
                        All_Nheights[i, yi, p] = neg_hts[p]
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
    Lt_havg = zeros(tlength)
    Lt_hmax = zeros(tlength)

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

    @info "Statistics Finished for dataset $m/$Lvals"
end


# plot statistics compared to delta values
δ = (0.05:.05:.55)./pm.Ñ
δ2 = repeat(δ, 2)

filescalename = apath * "DeltavAllScale.jld2"

jldsave(filescalename; setnames, 
    thorpe_havg_tavg, thorpe_hmax_tavg, 
    Cheight_havg_tavg, Cheight_hrms_trms,
    Nheight_havg_tavg, Nheight_hrms_trms,
    δ = δ2)