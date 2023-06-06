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
 
path = "/glade/scratch/whitleyv/NewAdvection/Parameters/Res/"
apath = path * "Analysis/"

sn = "U250N100Lz100g100"

include("../parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center, gausT_width)) * heaviside(y-ySlopeSameˢ)

####### DISSIPATION #########

nz = round(Int,pm.Lz/2),

YL_ep = 880*4
ylength_ep = round(Int, YL_ep/Δy)

# maximum average dissipation over waves 6-10
Δhoriz = [4,4,8,6,3,6,7,3,3.2]
Δvert = [2,4,4,3,3,6,3.5,2,1.6]

resSh = Δhoriz ./ 4
resSz = Δvert ./ 2

Lvals = length(Δhoriz)
resnames = []
for idx in 1:Lvals
    resnames = vcat(resnames, @sprintf("_Rz%0.0fRh%0.0f", resSz[idx]*100, resSh[idx]*100))
end

resnames = sn .* resnames

# data sets
eps_endMaxAvg = zeros(Lvals)
eps_endAvg = zeros(Lvals)
thorpe_hmax_tavg = zeros(Lvals)
Cheight_havg_tavg = zeros(Lvals)

cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns
isemp(A) = length(A) < 1
isNemp(A) = length(A) > 1
thresh(z) = z>10^(-6)

# rolling averages width
rolWidy = 20
rolWidz = 3

# number of individual results to log at each (y,t) location
numRes = 25

for (m, resname) in enumerate(resnames)

    @info "getting data from: " * resname

    name_prefix = "IntWave_" * resname
    filepath = path * name_prefix * ".jld2"

    Δh = Δhoriz[m]
    Δz = Δvert[m]

    zlength = round(Int, pm.Lz/Δz)
    xlength = round(Int, pm.Lx/Δh)
    ylength = round(Int, pm.Ly/Δh)

    @info "Dissip..."
    YL_ep = 880*4
    ylength_ep = round(Int, YL_ep/Δh)
    
    e_timeseries = FieldTimeSeries(filepath, "ϵ");
    # only include simulation in count if it ran long enough
        
    xe, ye, ze = nodes(e_timeseries) #CCC

    ei = interior(e_timeseries)[:,1:ylength_ep,:,:];

    YXgrid = reshape(repeat(yb[2:ylength_ep], zlength*xlength), ylength_ep, zlength, xlength)
    SlopeGrid = curvedslope.(YXgrid)
    boolZYX = (SlopeGrid .+ (2*Δz)) .<= zb[2:zlength]' # all the values greater than slope
    boolZ_perm = permutedims(boolZYX, [3, 1,2])

    # above slope values at each time step
    e_fluidvals = ei[boolZ_perm,:]
    e_xyzmax = maximum(e_fluidvals,dims =1)[1,:]


    e_xyzavg_cutoff = zeros(tlength)

    for i = 1:tlength
        tbool = e_fluidvals[:,i] .> 1e-7 # changing up a order to -4 orders less than expected avg for smallest avg
        if sum(tbool) > 0
            e_xyzavg_cutoff[i] = mean(e_fluidvals[tbool, i])
        end
    end

    tlength = length(e_timeseries.times)

    @info "Computing Rolling Wave Averages..." 
    include("../WaveValues.jl")
    wave_info=get_wave_indices(e_timeseries, pm, tlength)
    Wl = wave_info.Wl

    @info "Averaging over waves 6-10 and 2-5..."
    # you have completed 1 wave by T_Tσs[2]
    # this is the time index for full array
    Tσ6_idx = wave_info.T_Tσs[7] # completion of 6 waves
    Tσ10_idx = wave_info.T_Tσs[11] # completion of 10 waves

    eps_endMaxAvg[m] = mean(e_xyzmax[Tσ6_idx:Tσ10_idx])

    eps_endAvg[m] = mean(e_xyzavg_cutoff[Tσ6_idx:Tσ10_idx])

    @info "Thicknesses..."
    YL_th = 425*4
    ylength_th = round(Int, YL_th/Δh)

    # cutting off top and bottom to avoid boundaries:
    lastH = 450 #start at z = -50 (used to be -250)
    firstH = 50 #end at z = -450
    # indices to start and stop 
    z_st = round(Int, lastH/Δz) # start at z = -50
    z_en = round(Int, firstH/Δh) # end at z = -450 (smaller indices = deeper)
    y_st = round(Int, ((pm.Lz-lastH)/pm.Tanα)/4) # find corresponding y value on slope when z = 250
    y_en = round(Int, 2500/Δh) # just choosing this y value to include most of dye excursions

    zlength_sm = length(z_en:z_st)
    ylength_sm = length(y_st:y_en)

    c_timeseries = FieldTimeSeries(filepath, "Cs");
    xc, yc, zc = nodes(c_timeseries) #CCC

    b_timeseries = FieldTimeSeries(filepath,"b");
    xb, yb, zb = nodes(b_timeseries) #CCC

    @info "Calculating Wave Indices..."

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

    @info "Calculating Dye Heights at every time..."

    All_Cheights = zeros(cutWtlength,ylength_sm, numRes)
    for (j, i) in enumerate(wave_info.WavePeriods[W7length:W11length])

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
        dcdz = Δc /-Δz
        # finding all "roots"
        for (yi, y_md) in enumerate(y_st:y_en)
            roots = []
            for k = 2:size(dcdz)[2]-Δz
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
    All_Toverturns = zeros(LtcutWtlength, ylength_th-1, numRes)

    for (ni, i) in enumerate(wave_info.WavePeriods[W3length:W11length])

        Displace_Lt = zeros(ylength_th-1, zlength)
        #@info "time $i/$tlength..."
        b = b_timeseries[i]
        bi = interior(b)[:, 1:ylength_th, :]
    
        # average in x you get in terms of y and z
        b_xavg = mean(bi, dims = 1)[1,:,:]
    
        for (j,y) in enumerate(yb[2:ylength_th])

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

    @info "Calculating Dye Height stats..."
    Cht_havg = zeros(cutWtlength)

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
    end

    non0_idx = findall(cond_gt0, Cht_havg)
    # averaging in time to get one value per delta
    Cheight_havg_tavg[m] = isemp(non0_idx) ? 0 : mean(Cht_havg[non0_idx])

    @info "Calculating Thorpe Stats..."
    Lt_hmax = zeros(tlength)

    for i = 1:LtcutWtlength
        n0_ct_tim = count(cond_lte0, All_Toverturns[i,:,:], dims=1)
        Lt_std = sort(All_Toverturns[i,:,:],dims=1)
        Lt_time_i = []
        # collecting all the widths at each time step individually
        for q = 1:numRes
            Lt_time_i = [Lt_time_i; Lt_std[n0_ct_tim[q]+1:end,1]]
        end
        Lt_hmax[i] = isemp(Lt_time_i) ? 0 : maximum(Lt_time_i)
    end

    non0_idx = findall(cond_gt0, Lt_havg[wave_info.T_Tσs[3]:end]) .+ wave_info.T_Tσs[3] .- 1
    # averaging in time to get one value per delta
    # only averaging over past the 3rd wave
    thorpe_hmax_tavg[m] = isemp(non0_idx) ? 0 : mean(Lt_hmax[non0_idx])

end

Ozmidov = sqrt.(eps_endAvg./pm.Ñ^3)


filescalename = apath * "ConcergenceAnalysis.jld2"

jldsave(filescalename; resnames, 
eps_endMaxAvg, eps_endAvg,
thorpe_hmax_tavg, Cheight_havg_tavg, 
Ozmidov, Δvert, Δhoriz)

