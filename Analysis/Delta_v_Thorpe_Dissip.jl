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

#### ------------Each Simulation and Array Info----------###
@info "Setting up Save Arrays "
filesetnames = "SetList_mp.jld2"

scale_file = jldopen(filesetnames, "r+")

setnames = scale_file["setnames"]
sfiles = scale_file["setfilenames"]

cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns
cond_gt10(y) = y > 10 # cutoff for overturns
isemp(A) = length(A) < 1
isNemp(A) = length(A) >= 1

# number of individual results to log at each (y,t) location
numRes = 25
bnoise = 2*2*pm.Ñ^2
Lvals = length(setnames)

# rms over all values greater than 10
thorpe_rms = zeros(Lvals)

# maximum at every time 
# average by phase
# maximum of the phase
MaxAvg_overturn_byphase = zeros(Lvals)

# dissipation calculated using thorpe scale needs to be "insitu"
Lt_eps_610_Wavg = zeros(Lvals)
Lt_eps_full_Wavg =  zeros(Lvals)
Lt_eps_full_usingrms =  zeros(Lvals)

start_time = time_ns()

for (m, setname) in enumerate(setnames)
    #    setname = setnames[m]
    @info "$setname"
    name_prefix = sfiles[m] * setname

    filepath = path_name * name_prefix * ".jld2"

    b_timeseries = FieldTimeSeries(filepath,"b");
    xb, yb, zb = nodes(b_timeseries) #CCC

    tlength = length(b_timeseries.times)

    pm2 = getproperty(SimParams(), Symbol(setname))

    pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
                Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
                nz = round(Int,pm2.Lz/2),
                m = -π/pm2.Lz,
                l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
                Tf = 2*π/pm2.f, 
                Tσ = 2*π/pm2.σ))

    include("WaveValues.jl")
    wave_info=get_wave_indices(b_timeseries, pm2, tlength)
        
    # thorpe starts at 3Tσ
    W3length = wave_info.Wl * 3  + 1 
    LtcutWtlength = length(W3length:tlength)

    @info "Calculating Thorpe Displacements..."
    # calculating RMS value over only the over turn at every time and y position
    All_Toverturns = zeros(LtcutWtlength, ylength-1, numRes)
    Displace_Lt = zeros(ylength-1, zlength, LtcutWtlength)

    b = interior(b_timeseries)[:,1:ylength, :, W3length:tlength];
    b_xavg = mean(b, dims = 1)[1,:,:,:];

    for (ni, i) in enumerate(W3length:tlength)

        @info "time $ni/$LtcutWtlength..."
        b_xavgi = b_xavg[:,:, ni]

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

    @info "Calculating Thorpe Stats..."

    Missingwaves_idxs = (wave_info.WavePeriods[1,4] - 1)
    cut_wave_idxs = wave_info.WavePeriods[:,4:11] .- Missingwaves_idxs
    cut_Widxs = findall(cut_wave_idxs[:,end] .> LtcutWtlength)
    cut_wave_idxs[cut_Widxs, end] .= LtcutWtlength
    # organizing overturns by wave period
    WaveBasedOverturns = All_Toverturns[cut_wave_idxs, :,:]
    Max_waveperiod_overturn =  maximum(All_Toverturns, dims = (2,3))[cut_wave_idxs]
    MaxAvg_overturn_byphase[m] = maximum(mean(Max_waveperiod_overturn, dims = 2)[:,1])

    # only using overturns that are greater than 10 m (5 grid points)
    overturns_overthreshold = All_Toverturns[cond_gt4.(All_Toverturns)]
    thorpe_rms[m] = sqrt(mean(overturns_overthreshold.^2))

    Tσ6_idx = wave_info.T_Tσs[7] - Missingwaves_idxs# completion of 6 waves
    Tσ10_idx = wave_info.T_Tσs[11] - Missingwaves_idxs# completion of 10 waves

    All_Toverturns_610 = All_Toverturns[Tσ6_idx:Tσ10_idx,:,:]
    All_Toverturns_610_overthreshold = All_Toverturns_610[cond_gt4.(All_Toverturns_610)]
    # Lₜ² N³
    Lt_eps_610_Wavg[m] = pm2.Ñ^3 * mean((All_Toverturns_610_overthreshold.^2))
    Lt_eps_full_Wavg[m] = pm2.Ñ^3 * mean(overturns_overthreshold.^2)
    Lt_eps_full_usingrms[m] = pm2.Ñ^3 * thorpe_rms[m]^2

end

filescalename = apath * "DeltavThorpeDissip.jld2"

jldsave(filescalename; setnames, 
thorpe_rms, MaxAvg_overturn_byphase,
Lt_eps_610_Wavg,
Lt_eps_full_Wavg,
Lt_eps_full_usingrms)

