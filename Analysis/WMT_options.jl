using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

apath = path_name * "Analysis/"

include("parameters.jl")

setname = "U300N100Lz100g100"

# how many times were saved?
zlength = 250
ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)
xlength = 38
tlength = 161 

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

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

name1_prefix = "IntWave_" * setname
name2_prefix = "IntWave_smdt_" * setname
name3_prefix = "IntWave_postspin_" * setname

filepath1 = path_name * name1_prefix * ".jld2"
filepath2 = path_name * name2_prefix * ".jld2"
#filepath3 = path_name * name3_prefix * ".jld2"
filepath3 = path_name * "cgl_" * setname * ".jld2"

f3 = jldopen(filepath3)
CGli = f3["Cgli"];

@info "getting data from: " * setname

#########################
#                   BOUYANCY BINNING
#########################
zlength = 250 # only want slope values for those that include upp region
xlength = 38 
ylength = 1000 # out to 4000 m 

b_timeseries = FieldTimeSeries(filepath1,"b");
xb, yb, zb = nodes(b_timeseries) #CCC
# dissipation
e_timeseries = FieldTimeSeries(filepath1,"ϵ");
# gaussian at slope but with increased time resolution
Cg_timeseries = FieldTimeSeries(filepath2,"Cg");
Cs_timeseries = FieldTimeSeries(filepath2,"Cs");
Cgr_timeseries = FieldTimeSeries(filepath1,"Cgr");
#Cgl_timeseries = FieldTimeSeries(filepath3, "Cg");

CSi = interior(Cs_timeseries)[:, 1:ylength, 1:zlength,1:2:end];
CGi = interior(Cg_timeseries)[:, 1:ylength, 1:zlength,1:2:end];
CGri = interior(Cgr_timeseries)[:, 1:ylength, 1:zlength,:];
#Cgli = interior(Cgl_timeseries)[:, 1:ylength, 1:zlength,:];
ei = interior(e_timeseries)[:,1:ylength,1:zlength,:];
Bi = interior(b_timeseries)[:, 1:ylength, 1:zlength,:];
# Cgl2i = interior(Cg_timeseries[:, 1:1000, :,:];
@info "Computing Volumes..."
ini_Nisos = 40 # number of bins
Δb = -1.47*1e-4 #-500*pm.Ñ^2/ini_Nisos # "width" of bins

function VolumeBinning(ini_Nisos, tlength, Δb, Bi, CGi, CSi, ei)
    # buoynacy sum
    bsum = zeros(ini_Nisos, tlength)
    # total concentration and dissip
    contotG = zeros(ini_Nisos, tlength)
    contotS = zeros(ini_Nisos, tlength)
    etot = zeros(ini_Nisos, tlength)


    for j = 1:tlength
        @info "Time $j of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bj = Bi[:,:,:,j]
        cjG = CGi[:,:,:,j]
        cjS = CSi[:,:,:,j]
        ej = ei[:,:,:,j]

        # for each buoynacy class:
        for n = 1:ini_Nisos
            # (1) CCC locations in the density class
            boolB = (bj.< Δb*(n-1)) .& (bj.>= Δb *n)
            # finding the volume of that density class
            bsum[n,j]=sum(boolB)
            # dye in layer
            cS_inb = cjS[boolB]
            cG_inb = cjG[boolB]
            # dissip in layer
            e_inb = ej[boolB]
            # total concent in n
            contotG[n,j] = sum(cG_inb)
            contotS[n,j] = sum(cS_inb)
            # total dissip in n
            etot[n,j] = sum(e_inb)
            # average dye in n
        # cavg[n,j] = mean(cS_inb)
            # average dissip in n
            #eavg[n,j] = mean(e_inb)
        end
    end

    # change in volume from initial
    ΔVol = (bsum .- bsum[:,1]).* 16
    # change in concentration from initial
    ΔConG = (contotG .- contotG[:,1] ) .* 16
    ΔConS = (contotS .- contotS[:,1] ) .* 16

    reverse!(ΔVol, dims=1)
    reverse!(ΔConG, dims=1)
    reverse!(ΔConS, dims=1)
    reverse!(etot, dims=1)

    return ΔVol, ΔConG, ΔConS, etot
end

#####################
## Only including values within 1 delta of the slope
####################

function VolumeBinning(ini_Nisos, tlength, Δb, Bi, CGi, CSi, ei, δdep, ylength, zlength, xlength, yb, zb)

    YXgrid = reshape(repeat(yb[1:ylength], zlength*xlength), ylength, zlength, xlength)
    SlopeGridYX_pd = curvedslope.(YXgrid) .+ δdep
    boolYX = SlopeGridYX_pd .>= zb[1:zlength]';  # all the values less than slope + delta
    boolYX_perm = permutedims(boolYX, [3, 1, 2])

    # buoynacy sum
    bsum = zeros(ini_Nisos, tlength)
    # total concentration and dissip
    contotG = zeros(ini_Nisos, tlength)
    contotS = zeros(ini_Nisos, tlength)
    etot = zeros(ini_Nisos, tlength)

    for j = 1:tlength
        @info "Time $j of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bj = Bi[:,:,:,j]
        cjG = CGi[:,:,:,j]
        cjS = CSi[:,:,:,j]
        ej = ei[:,:,:,j]

        bj_slope = bj[boolYX_perm]
        cjG_slope = cjG[boolYX_perm]
        cjS_slope = cjS[boolYX_perm]
        ej_slope = ej[boolYX_perm]
    
        # for each buoynacy class:
        for n = 1:ini_Nisos
            # (1) CCC locations in the density class
            boolB = (bj_slope.< Δb*(n-1)) .& (bj_slope.>= Δb *n)
            # finding the volume of that density class
            bsum[n,j]=sum(boolB)
            # dye in layer
            cS_inb = cjS_slope[boolB]
            cG_inb = cjG_slope[boolB]
            # dissip in layer
            e_inb = ej_slope[boolB]
            # total concent in n
            contotG[n,j] = sum(cG_inb)
            contotS[n,j] = sum(cS_inb)
            # total dissip in n
            etot[n,j] = sum(e_inb)
        end    
    end

    # change in volume from initial
    ΔVol = (bsum .- bsum[:,1]).* 16
    # change in concentration from initial
    ΔConG = (contotG .- contotG[:,1] ) .* 16
    ΔConS = (contotS .- contotS[:,1] ) .* 16

    reverse!(ΔVol, dims=1)
    reverse!(ΔConG, dims=1)
    reverse!(ΔConS, dims=1)
    reverse!(etot, dims=1)

    return ΔVol, ΔConG, ΔConS, etot
end

(ΔVol_full, ΔConG_full, ΔConS_full, etot_full) = VolumeBinning(ini_Nisos, tlength, Δb, Bi, CGi, CSi, ei)
(ΔVol_1d, ΔConG_1d, ΔConS_1d, etot_1d) = VolumeBinning(ini_Nisos, tlength, Δb, Bi, CGi, CSi, ei,
                (pm.U₀/pm.Ñ), ylength, zlength, xlength, yb, zb)
(ΔVol_025d, ΔConG_025d, ΔConS_025d, etot_025d) = VolumeBinning(ini_Nisos, tlength, Δb, Bi, CGi, CSi, ei,
                0.25*(pm.U₀/pm.Ñ), ylength, zlength, xlength, yb, zb)

@info "Calculating Wave Indices..."
include("WaveValues.jl")

wave_info=get_wave_indices(b_timeseries, pm, tlength)

#########################
#                   Phase Averaging Buoyancy Bins
#########################

function PhaseAveraging(wave_info, ΔVol, ΔConG, ΔConS, etot, waveframe)
    #####################
    # Phase averaging 
    #####################
    end5waves = wave_info.WavePeriods[:,waveframe]

    ΔVol_Ph = ΔVol[:,end5waves];
    ΔConG_Ph = ΔConG[:,end5waves];
    ΔConS_Ph = ΔConS[:,end5waves];
    etot_Ph = etot[:, end5waves];

    ΔVol_Phavg = mean(ΔVol_Ph, dims = 5)[:,:,1];
    ΔConG_Phavg = mean(ΔConG_Ph, dims = 5)[:,:,1];
    ΔConS_Phavg = mean(ΔConS_Ph, dims = 5)[:,:,1];
    etot_Phavg = mean(etot_Ph, dims = 5)[:,:,1];

    return ΔVol_Phavg, ΔConG_Phavg, ΔConS_Phavg, etot_Phavg
end

(ΔVol_Phavg_full, ΔConG_Phavg_full, ΔConS_Phavg_full, etot_Phavg_full) = PhaseAveraging(wave_info, ΔVol_full, ΔConG_full, ΔConS_full, etot_full, 7:10)
(ΔVol_Phavg_full_b, ΔConG_Phavg_full_b, ΔConS_Phavg_full_b, etot_Phavg_full_b) = PhaseAveraging(wave_info, ΔVol_full, ΔConG_full, ΔConS_full, etot_full, 1:3)

(ΔVol_Phavg_1d, ΔConG_Phavg_1d, ΔConS_Phavg_1d, etot_Phavg_1d) = PhaseAveraging(wave_info, ΔVol_1d, ΔConG_1d, ΔConS_1d, etot_1d, 7:10)
(ΔVol_Phavg_1d_b, ΔConG_Phavg_1d_b, ΔConS_Phavg_1d_b, etot_Phavg_1d_b) = PhaseAveraging(wave_info, ΔVol_1d, ΔConG_1d, ΔConS_1d, etot_1d, 1:3)

(ΔVol_Phavg_025d, ΔConG_Phavg_025d, ΔConS_Phavg_025d, etot_Phavg_025d) = PhaseAveraging(wave_info, ΔVol_025d, ΔConG_025d, ΔConS_025d, etot_025d, 7:10)
(ΔVol_Phavg_b, ΔConG_Phavg_025d_b, ΔConS_Phavg_025d_b, etot_Phavg_025d_b) = PhaseAveraging(wave_info, ΔVol_025d, ΔConG_025d, ΔConS_025d, etot_025d, 1:3)

phase_times = b_timeseries.times[wave_info.WavePeriods[:,1]]/pm.Tσ

# buoynacy classes
# first index is the buoynacy found at lowest depth
# last index is buoyancy found at top ie. (Δb/2)
isos = Δb.*((ini_Nisos - .5):-1:0.5)

function phaseavg_volume_plot(phase_times, isos, ΔVol, ΔConG, ΔConS, etot, savename, apath,
    vmax, cmaxG, cmaxS, emax)
    f = Figure(resolution = (1600, 900), fontsize=26)
    ga = f[1, 1] = GridLayout() 

    # change in volume
    axvdif = Axis(ga[1, 1], ylabel = "Buoyancy [ms⁻²]", 
                title = "ΔVol = Vol - Vol₀")
    axesum = Axis(ga[1, 3],  title = "ε")
    axcdifG = Axis(ga[2, 1], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",
                title = "Gauss Δc = ∑c - ∑c₀")
    axcdifS = Axis(ga[2, 3], xlabel = "Wave Periods [Tσ]",
                title = "Tanh Δc = ∑c - ∑c₀")

    axcdifG.xticks =  .2:.2:1
    axcdifS.xticks =  .2:.2:1

    axvdif.yticks =   (-4e-3:1e-3:-2e-3,  ["-4×10⁻³", "-3×10⁻³", "-2×10⁻³"] )
    axcdifG.yticks =  (-4e-3:1e-3:-2e-3,  ["-4×10⁻³", "-3×10⁻³", "-2×10⁻³"] )

    # vol  []     e tot   e avg []]
    # c totG []]  c totS  c avg []

    limits!(axvdif, 0, 1, -5e-3, -1e-3)
    limits!(axcdifS, 0, 1, -5e-3, -1e-3)
    limits!(axcdifG, 0, 1, -5e-3, -1e-3)
    limits!(axesum, 0, 1, -5e-3, -1e-3)

    hidexdecorations!(axvdif)
    hidedecorations!(axesum)
    hideydecorations!(axcdifS)

    # just the slope values
    hvdif = heatmap!(axvdif, phase_times, isos, ΔVol', colormap = :balance, colorrange = (-vmax, vmax))
    hcGdif = heatmap!(axcdifG, phase_times, isos, ΔConG', colormap = :balance, colorrange = (-cmaxG, cmaxG))
    hcSdif = heatmap!(axcdifS, phase_times, isos, ΔConS', colormap = :balance, colorrange = (-cmaxS, cmaxS))
    hetot = heatmap!(axesum, phase_times, isos, etot', colormap = :thermal, colorrange = (0, emax))

    vtic = (vmax/2)*1e-5
    cticG = (cmaxG/2)*1e-4
    cticS = (cmaxS/2)*1e-4
    etic = (emax/2)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(ga[1,2], hvdif, ticks = (-vmax/2:vmax/2:vmax/2, [@sprintf("-%0.1f×10⁵", vtic), "0", @sprintf("%0.1f×10⁵", vtic)] ),
    size =25)
    cb2 = Colorbar(ga[2,2], hcGdif, ticks = (-cmaxG/2:cmaxG/2:cmaxG/2, [@sprintf("-%0.1f×10⁴", cticG), "0", @sprintf("%0.1f×10⁴", cticG)] ), 
    size =25)
    cb3 = Colorbar(ga[1,4], hetot, ticks = (0:etic:2*etic),
    size =25)
    cb4 = Colorbar(ga[2,4], hcSdif, ticks = (-cmaxS/2:cmaxS/2:cmaxS/2, [@sprintf("-%0.1f×10⁴", cticS), "0", @sprintf("%0.1f×10⁴", cticS)] ), 
    size =25)

    colsize!(ga, 2, Relative(0.01))
    colsize!(ga, 4, Relative(0.01))

    rowgap!(ga, 8)

    save(apath * savename * ".png", f)

end

phaseavg_volume_plot(phase_times, isos, ΔVol_Phavg_1d, ΔConG_Phavg_1d, ΔConS_Phavg_1d, etot_Phavg_1d,
 "bClasses_slopeend_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
 5e5, 1e4, 4e4, 0.05)
phaseavg_volume_plot(phase_times, isos, ΔVol_Phavg_1d_b, ΔConG_Phavg_1d_b, ΔConS_Phavg_1d_b, etot_Phavg_1d_b,
 "bClasses_slopebeg_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
 2e5, 5e3, 2e4, 0.05)

phaseavg_volume_plot(phase_times, isos, ΔVol_Phavg_full, ΔConG_Phavg_full, ΔConS_Phavg_full, etot_Phavg_full,
 "bClasses_end_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
 7e5, 5e4, 5e4, 0.05)
phaseavg_volume_plot(phase_times, isos, ΔVol_Phavg_full_b, ΔConG_Phavg_full_b, ΔConS_Phavg_full_b, etot_Phavg_full_b,
 "bClasses_beg_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
 4e5, 2e4, 2e4, 0.05)

phaseavg_volume_plot(phase_times, isos, ΔVol_Phavg_025d, ΔConG_Phavg_025d, ΔConS_Phavg_025d, etot_Phavg_025d,
 "bClasses_025d_end_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
 5e4, 5e3, 5e4, 0.05)
phaseavg_volume_plot(phase_times, isos, ΔVol_Phavg_b, ΔConG_Phavg_025d_b, ΔConS_Phavg_025d_b, etot_Phavg_025d_b,
 "bClasses_025d_beg_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
 1e4, 1e3, 5e3, 0.002)

#########################
#                   Wave averaging Volume Plots
#########################

function WaveAveraging(wave_info, ΔVol, ΔConG, ΔConS, etot)
    WL = wave_info.Wl
    nTσ = wave_info.nTσ

    Wtlength = length(WL+1:WL*nTσ-WL)

    # volume change since initial
    ΔVol_rWavg = zeros(ini_Nisos, Wtlength)
    # change in total concentration since initial
    ΔConG_rWavg = zeros(ini_Nisos, Wtlength)
    ΔConS_rWavg = zeros(ini_Nisos, Wtlength)
    # total dissipation
    etot_rWavg = zeros(ini_Nisos, Wtlength)

    for k = (WL+1):(WL*nTσ - WL)
        WaveIndices = wave_info.WavePeriods[k-WL:k+WL]
        # change in total volume since initial
        ΔVol_rWavg[:,k-WL] = mean(ΔVol[:,WaveIndices], dims=2)
        # change in total concentration since initial
        ΔConG_rWavg[:,k-WL] = mean(ΔConG[:,WaveIndices], dims=2)
        ΔConS_rWavg[:,k-WL] = mean(ΔConS[:,WaveIndices], dims=2)
        # total issipation
        etot_rWavg[:,k-WL] = mean(etot[:,WaveIndices], dims=2)

    end
    return ΔVol_rWavg, ΔConG_rWavg, ΔConS_rWavg, etot_rWavg
end

(ΔVol_rWavg_full, ΔConG_rWavg_full, ΔConS_rWavg_full, etot_rWavg_full) = WaveAveraging(wave_info, ΔVol_full, ΔConG_full, ΔConS_full, etot_full)
(ΔVol_rWavg_025d, ΔConG_rWavg_025d, ΔConS_rWavg_025d, etot_rWavg_025d) = WaveAveraging(wave_info, ΔVol_025d, ΔConG_025d, ΔConS_025d, etot_025d)

# wave period times for rolling average
rWtimes = b_timeseries.times[wave_info.WavePeriods[WL+1:WL*nTσ - WL]]/pm.Tσ

function volume_plot(times, isos, ΔVol, ΔConG, ΔConS, etot, savename, apath,
    vmax, cmaxG, cmaxS, emax)
    f = Figure(resolution = (1600, 900), fontsize=26)
    ga = f[1, 1] = GridLayout() 

    # change in volume
    axvdif = Axis(ga[1, 1], ylabel = "Buoyancy [ms⁻²]", 
                title = "ΔVol = Vol - Vol₀")
    axesum = Axis(ga[1, 3],  title = "ε")
    axcdifG = Axis(ga[2, 1], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",
                title = "Gauss Δc = ∑c - ∑c₀")
    axcdifS = Axis(ga[2, 3], xlabel = "Wave Periods [Tσ]",
                title = "Tanh Δc = ∑c - ∑c₀")

    axcdifG.xticks =  2:2:10
    axcdifS.xticks =  2:2:10

    axvdif.yticks =   (-4e-3:1e-3:-2e-3,  ["-4×10⁻³", "-3×10⁻³", "-2×10⁻³"] )
    axcdifG.yticks =  (-4e-3:1e-3:-2e-3,  ["-4×10⁻³", "-3×10⁻³", "-2×10⁻³"] )

    # vol  []     e tot   e avg []]
    # c totG []]  c totS  c avg []

    limits!(axvdif, 0, 11, -5e-3, -1e-3)
    limits!(axcdifS, 0, 11, -5e-3, -1e-3)
    limits!(axcdifG, 0, 11, -5e-3, -1e-3)
    limits!(axesum, 0, 11, -5e-3, -1e-3)

    hidexdecorations!(axvdif)
    hidedecorations!(axesum)
    hideydecorations!(axcdifS)

    # just the slope values
    hvdif = heatmap!(axvdif, times, isos, ΔVol', colormap = :balance, colorrange = (-vmax, vmax))
    hcGdif = heatmap!(axcdifG, times, isos, ΔConG', colormap = :balance, colorrange = (-cmaxG, cmaxG))
    hcSdif = heatmap!(axcdifS, times, isos, ΔConS', colormap = :balance, colorrange = (-cmaxS, cmaxS))
    hetot = heatmap!(axesum, times, isos, etot', colormap = :thermal, colorrange = (0, emax))

    vtic = (vmax/2)*1e-5
    cticG = (cmaxG/2)*1e-4
    cticS = (cmaxS/2)*1e-4
    etic = (emax/2)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(ga[1,2], hvdif, ticks = (-vmax/2:vmax/2:vmax/2, [@sprintf("-%0.1f×10⁵", vtic), "0", @sprintf("%0.1f×10⁵", vtic)] ),
    size =25)
    cb2 = Colorbar(ga[2,2], hcGdif, ticks = (-cmaxG/2:cmaxG/2:cmaxG/2, [@sprintf("-%0.1f×10⁴", cticG), "0", @sprintf("%0.1f×10⁴", cticG)] ), 
    size =25)
    cb3 = Colorbar(ga[1,4], hetot, ticks = (0:etic:2*etic),
    size =25)
    cb4 = Colorbar(ga[2,4], hcSdif, ticks = (-cmaxS/2:cmaxS/2:cmaxS/2, [@sprintf("-%0.1f×10⁴", cticS), "0", @sprintf("%0.1f×10⁴", cticS)] ), 
    size =25)

    colsize!(ga, 2, Relative(0.01))
    colsize!(ga, 4, Relative(0.01))

    rowgap!(ga, 8)

    save(apath * savename * ".png", f)
end

volume_plot(b_timeseries.times./pm.Tσ, isos, ΔVol_full, ΔConG_full, ΔConS_full, etot_full, 
"bClasses_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
1.5e6, 1.5e4, 7.5e4, 0.08)

volume_plot(b_timeseries.times./pm.Tσ, isos, ΔVol_025d, ΔConG_025d, ΔConS_025d, etot_025d, 
"bClasses_025d_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
5e5, 2e4, 5e4, 0.06)

function waveavg_volume_plot(rWtimes, isos, ΔVol_rWavg, ΔConG_rWavg, ΔConS_rWavg, etot_rWavg, savename, apath,
    vmax, cmaxG, cmaxS, emax)
    f = Figure(resolution = (1600, 900), fontsize=26)
    ga = f[1, 1] = GridLayout() 

    # change in volume
    axvdif = Axis(ga[1, 1], ylabel = "Buoyancy [ms⁻²]", 
                title = "ΔVol = Vol - Vol₀")
    axesum = Axis(ga[1, 3], ylabel = "Buoyancy [ms⁻²]", 
                title = "ε")
    axcdifG = Axis(ga[2, 1], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",
                title = "Gauss Δc = ∑c - ∑c₀")
    axcdifS = Axis(ga[2, 3], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",
                title = "Tanh Δc = ∑c - ∑c₀")

    axcdifG.xticks =  2:2:10
    axcdifS.xticks =  2:2:10

    axvdif.yticks =   (-4e-3:1e-3:-2e-3,  ["-4×10⁻³", "-3×10⁻³", "-2×10⁻³"] )
    axcdifG.yticks =  (-4e-3:1e-3:-2e-3,  ["-4×10⁻³", "-3×10⁻³", "-2×10⁻³"] )

    # vol  []     e tot   e avg []]
    # c totG []]  c totS  c avg []

    limits!(axvdif, 1, 10, -5e-3, -1e-3)
    limits!(axcdifS, 1, 10, -5e-3, -1e-3)
    limits!(axcdifG, 1, 10, -5e-3, -1e-3)
    limits!(axesum, 1, 10, -5e-3, -1e-3)

    hidexdecorations!(axvdif)
    hidedecorations!(axesum)
    hideydecorations!(axcdifS)

    hvdif = heatmap!(axvdif, rWtimes, isos, ΔVol_rWavg', colormap = :balance, colorrange = (-vmax, vmax))
    hcGdif = heatmap!(axcdifG, rWtimes, isos, ΔConG_rWavg', colormap = :balance, colorrange = (-cmaxG, cmaxG))
    hcSdif = heatmap!(axcdifS, rWtimes, isos, ΔConS_rWavg', colormap = :balance, colorrange = (-cmaxS, cmaxS))
    hetot = heatmap!(axesum, rWtimes, isos, etot_rWavg', colormap = :thermal, colorrange = (0, emax))

    vtic = (vmax/2)*1e-5
    cticG = (cmaxG/2)*1e-4
    cticS = (cmaxS/2)*1e-4
    etic = (emax/2)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(ga[1,2], hvdif, ticks = (-vmax/2:vmax/2:vmax/2, [@sprintf("-%0.1f×10⁵", vtic), "0", @sprintf("%0.1f×10⁵", vtic)] ),
    size =25)
    cb2 = Colorbar(ga[2,2], hcGdif, ticks = (-cmaxG/2:cmaxG/2:cmaxG/2, [@sprintf("-%0.1f×10⁴", cticG), "0", @sprintf("%0.1f×10⁴", cticG)] ), 
    size =25)
    cb3 = Colorbar(ga[1,4], hetot, ticks = (0:etic:2*etic),
    size =25)
    cb4 = Colorbar(ga[2,4], hcSdif, ticks = (-cmaxS/2:cmaxS/2:cmaxS/2, [@sprintf("-%0.1f×10⁴", cticS), "0", @sprintf("%0.1f×10⁴", cticS)] ), 
    size =25)

    colsize!(ga, 2, Relative(0.01))
    colsize!(ga, 4, Relative(0.01))

    rowgap!(ga, 8)

    save(apath * savename * ".png", f)
end

waveavg_volume_plot(b_timeseries.times./pm.Tσ, isos, ΔVol_rWavg, ΔConG_rWavg, ΔConS_rWavg, etot_rWavg, 
"bClasses_rWavg_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
1e6, 1.5e4, 5.5e4, 0.04)

waveavg_volume_plot(b_timeseries.times./pm.Tσ, isos, ΔVol_rWavg_025d, ΔConG_rWavg_025d, ΔConS_rWavg_025d, etot_rWavg_025d, 
"bClasses_rWavg_025d_" * @sprintf("n%d_", ini_Nisos) * setname, apath,
5e5, 2e4, 5.5e4, 0.04)

#########################
#                   BOUYANCY dM/dB Calculation
#########################

ini_Nisos = 250
Δb = -500*pm.Ñ^2/ini_Nisos

spinlength = size(CGli)[4]

function dMdb_binning(ini_Nisos, nrange, tlength, spinlength, Bi, CGi, CGri, CGli, CSi)
   
    @info "Calculating isopycnal volume..."
    M_btG = zeros(ini_Nisos, tlength)
    M_btS = zeros(ini_Nisos, tlength)
    M_btGr = zeros(ini_Nisos, tlength)
    M_btGl = zeros(ini_Nisos, spinlength)

    prelength = tlength - spinlength

    for i = 1:prelength
        @info "Time $i of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bi = Bi[:,:, :, i];
        Cgi = CGi[:,:,:, i];
        Cgri = CGri[:,:,:, i];
        Csi = CSi[:,:,:, i];

        # for each buoyancy class:
        for n = nrange
            # (1) CCC locations in the density class
            # starting with b < 0 should be almost the whole domain
            boolB = (bi.< Δb*(n-1))
            # dye within the isopycnal layer
            cG_inb = Cgi[boolB]
            cGr_inb = Cgri[boolB]
            cS_inb = Csi[boolB]
            # volume integrated dye concentration:
            M_btG[n,i] = sum(cG_inb)*16
            M_btGr[n,i] = sum(cGr_inb)*16
            M_btS[n,i] = sum(cS_inb)*16
        end
    end

    for i = prelength+1:tlength
        @info "Time $i of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bi = Bi[:,:, :, i];
        Cgi = CGi[:,:,:, i];
        Cgri = CGri[:,:,:, i];
        #Csi = CSi[:,:,:, i];
        Cgli = CGli[:,:,:, i - prelength];

        # for each buoyancy class:
        for n = nrange
            # (1) CCC locations in the density class
            # starting with b < 0 should be almost the whole domain
            boolB = (bi.< Δb*(n-1))
            # dye within the isopycnal layer
            cG_inb = Cgi[boolB]
            cGr_inb = Cgri[boolB]
            cGl_inb = Cgli[boolB]
            #cS_inb = Csi[boolB]
            # volume integrated dye concentration:
            M_btG[n,i] = sum(cG_inb)*16
            M_btGr[n,i] = sum(cGr_inb)*16
            M_btGl[n,i - prelength] = sum(cGl_inb)*16
            #M_btS[n,i] = sum(cS_inb)*16
        end
    end

    return M_btG, M_btGr, M_btGl, M_btS
end

function dMdb_binning(ini_Nisos, nrange, tlength, spinlength, Bi, CGi, CGri, CGli)
   
    @info "Calculating isopycnal volume..."
    M_btG = zeros(ini_Nisos, tlength)
    M_btGr = zeros(ini_Nisos, tlength)
    M_btGl = zeros(ini_Nisos, spinlength)

    prelength = tlength - spinlength

    for i = 1:prelength
        @info "Time $i of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bi = Bi[:,:, :, i];
        Cgi = CGi[:,:,:, i];
        Cgri = CGri[:,:,:, i];

        # for each buoyancy class:
        for n = nrange
            # (1) CCC locations in the density class
            # starting with b < 0 should be almost the whole domain
            boolB = (bi.< Δb*(n-1))
            # dye within the isopycnal layer
            cG_inb = Cgi[boolB]
            cGr_inb = Cgri[boolB]
            # volume integrated dye concentration:
            M_btG[n,i] = sum(cG_inb)*16
            M_btGr[n,i] = sum(cGr_inb)*16
        end
    end

    for i = prelength+1:tlength
        @info "Time $i of $tlength..."
        # at each time step get the data at (x,ycut,z)
        bi = Bi[:,:, :, i];
        Cgi = CGi[:,:,:, i];
        Cgri = CGri[:,:,:, i];
        Cgli = CGli[:,:,:, i - prelength];

        # for each buoyancy class:
        for n = nrange
            # (1) CCC locations in the density class
            # starting with b < 0 should be almost the whole domain
            boolB = (bi.< Δb*(n-1))
            # dye within the isopycnal layer
            cG_inb = Cgi[boolB]
            cGr_inb = Cgri[boolB]
            cGl_inb = Cgli[boolB]
            # volume integrated dye concentration:
            M_btG[n,i] = sum(cG_inb)*16
            M_btGr[n,i] = sum(cGr_inb)*16
            M_btGl[n,i - prelength] = sum(cGl_inb)*16
        end
    end

    return M_btG, M_btGr, M_btGl
end

(M_btG, M_btGr, M_btGl) = dMdb_binning(ini_Nisos, 1:ini_Nisos, tlength, spinlength, Bi, CGi, CGri, CGli)

Δt = b_timeseries.times[2] - b_timeseries.times[1]

function dMdb_Concetration_Calc(M_bt, Δb, Δt)
    # rate of change of concentration in buoyancy class
    ∂M_bt = M_bt[2:end,:] .- M_bt[1:end-1,:]
    ∂M∂b = ∂M_bt ./ Δb

    # change from the initial amount of concentration in this buoyancy class
    Δ∂M∂b = ∂M∂b .-  ∂M∂b[:,1]

    # calculating the rate of change of dMdb
    ∂_∂M∂b = (∂M∂b[:,2:end] .- ∂M∂b[:,1:end-1] )
    ∂_∂M∂b_∂t = ∂_∂M∂b ./ Δt

    return reverse!(M_bt, dims =1), reverse!(∂_∂M∂b_∂t, dims =1), reverse!(∂M∂b, dims=1), reverse!(Δ∂M∂b, dims =1)
end

(M_btG, ∂_∂M∂b_∂tG, ∂M∂bG, Δ∂M∂bG) = dMdb_Concetration_Calc(M_btG,  Δb, Δt)
(M_btGr, ∂_∂M∂b_∂tGr, ∂M∂bGr, Δ∂M∂bGr) = dMdb_Concetration_Calc(M_btGr,  Δb, Δt)
(M_btGl, ∂_∂M∂b_∂tGl, ∂M∂bGl, Δ∂M∂bGl) = dMdb_Concetration_Calc(M_btGl,  Δb, Δt)
(M_btS, ∂_∂M∂b_∂tS, ∂M∂bS, Δ∂M∂bS) = dMdb_Concetration_Calc(M_btS,  Δb, Δt)

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)

function WaveAveraging(wave_info, ini_Nisos, ∂_∂M∂b_∂tS, ∂M∂bS, Δ∂M∂bS)

    # total number of data points
    Wtlength = wave_info.Wl * wave_info.nTσ
    # rolling wave avg
    hWl = floor(Int, wave_info.Wl/2)

    ∂_∂M∂b_∂tS_rW = zeros(ini_Nisos-1, Wtlength-2*hWl)
    ∂M∂bS_rW = zeros(ini_Nisos-1, Wtlength-2*hWl)
    Δ∂M∂bS_rW = zeros(ini_Nisos-1, Wtlength-2*hWl)
    wtims = b_timeseries.times[wave_info.WavePeriods[hWl+1:Wtlength-hWl]]

    for (i,l) in enumerate(hWl+1:Wtlength-hWl-1)
        wdxs = wave_info.WavePeriods[l-hWl:l+hWl]
        ∂_∂M∂b_∂tS_rW[:,i] = mean(∂_∂M∂b_∂tS[:,wdxs], dims=2)[:,1]
        ∂M∂bS_rW[:,i] = mean(∂M∂bS[:,wdxs], dims=2)[:,1]
        Δ∂M∂bS_rW[:,i] = mean(Δ∂M∂bS[:,wdxs], dims=2)[:,1]
    end

    return ∂_∂M∂b_∂tS_rW, ∂M∂bS_rW, Δ∂M∂bS_rW, wtims
end

(∂_∂M∂b_∂tG_rW, ∂M∂bG_rW, Δ∂M∂bG_rW, wtims) = WaveAveraging(wave_info, ini_Nisos, ∂_∂M∂b_∂tG, ∂M∂bG, Δ∂M∂bG)
(∂_∂M∂b_∂tGr_rW, ∂M∂bGr_rW, Δ∂M∂bGr_rW, _) = WaveAveraging(wave_info, ini_Nisos, ∂_∂M∂b_∂tGr, ∂M∂bGr, Δ∂M∂bGr)
(∂_∂M∂b_∂tS_rW, ∂M∂bS_rW, Δ∂M∂bS_rW, _) = WaveAveraging(wave_info, ini_Nisos, ∂_∂M∂b_∂tS, ∂M∂bS, Δ∂M∂bS)

isos = Δb.*((ini_Nisos - 0.5):-1:0.5)
isosdb = Δb.*((ini_Nisos - 1.5):-1:0.5)

function dMdb_SlopeGauss_plot(times, isos, ∂M∂bG, ∂_∂M∂b_∂tG, Δ∂M∂bG, ∂M∂bS, ∂_∂M∂b_∂tS, Δ∂M∂bS,
    ∂M∂bG_bounds, ∂_∂M∂b_∂tG_bounds, Δ∂M∂bG_bounds, ∂M∂bS_bounds, ∂_∂M∂b_∂tS_bounds, Δ∂M∂bS_bounds,
    ∂M∂bG_ticks, ∂_∂M∂b_∂tG_ticks, Δ∂M∂bG_ticks, ∂M∂bS_ticks, ∂_∂M∂b_∂tS_ticks, Δ∂M∂bS_ticks,
    savename, apath)
    f = Figure(resolution = (1000, 1100), fontsize=26)
    ga = f[1, 1] = GridLayout() 
    
    # change in volume
    axMdbpG = Axis(ga[1, 1], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",)
    axMdbpS = Axis(ga[1, 2], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",)
    axdtMpG = Axis(ga[2, 1], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",)
    axdtMpS = Axis(ga[2, 2], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",)
    axdMpG = Axis(ga[3, 1], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",)
    axdMpS = Axis(ga[3, 2], ylabel = "Buoyancy [ms⁻²]", 
                xlabel = "Wave Periods [Tσ]",)
    
    Label(ga[1, 1:2, Top()],"Tracer Distribution ∂M/∂b", valign = :bottom,
    font = :bold, fontsize = 25,
    padding = (0, 0, 5, 0))
    
    Label(ga[2, 1:2, Top()],"Diapycnal Transport ∂(∂M/∂b)/∂t", valign = :bottom,
        font = :bold, fontsize = 25,
        padding = (0, 0, 5, 0))
    
    Label(ga[3, 1:2, Top()],"Distribution Changes ∂M/∂b - ∂M₀/∂b", valign = :bottom,
    font = :bold, fontsize = 25,
    padding = (0, 0, 5, 0))
    axdMpG.xticks =  2:2:10
    axdMpS.xticks =  2:2:10
    
    axMdbpG.yticks =  -6*1e-3:2e-3:0
    axdtMpG.yticks =  -6*1e-3:2e-3:0
    axdMpG.yticks =  -6*1e-3:2e-3:0
    
    limits!(axMdbpG, 1, 10, -5e-3, -1e-3)
    limits!(axMdbpS, 1, 10, -5e-3, -1e-3)
    limits!(axdtMpG, 1, 10, -5e-3, -1e-3)
    limits!(axdtMpS, 1, 10, -5e-3, -1e-3)
    limits!(axdMpG, 1, 10, -5e-3, -1e-3)
    limits!(axdMpS, 1, 10, -5e-3, -1e-3)
    
    hidexdecorations!(axMdbpG)
    hidexdecorations!(axdtMpG)
    hidedecorations!(axMdbpS)
    hidedecorations!(axdtMpS)
    hideydecorations!(axdMpS)
        
    hMG = heatmap!(axMdbpG, times, isos, ∂M∂bG', colormap = :tempo, colorrange = ∂M∂bG_bounds)
    hMS = heatmap!(axMdbpS, times, isos, ∂M∂bS', colormap = :tempo, colorrange = ∂M∂bS_bounds)
    hdtMG = heatmap!(axdtMpG, times, isos, ∂_∂M∂b_∂tG', colormap = :balance, colorrange = ∂_∂M∂b_∂tG_bounds)
    hdtMS = heatmap!(axdtMpS, times, isos, ∂_∂M∂b_∂tS', colormap = :balance, colorrange = ∂_∂M∂b_∂tS_bounds)
    hdMG = heatmap!(axdMpG, times, isos, Δ∂M∂bG', colormap = :balance, colorrange = Δ∂M∂bG_bounds)
    hdMS = heatmap!(axdMpS, times, isos, Δ∂M∂bS', colormap = :balance, colorrange = Δ∂M∂bS_bounds)

    Label(ga[1, 3, Top()],"Gaussian", valign = :bottom,
    font = :bold, fontsize = 25, halign = :right,
    padding = (-5, 0, 5, 0))

    Label(ga[1, 3, Top()],"Hyperbolic\nTangent", valign = :bottom,
    font = :bold, fontsize = 25, halign = :left,
    padding = (5, 0, 5, 0))

    # create colorbars the size of the whole data set
    cb1 = Colorbar(ga[1,3], hMG, ticks = ∂M∂bG_ticks,
    size =25, flipaxis=false)
    cb2 = Colorbar(ga[1,3], hMS, ticks = ∂M∂bS_ticks,
    size =25)

    cb3 = Colorbar(ga[2,3], hdtMG, ticks = ∂_∂M∂b_∂tG_ticks, 
    size =25, flipaxis=false)
    cb4 = Colorbar(ga[2,3], hdtMS, ticks = ∂_∂M∂b_∂tS_ticks, 
    size =25)

    cb5 = Colorbar(ga[3,3], hdMG, ticks = Δ∂M∂bG_ticks,
    size =25,  flipaxis=false)
    cb6 = Colorbar(ga[3,3], hdMS, ticks = Δ∂M∂bS_ticks,
    size =25)

    colsize!(ga, 3, Relative(0.03))

    save(apath * savename * ".png", f, px_per_unit = 2)
end

#########################
#                   FULL BOUYANCY dM/dB PLOT
#########################
# dMdb is in decsending order by class number but backwards buoynacy order

∂M∂bG_bounds = (0, 1.5e8)
∂M∂bS_bounds = (0, 1e9)
∂_∂M∂b_∂tG_bounds = (-7e4, 7e4)
∂_∂M∂b_∂tS_bounds = (-5e5, 5e5)
Δ∂M∂bG = (-2e8, 2e8)
Δ∂M∂bS = (-4e8, 4e8)

# create colorbars the size of the whole data set
∂M∂bG_ticks = (0:5e7:1e8, ["0", "5×10⁷","1×10⁸"] )
∂M∂bS_ticks = (0:5e8:1e9, ["0", "5×10⁸","1×10⁹"] )
∂_∂M∂b_∂tG_ticks = (-3e4:3e4:3e4, ["-3×10⁴", "0", "3×10⁴"] )
∂_∂M∂b_∂tS_ticks = (-3e5:3e5:3e5,  ["-3×10⁵", "0", "3×10⁵"] )
Δ∂M∂bG_ticks = (-1e8:1e8:1e8,  ["-10⁸", "0", "10⁸"] )
Δ∂M∂bS_ticks = (-2e8:2e8:2e8,  ["-2×10⁸", "0", "2×10⁸"] )

dMdb_SlopeGauss_plot(wtims./pm.Tσ, isosdb, ∂M∂bG_rW, ∂_∂M∂b_∂tG_rW, Δ∂M∂bG, ∂M∂bS_rW, ∂_∂M∂b_∂tS_rW, Δ∂M∂bS_rW,
    ∂M∂bG_bounds, ∂_∂M∂b_∂tG_bounds, Δ∂M∂bG_bounds, ∂M∂bS_bounds, ∂_∂M∂b_∂tS_bounds, Δ∂M∂bS_bounds,
    ∂M∂bG_ticks, ∂_∂M∂b_∂tG_ticks, Δ∂M∂bG_ticks, ∂M∂bS_ticks, ∂_∂M∂b_∂tS_ticks, Δ∂M∂bS_ticks,
    "MbClasses_rWavg_" * @sprintf("n%d_", ini_Nisos) * setname, apath)

#########################
#                   CLOSE UP BOUYANCY dM/dB Calculation
#########################

cut_ini_Nisos = 46
Δb = -500*pm.Ñ^2/ini_Nisos

(M_btG_cut, M_btGr_cut, M_btS_cut) = dMdb_Concetration_Calc(cut_ini_Nisos, 120:165, tlength, Bi, CGi, CGri, CSi)

(M_btG_cut, ∂_∂M∂b_∂tG_cut, ∂M∂bG_cut, Δ∂M∂bG_cut) = dMdb_Concetration_Calc(M_btG_cut,  Δb, Δt)
(M_btS_cut, ∂_∂M∂b_∂tS_cut, ∂M∂bS_cut, Δ∂M∂bS_cut) = dMdb_Concetration_Calc(M_btS_cut,  Δb, Δt)

isos_cut = Δb.*((165 - 0.5):-1:(120-0.5))
isosdb_cut = Δb.*((165 - 1.5):-1:(120-0.5))

∂M∂bG_cutbounds = (2e7, 8e7)
∂M∂bS_cutbounds = (2e8, 8e8)
∂_∂M∂b_∂tG_cutbounds = (-2e4, 2e4)
∂_∂M∂b_∂tS_cutbounds = (-2e5, 2e5)
Δ∂M∂bG_cutbounds = (-1.5e8, 1.5e8)
Δ∂M∂bS_cutbounds = (-1.5e8, 1.5e8)

∂M∂bG_cutticks = (3e7:3e7:6e7, ["3×10⁷","6×10⁷"] )
∂M∂bS_cutticks = (3e8:3e8:6e8, ["3×10⁸","6×10⁸"] )
∂_∂M∂b_∂tG_cutticks = (-1e4:1e4:1e4, ["-1×10⁴", "0", "1×10⁴"] )
∂_∂M∂b_∂tS_cutticks = (-1e5:1e5:1e5,  ["-1×10⁵", "0", "1×10⁵"] ) 
Δ∂M∂bG_cutticks = (-1e8:1e8:1e8,  ["-1×10⁸", "0", "1×10⁸"] )
Δ∂M∂bG_cutticks = (-1e8:1e8:1e8,  ["-1×10⁸", "0", "1×10⁸"] )

dMdb_SlopeGauss_plot(b_timeseries.times./pm.Tσ, isosdb_cut, ∂M∂bG_cut, ∂_∂M∂b_∂tG_cut, Δ∂M∂bG, ∂M∂bS_cut, ∂_∂M∂b_∂tS_cut, Δ∂M∂bS_cut,
    ∂M∂bG_cutbounds, ∂_∂M∂b_∂tG_cutbounds, Δ∂M∂bG_cutbounds, ∂M∂bS_cutbounds, ∂_∂M∂b_∂tS_cutbounds, Δ∂M∂bS_cutbounds,
    ∂M∂bG_cutticks, ∂_∂M∂b_∂tG_cutticks, Δ∂M∂bG_cutticks, ∂M∂bS_cutticks, ∂_∂M∂b_∂tS_cutticks, Δ∂M∂bS_cutticks,
    "MbClasses_cut_" * @sprintf("n%d_", cut_ini_Nisos) * setname, apath)

#########################
#                   JUST SLOPE BOUYANCY dM/dB PLOT
#########################

f = Figure(resolution = (700, 900), fontsize=26)
ga = f[1, 1] = GridLayout() 

# change in volume
axMdbpS = Axis(ga[1, 1], ylabel = "Buoyancy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdMpS = Axis(ga[2, 1], ylabel = "Buoyancy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdtMpS = Axis(ga[3, 1], ylabel = "Buoyancy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)

Label(ga[1, 1, Top()],"Tracer Distribution", valign = :bottom,
font = :bold, fontsize = 25,
padding = (0, 0, 5, 0))

Label(ga[2, 1, Top()],"Change in Distribution", valign = :bottom,
font = :bold, fontsize = 25,
padding = (0, 0, 5, 0))

Label(ga[3, 1, Top()],"Diapycnal Transport", valign = :bottom,
    font = :bold, fontsize = 25,
    padding = (0, 0, 5, 0))

axdtMpS.xticks =  4:2:10

axMdbpS.yticks =  -6*1e-3:2e-3:0
axdtMpS.yticks =  -6*1e-3:2e-3:0
axdMpS.yticks =  -6*1e-3:2e-3:0

limits!(axMdbpS, 3, 10, -5e-3, -1e-3)
limits!(axdtMpS, 3, 10, -5e-3, -1e-3)
limits!(axdMpS, 3, 10, -5e-3, -1e-3)

hidexdecorations!(axMdbpS)
hidexdecorations!(axdMpS)

hMS = heatmap!(axMdbpS, wtims./pm.Tσ, isosdb, ∂M∂bS_rW', colormap = :tempo, colorrange = (0, 6e8))
hdtMS = heatmap!(axdtMpS, wtims./pm.Tσ, isosdb, ∂_∂M∂b_∂tS_rW', colormap = :balance, colorrange = (-4e4, 4e4))
hdMS = heatmap!(axdMpS, wtims./pm.Tσ, isosdb, Δ∂M∂bS_rW', colormap = :balance, colorrange = (-2.5e8, 2.5e8))

# create colorbars the size of the whole data set
cb2 = Colorbar(ga[1,2], hMS, ticks = (0:3e8:6e8, ["0", "3×10⁸","6×10⁸"] ),
size =25)

cb6 = Colorbar(ga[2,2], hdMS, ticks = (-1e8:1e8:1e8, ["-1×10⁸", "0", "1×10⁸"] ),
 size =25)

cb4 = Colorbar(ga[3,2], hdtMS, ticks = (-3e4:3e4:3e4,  ["-3×10⁴", "0", "3×10⁴"] ), 
size =25)

colsize!(ga, 2, Relative(0.03))

savename = "MbSlopeClasses_" * @sprintf("n%d_", ini_Nisos) * setname

save(apath * savename * ".png", f)

#########################
#                   TRACER WEIGHTING Calculation
#########################

function calculate_tracer_weight(CGi, Bi)
    
    @info "Calculating Denominator of tracer wtd average..."
    # ∫ c dV
    csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

    @info "Calculating first moment..."
    # ∫ cb dV
    cb = CGi .* Bi;
    cbsum = sum(cb, dims=(1,2,3))[1,1,1,:];

    # b̄ =  ∫ cb dV / ∫ c dV
    c_weighted_bavg = cbsum ./ csum;

    return c_weighted_bavg, csum
end

function calculate_cb_velocities(CGi, CGli, CGri, Bi, Δt)
    @info "Calculating Tracer Velocity..."
    ∂Cg_∂t = (CGi[:,:,:,2:end] .- CGi[:,:,:,1:end-1])./ Δt
    ∂Cgl_∂t = (CGli[:,:,:,2:end] .- CGli[:,:,:,1:end-1])./ Δt
    ∂Cgr_∂t = (CGri[:,:,:,2:end] .- CGri[:,:,:,1:end-1])./ Δt
    @info "Calculating Buoynacy Velocity..."
    ∂b_∂t = (Bi[:,:,:,2:end] .- Bi[:,:,:,1:end-1])./ Δt

    return ∂Cg_∂t, ∂Cgl_∂t, ∂Cgr_∂t, ∂b_∂t
end

function calculate_tracer_weight_velocity_terms(b̄_Cg, csum, ∂b_∂t, ∂C_∂t, CGi, Bi, Δt)
    
    @info "Calculating Tracer Wtd Velocity"
    ∂b̄_∂t =  (b̄_Cg[2:end] .- b̄_Cg[1:end-1])./ Δt

    @info "Calculating total two terms..."
    b∂C_∂t = ∂C_∂t .*Bi[:,:,:,2:end]
    C∂b_∂t = ∂b_∂t .*CGi[:,:,:,2:end]

    @info "Integrating..."
    ∫b∂C_∂t = sum(b∂C_∂t, dims = (1,2,3))[1,1,1,:]
    ∫C∂b_∂t = sum(C∂b_∂t, dims = (1,2,3))[1,1,1,:]

    ∫b∂C_∂t_csum = ∫b∂C_∂t ./ csum[2:end]
    ∫C∂b_∂t_csum = ∫C∂b_∂t ./ csum[2:end]

    return ∫b∂C_∂t_csum, ∫C∂b_∂t_csum, ∂b̄_∂t
end

(b̄_Cg, Cg_sum) = calculate_tracer_weight(CGi, Bi)
(b̄_Cgr, Cgr_sum) = calculate_tracer_weight(CGri, Bi)
(b̄_Cgl, Cgl_sum) = calculate_tracer_weight(CGli, Bi[:,:,:,(tlength-spinlength+1):end])

(∂Cg_∂t, ∂Cgl_∂t, ∂Cgr_∂t, ∂b_∂t) = calculate_cb_velocities(CGi, CGli, CGri, Bi, Δt)

(∫b∂Cg_∂t_Cg_sum, ∫Cg∂b_∂t_Cgsum, ∂b̄_∂t_Cg) = calculate_tracer_weight_velocity_terms(b̄_Cg, Cg_sum, ∂b_∂t, ∂Cg_∂t, CGi, Bi, 600.0)
(∫b∂Cgr_∂t_Cgr_sum, ∫Cgr∂b_∂t_Cgrsum, ∂b̄_∂t_Cgr) = calculate_tracer_weight_velocity_terms(b̄_Cgr, Cgr_sum, ∂b_∂t, ∂Cgr_∂t, CGri, Bi, 600.0)
(∫b∂Cgl_∂t_Cgl_sum, ∫Cgl∂b_∂t_Cglsum, ∂b̄_∂t_Cgl) = calculate_tracer_weight_velocity_terms(b̄_Cgl, Cgl_sum, ∂b_∂t[:,:,:,(tlength-spinlength+1):end], ∂Cgl_∂t, CGli, Bi[:,:,:,(tlength-spinlength+1):end], 600.0)

######################
#                  GAUSSIAN, RAISED, and LATE RELEASE WTD VELOCITY TERMS PLOT
#######################

f = Figure(resolution = (1200, 1000), fontsize=26)
ga = f[1, 1] = GridLayout() 

ax1 = Axis(ga[1, 1], ylabel = "∂b̄/∂t [ms⁻³]")
ax2 = Axis(ga[2, 1], ylabel = "∂b̄/∂t [ms⁻³]")
ax3 = Axis(ga[3, 1], ylabel = "∂b̄/∂t [ms⁻³]", xlabel = "Wave Periods [Tσ]")
ax3.xticks =  2:2:10

ax1.yticks =  (-5e-7:5e-7:5e-7, ["-5×10⁻⁷", "0", "5×10⁻⁷"])
ax2.yticks =   (-5e-7:5e-7:5e-7, ["-5×10⁻⁷", "0", "5×10⁻⁷"])
ax3.yticks =  (-5e-7:5e-7:5e-7, ["-5×10⁻⁷", "0", "5×10⁻⁷"])

limits!(ax1, 1, 10, -1e-6,   1e-6)
limits!(ax2, 1, 10, -1e-6,   1e-6)
limits!(ax3, 1, 10, -1e-6,   1e-6)

hidexdecorations!(ax1)
hidexdecorations!(ax2)

dt_times = b_timeseries.times[2:end]./pm.Tσ
dt_times_spin = b_timeseries.times[tlength-spinlength+2:end]./pm.Tσ

dbb = lines!(ax1, dt_times, ∂b̄_∂t_Cg, color = :gray32, linewidth = 5)
cdb = lines!(ax1, dt_times, ∫Cg∂b_∂t_Cgsum, color = :firebrick2, linewidth = 5)
bdc = lines!(ax1, dt_times, ∫b∂Cg_∂t_Cg_sum, color = :dodgerblue2, linewidth = 5)

lines!(ax2, dt_times_spin, ∂b̄_∂t_Cgl, color = :gray32, linewidth = 5)
lines!(ax2, dt_times_spin, ∫Cgl∂b_∂t_Cglsum, color = :firebrick2, linewidth = 5)
lines!(ax2, dt_times_spin, ∫b∂Cgl_∂t_Cgl_sum, color = :dodgerblue2, linewidth = 5)

lines!(ax3, dt_times, ∂b̄_∂t_Cgr, color = :gray32, linewidth = 5)
lines!(ax3, dt_times, ∫Cgr∂b_∂t_Cgrsum, color = :firebrick2, linewidth = 5)
lines!(ax3, dt_times, ∫b∂Cgr_∂t_Cgr_sum, color = :dodgerblue2, linewidth = 5)

Legend(ga[1, 1, Top()], [dbb, cdb, bdc], ["∂b̄/∂t", "∫c ∂b/∂t dV", "∫b ∂c/∂t dV"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 30), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :bottom, orientation = :horizontal)

Label(ga[1, 1, Top()],"On Slope, Release t = 0 Tσ", valign = :top,
                font = :bold, fontsize = 25,
                padding = (0, 0, 10, 0))
                
Label(ga[2, 1, Top()],"On Slope, Release t = 4 Tσ", valign = :top,
                font = :bold, fontsize = 25,
                padding = (0, 0, 10, 0))
    
Label(ga[3, 1, Top()], "Off Slope, Release t = 0 Tσ", valign = :top,
                font = :bold, fontsize = 25,
                padding = (0, 0, 10, 0))

savename = "Gauss_cwtd_dbdt_" * setname

save(apath * savename * ".png", f)

#########################
#                   GAUSSIAN, RAISED, and LATE RELEASE BOUYANCY dM/dB PLOT With TRACER WTD
#########################
isosdb = Δb.*((ini_Nisos - 1.5):-1:0.5)

braised₀ = b̄_Cgr[1]
bonslope₀ = b̄_Cg[1]
blate₀ = b̄_Cgl[1]

pre_∂M∂bGl = zeros(ini_Nisos-1, tlength - spinlength)
full_∂M∂bGl = hcat(pre_∂M∂bGl, ∂M∂bGl)
pre_Δ∂M∂bGl = zeros(ini_Nisos-1, tlength - spinlength)
full_Δ∂M∂bGl = hcat(pre_Δ∂M∂bGl, Δ∂M∂bGl)

f = Figure(resolution = (2000, 1400), fontsize=26)
ga = f[1, 1] = GridLayout() 

# change in volume
axMdbpG = Axis(ga[1, 1], ylabel = "Buoyancy [ms⁻²]")
axMdbpGl = Axis(ga[1, 2])
axMdbpGr = Axis(ga[1, 3])
axdMpG = Axis(ga[2, 1], ylabel = "Buoyancy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdMpGl = Axis(ga[2, 2], xlabel = "Wave Periods [Tσ]",)
axdMpGr = Axis(ga[2, 3], xlabel = "Wave Periods [Tσ]",)

Label(ga[1, 1:3, Top()],"Tracer Distribution ∂M/∂b", valign = :top,
font = :bold, fontsize = 25,
padding = (0, 0, 10, 0))

Label(ga[2, 1:3, Top()],"Distribution Changes ∂M/∂b - ∂M₀/∂b", valign = :top,
font = :bold, fontsize = 25,
padding = (0, 0, 10, 0))

axdMpG.xticks =  2:2:10
axdMpGr.xticks =  2:2:10
axdMpGl.xticks =  2:2:10

axMdbpG.yticks =  ((bonslope₀ - 5e-4):5e-4:(bonslope₀ + 6e-4), ["b₀ - 5×10⁻⁴", "b₀", "b₀ + 5×10⁻⁴"])
axdMpG.yticks =  ((bonslope₀ - 5e-4):5e-4:(bonslope₀ + 6e-4), ["b₀ - 5×10⁻⁴", "b₀", "b₀ + 5×10⁻⁴"])

limits!(axMdbpG, 1, 10,bonslope₀ - 7e-4,bonslope₀ + 7e-4)
limits!(axMdbpGr, 1, 10, braised₀ - 7e-4 , braised₀ + 7e-4)
limits!(axMdbpGl, 1, 10, blate₀ - 7e-4 , blate₀ + 7e-4)
limits!(axdMpG, 1, 10, bonslope₀ - 7e-4, bonslope₀ + 7e-4)
limits!(axdMpGr, 1, 10, braised₀ - 7e-4 , braised₀ + 7e-4)
limits!(axdMpGl, 1, 10, blate₀ - 7e-4 , blate₀ + 7e-4)

hidexdecorations!(axMdbpG)
hidedecorations!(axMdbpGr)
hidedecorations!(axMdbpGl)
hideydecorations!(axdMpGr)
hideydecorations!(axdMpGl)

hMG = heatmap!(axMdbpG, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bG', colormap = :tempo, colorrange = (0, 2e8))
text!(axMdbpG, Point.(2, bonslope₀ - 0.6*1e-3), text = "On Slope: (-250 m, 703 m)", align = (:left, :center), 
        color = :black, fontsize = 20, font = :bold)
lines!(axMdbpG, b_timeseries.times./pm.Tσ, b̄_Cg, color = :dodgerblue2, linewidth = 5)

hMGl = heatmap!(axMdbpGl, b_timeseries.times./pm.Tσ, isosdb, full_∂M∂bGl', colormap = :tempo, colorrange = (0, 2e8))
text!(axMdbpGl, Point.(2, blate₀ - 0.6*1e-3), text = "On Slope: (-250 m, 703 m)", align = (:left, :center), 
        color = :black, fontsize = 20, font = :bold)
lines!(axMdbpGl, b_timeseries.times[tlength-spinlength+1:end]./pm.Tσ, b̄_Cgl, color = :dodgerblue2, linewidth = 5)

hMGr = heatmap!(axMdbpGr, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bGr', colormap = :tempo, colorrange = (0, 2e8))
text!(axMdbpGr, Point.(2, braised₀ - 0.6*1e-3), text = "Up 100m: (-136 m , 743 m)", align = (:left, :center), 
        color = :black, fontsize = 20, font = :bold)
lines!(axMdbpGr, b_timeseries.times./pm.Tσ, b̄_Cgr, color = :dodgerblue2, linewidth = 5)

hdMG = heatmap!(axdMpG, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bG', colormap = :balance, colorrange = (-2e8, 2e8))
hdMGl = heatmap!(axdMpGl, b_timeseries.times./pm.Tσ, isosdb, full_Δ∂M∂bGl', colormap = :balance, colorrange = (-2e8, 2e8))
hdMGr = heatmap!(axdMpGr, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bGr', colormap = :balance, colorrange = (-2e8, 2e8))

# create colorbars the size of the whole data set

cb1 = Colorbar(ga[1,4], hMG, ticks = (0:1e8:2e8, ["0", "1×10⁸", "2×10⁸"] ),
size =25)
cb2 = Colorbar(ga[2,4], hdMG, ticks = (-1e8:1e8:1e8,  ["-1×10⁸", "0", "1×10⁸"] ),
 size =25)

colsize!(ga, 4, Relative(0.05))

savename = "MbGaussClasses_cwtd_" * @sprintf("n%d_", ini_Nisos) * setname

save(apath * savename * ".png", f)