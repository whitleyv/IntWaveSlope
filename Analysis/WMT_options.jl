using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
apath = path_name * "Analysis/"

include("../parameters.jl")

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

filepath1 = path_name * name1_prefix * ".jld2"
filepath2 = path_name * name2_prefix * ".jld2"

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

CSi = interior(Cs_timeseries)[:, 1:ylength, 1:zlength,1:2:end];
CGi = interior(Cg_timeseries)[:, 1:ylength, 1:zlength,1:2:end];
CGri = interior(Cgr_timeseries)[:, 1:ylength, 1:zlength,:];
ei = interior(e_timeseries)[:,1:ylength,1:zlength,:];
Bi = interior(b_timeseries)[:, 1:ylength, 1:zlength,:];

@info "Computing Volumes..."
ini_Nisos = 35 # number of bins
Δb = -500*pm.Ñ^2/ini_Nisos # "width" of bins
#Nisos = ini_Nisos + 3 

# buoynacy sum
bsum = zeros(ini_Nisos, tlength)
# total concentration and dissip
contotG = zeros(ini_Nisos, tlength)
contotS = zeros(ini_Nisos, tlength)
etot = zeros(ini_Nisos, tlength)
# average concentration and dissip
cavg = zeros(ini_Nisos, tlength)
eavg = zeros(ini_Nisos, tlength)

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
        cavg[n,j] = mean(cS_inb)
        # average dissip in n
        eavg[n,j] = mean(e_inb)
    end
end

# change in volume from initial
ΔVol = (bsum .- bsum[:,1]).* 16
# change in concentration from initial
ΔConG = (contotG .- contotG[:,1] ) .* 16
ΔConS = (contotS .- contotS[:,1] ) .* 16

@info "Calculating Wave Indices..."
include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
# rolling wave average over two waves
WL = wave_info.Wl
nTσ = wave_info.nTσ
Wtlength = length(WL+1:WL*nTσ-WL)

# volume change since initial
ΔVol_rWavg = zeros(ini_Nisos, Wtlength)
# change in total concentration since initial
ΔConG_rWavg = zeros(ini_Nisos, Wtlength)
ΔConS_rWavg = zeros(ini_Nisos, Wtlength)
# average dye in n
Conavg_rWavg = zeros(ini_Nisos, Wtlength)
# average dissip in n
eavg_rWavg = zeros(ini_Nisos, Wtlength)
# total issipation
etot_rWavg = zeros(ini_Nisos, Wtlength)

for k = (WL+1):(WL*nTσ - WL)
    WaveIndices = wave_info.WavePeriods[k-WL:k+WL]
    # change in total volume since initial
    ΔVol_rWavg[:,k-WL] = mean(ΔVol[:,WaveIndices], dims=2)
    # change in total concentration since initial
    ΔConG_rWavg[:,k-WL] = mean(ΔConG[:,WaveIndices], dims=2)
    ΔConS_rWavg[:,k-WL] = mean(ΔConS[:,WaveIndices], dims=2)
    # average dye in n
    Conavg_rWavg[:,k-WL] = mean(cavg[:,WaveIndices], dims=2)
    # average dissip in n
    eavg_rWavg[:,k-WL] = mean(eavg[:,WaveIndices], dims=2)
    # total issipation
    etot_rWavg[:,k-WL] = mean(etot[:,WaveIndices], dims=2)

end

# n = 1 is really at the top of the domain where b = 0
# switch the ordering so n=1 is the deepest?
reverse!(ΔVol_rWavg, dims=1)
reverse!(ΔConG_rWavg, dims=1)
reverse!(ΔConS_rWavg, dims=1)
reverse!(Conavg_rWavg, dims=1)
reverse!(eavg_rWavg, dims=1)
reverse!(etot_rWavg, dims=1)

# wave period times for rolling average
rWtimes = b_timeseries.times[wave_info.WavePeriods[WL+1:WL*nTσ - WL]]/pm.Tσ

# buoynacy classes
# first index is the buoynacy found at lowest depth
# last index is buoyancy found at top ie. (Δb/2)
isos = Δb.*((ini_Nisos - .5):-1:0.5)

eavg_rWavg_log = log10.(clamp.(eavg_rWavg, 1e-14, 1e-1))
Conavg_rWavg_log = log10.(clamp.(Conavg_rWavg, 1e-14, 1e-1))
etot_rWavg_log = log10.(clamp.(etot_rWavg, 1e-14, 1e-1))

f = Figure(resolution = (1500, 900), fontsize=26)
ga = f[1, 1] = GridLayout() 

# change in volume
axvdif = Axis(ga[1, 1], ylabel = "nΔb [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",
            title = "ΔVol = Vol - Vol₀")
axesum = Axis(ga[1, 3], ylabel = "nΔb [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",
            title = "∑ε")
axeavg = Axis(ga[1, 4], ylabel = "nΔb [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",
            title = "⟨ε⟩")

axcdifG = Axis(ga[2, 1], ylabel = "nΔb [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",
            title = "Gauss Δc = ∑c - ∑c₀")
axcdifS = Axis(ga[2, 3], ylabel = "nΔb [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",
            title = "Tanh Δc = ∑c - ∑c₀")
axcavg = Axis(ga[2, 4], ylabel = "nΔb [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",
            title = "Tanh ⟨c⟩")

axcdifG.xticks =  1:2:10
axcdifS.xticks =  1:2:10
axcavg.xticks =  1:2:10

axvdif.yticks =   (-5e-3:2e-3:-1e-3,  ["5×10⁻³", "3×10⁻³", "1×10⁻³"] )
axcdifG.yticks =  (-5e-3:2e-3:-1e-3,  ["5×10⁻³", "3×10⁻³", "1×10⁻³"] )

# vol  []     e tot   e avg []]
# c totG []]  c totS  c avg []

limits!(axvdif, 1, 10, isos[1], isos[end])
limits!(axcdifS, 1, 10, isos[1], isos[end])
limits!(axcdifG, 1, 10, isos[1], isos[end])
limits!(axcavg, 1, 10, isos[1], isos[end])
limits!(axesum, 1, 10, isos[1], isos[end])
limits!(axeavg, 1, 10, isos[1], isos[end])

hidexdecorations!(axvdif)
hidedecorations!(axesum)
hidedecorations!(axeavg)
hideydecorations!(axcdifS)
hideydecorations!(axcavg)

hvdif = heatmap!(axvdif, rWtimes, isos, ΔVol_rWavg', colormap = :balance, colorrange = (-1.2e6, 1.2e6))
hcGdif = heatmap!(axcdifG, rWtimes, isos, ΔConG_rWavg', colormap = :balance, colorrange = (-1.7e4, 1.7e4))
hcSdif = heatmap!(axcdifS, rWtimes, isos, ΔConS_rWavg', colormap = :balance, colorrange = (-5.5e4, 5.5e4))
hcavg = heatmap!(axcavg, rWtimes, isos, Conavg_rWavg', colormap = :thermal, colorrange = (0, 0.03))
hetot = heatmap!(axesum, rWtimes, isos, etot_rWavg', colormap = :thermal, colorrange = (0, 0.08))
heavg = heatmap!(axeavg, rWtimes, isos, eavg_rWavg', colormap = :thermal, colorrange = (0, 3e-7))

# create colorbars the size of the whole data set
cb1 = Colorbar(ga[1,2], hvdif, ticks = (-1e6:1e6:1e6, ["-10⁶", "0", "10⁶"] ),
size =25, flipaxis=false, label = "ΔVol")

cb5 = Colorbar(ga[1,5], hetot, ticks = (0:0.02:0.06),
 size =25, label = "∑ε", flipaxis=false)
cb6 = Colorbar(ga[1,5], heavg, ticks = (0:1e-7:3e-7,  ["0", "1×10⁻⁷", "2×10⁻⁷","3×10⁻⁷"]),
 size =25, label = "⟨ε⟩")

cb2 = Colorbar(ga[2,2], hcGdif, ticks = (-1e4:1e4:1e4, ["-10⁴", "0", "10⁴"] ), 
size =25, label = "Gauss Δc", flipaxis=false)
cb3 = Colorbar(ga[2,2], hcSdif, ticks = (-5e4:5e4:5e4, ["-5×10⁴", "0", "5×10⁴"] ), 
size =25, label = "Δc × 10⁻⁵")

cb4 = Colorbar(ga[2,5], hcavg, ticks = (0:.01:.03 ), 
size =25, label = "⟨c⟩")

colsize!(ga, 2, Relative(0.01))
colsize!(ga, 5, Relative(0.01))

bigtitle = @sprintf("Buoyancy Space Analysis, U₀=%0.2f, N=%0.2f ×10⁻³, δ=%0.1f", pm.U₀, pm.Ñ*1e3, pm.U₀/pm.Ñ)

Label(f[1, 1, Top()],bigtitle, valign = :bottom,
    font = :bold, fontsize = 30,
    padding = (0, 0, 50, 0))

rowgap!(ga, 10)

savename = "bClasses_rWavg_" * @sprintf("n%d_", ini_Nisos) * setname

save(apath * savename * ".png", f)


#########################
#                   BOUYANCY dM/dB Calculation
#########################

ini_Nisos = 250
Δb = -500*pm.Ñ^2/ini_Nisos

@info "Calculating isopycnal volume..."
M_btG = zeros(ini_Nisos, tlength)
M_btS = zeros(ini_Nisos, tlength)
M_btGr = zeros(ini_Nisos, tlength)
for i = 1:tlength
    @info "Time $i of $tlength..."
    # at each time step get the data at (x,ycut,z)
    bi = Bi[:,:, :, i];
    #Cgi = CGi[:,:,:, i];
    #Cgri = CGri[:,:,:, i];
    Csi = CSi[:,:,:, i];

    # for each buoyancy class:
    for n = 1:ini_Nisos
        # (1) CCC locations in the density class
        # starting with b < 0 should be almost the whole domain
        boolB = (bi.< Δb*(n-1))
        # dye within the isopycnal layer
        #cG_inb = Cgi[boolB]
        #cGr_inb = Cgri[boolB]
        cS_inb = Csi[boolB]
        # volume integrated dye concentration:
        #M_btG[n,i] = sum(cG_inb)*16
        #M_btGr[n,i] = sum(cGr_inb)*16
        M_btS[n,i] = sum(cS_inb)*16
    end
end

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
(M_btS, ∂_∂M∂b_∂tS, ∂M∂bS, Δ∂M∂bS) = dMdb_Concetration_Calc(M_btS,  Δb, Δt)

include("WaveValues.jl")
wave_info=get_wave_indices(b_timeseries, pm, tlength)
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

isos = Δb.*((ini_Nisos - 0.5):-1:0.5)
isosdb = Δb.*((ini_Nisos - 1.5):-1:0.5)

#########################
#                   CLOSE UP BOUYANCY dM/dB Calculation
#########################

cut_ini_Nisos = 46
Δb = -500*pm.Ñ^2/ini_Nisos

@info "Calculating isopycnal volume..."
M_btG_cut = zeros(cut_ini_Nisos, tlength)
M_btS_cut = zeros(cut_ini_Nisos, tlength)

for i = 1:tlength
    @info "Time $i of $tlength..."
    # at each time step get the data at (x,ycut,z)
    bi = Bi[:,:, :, i];
    Cgi = CGi[:,:,:, i];
    Csi = CSi[:,:,:, i];

    # for each buoyancy class:
    for n = 120:165
        # (1) CCC locations in the density class
        # starting with b < 0 should be almost the whole domain
        boolB = (bi.< Δb*(n-1))
        # dye within the isopycnal layer
        cG_inb = Cgi[boolB]
        cS_inb = Csi[boolB]
        # volume integrated dye concentration:
        M_btG_cut[n-119,i] = sum(cG_inb)*16
        M_btS_cut[n-119,i] = sum(cS_inb)*16
    end
end

(M_btG_cut, ∂_∂M∂b_∂tG_cut, ∂M∂bG_cut, Δ∂M∂bG_cut) = dMdb_Concetration_Calc(M_btG_cut,  Δb, Δt)
(M_btS_cut, ∂_∂M∂b_∂tS_cut, ∂M∂bS_cut, Δ∂M∂bS_cut) = dMdb_Concetration_Calc(M_btS_cut,  Δb, Δt)

isos_cut = Δb.*((165 - 0.5):-1:(120-0.5))
isosdb_cut = Δb.*((165 - 1.5):-1:(120-0.5))

f = Figure(resolution = (1100, 900), fontsize=26)
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
axdMpG.xticks =  4:2:10
axdMpS.xticks =  4:2:10

axMdbpG.yticks =  (-3.75*1e-3:0.5*1e-3:-3e-3, ["-3.75×10⁻³", "-3.25×10⁻³"])
axdtMpG.yticks =  (-3.75*1e-3:0.5*1e-3:-3e-3, ["-3.75×10⁻³", "-3.25×10⁻³"])
axdMpG.yticks =  (-3.75*1e-3:0.5*1e-3:-3e-3, ["-3.75×10⁻³", "-3.25×10⁻³"])

limits!(axMdbpG, 3, 10, -4e-3, -3e-3)
limits!(axMdbpS, 3, 10, -4e-3, -3e-3)
limits!(axdtMpG, 3, 10, -4e-3, -3e-3)
limits!(axdtMpS, 3, 10, -4e-3, -3e-3)
limits!(axdMpG, 3, 10, -4e-3, -3e-3)
limits!(axdMpS, 3, 10, -4e-3, -3e-3)

hidexdecorations!(axMdbpG)
hidexdecorations!(axdtMpG)
hidedecorations!(axMdbpS)
hidedecorations!(axdtMpS)
hideydecorations!(axdMpS)

hMG = heatmap!(axMdbpG, b_timeseries.times./pm.Tσ, isosdb_cut, ∂M∂bG_cut', colormap = :tempo, colorrange = (2e7, 8e7))
hMS = heatmap!(axMdbpS, b_timeseries.times./pm.Tσ, isosdb_cut, ∂M∂bS_cut', colormap = :tempo, colorrange = (2e8, 8e8))
hdtMG = heatmap!(axdtMpG , b_timeseries.times./pm.Tσ, isosdb_cut, ∂_∂M∂b_∂tG_cut', colormap = :balance, colorrange = (-2e4, 2e4))
hdtMS = heatmap!(axdtMpS, b_timeseries.times./pm.Tσ, isosdb_cut, ∂_∂M∂b_∂tS_cut', colormap = :balance, colorrange = (-2e5, 2e5))
hdMG = heatmap!(axdMpG, b_timeseries.times./pm.Tσ, isosdb_cut, Δ∂M∂bG_cut', colormap = :balance, colorrange = (-1.5e8, 1.5e8))
hdMS = heatmap!(axdMpS, b_timeseries.times./pm.Tσ, isosdb_cut, Δ∂M∂bS_cut', colormap = :balance, colorrange = (-1.5e8, 1.5e8))

Label(ga[1, 3, Top()],"Gaussian", valign = :bottom,
font = :bold, fontsize = 25, halign = :right,
padding = (-5, 0, 5, 0))

Label(ga[1, 3, Top()],"Hyperbolic\nTangent", valign = :bottom,
font = :bold, fontsize = 25, halign = :left,
padding = (5, 0, 5, 0))

# create colorbars the size of the whole data set
cb1 = Colorbar(ga[1,3], hMG, ticks = (3e7:3e7:6e7, ["3×10⁷","6×10⁷"] ),
size =25, flipaxis=false)
cb2 = Colorbar(ga[1,3], hMS, ticks = (3e8:3e8:6e8, ["3×10⁸","6×10⁸"] ),
size =25)

cb3 = Colorbar(ga[2,3], hdtMG, ticks = (-1e4:1e4:1e4, ["-1×10⁴", "0", "1×10⁴"] ), 
size =25, flipaxis=false)
cb4 = Colorbar(ga[2,3], hdtMS, ticks = (-1e5:1e5:1e5,  ["-1×10⁵", "0", "1×10⁵"] ), 
size =25)

cb5 = Colorbar(ga[3,3], hdMG, ticks = (-1e8:1e8:1e8,  ["-1×10⁸", "0", "1×10⁸"] ),
 size =25,  flipaxis=false)
cb6 = Colorbar(ga[3,3], hdMS, ticks = (-1e8:1e8:1e8,  ["-1×10⁸", "0", "1×10⁸"] ),
 size =25)

colsize!(ga, 3, Relative(0.03))

savename = "MbClasses_cut_" * @sprintf("n%d_", cut_ini_Nisos) * setname

save(apath * savename * ".png", f, px_per_unit = 2)

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
#                   FULL BOUYANCY dM/dB PLOT
#########################
# dMdb is in decsending order by class number but backwards buoynacy order

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

hMG = heatmap!(axMdbpG, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bG', colormap = :tempo, colorrange = (0, 1.5e8))
hMS = heatmap!(axMdbpS, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bS', colormap = :tempo, colorrange = (0, 1e9))
hdtMG = heatmap!(axdtMpG , b_timeseries.times./pm.Tσ, isosdb, ∂_∂M∂b_∂tG', colormap = :balance, colorrange = (-7e4, 7e4))
hdtMS = heatmap!(axdtMpS, b_timeseries.times./pm.Tσ, isosdb, ∂_∂M∂b_∂tS', colormap = :balance, colorrange = (-5e5, 5e5))
hdMG = heatmap!(axdMpG, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bG', colormap = :balance, colorrange = (-2e8, 2e8))
hdMS = heatmap!(axdMpS, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bS', colormap = :balance, colorrange = (-4e8, 4e8))


# create colorbars the size of the whole data set
cb1 = Colorbar(ga[1,3], hMG, ticks = (0:5e7:1e8, ["0", "5×10⁷","1×10⁸"] ),
size =25, flipaxis=false, label = "Gauss")
cb2 = Colorbar(ga[1,3], hMS, ticks = (0:5e8:1e9, ["0", "5×10⁸","1×10⁹"] ),
size =25, label = "Tanh")

cb3 = Colorbar(ga[2,3], hdtMG, ticks = (-3e4:3e4:3e4, ["-3×10⁴", "0", "3×10⁴"] ), 
size =25, label = "Gauss", flipaxis=false)
cb4 = Colorbar(ga[2,3], hdtMS, ticks = (-3e5:3e5:3e5,  ["-3×10⁵", "0", "3×10⁵"] ), 
size =25, label = "Tanh")

cb5 = Colorbar(ga[3,3], hdMG, ticks = (-1e8:1e8:1e8,  ["-10⁸", "0", "10⁸"] ),
 size =25, label = "Gauss", flipaxis=false)
cb6 = Colorbar(ga[3,3], hdMS, ticks = (-2e8:2e8:2e8,  ["-2×10⁸", "0", "2×10⁸"] ),
 size =25, label = "Tanh")

colsize!(ga, 3, Relative(0.05))

bigtitle = @sprintf("Tracers in Buoyancy Space, U₀=%0.2f, N=%0.2f ×10⁻³, δ=%0.1f", pm.U₀, pm.Ñ*1e3, pm.U₀/pm.Ñ)

Label(f[1, 1, Top()],bigtitle, valign = :bottom,
    font = :bold, fontsize = 30,
    padding = (0, 0, 50, 0))

savename = "MbClasses_rWavg_" * @sprintf("n%d_", ini_Nisos) * setname

save(apath * savename * ".png", f)

#########################
#                   GAUSSIAN AND RAISED BOUYANCY dM/dB PLOT
#########################

f = Figure(resolution = (1000, 1100), fontsize=26)
ga = f[1, 1] = GridLayout() 

# change in volume
axMdbpG = Axis(ga[1, 1], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axMdbpS = Axis(ga[1, 2], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdtMpG = Axis(ga[2, 1], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdtMpS = Axis(ga[2, 2], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdMpG = Axis(ga[3, 1], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdMpS = Axis(ga[3, 2], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)

Label(ga[1, 1:2, Top()],"Tracer Distribution ∂M/∂b", valign = :top,
font = :bold, fontsize = 25,
padding = (0, 0, 15, 0))

Label(ga[2, 1:2, Top()],"Diapycnal Transport ∂(∂M/∂b)/∂t", valign = :top,
    font = :bold, fontsize = 25,
    padding = (0, 0, 5, 0))

Label(ga[3, 1:2, Top()],"Distribution Changes ∂M/∂b - ∂M₀/∂b", valign = :top,
font = :bold, fontsize = 25,
padding = (0, 0, 5, 0))

#Label(ga[1, 1, Top()],"Gaussian", valign = :bottom,
#font = :bold, fontsize = 25,
#padding = (0, 0, 5, 0))

#Label(ga[1, 2, Top()],"Raised Gaussian", valign = :bottom,
#font = :bold, fontsize = 25,
#padding = (0, 0, 0, 0))

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

hMG = heatmap!(axMdbpG, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bG', colormap = :tempo, colorrange = (0, 1.5e8))
text!(axMdbpG, Point.(6, -0.0045), text = "(-250 m, 703 m)", align = (:left, :center), 
        color = :black, fontsize = 20)

hMS = heatmap!(axMdbpS, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bGr', colormap = :tempo, colorrange = (0, 1.5e8))
text!(axMdbpS, Point.(6, -0.0045), text = "(-136 m , 743 m)", align = (:left, :center), 
        color = :black, fontsize = 20)

hdtMG = heatmap!(axdtMpG , b_timeseries.times./pm.Tσ, isosdb, ∂_∂M∂b_∂tG', colormap = :balance, colorrange = (-7e4, 7e4))
hdtMS = heatmap!(axdtMpS, b_timeseries.times./pm.Tσ, isosdb, ∂_∂M∂b_∂tGr', colormap = :balance, colorrange = (-7e4, 7e4))
hdMG = heatmap!(axdMpG, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bG', colormap = :balance, colorrange = (-2e8, 2e8))
hdMS = heatmap!(axdMpS, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bGr', colormap = :balance, colorrange = (-2e8, 2e8))

# create colorbars the size of the whole data set
cb2 = Colorbar(ga[1,3], hMS, ticks = (0:5e7:1e8, ["0", "5×10⁷","1×10⁸"] ),
size =25)
cb4 = Colorbar(ga[2,3], hdtMS, ticks = (-3e4:3e4:3e4, ["-3×10⁴", "0", "3×10⁴"] ), 
size =25)
cb6 = Colorbar(ga[3,3], hdMS, ticks = (-1e8:1e8:1e8,  ["-10⁸", "0", "10⁸"] ),
 size =25)

colsize!(ga, 3, Relative(0.05))

dye_centy = dye_centz/-pm.Tanα
dye_centz_raised = dye_centz + 120*cos(al)
dye_centy_raised = dye_centy + 120*sin(al)

bigtitle = @sprintf("Tracers in Buoyancy Space, U₀=%0.2f, N=%0.2f ×10⁻³, δ=%0.1f", pm.U₀, pm.Ñ*1e3, pm.U₀/pm.Ñ)

Label(f[1, 1, Top()],bigtitle, valign = :bottom,
    font = :bold, fontsize = 30,
    padding = (0, 0, 50, 0))

savename = "MbGaussClasses_rWavg_" * @sprintf("n%d_", ini_Nisos) * setname

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

    return c_weighted_bavg
end

(b̄_Cg) = calculate_tracer_weight(CGi, Bi)
(b̄_Cgr) = calculate_tracer_weight(CGri, Bi)

#########################
#                   BOUYANCY dM/dB PLOT w/ TRCAER WTD
#########################

f = Figure(resolution = (1000, 1100), fontsize=26)
ga = f[1, 1] = GridLayout() 

# change in volume
axMdbpG = Axis(ga[1, 1], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axMdbpS = Axis(ga[1, 2], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdtMpG = Axis(ga[2, 1], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdtMpS = Axis(ga[2, 2], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdMpG = Axis(ga[3, 1], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)
axdMpS = Axis(ga[3, 2], ylabel = "Buoynacy [ms⁻²]", 
            xlabel = "Wave Periods [Tσ]",)

Label(ga[1, 1:2, Top()],"Tracer Distribution ∂M/∂b", valign = :top,
font = :bold, fontsize = 25,
padding = (0, 0, 15, 0))

Label(ga[2, 1:2, Top()],"Diapycnal Transport ∂(∂M/∂b)/∂t", valign = :top,
    font = :bold, fontsize = 25,
    padding = (0, 0, 5, 0))

Label(ga[3, 1:2, Top()],"Distribution Changes ∂M/∂b - ∂M₀/∂b", valign = :top,
font = :bold, fontsize = 25,
padding = (0, 0, 5, 0))

axdMpG.xticks =  2:2:10
axdMpS.xticks =  2:2:10

axMdbpG.yticks =  (bonslope₀ - 5e-4:5e-4:bonslope₀ + 5e-4, ["b₀ - 5×10⁻⁴", "b₀", "b₀ + 5×10⁻⁴"])
axdtMpG.yticks =  (bonslope₀- 5e-4:5e-4:bonslope₀ + 5e-4, ["b₀ - 5×10⁻⁴", "b₀", "b₀ + 5×10⁻⁴"])
axdMpG.yticks =  (bonslope₀ - 5e-4:5e-4:bonslope₀ + 5e-4, ["b₀ - 5×10⁻⁴", "b₀", "b₀ + 5×10⁻⁴"])

limits!(axMdbpG, 1, 10,bonslope₀ - 7e-4,bonslope₀ + 7e-4)
limits!(axMdbpS, 1, 10, braised₀ - 7e-4 , braised₀ + 7e-4)
limits!(axdtMpG, 1, 10, bonslope₀ -7e-4, bonslope₀ + 7e-4)
limits!(axdtMpS, 1, 10, braised₀ - 7e-4 , braised₀ + 7e-4)
limits!(axdMpG, 1, 10, bonslope₀ - 7e-4, bonslope₀ + 7e-4)
limits!(axdMpS, 1, 10, braised₀ - 7e-4 , braised₀ + 7e-4)

hidexdecorations!(axMdbpG)
hidexdecorations!(axdtMpG)
hidedecorations!(axMdbpS)
hidedecorations!(axdtMpS)
hideydecorations!(axdMpS)

hMG = heatmap!(axMdbpG, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bG', colormap = :tempo, colorrange = (0, 8e7))
text!(axMdbpG, Point.(3, bonslope₀ - 0.6*1e-3), text = "On Slope: (-250 m, 703 m)", align = (:left, :center), 
        color = :black, fontsize = 20, font = :bold)
lines!(axMdbpG, b_timeseries.times./pm.Tσ, b̄_Cg, color = :dodgerblue2, linewidth = 5)

hMS = heatmap!(axMdbpS, b_timeseries.times./pm.Tσ, isosdb, ∂M∂bGr', colormap = :tempo, colorrange = (0, 2.2e8))
text!(axMdbpS, Point.(3, braised₀ - 0.6*1e-3), text = "Up 100m: (-136 m , 743 m)", align = (:left, :center), 
        color = :black, fontsize = 20, font = :bold)
lines!(axMdbpS, b_timeseries.times./pm.Tσ, b̄_Cgr, color = :dodgerblue2, linewidth = 5)

hdtMG = heatmap!(axdtMpG , b_timeseries.times./pm.Tσ, isosdb, ∂_∂M∂b_∂tG', colormap = :balance, colorrange = (-7e4, 7e4))
hdtMS = heatmap!(axdtMpS, b_timeseries.times./pm.Tσ, isosdb, ∂_∂M∂b_∂tGr', colormap = :balance, colorrange = (-7e4, 7e4))
hdMG = heatmap!(axdMpG, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bG', colormap = :balance, colorrange = (-2e8, 2e8))
hdMS = heatmap!(axdMpS, b_timeseries.times./pm.Tσ, isosdb, Δ∂M∂bGr', colormap = :balance, colorrange = (-2e8, 2e8))

# create colorbars the size of the whole data set
cb2 = Colorbar(ga[1,3], hMG, ticks = (0:5e7:1e8, ["0", "5×10⁷","1×10⁸"] ),
size =25, flipaxis = false)
cb2 = Colorbar(ga[1,3], hMS, ticks = (0:1e8:2e8, ["0", "1×10⁸", "2×10⁸"] ),
size =25)
cb4 = Colorbar(ga[2,3], hdtMS, ticks = (-3e4:3e4:3e4, ["-3×10⁴", "0", "3×10⁴"] ), 
size =25)
cb6 = Colorbar(ga[3,3], hdMS, ticks = (-1e8:1e8:1e8,  ["-10⁸", "0", "10⁸"] ),
 size =25)

colsize!(ga, 3, Relative(0.05))

bigtitle = @sprintf("Tracers in Buoyancy Space, U₀=%0.2f, N=%0.2f ×10⁻³, δ=%0.1f", pm.U₀, pm.Ñ*1e3, pm.U₀/pm.Ñ)

Label(f[1, 1, Top()],bigtitle, valign = :bottom,
    font = :bold, fontsize = 30,
    padding = (0, 0, 50, 0))

savename = "MbGaussClasses_wcwtdsc_" * @sprintf("n%d_", ini_Nisos) * setname

save(apath * savename * ".png", f)

dye_centy = dye_centz/-pm.Tanα
dye_centz_raised = dye_centz + 120*cos(al)
dye_centy_raised = dye_centy + 120*sin(al)
braised₀ = dye_centz_raised.* pm.Ñ^2
bonslope₀ = dye_centz.* pm.Ñ^2