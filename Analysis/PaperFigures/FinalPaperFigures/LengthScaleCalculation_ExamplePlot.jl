using Measures
using Printf
using JLD2
using Oceananigans
using Statistics
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"

setname = "U250N100Lz100g100"
include("parameters.jl")

pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSame = 1336.6                           # point where planar and curved corner math up the best

ΔySlopeSame = 0

@inline heaviside(X) = ifelse(X < 0, 0., 1.) # heaviside returns 1 if x >= 0
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp( - (y - ystar)^2 / (2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα * y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = @. linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

@info "getting data from: " * setname

name_prefix = "IntWave_" * setname
filepath = path_name * name_prefix * ".jld2"
@info "getting data from: " * setname

c_timeseries = FieldTimeSeries(filepath, "Cs");
xc, yc, zc = nodes(c_timeseries) #CCC

N_timeseries = FieldTimeSeries(filepath, "N2");
xn, yn, zn = nodes(N_timeseries) #CCC

land = curvedslope.(yc) 

lastH = 450 #start at z = -50 (used to be -250)
firstH = 50 #end at z = -450

# indices to start and stop 
z_st = round(Int, lastH/2) # start at z = -50
z_en = round(Int, firstH/2) # end at z = -450 (smaller indices = deeper)
y_st = round(Int, ((pm.Lz-lastH)/pm.Tanα)/4) # find corresponding y value on slope when z = -50
y_en = round(Int, 2500/4) # just choosing this y value to include most of dye excursions

zlength_sm = length(z_en:z_st)
ylength_sm = length(y_st:y_en)

isemp(A) = length(A) < 1
isNemp(A) = length(A) > 1
thresh(z) = z>=10^(-6)
cond_gt0(y) = y > 0 
cond_lt0(y) = y < 0
cond_lte0(y) = y <= 0
cond_gt1(y) = y > 1 # find next overturn
cond_gt4(y) = y > 4 # cutoff for overturns

# rolling averages width
rolWidy = 20
rolWidz = 3

# number of individual results to log at each (y,t) location
numRes = 25

#choosing a good time index
tdx = 160

@info "Pulling Dye values at given time..."

# pulling c values at time step
c = c_timeseries[tdx]
# averaging in x and cutting z range
ci = mean(interior(c)[:,:, z_en:z_st], dims =1)[1,:,:]
# averaging in x and not cutting z range
cifull = mean(interior(c), dims =1)[1,:,:]

zlength = length(zc)

# averaging rolling window

# yi = y index in cut range
# y_md y index from total range
yi= 560
y_md = 560 + y_st - 1

# rolling average at the one y val we wanted in the cut z range
cmi = mean(ci[y_md-rolWidy:rolWidy+y_md, :], dims = 1)

Δc = cmi[1:end-1]-cmi[2:end]
dcdz = Δc /-2


@info "Finidng all tracer minima..."

# finding all "roots"
roots = []
for k = 2:length(dcdz)-2
    if dcdz[k] <= 0 && dcdz[k+1] >0
        roots = [roots;k]
    end
end

top_cond(z) = z>=curvedslope(yc[y_md])
# first index above the slope, fy : end in the fluid
fy = findfirst(top_cond, zc)
# cutting off any roots that are within 4 indices of bottom (6m)
cutbot = sum( roots .< fy + 3)
roots = roots[cutbot+1:end]

rt_widths = []
new_roots_bot = []
new_roots_top = []

for p = 2:length(roots)
    maxc = maximum(cmi[roots[p-1]:roots[p]])
    minr = maxc*.5 
    if maxc > 10^(-4) && cmi[roots[p]] <= minr && cmi[roots[p-1]] <= minr
        rng = roots[p-1]:roots[p]
        fr = findfirst(thresh, cmi[rng])
        lr = findfirst(thresh, reverse(cmi[rng]))
        newr1 = roots[p-1] + fr -1
        newr2 = roots[p] - lr
        z1 = zc[newr1 + z_en]
        z2 = zc[newr2 + z_en]
        new_wid = z2-z1
        new_roots_bot = [new_roots_bot; newr1]
        new_roots_top = [new_roots_top; newr2]
        rt_widths=[rt_widths; new_wid]
    end
end

@info "startification "

ñ2i = mean(interior(N_timeseries), dims=1)[1,:,:,:];
tlength = length(N_timeseries.times)
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

include("../WaveValues.jl")
wave_info=get_wave_indices(N_timeseries, pm, 161)

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

if Wtlength > (W11length + hWl)
    # if room to go to end of "last" wave, then go!
    rWtlength = cutWtlength
else 
    rWtlength = cutWtlength - hWl
end

# rolling wave avg
ñ2i_pert_ryrzrW = zeros(ylength_sm, zlength_sm, rWtlength)
for (i,l) in enumerate(W7length:W7length+rWtlength-1)
    ñ2i_pert_ryrzrW[:,:,i] = mean(ñ2i_pert_ryrz[:,:,wave_info.WavePeriods[l-hWl:l+hWl]], dims=3)[:,:,1]
end

@info "Calculating Unstable Heights..."

All_Nheights = zeros(5)
negs = []
negs = findall(cond_lt0, ñ2i_pert_ryrzrW[yi,:,end])

# assuming there was at least one value that was negative
top_cond(z) = z>=curvedslope(yn[y_md])
# first index above the slope, fy : end in the fluid
fy = findfirst(top_cond, zn[z_en:z_st])
# cutting off any negative values that are within 4 indices of bottom (6m)
cutbot = sum( negs .< fy + 3)
negs = negs[cutbot+1:end]
# difference in anomalies to find sign changes
Nanom_dif = ñ2i_pert_ryrzrW[yi,2:end,end] .- ñ2i_pert_ryrzrW[yi,1:end-1,end] 
# different in negative anomalies
dnegs = negs[2:end] .- negs[1:end-1]
# find all the places there was a jump in negative values
skip_negs = findall(cond_gt1, dnegs)
# the starting points will be the 1st neg + all the values after the jump
st_negs = [negs[1] ; negs[skip_negs.+1]]
en_negs = [negs[skip_negs] ; negs[end]]

new_negs_bot = []
new_negs_top = []

for p = 1:length(st_negs)
    if (en_negs[p] < zlength_sm) & (st_negs[p] > 1)

        anom_en = findfirst(cond_lte0, Nanom_dif[en_negs[p]+1:end])
        anom_st = findlast(cond_gt0,Nanom_dif[1:st_negs[p]-1])
        if (typeof(anom_st) == Int64) & (typeof(anom_en) == Int64) 
            All_Nheights[p] = zn[anom_en+en_negs[p]] .- zn[anom_st]
            new_negs_bot = [new_negs_bot; anom_st]
            new_negs_top = [new_negs_top; anom_en+en_negs[p]]
        end

    end
     
end 


@info "Plotting Appendix Plot..."

# cmifull is length(y_st:100) and full zlength but with y rolling average
# ci full is just x averaged with full y and z1
# cmi is length(y_st:y_en) ie ylength_sm and zlength_sm

cifull_log = log10.(clamp.(cifull, 1e-15, 1e-0))
cmi_logy = log10.(clamp.(cmi, 1e-15, 1e-0))

fullbot_idxs = new_negs_bot.+ z_en
fulltop_idxs = new_negs_top .+ z_en
for i = 2:4
    if fulltop_idxs[i-1] > fullbot_idxs[i]
        fullbot_idxs[i] = new_negs_top[i-1] .+ z_en
        fulltop_idxs[i-1] = new_negs_bot[i] .+ z_en
    end
end

f = Figure(resolution = (1500, 700), fontsize=26)
    ga = f[1, 1] = GridLayout()
    gb = f[1, 2] = GridLayout()

    axtop = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y[m]",
    title = "Tracer at t = 26.7 hrs = 11 Tσ")
    axtop.xticks = 1000:1000:3000
    axtop.yticks = -400:200:-200
    limits!(axtop, 138, 3400, -450, -50)

    δ = pm.U₀/pm.Ñ

    axrt = Axis(gb[1, 1], xlabel = "Concentration", ylabel = "z [m]", 
    title ="Tracer at y = 2374 m")
    axrtN = Axis(gb[1, 2], xlabel = "N² Anomaly", 
    title ="N²' at y = 2374 m")

    # we want a log y axis with ticks at -8 through 0
    axrt.xticks = (-10:4:-2, ["10⁻¹⁰", "10⁻⁶","10⁻²"])
    axrt.yticks = -400:100:-100
    limits!(axrt, -15, 0, -450, -50)

    hideydecorations!(axrtN, grid = false)
    axrtN.xticks = (-1e-5:1e-5:1e-5, ["10⁻⁵", "0","10⁻⁵"])
    limits!(axrtN, -2e-5, 2e-5, -450, -50)

    hmc = heatmap!(axtop, yc, zc, cifull_log, colormap= :thermal, colorrange = (-6, 0),)
    lines!(axtop, yc, land, color=:black, lw = 4)
    ploc = lines!(axtop, yc[y_md].*ones(zlength_sm), zc[z_en:z_st], linewidth = 6, color = :dodgerblue2)

    del = hspan!(axrt, zc[new_roots_bot .+ z_en],zc[new_roots_top .+ z_en], color = (:firebrick2, 0.2),)
    pf = lines!(axrt, vec(cmi_logy), zc[z_en:z_st], color = :gray, linewidth = 7,)
    thp = lines!(axrt, -6 .* ones(zlength_sm), zc[z_en:z_st], color = :black, linewidth = 2, linestyle = :dash)
    #text!(axrt, Point.(-6, -130), text = " Threshold", align = (:left, :center), color = :black,
    #                fontsize = 26)
    minp = scatter!(axrt, cmi_logy[new_roots_bot], zc[new_roots_bot .+ z_en], 
        marker=:rect, markersize = 20, color = :dodgerblue2,
        strokecolor = :black, strokewidth = 1)
    scatter!(axrt, cmi_logy[new_roots_top], zc[new_roots_top .+ z_en], 
        marker=:rect, markersize = 20, color = :dodgerblue2,
        strokecolor = :black, strokewidth = 1)
    for (i, wd) in enumerate(rt_widths)
        zval = zc[new_roots_bot .+ z_en][i] .+ 7
        dellabel =  rich(@sprintf("%0.0f m = %0.1f ", wd, wd/(pm.U₀/pm.Ñ)) , "h", subscript("w"))
        text!(axrt, Point.(-6.4, zval), text = dellabel, align = (:right, :center), color = :black,
                    fontsize = 26)
    end

    del = hspan!(axrtN, zn[fullbot_idxs], zn[fulltop_idxs], color = (:firebrick2, 0.2),)
    pf = lines!(axrtN, ñ2i_pert_ryrzrW[yi,:,end], zn[z_en:z_st], color = :gray, linewidth = 7,)
    minp = scatter!(axrtN, ñ2i_pert_ryrzrW[yi,fullbot_idxs .- z_en,end], zn[fullbot_idxs], 
        marker=:rect, markersize = 20, color = :dodgerblue2,
        strokecolor = :black, strokewidth = 1)
    scatter!(axrtN, ñ2i_pert_ryrzrW[yi,fulltop_idxs .- z_en, end], zn[fulltop_idxs], 
        marker=:rect, markersize = 20, color = :dodgerblue2,
        strokecolor = :black, strokewidth = 1)
    for (i, wd) in enumerate(new_negs_top) 
        zval = zn[fullbot_idxs][i] .+ 7
        dellabel =  rich(@sprintf("%0.0f m = %0.1f ", wd, wd/(pm.U₀/pm.Ñ)) , "h", subscript("w"))
        text!(axrtN, Point.(-1e-6, zval), text = dellabel, align = (:right, :center), color = :black,
                    fontsize = 26)
    end

    leg = Legend(ga[1, 1], 
    [ploc, pf, thp, minp, del], tellheight = false,
    ["Profile Location", "Profile", "Measuring Threshold", "Minima After Filtering", "Measured Thickness Between Minima"])

    #colsize!(f.layout, 2, Auto(.8))
    rowsize!(ga, 1, Auto(.7))
    rowgap!(ga, 15)

    Label(ga[2, 1, TopLeft()], "a",
    fontsize = 30,
    font = :bold,
    padding = (5, 5, 5, 5),
    halign = :left)
    Label(gb[1, 1, TopLeft()], "b",
        fontsize = 30,
        font = :bold,
        padding = (5, 5, 5, 5),
        halign = :left)
    Label(gb[1, 2, TopLeft()], "c",
        fontsize = 30,
        font = :bold,
        padding = (5, 5, 5, 5),
        halign = :left)

    apath = path_name * "Analysis/"

savename = apath * "Paper_profileLth_" * setname 
save(savename * ".png", f, px_per_unit = 2)

