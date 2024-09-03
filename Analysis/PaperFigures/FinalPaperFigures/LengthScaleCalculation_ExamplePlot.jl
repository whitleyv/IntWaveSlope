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

function Find_Tracer_Intrusions(cutWtlength, W7length, W11length, numRes, c_timeseries, yc, zc, wave_info, rolWidy)
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
 
    return All_Cheights
end

tlength = length(c_timeseries.times)

include("WaveValues.jl")
wave_info=get_wave_indices(c_timeseries, pm, tlength)
# total number of data points
Wtlength = wave_info.Wl * wave_info.nTσ
W7length = wave_info.Wl * 7  + 1 
# 4 waves long
cutWtlength = 4*wave_info.Wl
# ends at 10.9Tσ
W11length = wave_info.Wl * 11

All_Cheights = Find_Tracer_Intrusions(cutWtlength, W7length, W11length, numRes, c_timeseries, yc, zc, wave_info, rolWidy)
All_non0Cheights = All_Cheights[All_Cheights .> 2]
All_Cheights_Phavg = zeros(15)
All_non0Cheights_Avg = mean(All_non0Cheights)
for i = 1:15
    Phase_idxs = [i, i+wave_info.Wl, i+ (2*wave_info.Wl), i+ (3*wave_info.Wl)]
    Phase_intrusions = All_Cheights[Phase_idxs, :,:]
    All_Cheights_Phavg[i] = mean(Phase_intrusions[Phase_intrusions .> 2])
end

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


@info "Finidng all tracer minima at time step..."

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

#include("../WaveValues.jl")
#wave_info=get_wave_indices(N_timeseries, pm, 161)

# gather the last 4 waves in a usual sim
# this is 7Tσ - 10.9Tσ or waves 7-10
# in the array this is columns 8:11

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

                end

            end   
        end 
        
    end
end
All_non0Nheights = All_Nheights[cond_gt0.(All_Nheights)]
All_non0Nheights_avg = mean(All_non0Nheights)

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

f = Figure(resolution = (1500, 1000), fontsize=26)
    ga = f[1, 1] = GridLayout()
    gb = f[1, 2] = GridLayout()
    #gc = f[2, 1:2] = GridLayout()
    #gd = f[2, 2] = GridLayout()

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

    axbotlt = Axis(ga[3, 1], ylabel = rich("L", subscript("tr"),"/h",subscript("w")), xlabel = "Tσ",
    )#xminorticksvisible = true, xminorgridvisible = true, yminorticksvisible = true, yminorgridvisible = true)#title = "Phase-Averaged Lₜᵣ")
    axbotlt.xticks = (0:0.25:1, ["0", "0.25", "0.5", "0.75", "1"])
    #axbotlt.xminorticks = IntervalsBetween(2)
    #axbotlt.yminorticks = IntervalsBetween(2)
    axbotlt.yticks = 1:0.2:1.4
    limits!(axbotlt, 1/15, 1, 1, 1.4)

    axbotm = Axis(gb[2, 1], ylabel = "Probability", xlabel = rich("L", subscript("tr"),"/h",subscript("w")),
    title = rich("L", subscript("tr"), " Intrusion Distribution"))
    axbotm.xticks = 0:1:3
    #axbotm.yticks = (1000:1000:2000 , ["1×10³", "2×10³"])
    #axbotm.yminorticks = IntervalsBetween(2)
    limits!(axbotm, 0, 3.5, 0, 0.2)

    axbotrt = Axis(gb[2, 2],  xlabel = rich("L", subscript("N", superscript("2")), "/h",subscript("w")),
    title = rich("L", subscript("N", superscript("2")), " Intrusion Distribution"))
    axbotrt.xticks = 0:1:3
    #axbotrt.yticks = (5000:5000:10000, ["5×10³", "1×10⁴"])
    #axbotrt.yminorticks = IntervalsBetween(2)
    limits!(axbotrt, 0, 3.5, 0, 0.2)
    hideydecorations!(axbotrt, grid = false)

        # phase averaged plot and heatmap have same y axis labels size
        yspace = maximum(tight_yticklabel_spacing!, [axbotlt, axtop])
        axbotlt.yticklabelspace = yspace
        axtop.yticklabelspace = yspace
    
        # phase averaged plot and heatmap have same y axis labels size
        yspace = maximum(tight_yticklabel_spacing!, [axbotrt, axrt])
        axbotrt.yticklabelspace = yspace
        axrt.yticklabelspace = yspace
    
        # heatmap and profiles have same x axis labels size
        xspace = maximum(tight_xticklabel_spacing!, [axtop, axrt])
        axtop.xticklabelspace = xspace
        axrt.xticklabelspace = xspace
    
    
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
    
    avgc = lines!(axbotlt, 1/15:1/15:1, All_non0Cheights_Avg./δ .* ones(15), color = :gray30, linestyle = :dash, linewidth = 5)
    avgphc = lines!(axbotlt, 1/15:1/15:1, All_Cheights_Phavg./δ, color = :black, linewidth = 7)

    h1 = hist!(axbotm, All_non0Cheights./δ, bins = 20, strokewidth = 1, strokecolor = :black, color = :dodgerblue2, normalization = :probability)
    lines!(axbotm, All_non0Cheights_Avg./δ .* ones(21), 0:1e3:2e4, color = :gray30, linestyle = :dash, linewidth = 5)

    #lines!(axtop, yc[y_md].*ones(zlength_sm), zc[z_en:z_st], linewidth = 6, color = :dodgerblue2)

    h2 = hist!(axbotrt, All_non0Nheights./δ, bins = 20, strokewidth = 1, strokecolor = :black, color = :dodgerblue2, normalization = :probability)
    lines!(axbotrt, All_non0Nheights_avg./δ .* ones(21), 0:1e3:2e4, color = :gray30, linestyle = :dash, linewidth = 5)

    leg = Legend(ga[1, 1], 
    [ploc, pf, thp, minp, del, avgc, avgphc], tellheight = false,
    ["Profile Location", "Profile", "Measuring Threshold", 
    "Minima After Filtering", "Measured Thickness Between Minima",
    "Average Intrusion Thickness", "Phase-averaged Intrusion Thickness"])

    #colsize!(f.layout, 2, Auto(.8))
    rowsize!(ga, 1, Auto(.8)) # the legend proportional to the rest of ga takes up 0.7
    #rowgap!(ga, 20)
    #rowsize!(f.layout, 2, Auto(0.4))
    rowsize!(ga, 3, Relative(0.3))
    rowsize!(gb, 2, Relative(0.3))
    #rowgap!(f.layout, 5)
    #colsize!(gc, 1, Relative(0.4))
    #colgap!(gc, 15)
    #rowgap!(gd, 15)

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

    Label(ga[3, 1, TopLeft()], "d",
        fontsize = 30,
        font = :bold,
        padding = (5, 5, 5, 5),
        halign = :left)
    Label(gb[2, 1, TopLeft()], "e",
        fontsize = 30,
        font = :bold,
        padding = (5, 5, 5, 5),
        halign = :left)
    Label(gb[2, 2, TopLeft()], "f",
        fontsize = 30,
        font = :bold,
        padding = (5, 5, 5, 5),
        halign = :left)

    apath = path_name * "Analysis/"

savename = apath * "Paper_profileLth_whist_" * setname 
save(savename * ".png", f, px_per_unit = 2)

