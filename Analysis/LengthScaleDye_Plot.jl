using Measures
using Printf
using JLD2
using Oceananigans
using Statistics

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

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

N2_timeseries = FieldTimeSeries(filepath, "N2");
xn, yn, zn = nodes(N2_timeseries) #CCC

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

# pulling c values within range
c = c_timeseries[tdx]
ci = mean(interior(c)[:,:, z_en:z_st], dims =1)[1,:,:]
cifull = mean(interior(c), dims =1)[1,:,:]

zlength = length(zc)

# averaging rolling window
cmi = zeros(ylength_sm, zlength_sm)
cmifull = zeros(length(y_st:1000), zlength)

# yi = y index in cut range
# y_md y index from total range
for (yi, y_md) in enumerate(y_st:y_en)
    cmi[yi,:] = mean(ci[y_md-rolWidy:rolWidy+y_md, :], dims = 1)
end

for (yi, y_md) in enumerate(y_st:1000)
    cmifull[yi,:] = mean(cifull[y_md-rolWidy:rolWidy+y_md, :], dims = 1)
end

Δc = cmi[:,1:end-1]-cmi[:,2:end]
dcdz = Δc /-2

yi= 560
y_md = 560 + y_st - 1

@info "Finidng all tracer minima..."

# finding all "roots"
roots = []
for k = 2:size(dcdz)[2]-2
    if dcdz[yi,k] <= 0 && dcdz[yi,k+1] >0
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
        new_roots_bot = [new_roots_bot; newr1]
        new_roots_top = [new_roots_top; newr2]
        rt_widths=[rt_widths; new_wid]
    end
end

All_Cheights = zeros(ylength_sm, numRes)

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

        rt_widths2 = []

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
                rt_widths2=[rt_widths2; new_wid]
                All_Cheights[yi, w] = new_wid
                w +=1
            end
        end
    end    
end         

@info "Averaging over each vertical profile..."

avgthickness = zeros(length(y_st:1000))
for j = 1:ylength_sm
    yth = All_Cheights[j,:]
    if sum(cond_gt0.(yth))>0
        avgthickness[j] = mean(yth[cond_gt0.(yth)])
    end
end

@info "Calculating Wave Indices..."
tlength = 161
include("WaveValues.jl")
wave_info=get_wave_indices(N2_timeseries, pm, tlength)

hWl = round(Int, wave_info.Wl/2)

@info "Pulling N2 values..."
tdx = 153
# pulling c values within range
N20 = N2_timeseries[1]
N2i = mean(interior(N2_timeseries), dims =1)[1,:,:,tdx-hWl:tdx]
N20i = mean(interior(N20), dims =1)[1,:,:]

N2i_pert = N2i .- N20i

@info "Smoothing N2 values..."

# rolling y avg
ñ2i_pert_ry = zeros(ylength_sm, zlength, 9)
for (yi, y_md) in enumerate(y_st:y_en)
    ñ2i_pert_ry[yi,:, :] = mean(N2i_pert[y_md-rolWidy:rolWidy+y_md, :, :], dims = 1)
end

# rolling z avg
ñ2i_pert_ryrz = zeros(ylength_sm, zlength_sm, 9)
# rolling vertical average over 7 grid points, to smooth things out
for (ki, zk) in enumerate(z_en:z_st)
    ñ2i_pert_ryrz[:,ki,:] = mean(ñ2i_pert_ry[:,zk-rolWidz:zk+rolWidz,:], dims =2)[:,1,:]
end

ñ2i_pert_ryrzrW = mean(ñ2i_pert_ryrz, dims = 3)[:,:,1]

@info "Finding Thicknesses at a y value..."

yi= 560
y_md = 560 + y_st - 1
negs = []
negs = findall(cond_lt0, ñ2i_pert_ryrzrW[yi,:])
top_cond(z) = z>=curvedslope(yn[y_md])
# first index above the slope, fy : end in the fluid
fy = findfirst(top_cond, zn)
# cutting off any negative values that are within 4 indices of bottom (6m)
cutbot = sum( negs .< fy + 3)
negs = negs[cutbot+1:end]

dnegs = negs[2:end] .- negs[1:end-1]
skip_negs = findall(cond_gt1, dnegs)
st_negs = [negs[1] ; negs[skip_negs.+1]]
en_negs = [negs[skip_negs] ; negs[end]]
neg_hts = zn[en_negs] .- zn[st_negs]

@info "Finding results for all time to plot..."


@info "Calculating Wave Indices..."

All_Cheights = zeros(tlength,ylength_sm, numRes)

@info "Calculating Dye Heights at every time..."
for i in 1:tlength

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
                    All_Cheights[i, yi, w] = new_wid
                    w +=1
                end
            end
        end    
        
    end

end

include("WaveValues.jl")
wave_info=get_wave_indices(N_timeseries, pm, tlength)

# total number of data points
Wtlength = wave_info.Wl * wave_info.nTσ
#for rolling wave avg, can't use the last half of wave if no more waves after that...
hWl = floor(Int64, wave_info.Wl/2)

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
ñ2i_pert_ryrzrW = zeros(ylength_sm, zlength_sm, Wtlength-2*hWl)
for (i,l) in enumerate(hWl+1:Wtlength-hWl)
    wdxs = wave_info.WavePeriods[l-hWl:l+hWl]
    ñ2i_pert_ryrzrW[:,:,i] = mean(ñ2i_pert_ryrz[:,:,wdxs], dims=3)[:,:,1]
end

@info "Calculating Unstable Heights..."
All_Nheights = zeros(Wtlength-2*hWl, ylength_sm, numRes)

for i = 1:Wtlength-2*hWl
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

@info "Calculating Dye Height stats..."
Cht_havg = zeros(tlength)

for i = 1:161
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

Nht_havg = zeros(Wtlength-2*hWl)

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
end



#################################################################
############# Plotting
##################################################################

@info "Plotting Both Profiles..."

using CairoMakie

δ = pm.U₀/pm.Ñ
ñ2i_yval = ñ2i_pert_ryrzrW[yi,:]
cmi_log = log10.(clamp.(cmi, 1e-12, 1e-0))
# the one profiel we actually want
cmi_logy = cmi_log[yi,:]

posidx = findall( ñ2i_yval .>= 0)
negidx = findall( ñ2i_yval .< 0)

n2_log = zeros(length(ñ2i_yval))
n2_log[posidx] = log10.(ñ2i_yval[posidx])
n2_log[negidx] = log10.(abs.(ñ2i_yval[negidx]))

f = Figure(resolution = (1100, 900), fontsize=26)
ga = f[1, 1] = GridLayout()

axlt = Axis(ga[1, 1], xlabel = "Tracer Concentration", ylabel = "z [m]")
# we want a log y axis with ticks at -8 through 0
axlt.xticks = (-10:4:-2, ["10⁻¹⁰", "10⁻⁶","10⁻²"])
axlt.yticks = -400:100:-100
limits!(axlt, -12.1, 0, -450, -50)

axrt = Axis(ga[1, 2], xlabel = "Stratification Anomaly",)
#axrt.xscale = Makie.Symlog10(10.0)
#axrt.xtickformat = :Scientific
# we want a log y axis with ticks at -8 through 0
axrt.xticks = (-1e-5:1e-5:1e-5, ["-10⁻⁵", "0","10⁻⁵"])
axrt.yticks = -400:100:-100
limits!(axrt, -2e-5, 2e-5, -450, -50)

hideydecorations!(axrt)

del = hspan!(axlt, zc[new_roots_bot .+ z_en],zc[new_roots_top .+ z_en], color = (:dodgerblue2, 0.2),)
pf = lines!(axlt, cmi_logy, zc[z_en:z_st], color = :gray, linewidth = 7,)
#lines!(axrt, -6 .* ones(zlength_sm), zc[z_en:z_st], color = :black, linewidth = 4)
min = scatter!(axlt, cmi_logy[new_roots_bot], zc[new_roots_bot .+ z_en], marker = :hline, color = :dodgerblue2, 
            markersize = 20, strokewidth = 5)
scatter!(axlt, cmi_logy[new_roots_top], zc[new_roots_top .+ z_en], marker = :hline, color = :dodgerblue2, 
            markersize = 20, strokewidth = 5)
for (i, wd) in enumerate(rt_widths)
    zval = zc[new_roots_bot .+ z_en][i] .+ 7
    dellabel =  @sprintf("%0.0f m = %0.1f δ", wd, wd/δ) 
    text!(axlt, Point.(-6, zval), text = dellabel, align = (:right, :center), color = :black,
                fontsize = 26)
end

del = hspan!(axrt, zn[st_negs .+ z_en],zn[en_negs .+ z_en], color = (:firebrick2, 0.2),)
pf = lines!(axrt, ñ2i_yval, zn[z_en:z_st], color = :gray, linewidth = 7,)
min = scatter!(axrt, ñ2i_yval[st_negs], zn[st_negs .+ z_en], marker = :hline, color = :firebrick2, 
            markersize = 20, strokewidth = 5)
scatter!(axrt, ñ2i_yval[en_negs], zn[en_negs .+ z_en], marker = :hline, color = :firebrick2, 
            markersize = 20, strokewidth = 5)
for (i, wd) in enumerate(neg_hts)
    zval = zn[st_negs .+ z_en][i] .+ 7
    dellabel =  @sprintf("%0.0f m = %0.1f δ", wd, wd/δ) 
    text!(axrt, Point.(-6e-6, zval), text = dellabel, align = (:right, :center), color = :black,
                fontsize = 26)
end

apath = path_name * "Analysis/"

savename = apath * "ProfilesThick_" * setname 
save(savename * ".png", f)




@info "Plotting Both Profiles with ful plot in time too..."

f = Figure(resolution = (1500, 900), fontsize=26)
ga = f[1, 1] = GridLayout()

axlt = Axis(ga[1, 1], xlabel = "Tracer Concentration", ylabel = "z [m]")
# we want a log y axis with ticks at -8 through 0
axlt.xticks = (-10:4:-2, ["10⁻¹⁰", "10⁻⁶","10⁻²"])
axlt.yticks = -400:100:-100
limits!(axlt, -12.1, 0, -450, -50)

axrt = Axis(ga[1, 2], xlabel = "Stratification Anomaly",)
#axrt.xscale = Makie.Symlog10(10.0)
#axrt.xtickformat = :Scientific
# we want a log y axis with ticks at -8 through 0
axrt.xticks = (-1e-5:1e-5:1e-5, ["-10⁻⁵", "0","10⁻⁵"])
axrt.yticks = -400:100:-100
limits!(axrt, -2e-5, 2e-5, -450, -50)

axbt = Axis(ga[2, 1:2], xlabel = "Wave Periods, Tσ", ylabel = "Intrusion Thickness/δ")
# we want a log y axis with ticks at -8 through 0
axbt.xticks = 2:2:10
axbt.yticks = 0:0.5:1.5
limits!(axbt, 0.5, 10.5, 0, 2)

hideydecorations!(axrt)

del = hspan!(axlt, zc[new_roots_bot .+ z_en],zc[new_roots_top .+ z_en], color = (:dodgerblue2, 0.2),)
pf = lines!(axlt, cmi_logy, zc[z_en:z_st], color = :gray, linewidth = 7,)
#lines!(axrt, -6 .* ones(zlength_sm), zc[z_en:z_st], color = :black, linewidth = 4)
min = scatter!(axlt, cmi_logy[new_roots_bot], zc[new_roots_bot .+ z_en], marker = :hline, color = :dodgerblue2, 
            markersize = 20, strokewidth = 5)
scatter!(axlt, cmi_logy[new_roots_top], zc[new_roots_top .+ z_en], marker = :hline, color = :dodgerblue2, 
            markersize = 20, strokewidth = 5)
for (i, wd) in enumerate(rt_widths)
    zval = zc[new_roots_bot .+ z_en][i] .+ 7
    dellabel =  @sprintf("%0.0f m = %0.1f δ", wd, wd/δ) 
    text!(axlt, Point.(-6, zval), text = dellabel, align = (:right, :center), color = :black,
                fontsize = 26)
end

del = hspan!(axrt, zn[st_negs .+ z_en],zn[en_negs .+ z_en], color = (:firebrick2, 0.2),)
pf = lines!(axrt, ñ2i_yval, zn[z_en:z_st], color = :gray, linewidth = 7,)
min = scatter!(axrt, ñ2i_yval[st_negs], zn[st_negs .+ z_en], marker = :hline, color = :firebrick2, 
            markersize = 20, strokewidth = 5)
scatter!(axrt, ñ2i_yval[en_negs], zn[en_negs .+ z_en], marker = :hline, color = :firebrick2, 
            markersize = 20, strokewidth = 5)
for (i, wd) in enumerate(neg_hts)
    zval = zn[st_negs .+ z_en][i] .+ 7
    dellabel =  @sprintf("%0.0f m = %0.1f δ", wd, wd/δ) 
    text!(axrt, Point.(-6e-6, zval), text = dellabel, align = (:right, :center), color = :black,
                fontsize = 26)
end

lines!(axbt, c_timeseries.times./pm.Tσ, Cht_havg, color = :dodgerblue2, linewidth = 7,)
lines!(axbt, N2_timeseries.times[wave_info.WavePeriods[hWl+1:Wtlength-hWl]]./pm.Tσ, Nht_havg, 
            color = :firebrick2, linewidth = 7,)

apath = path_name * "Analysis/"

savename = apath * "ProfilesThick_wTimes_" * setname 
save(savename * ".png", f)




@info "Plotting Appendix Plot..."

# cmifull is length(y_st:100) and full zlength but with y rolling average
# ci full is just x averaged with full y and z1
# cmi is length(y_st:y_en) ie ylength_sm and zlength_sm

cifull_log = log10.(clamp.(cifull, 1e-12, 1e-0))
cmifull_log = log10.(clamp.(cmifull, 1e-12, 1e-0))
cmi_log = log10.(clamp.(cmi, 1e-12, 1e-0))

# the one profiel we actually want
cmi_logy = cmi_log[yi,:]

f = Figure(resolution = (1100, 900), fontsize=26)
ga = f[1, 1] = GridLayout()
gb = f[1, 2] = GridLayout()

axtop = Axis(ga[2, 1], ylabel = "z [m]", xlabel = "y[m]",
title = "Tracer at t = 26.7 hrs = 11 Tσ")
axtop.xticks = 1000:1000:3000
axtop.yticks = -400:200:-200
limits!(axtop, 138, 3400, -450, -50)
hidexdecorations!(axtop)

δ = pm.U₀/pm.Ñ

axbot = Axis(ga[3, 1], ylabel = "Lₜₕ/δ", xlabel = "y[m]", 
title = "Profile-Averaged Thickness")
axbot.xticks = 1000:1000:3000
axbot.yticks = 0:0.4:1.5
limits!(axbot, 138, 3400, 0, 1.5)

axrt = Axis(gb[1, 1], xlabel = "Concentration", ylabel = "z [m]", 
title ="Tracer at y = 2374 m")
# we want a log y axis with ticks at -8 through 0
axrt.xticks = (-10:4:-2, ["10⁻¹⁰", "10⁻⁶","10⁻²"])
axrt.yticks = -400:100:-100
limits!(axrt, -12.1, 0, -450, -50)

hmc = heatmap!(axtop, yc, zc, cifull_log, colormap= :thermal, colorrange = (-6, 0),)
lines!(axtop, yc, land, color=:black, lw = 4)
ploc = lines!(axtop, yc[y_md].*ones(zlength_sm), zc[z_en:z_st], linewidth = 6, color = :dodgerblue2)

spanst = findfirst(cond_gt0, avgthickness) -1 + y_st

barplot!(axbot,yc[y_st:1000],  avgthickness./δ, color = (:firebrick2, 1.0), transparency = false)
barplot!(axbot,yc[y_md],  avgthickness[y_md]./δ, color = :dodgerblue2, width = 1)
#vspan!(axbot, yc[y_st], yc[spanst], color = (:gray, 0.2))
vspan!(axbot, yc[y_en], 3400, color = (:gray, 0.2))
#text!(axbot, Point.(150, 0.4), text = "Does Not Meet Thresholds", align = (:left, :center), color = :black,
#                fontsize = 20)
text!(axbot, Point.(2700, 0.8), text = "Excluded\nDomain", align = (:left, :center), color = :black,
                fontsize = 20)


del = hspan!(axrt, zc[new_roots_bot .+ z_en],zc[new_roots_top .+ z_en], color = (:firebrick2, 0.2),)
pf = lines!(axrt, cmi_logy, zc[z_en:z_st], color = :gray, linewidth = 7,)
#lines!(axrt, -6 .* ones(zlength_sm), zc[z_en:z_st], color = :black, linewidth = 4)
min = scatter!(axrt, cmi_logy[roots], zc[roots .+ z_en], marker = :hline, color = :black, 
            markersize = 20, strokewidth = 5)
minp = scatter!(axrt, cmi_logy[new_roots_bot], zc[new_roots_bot .+ z_en], 
        marker=:rect, markersize = 30, color = :dodgerblue2,
        strokecolor = :black, strokewidth = 1)
scatter!(axrt, cmi_logy[new_roots_top], zc[new_roots_top .+ z_en], 
        marker=:rect, markersize = 30, color = :dodgerblue2,
        strokecolor = :black, strokewidth = 1)
for (i, wd) in enumerate(rt_widths)
    zval = zc[new_roots_bot .+ z_en][i] .+ 7
    dellabel =  @sprintf("%0.0f m = %0.1f δ", wd, wd/(pm.U₀/pm.Ñ)) 
    text!(axrt, Point.(-6, zval), text = dellabel, align = (:right, :center), color = :black,
                fontsize = 26)
end

leg = Legend(ga[1, 1], #patchlabelgap = 10, rowgap = 5,
[ploc, pf, min, minp, del], tellheight = false, #margin = (10.0,0.0, 20.0,10.0),
["Profile Location", "Profile", "All Discovered Minima", "Minima After Filtering Threshold", "Measured Thickness Between Minima"])

yspace = maximum(tight_yticklabel_spacing!, [axtop, axbot])
axtop.yticklabelspace = yspace
axbot.yticklabelspace = yspace

colsize!(f.layout, 2, Auto(.3))
rowsize!(ga, 1, Auto(.7))
rowsize!(ga, 3, Auto(.5))
rowgap!(ga, 15)

apath = path_name * "Analysis/"

savename = apath * "profileLth_" * setname 
save(savename * ".png", f)
