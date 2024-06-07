using CairoMakie
using Oceananigans
using Printf

path_name = "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"
setname = "U250N100Lz100g100"

@info "getting data from: " * setname

name_prefix = "IntWave_" * setname
filepath = path_name * name_prefix * ".jld2"

v_timeseries = FieldTimeSeries(filepath, "v");
b_timeseries = FieldTimeSeries(filepath, "b");
c_timeseries = FieldTimeSeries(filepath, "Cs");

xc, yc, zc = nodes(c_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# last index on the slope
xlocat = 19

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSame = 1336.6                           # point where planar and curved corner math up the best

ΔySlopeSame = 0

@inline heaviside(X) = ifelse(X < 0, 0., 1.) # heaviside returns 1 if x >= 0
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lz + pm.Lz * exp( - (y - ystar)^2 / (2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα * y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

land = curvedslope.(yc) 

interest_idx3 = [round(Int,(pm.Tσ*3.0)/600), round(Int,(pm.Tσ*3.2)/600)]
interest_idx4 = [round(Int,(pm.Tσ*3.5)/600), round(Int,(pm.Tσ*3.7)/600)]

interest_idxb = [interest_idx3, interest_idx4]

tims = interest_idxb.* 600.0
wavs = tims ./ pm.Tσ

###########
#   UPSLOPE
###########
 
# figure of full progression
f = Figure(resolution = (2000, 700), fontsize=35)
    g = f[1, 1] = GridLayout()

    gcb1 = f[:, 2] = GridLayout()
    gcb2 = f[:, 3] = GridLayout()        

    ##################################
    #
    # z  v1(1,1)    v2(1,2)  (bar)
    # z  c1(2,1)    c2(2,2)  (bar)
    #      y           y           y
    ##################################

    # set up axes and axes labels
    axv1 = Axis(g[1, 1], ylabel = "z [m]")
    axc1 = Axis(g[2, 1], ylabel = "z [m]", xlabel = "y [m]")

    axv2 = Axis(g[1, 2])
    axc2 = Axis(g[2, 2], xlabel = "y [m]")

    axc1.xticks = 500:1000:1500
    axc2.xticks = 500:1000:1500

    # cut limits of axes to specific size

    limits!(axv1, 0, 2000, -450, 0)
    limits!(axc1, 0, 2000, -450, 0)
    limits!(axv2, 0, 2000, -450, 0)
    limits!(axc2, 0, 2000, -450, 0)

    axv = [axv1, axv2]
    axc = [axc1, axc2]

    axv1.yticks = [-250, 0]
    axc1.yticks = [-250, 0]

    time_pre = "t = "
    time_post = " Tσ"

    # for each time
    for (p,j) in enumerate(interest_idx3)

        # get the data
        v = interior(v_timeseries[j], xlocat, :, :)
        b = interior(b_timeseries[j], xlocat, :, :)
        c = log10.(clamp.(interior(c_timeseries[j], xlocat, :, :),1e-8,1))

        phaselabel = time_pre * @sprintf("%0.1f", wavs[1][p]) * time_post

        # create the heatmaps
        global hmv = heatmap!(axv[p], yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
        contour!(axv[p], yc, zc, b, color = :black, linewidth = 3, levels = -0.006:0.001:0, alpha = 0.5)
        band!(axv[p], yc, land, -500, color=:black)
        text!(axv[p],Point.(50, -400), text = phaselabel, align = (:left, :center), color = :white, 
        font = :bold, fontsize = 35)

        global hmc = heatmap!(axc[p], yc, zc, c, colormap= :thermal, colorrange = (-4, 0),)
        contour!(axc[p], yc, zc, b, color = :white, linewidth = 3, levels = -0.006:0.001:0, alpha = 0.5)
        band!(axc[p], yc, land, -500, color=:black)
    end

    # get rid of the inner values
    hidedecorations!(axv2)
    hidexdecorations!(axv1)
    hideydecorations!(axc2)

    colgap!(g, 15)
    rowgap!(g, 5)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hmv, ticks = (-0.2:.1:0.2), size =35, label = "Velocity, v [ms⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hmc, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), 
            size = 35, label = "Tracer Concentration")

# making colorbars take up less space
colsize!(f.layout, 2, Relative(0.05))
colsize!(f.layout, 3, Relative(0.05))    

savename = "waveseries_OSM_upslope_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f, px_per_unit = 2)



###########
#   DOWNSLOPE
###########
 
# figure of full progression
f = Figure(resolution = (2000, 700), fontsize=35)
    g = f[1, 1] = GridLayout()

    gcb1 = f[:, 2] = GridLayout()
    gcb2 = f[:, 3] = GridLayout()        

    ##################################
    #
    # z  v1(1,1)    v2(1,2)  (bar)
    # z  c1(2,1)    c2(2,2)  (bar)
    #      y           y           y
    ##################################

    # set up axes and axes labels
    axv1 = Axis(g[1, 1], ylabel = "z [m]")
    axc1 = Axis(g[2, 1], ylabel = "z [m]", xlabel = "y [m]")

    axv2 = Axis(g[1, 2])
    axc2 = Axis(g[2, 2], xlabel = "y [m]")

    axc1.xticks = 500:1000:1500
    axc2.xticks = 500:1000:1500

    # cut limits of axes to specific size

    limits!(axv1, 0, 2000, -450, 0)
    limits!(axc1, 0, 2000, -450, 0)
    limits!(axv2, 0, 2000, -450, 0)
    limits!(axc2, 0, 2000, -450, 0)

    axv = [axv1, axv2]
    axc = [axc1, axc2]

    axv1.yticks = [-250, 0]
    axc1.yticks = [-250, 0]

    time_pre = "t = "
    time_post = " Tσ"

    # for each time
    for (p,j) in enumerate(interest_idx4)

        # get the data
        v = interior(v_timeseries[j], xlocat, :, :)
        b = interior(b_timeseries[j], xlocat, :, :)
        c = log10.(clamp.(interior(c_timeseries[j], xlocat, :, :),1e-8,1))

        phaselabel = time_pre * @sprintf("%0.1f", wavs[2][p]) * time_post

        # create the heatmaps
        global hmv = heatmap!(axv[p], yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
        contour!(axv[p], yc, zc, b, color = :black, linewidth = 3, levels = -0.006:0.001:0, alpha = 0.5)
        band!(axv[p], yc, land, -500, color=:black)
        text!(axv[p],Point.(50, -400), text = phaselabel, align = (:left, :center), color = :white, 
        font = :bold, fontsize = 35)

        global hmc = heatmap!(axc[p], yc, zc, c, colormap= :thermal, colorrange = (-4, 0),)
        contour!(axc[p], yc, zc, b, color = :white, linewidth = 3, levels = -0.006:0.001:0, alpha = 0.5)
        band!(axc[p], yc, land, -500, color=:black)
    end

    # get rid of the inner values
    hidedecorations!(axv2)
    hidexdecorations!(axv1)
    hideydecorations!(axc2)

    colgap!(g, 15)
    rowgap!(g, 5)

    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hmv, ticks = (-0.2:.1:0.2), size =35, label = "Velocity, v [ms⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hmc, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), 
            size = 35, label = "Tracer Concentration")

# making colorbars take up less space
colsize!(f.layout, 2, Relative(0.05))
colsize!(f.layout, 3, Relative(0.05))    

savename = "waveseries_OSM_downslope_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f, px_per_unit = 2)

