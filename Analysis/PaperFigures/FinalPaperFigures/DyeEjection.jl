using CairoMakie
using Oceananigans
using Printf

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"

setname = "U250N100Lz100g100"

@info "getting data from: " * setname

name_prefix = "IntWave_" * setname
filepath = path_name * name_prefix * ".jld2"

v_timeseries = FieldTimeSeries(filepath, "v");
b_timeseries = FieldTimeSeries(filepath, "b");
c_timeseries = FieldTimeSeries(filepath, "Cs");
e_timeseries = FieldTimeSeries(filepath, "ϵ");

xc, yc, zc = nodes(c_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC

include("../parameters.jl")
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

interest_idx1 = [round(Int,(pm.Tσ*4.67)/600), round(Int,(pm.Tσ*4.95)/600), round(Int,(pm.Tσ*5.1)/600)]

tims = interest_idx1 .* 600.0
wavs = tims ./ pm.Tσ

wavsa = interest_idx1 .* (600/pm.Tσ)

# figure of singular initial progression # resolution (x,y)

f1 = Figure(resolution = (1800, 800), fontsize=26)
    ga = f1[1, 1] = GridLayout() # velocity
    gb = f1[2, 1] = GridLayout() # dye
    gc = f1[3, 1] = GridLayout() # dissipation  

    gcb1 = f1[1,2] = GridLayout()
    gcb2 = f1[2,2] = GridLayout()
    gcb3 = f1[3,2] = GridLayout()

    axv1 = Axis(ga[1, 1])
    axv2 = Axis(ga[1, 2])
    axv3 = Axis(ga[1, 3])

    axc1 = Axis(gb[1, 1])
    axc2 = Axis(gb[1, 2])
    axc3 = Axis(gb[1, 3])

    axe1 = Axis(gc[1, 1], ylabel = "z [m]", xlabel = "y [m]")
    axe2 = Axis(gc[1, 2])
    axe3 = Axis(gc[1, 3])

    axs = [axv1, axv2, axv3, axc1, axc2, axc3, axe1, axe2, axe3]
    axv = [axv1, axv2, axv3]
    axc = [axc1, axc2, axc3]
    axe = [axe1, axe2, axe3]

    for ax in axs
        limits!(ax, 285, 1800, -450, -100)
    end

    for ax in axe
        ax.xticks = 500:500:1500
    end

    axv1.yticks = [-350, -150]
    axc1.yticks = [-350, -150]
    axe1.yticks = [-350, -150]

    hidedecorations!(axv2)
    hidedecorations!(axv3)
    hidedecorations!(axc2)
    hidedecorations!(axc3)
    hidexdecorations!(axv1)
    hidexdecorations!(axc1)
    hideydecorations!(axe2)
    hideydecorations!(axe3)

    time_pre = "t = "
    time_post = " Tσ"

    # for each time
    for (p,j) in enumerate(interest_idx1)

        # get the data
        v = interior(v_timeseries[j], xlocat, :, :)
        b = interior(b_timeseries[j], xlocat, :, :)
        c = log10.(clamp.(interior(c_timeseries[j], xlocat, :, :),1e-8,1))
        e = log10.(clamp.(interior(e_timeseries[j], xlocat, :, :),1e-9,1))

        phaselabel = time_pre * @sprintf("%0.1f", wavsa[p]) * time_post

        # create the heatmaps
        global hmv = heatmap!(axv[p], yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
        contour!(axv[p], yc, zc, b, color = :black, lw = 6, levels = -0.006:0.0005:0, alpha = 0.5)
        lines!(axv[p], yc, land, color=:black, lw = 4)
        text!(axv[p], Point.(300, -400), text = phaselabel, align = (:left, :center), color = :black, 
        font = :bold, fontsize = 30)
        #text!(axv[p], Point.(1350, -140), text = phaselabel, align = (:left, :center), color = :black, 
        #font = :bold, fontsize = 30)


        global hmc = heatmap!(axc[p], yc, zc, c, colormap= :thermal, colorrange = (-4, 0),)
        contour!(axc[p], yc, zc, b, color = :white, lw = 6, levels = -0.006:0.0005:0, alpha = 0.5)
        lines!(axc[p], yc, land, color=:black, lw = 4)

        global hme = heatmap!(axe[p], yc, zc, e, colormap= :thermal, colorrange = (-9, -6),)
        contour!(axe[p], yc, zc, b, color = :white, lw = 6, levels = -0.006:0.0005:0, alpha = 0.5)
        lines!(axe[p], yc, land, color=:black, lw = 4)
    end

    δ = pm.U₀/pm.Ñ
    α = atan(pm.Tanα)
    δy = δ*sin(α)
    δz = δ*cos(α)

    arrows!(axc[1], [(250-10)/pm.Tanα], [-250], [0], [-δ], arrowsize = 15, linewidth = 5,  
        arrowcolor = :firebrick2, linecolor = :firebrick2)
    arrows!(axc[1], [(250-10)/pm.Tanα], [-250 - δ], [0], [δ], arrowsize = 15, linewidth = 5,  
        arrowcolor = :firebrick2, linecolor = :firebrick2)
    text!(axc[1], [250/pm.Tanα - 130], [-250 - δ/2], text=rich("h", subscript("w")), align = (:left, :center), 
        color = :firebrick2, font = :bold, fontsize = 25)

    Label(ga[1, 1, Top()], "Horizontal Velocity, v [ms⁻¹]",
        fontsize = 30, font = :bold,
        padding = (5, 5, 5, 10),
        halign = :left)
    Label(gb[1, 1, Top()], "Tracer Concentration, c",
        fontsize = 30, font = :bold,
        padding = (5, 5, 5, 10),
        halign = :left)
    Label(gc[1, 1, Top()], "Dissipation Rate, ε [m²s⁻³]",
        fontsize = 30,font = :bold,
        padding = (5, 5, 5, 10),
        halign = :left)
#=
    text!(axv[1, 1], Point.(300, -400), text = "v [ms⁻¹]", align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26)
    text!(axc[1, 1], Point.(300, -400), text = "c", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)
    text!(axe[1, 1], Point.(300, -400), text = "ε [m²s⁻³]", align = (:left, :center), color = :white, 
        font = :bold, fontsize = 26)
=#
    # create colorbars the size of the whole data set
    cb1 = Colorbar(gcb1[1,1], hmv, ticks = (-0.2:.2:0.2), size =35) #, label = "Horizontal Velocity, v [ms⁻¹]")
    cb2 = Colorbar(gcb2[1,1], hmc, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), 
            size = 35,)# label = "Tracer Concentration, c")
    cb3 = Colorbar(gcb3[1,1], hme, ticks = (-9:1:-7, ["10⁻⁹", "10⁻⁸", "10⁻⁷"] ), 
            size = 35,)# label = "Dissipation, ε [m²s⁻³]")

    # making colorbars take up less space
    colsize!(f1.layout, 2, Relative(0.02))
    #colsize!(f1.layout, 3, Relative(0.05))
    #colsize!(f1.layout, 4, Relative(0.05))
    rowgap!(f1.layout, 7)
savename = "paper_dyeejections_fix_" * setname
apath  = path_name * "Analysis/"

save(apath * savename * ".png", f1, px_per_unit = 2)

