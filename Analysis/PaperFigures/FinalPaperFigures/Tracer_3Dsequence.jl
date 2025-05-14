using GLMakie
using Measures
using Printf
using JLD2

#setname = "U350N100Lz130g100Lx300"
setname = "U350N100Lz100g100"
vsn = "U250N100Lz100g100"

include("../../parameters.jl")

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
path_name = "Data/"

filepath = path_name * "Cs_mp_noS_" * setname * ".jld2"

#filepath = path_name * "DyesOnly_" * setname * ".jld2"

f1 = jldopen(filepath)

ctimes = f1["times"];
wtimes = ctimes./ pm.Tσ
#time_idxs = [45, 50, 56]
#time_idxs = [98, 100, 103]
time_idxs = [69, 74, 78]
ycut = 875
xc = f1["xc"];
yc = f1["yc"][1:ycut];
zc = f1["zc"];

c_timeseries1 = f1["ci"][:,1:ycut,:,time_idxs];
(cLx, cLy, cLz, tlength) = size(c_timeseries1)


LxCex = 152
land = [curvedslope(y) for x in -5:4:LxCex, y in yc];
      #=
# figure of full progression
f = Figure(resolution = (1200, 1300), fontsize=30)
    gc = f[1, 1] = GridLayout()

    gcb1 = f[1, 2] = GridLayout()

    zt = [-500, -250, 0]
    yt = [0, 1000, 2000, 3000]
    xt = [0, 150]
    # set up axes and axes labels
    axc1 = Axis3(gc[1, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, 
        zticks = (zt), yticks = (yt), ylabeloffset = 25, zlabeloffset = 100,
        ylabel = "y [m]", zlabel = "z [m]", xticks = (xt), xticklabelpad = -25, xlabeloffset = 35,
        xlabel="x [m]")
    axc2 = Axis3(gc[2, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, 
        zticks = (zt), yticks = (yt), ylabeloffset = 25,  zlabeloffset = 100,
        ylabel = "y [m]", zlabel = "z [m]", xticks = (xt), xticklabelpad = -25, xlabeloffset = 35,
        xlabel="x [m]")
    axc3 = Axis3(gc[3, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, 
    zticks = (zt), yticks = (yt), ylabeloffset = 25, zlabeloffset = 100,
        xticks = (xt), xticklabelpad = -25, xlabeloffset = 35,
        xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]" )
        #left, right, bottom, top

    #colsize!(f.layout, 2, Relative(0.05))

    axes_c = [axc1, axc2, axc3]

    lower = Point3f.(LxCex, yc, -500);
    upper = Point3f.(LxCex, yc, curvedslope.(yc));

    for (p,ax) in enumerate(axes_c)
        ax.limits = ((0,LxCex), (0,3500), (-500,0))
        surface!(ax, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100))
        band!(ax, lower, upper, color = (:black, 0.9), shading = true)

    end

    # get rid of the inner values
    hideydecorations!(axc1, grid = false)
    hideydecorations!(axc2, grid = false)
    hidexdecorations!(axc1, grid = false)
    hidexdecorations!(axc2, grid = false)

    #colgap!(gc, -100)
    cb = Colorbar(gcb1[1, 1],  limits = (-4,-.5), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
    label = "Tracer Concentration",lowclip = :white, height = Relative(3/4), width = 40,
    labelpadding = 10.0)
    
    colgap!(f.layout, 30)
    rowgap!(gc, -250)

    time_pre = "t = "
    time_post = " Tσ"

    # for each time
    for (p,j) in enumerate(time_idxs)
        c = log10.(clamp.(c_timeseries1[:,:,:,p], 7e-5,10^(-1/2)));
        c_flat = log10.(clamp.(c_timeseries1[cLx,:,:,p], 7e-5,10^(-1/2)));

        phaselabel = time_pre * @sprintf("%0.1f", wtimes[j]) * time_post
        clevels = -4:0.02:-.5
        # create the heatmaps
        cg = contour!(axes_c[p], xc, yc, zc, c, levels = clevels, colormap = (:thermal), alpha = 0.10)
        contourf!(axes_c[p], yc, zc, c_flat; levels = clevels, colormap = :thermal,
            transformation=(:yz, LxCex))
        text!(axes_c[p], Point.(2, 2800, -70), text = phaselabel, align = (:left, :center), color = :black, 
        font = :bold, fontsize = 30, rotation = -π/60)

    end

    for ax in axes_c
        ax.protrusions = (150,0,0,0)
    end
savename = "tracer_3dwaveseries_" * setname
apath  =  "Analysis/Plots/"

display(f)
save(apath * savename * ".png", f, px_per_unit = 2)
=#

c1 = log10.(clamp.(c_timeseries1[:,:,:,1], 7e-5,10^(-1/2)));
c2 = log10.(clamp.(c_timeseries1[:,:,:,2], 7e-5,10^(-1/2)));
c3 = log10.(clamp.(c_timeseries1[:,:,:,3], 7e-5,10^(-1/2)));
c_flat1 = log10.(clamp.(c_timeseries1[cLx,:,:,1], 7e-5,10^(-1/2)));
c_flat2 = log10.(clamp.(c_timeseries1[cLx,:,:,2], 7e-5,10^(-1/2)));
c_flat3 = log10.(clamp.(c_timeseries1[cLx,:,:,3], 7e-5,10^(-1/2)));

# higher resolution:
f = Figure(resolution = (2240, 2122), fontsize=45);
    gc = f[1, 1] = GridLayout()

    #gcb1 = f[1, 2] = GridLayout()

    zt = [-500, -250, 0]
    yt = [0, 1000, 2000, 3000]
    xt = [0, 150]
    # set up axes and axes labels
    axc1 = Axis3(gc[1, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), 
        zticks = (zt), zlabeloffset = 160, zlabel = "z [m]")
    axc2 = Axis3(gc[2, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1),
        zticks = (zt), zlabeloffset = 160, zlabel = "z [m]")
    axc3 = Axis3(gc[3, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, 
        zticks = (zt), yticks = (yt), ylabeloffset = 55, zlabeloffset = 160,
        xticks = (xt), xticklabelpad = -42, xlabeloffset = 50,
        xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]" )
        #left, right, bottom, top

    hideydecorations!(axc1, grid = false)
    hideydecorations!(axc2, grid = false)
    hidexdecorations!(axc1, grid = false)
    hidexdecorations!(axc2, grid = false)

    #colsize!(f.layout, 2, Relative(0.05))

    axes_c = [axc1, axc2, axc3]

    lower = Point3f.(LxCex, yc, -500);
    upper = Point3f.(LxCex, yc, curvedslope.(yc));

    for (p,ax) in enumerate(axes_c)
        ax.limits = ((0,LxCex), (0,3500), (-500,0))
        surface!(ax, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100))
        band!(ax, lower, upper, color = (:black, 0.9), shading = true)

    end


    #colgap!(gc, -100)
    #=
    cb = Colorbar(gcb1[1, 1],  limits = (-4,-.5), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
    label = "Tracer Concentration",lowclip = :white, height = Relative(0.9), width = 80,
    labelpadding = 10.0)
    =#
    cb = Colorbar(gc[1:3, 2],  limits = (-4,-.5), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
    label = "Tracer Concentration",lowclip = :white, height = Relative(0.9), width = 80,
    labelpadding = 10.0)

    colsize!(gc, 2, Relative(0.05))
    rowgap!(gc, -100)

    time_pre = "t = "
    time_post = " Tσ"

    # for each time
    # create the heatmaps
    clevels = -4:0.02:-.5

    cg = contour!(axes_c[1], xc, yc, zc, c1, levels = clevels, colormap = (:thermal), alpha = 0.10)
    contourf!(axes_c[1], yc, zc, c_flat1; levels = clevels, colormap = :thermal,
        transformation=(:yz, LxCex))
    contour!(axes_c[2], xc, yc, zc, c2, levels = clevels, colormap = (:thermal), alpha = 0.10)
    contourf!(axes_c[2], yc, zc, c_flat2; levels = clevels, colormap = :thermal,
        transformation=(:yz, LxCex))
    contour!(axes_c[3], xc, yc, zc, c3, levels = clevels, colormap = (:thermal), alpha = 0.10)
    contourf!(axes_c[3], yc, zc, c_flat3; levels = clevels, colormap = :thermal,
        transformation=(:yz, LxCex))


    for (p,j) in enumerate(time_idxs)
        phaselabel = time_pre * @sprintf("%0.1f", wtimes[j]) * time_post
        text!(axes_c[p], Point.(2, 2800, -70), text = phaselabel, align = (:left, :center), color = :black, 
        font = :bold, fontsize = 50, rotation = -π/60)

    end

    for ax in axes_c
        ax.protrusions = (100,0,100,0)
    end
    colgap!(gc, -50)

    #colgap!(f.layout, -400)

    #resize!(f.scene, (2240, 2122))     # resize the actual framebuffer

    #img = Makie.colorbuffer(f.scene)

    #save("highres_3d_output.png", img)


savename = "tracer_3dwaveseries_betterres_" * setname
apath  =  "Analysis/PaperFigures/FinalPaperFigures/FinalFiguresUsed/"

#display(f)
save(apath * savename * ".png", f)

