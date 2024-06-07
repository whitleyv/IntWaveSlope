using GLMakie
using Measures
using Printf
using JLD2

setname = "U250N100Lz100g100"

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

filepath = path_name * "vorticity_terms_" * setname * ".jld2"

f1 = jldopen(filepath)

ctimes = f1["times"];
wtimes = ctimes./ pm.Tσ
time_idxs = [45, 50, 56]
ycut = 525
xc = f1["xc"];
yc = f1["yc"][1:ycut];
zc = f1["zc"];

ω_x = f1["ω_x"][:,1:ycut,:,time_idxs];
ω_y = f1["ω_y"][:,1:ycut,:,time_idxs];

(cLx, cLy, cLz, tlength) = size(ω_x)
LxCex = 152
land = [curvedslope(y) for x in -5:4:LxCex, y in yc];

#title_prefix = "Vertical Vorticity (ω_z)"


#title = @lift @sprintf( "Vertical Vorticity (ω_z), t = %.2f hrs, Tσ = %.2f", ctimes[$n]/3600, ctimes[$n]/pm.Tσ)

##################################
#
# z  x1(1,1)    x2(1,2)    x3(1,3)  x4(1,4) (bar)
# z  y1(2,1)    y2(2,2)    y3(2,3)  y4(2,4) (bar)
#      y           y           y
##################################

    # creating a custom color map to cover the values I want and show white in the middle
    # there are 256 colors in a names colormap

color_ranges_from =  -.02:(.04/255):.02 # full range of values that we will be covering in the bar
n005_idx = sum(color_ranges_from .< -0.005) + 1
p005_idx = sum(color_ranges_from .< 0.005) + 1
white_num = p005_idx - n005_idx
fraction_ofwhite = white_num/length(color_ranges_from)
white_grad_length = round(Int, fraction_ofwhite*(2*256)/(1-fraction_ofwhite))
newcolorbar = vcat(vcat(Makie.to_colormap(:curl)[1:256],
                        Makie.to_colormap(cgrad([:white, :white]))[1:white_grad_length]), 
                        Makie.to_colormap(:curl)[256:end])

plevels = .005:1e-4:.02
nlevels = -.02:1e-4:-.005
      
#=
# figure of full progression
f = Figure(resolution = (2300, 700), fontsize=26)
    gc = f[1, 1] = GridLayout()

    gcb1 = f[1, 2] = GridLayout()

    # set up axes and axes labels
    axωx1 = Axis3(gc[1, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), 
        zticks = ([-500, -250, 0]), yticks = ([1000, 2000]), ylabeloffset = 25, zlabeloffset = 100,
        ylabel = "y [m]", zlabel = "z [m]")
    axωx2 = Axis3(gc[1, 2], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), zticks = ([-500, -250, 0]),
        yticks = ([1000, 2000]), ylabeloffset = 25,  ylabel = "y [m]")
    axωx3 = Axis3(gc[1, 3], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, zticks = ([-500, -250, 0]),
        xticks = ([0, 150]), xticklabelpad = -25, xlabeloffset = 30,
        yticks = ([1000, 2000]), ylabeloffset = 25,
        xlabel="x [m]", ylabel = "y [m]")
        #left, right, bottom, top

    axωy1 = Axis3(gc[2, 1],  azimuth = π/8, elevation = 0.15, aspect =(1,3,1), 
        zticks = ([-500, -250, 0]), yticks = ([1000, 2000]), ylabeloffset = 25, zlabeloffset = 100,
        ylabel = "y [m]", zlabel = "z [m]")
    axωy2 = Axis3(gc[2, 2],  azimuth = π/8, elevation = 0.15, aspect =(1,3,1), zticks = ([-500, -250, 0]),
        yticks = ([1000, 2000]), ylabeloffset = 25,  ylabel = "y [m]")
    axωy3 = Axis3(gc[2, 3],  azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0,zticks = ([-500, -250, 0]),
        xticks = ([0, 150]), xticklabelpad = -25, xlabeloffset = 30,
        yticks = ([1000, 2000]), ylabeloffset = 25,
        xlabel="x [m]", ylabel = "y [m]", )

    colsize!(f.layout, 2, Relative(0.05))

    axes_ωx = [axωx1, axωx2, axωx3]
    axes_ωy = [axωy1, axωy2, axωy3]

    lower = Point3f.(LxCex, yc, -500);
    upper = Point3f.(LxCex, yc, curvedslope.(yc));

    for (p,ax) in enumerate(axes_ωx)
        ax.limits = ((0,LxCex), (0,2100), (-500,0))
        surface!(ax, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100))
        band!(ax, lower, upper, color = (:black, 0.9), shading = true)

    end
    for ax in axes_ωy
        ax.limits = ((0,LxCex), (0,2100), (-500,0))
        surface!(ax, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100))
        band!(ax, lower, upper, color = (:black, 0.9), shading = true)
    end

    # get rid of the inner values
    hidexdecorations!(axωx1, grid = false)
    hidexdecorations!(axωx2, grid = false)

    hidezdecorations!(axωx2, grid = false)
    hidezdecorations!(axωx3, grid = false)

    hidezdecorations!(axωy2, grid = false)
    hidezdecorations!(axωy3, grid = false)

    hidexdecorations!(axωy1, grid = false)
    hidexdecorations!(axωy2, grid = false)

    cb = Colorbar(gcb1[1:2, 1], limits = (-.02, .02), colormap = newcolorbar, ticks = -.02:.005:.02,
    size = 25,label = "ω [s⁻¹]")

    colgap!(gc, -100)
    rowgap!(gc, 20)
    colgap!(f.layout, -20)
    
    time_pre = "t = "
    time_post = " Tσ"

    # for each time
    for (p,j) in enumerate(time_idxs)
        c_pos_ωx = clamp.(ω_x[:,:,:,p], .002, .02);
        c_pos_ωy = clamp.(ω_y[:,:,:,p], .002, .02);

        c_neg_ωx = clamp.(ω_x[:,:,:,p], -.02, -.002);
        c_neg_ωy = clamp.(ω_y[:,:,:,p], -.02, -.002);

        phaselabel = time_pre * @sprintf("%0.1f", wtimes[j]) * time_post

        # create the heatmaps
        cp_ωx = contour!(axes_ωx[p], xc, yc, zc, c_pos_ωx, levels = plevels, colormap = (cgrad(Makie.to_colormap(:curl)[256:end])), alpha=0.5)
        cn_ωx = contour!(axes_ωx[p], xc, yc, zc, c_neg_ωx, levels = nlevels, colormap = (Reverse(cgrad(Makie.to_colormap(:curl)[1:256]))), alpha = 0.5)
        text!(axes_ωx[p], Point.(2, 1500, -70), text = phaselabel, align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26, rotation = -π/60)

        contour!(axes_ωy[p], xc, yc, zc, c_pos_ωy, levels = plevels, colormap = (cgrad(Makie.to_colormap(:curl)[256:end])), alpha = 0.5)
        contour!(axes_ωy[p], xc, yc, zc, c_neg_ωy, levels = nlevels, colormap = (Reverse(cgrad(Makie.to_colormap(:curl)[1:256]))), alpha = 0.5)
        text!(axes_ωy[p], Point.(2, 1500, -70), text = phaselabel, align = (:left, :center), color = :black, 
        font = :bold, fontsize = 26, rotation = -π/60)

        #Label(gc[1, p, Top()], phaselabel, fontsize = 30, font = :bold, halign = :center, padding = (0,0,70,0))
    end

    axωx1.protrusions = (100,0,0,0)

    Label(gc[1, 1, Top()], "Spanwise Vorticity", fontsize = 30, font = :bold, halign = :left )
    Label(gc[2, 1, Top()], "Streamwise Vorticity", fontsize = 30, font = :bold, halign = :left)

savename = "vorticity_wxy_3dwaveseries_" * setname
apath  =  "Analysis/Plots/"

save(apath * savename * ".png", f, px_per_unit = 2)

=#

##################################
#             
# z  x1(1,1) ||
# z  x2(2,1) ||
# z  x3(3,1) ||
#     y       
##################################

      
# figure of full progression
f = Figure(resolution = (1200, 1300), fontsize=30)
    gc = f[1, 1] = GridLayout()

    gcb1 = f[1, 2] = GridLayout()

    # set up axes and axes labels
    axωx1 = Axis3(gc[1, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, 
        zticks = ([-500, -250, 0]), yticks = ([0, 1000, 2000]), ylabeloffset = 25, zlabeloffset = 100,
        ylabel = "y [m]", zlabel = "z [m]", xticks = ([0, 150]), xticklabelpad = -25, xlabeloffset = 35,
        xlabel="x [m]")
    axωx2 = Axis3(gc[2, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, 
        zticks = ([-500, -250, 0]), yticks = ([0, 1000, 2000]), ylabeloffset = 25,  zlabeloffset = 100,
        ylabel = "y [m]", zlabel = "z [m]", xticks = ([0, 150]), xticklabelpad = -25, xlabeloffset = 35,
        xlabel="x [m]")
    axωx3 = Axis3(gc[3, 1], azimuth = π/8, elevation = 0.15, aspect =(1,3,1), xtickwidth = 0, 
        zticks = ([-500, -250, 0]), yticks = ([0, 1000, 2000]), ylabeloffset = 25, zlabeloffset = 100,
        xticks = ([0, 150]), xticklabelpad = -25, xlabeloffset = 35,
        xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]" )
        #left, right, bottom, top

    #colsize!(f.layout, 2, Relative(0.05))

    axes_ωx = [axωx1, axωx2, axωx3]

    lower = Point3f.(LxCex, yc, -500);
    upper = Point3f.(LxCex, yc, curvedslope.(yc));

    for (p,ax) in enumerate(axes_ωx)
        ax.limits = ((0,LxCex), (0,2100), (-500,0))
        surface!(ax, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100))
        band!(ax, lower, upper, color = (:black, 0.9), shading = true)

    end

    # get rid of the inner values
    hideydecorations!(axωx1, grid = false)
    hideydecorations!(axωx2, grid = false)
    hidexdecorations!(axωx1, grid = false)
    hidexdecorations!(axωx2, grid = false)

    #colgap!(gc, -100)

    cb = Colorbar(gcb1[1, 1], limits = (-.02, .02), colormap = newcolorbar, ticks = -.02:.005:.02,
    label = "Spanwise Vorticity, ωₓ [s⁻¹]", height = Relative(3/4), width = 40,
    labelpadding = 10.0, highclip = newcolorbar[end], lowclip = newcolorbar[1])
    
    colgap!(f.layout, 30)
    rowgap!(gc, -250)

    time_pre = "t = "
    time_post = " Tσ"

    # for each time
    for (p,j) in enumerate(time_idxs)
        c_pos_ωx = clamp.(ω_x[:,:,:,p], .002, .02);
        c_neg_ωx = clamp.(ω_x[:,:,:,p], -.02, -.002);

        phaselabel = time_pre * @sprintf("%0.1f", wtimes[j]) * time_post

        # create the heatmaps
        cp_ωx = contour!(axes_ωx[p], xc, yc, zc, c_pos_ωx, levels = plevels, colormap = (cgrad(Makie.to_colormap(:curl)[256:end])), alpha=0.5)
        cn_ωx = contour!(axes_ωx[p], xc, yc, zc, c_neg_ωx, levels = nlevels, colormap = (Reverse(cgrad(Makie.to_colormap(:curl)[1:256]))), alpha = 0.5)
        text!(axes_ωx[p], Point.(2, 1500, -70), text = phaselabel, align = (:left, :center), color = :black, 
        font = :bold, fontsize = 30, rotation = -π/60)

    end

    axωx1.protrusions = (150,0,0,0)
    axωx2.protrusions = (150,0,0,0)
    axωx3.protrusions = (150,0,0,0)

savename = "vorticity_wx_3dwaveseries_" * setname
apath  =  "Analysis/Plots/"

display(f)
save(apath * savename * ".png", f, px_per_unit = 2)
