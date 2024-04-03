using GLMakie
using Measures
using Printf
using JLD2

setname = "U250N100Lz100g100"

include("../parameters.jl")

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

ycut = 525


f1 = jldopen(filepath)

ω_x = f1["ω_x"][:,1:ycut,:,:];
#ω_y = f1["ω_y"];
#ω_z = f1["ω_z"];

xc = f1["xc"];
yc = f1["yc"][1:ycut];
zc = f1["zc"];
ctimes = f1["times"];

(cLx, cLy, cLz, tlength) = size(ω_x)
LxCex = 152
land = [curvedslope(y) for x in -5:4:LxCex, y in yc];
title_prefix = "Vertical Vorticity (ω_z)"

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
@info "Starting observable plot..."
n = Observable(1)

title = @lift @sprintf( "Vertical Vorticity (ω_z), t = %.2f hrs, Tσ = %.2f", ctimes[$n]/3600, ctimes[$n]/pm.Tσ)


f = Figure(resolution = (1600, 600),fontsize=26) 
    ga = f[1, 1] = GridLayout()

    ax1 = Axis3(ga[1, 1], azimuth = π/8, # rotation of plot
                elevation = 0.15, 
                aspect =(1,3,1), xtickwidth = 0,
                xticks = ([2,LxCex], [string.(cLx*4), "0"]), xticklabelpad = -25, 
                zticks = ([-500, -250, 0]),
                yticks = ([1000, 2000, 3000]),
                xlabeloffset = 10, ylabeloffset = 25, zlabeloffset = 100,
                xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",
                title = title)
    limits!((0,LxCex), (0,2500), (-500,0))
    surface!(ax1, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100),)

    # topography yz cut lower and uppbounds
    lower = Point3f.(LxCex, yc, -500);
    upper = Point3f.(LxCex, yc, curvedslope.(yc));
    band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

    c_pos = @lift clamp.(ω_z[:,:,:,$n], .002, .02);
    c_neg = @lift clamp.(ω_z[:,:,:,$n], -.02, -.002);

    #c_flat = @lift ω[cLx,:,:,$n];

    cg = contour!(ax1, xc, yc, zc, c_pos, levels = (.005:1e-4:.02), colormap = (cgrad(Makie.to_colormap(:curl)[256:end])), alpha = 0.5)
    contour!(ax1, xc, yc, zc, c_neg, levels = (-.02:1e-4:-.005), colormap = (Reverse(cgrad(Makie.to_colormap(:curl)[1:256]))), alpha = 0.5)

    #contourf!(ax1, yc, zc, c_flat; levels = (-4:0.05:0), colormap = :balance,
    #    transformation=(:yz, LxCex))

    cb = Colorbar(ga[1, 2], limits = (-.02, -.005), colormap = (Reverse(:tempo)), ticks = [-.02,-.015,-.01],
        size = 25, flipaxis=false)
    cb = Colorbar(ga[1, 3], limits = (.005, .02), colormap = (:matter), ticks = [0.01,0.015,0.02],
        label = title_prefix, size = 25)
    colgap!(ga, 5)

    savename = "vorticity_wz_3dwaveseries_" * setname
    #apath  = path_name * "Analysis/"
    apath  =  "Analysis/Plots/"

frames = 1:tlength
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end


@info "Starting observable plot..."
n = Observable(1)

title = @lift @sprintf( "Stream-wise Vorticity (ω_y), t = %.2f hrs, Tσ = %.2f", ctimes[$n]/3600, ctimes[$n]/pm.Tσ)
title_prefix = "Streamwise Vorticity (ω_y)"


f = Figure(resolution = (1600, 600),fontsize=26) 
    ga = f[1, 1] = GridLayout()

    ax1 = Axis3(ga[1, 1], azimuth = π/8, # rotation of plot
                elevation = 0.15, 
                aspect =(1,3,1), xtickwidth = 0,
                xticks = ([2,LxCex], ["0", string.(cLx*4)]), xticklabelpad = -25, 
                zticks = ([-500, -250, 0]),
                yticks = ([1000, 2000, 3000]),
                xlabeloffset = 10, ylabeloffset = 25, zlabeloffset = 100,
                xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",
                title = title)
    limits!((0,LxCex), (0,2500), (-500,0))
    surface!(ax1, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100),)

    # topography yz cut lower and uppbounds
    lower = Point3f.(LxCex, yc[1:625], -500);
    upper = Point3f.(LxCex, yc[1:625], curvedslope.(yc[1:625]));
    band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

    c_pos = @lift clamp.(ω_y[:,1:625,:,$n], .002, .02);
    c_neg = @lift clamp.(ω_y[:,1:625,:,$n], -.02, -.002);

    #c_flat = @lift ω[cLx,:,:,$n];

    cg = contour!(ax1, xc, yc[1:625], zc, c_pos, levels = (.005:1e-4:.02), colormap = (cgrad(Makie.to_colormap(:curl)[256:end])), alpha = 0.5)
    contour!(ax1, xc, yc[1:625], zc, c_neg, levels = (-.02:1e-4:-.005), colormap = (Reverse(cgrad(Makie.to_colormap(:curl)[1:256]))), alpha = 0.5)

    #contourf!(ax1, yc, zc, c_flat; levels = (-4:0.05:0), colormap = :balance,
    #    transformation=(:yz, LxCex))

    cb = Colorbar(ga[1, 2], limits = (-.02, -.005), colormap = (Reverse(:tempo)), ticks = [-.02,-.015,-.01],
        size = 25, flipaxis=false)
    cb = Colorbar(ga[1, 3], limits = (.005, .02), colormap = (:matter), ticks = [0.01,0.015,0.02],
        label = title_prefix, size = 25)
    colgap!(ga, 5)

    savename = "vorticity_wy_3dwaveseries_" * setname
    #apath  = path_name * "Analysis/"
    apath  =  "Analysis/Plots/"

frames = 1:tlength
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end

=#

@info "Starting observable plot..."
n = Observable(1)

title = @lift @sprintf( "Spanwise Vorticity (ωₓ), t = %.2f hrs, Tσ = %.2f", ctimes[$n]/3600, ctimes[$n]/pm.Tσ)
#title_prefix = "Spanwise Vorticity (ω_x)"


f = Figure(resolution = (1650, 600),fontsize=26) 
    ga = f[1, 1] = GridLayout()

    ax1 = Axis3(ga[1, 1], azimuth = π/8, # rotation of plot
                elevation = 0.15, 
                aspect =(1,3,1), xtickwidth = 0,
                xticks = ([2,LxCex], ["0", string.(cLx*4)]), xticklabelpad = -25, 
                zticks = ([-500, -250, 0]),
                yticks = 0:500:2000,
                xlabeloffset = 5, ylabeloffset = 25, zlabeloffset = 100,
                xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",
                title = title)
    limits!((0,LxCex), (0,2100), (-500,0))
    surface!(ax1, -5:4:LxCex, yc, land, shading = true, color=fill((:black, .7),100,100),)

    # topography yz cut lower and uppbounds
    lower = Point3f.(LxCex, yc, -500);
    upper = Point3f.(LxCex, yc, curvedslope.(yc));
    band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

    c_pos = @lift clamp.(ω_x[:,:,:,$n], .002, .02);
    c_neg = @lift clamp.(ω_x[:,:,:,$n], -.02, -.002);

    cg = contour!(ax1, xc, yc, zc, c_pos, levels = plevels, colormap = (cgrad(Makie.to_colormap(:curl)[256:end])), alpha = 0.5)
    contour!(ax1, xc, yc, zc, c_neg, levels = nlevels, colormap = (Reverse(cgrad(Makie.to_colormap(:curl)[1:256]))), alpha = 0.5)

    #contourf!(ax1, yc, zc, c_flat; levels = (-4:0.05:0), colormap = :balance,
    #    transformation=(:yz, LxCex))
    cb = Colorbar(ga[1, 2], limits = (-.02, .02), colormap = newcolorbar, ticks = -.02:.005:.02,
    label = "Spanwise Vorticity, ωₓ [s⁻¹]", size = 25,)
    colgap!(ga, -50)

    savename = "vorticity_wx_3dwaveseries_" * setname
    #apath  = path_name * "Analysis/"
    apath  =  "Analysis/PaperFigures/FinalPaperFigures/FinalFiguresUsed/"

#    ax1.protrusions = (115,0,0,-100)
display(f)


frames = 1:tlength
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end
