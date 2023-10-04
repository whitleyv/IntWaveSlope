using GLMakie
using LinearAlgebra
using Measures

yy = 0:4:5500
zz = -500:2:0
xx = 0:4:152

# curved topo parameters
gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
gausT_width = 180
ySlopeSame = 1332.22                           # point where planar and curved corner math up the best

# forcing parameters
gausW_width = 500/3
gausW_center = 4500 

Lz = 500
Lx = 152
Ly = 5500

# sponge params
Sp_Region_right = 500                               # size of sponge region on RHS
Sp_Region_left = 500
Sp_Center = 20
Sp_Region_z = 50

# parameters
Ñ = 3.5 * 10^(-3)              # buoyancy frequency
f = Ñ/10.7                    # inertial frequency
σ = 2.2*f                     # tidal frequency

# topography parameters
Tanθ = sqrt((σ^2 - f^2)/(Ñ^2-σ^2))        # slope of internal wave energy propogation
γᶜ = 1.9                                  # when bulk slope is supercritical
Tanα = γᶜ * Tanθ                          # bulk topographic slope
m = -π/Lz
l = sqrt(((π/Lz)^2 * (f^2 - σ^2)) / (σ^2 - Ñ^2))

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -Lz + Lz * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center, gausT_width)) * heaviside(y-ySlopeSame)

land = [curvedslope(y) for x in -5:4:155, y in yy];
land_cut = [curvedslope(y) for x in -5:4:155, y in 0:4:2000];

Fv_wave(y,z, t) = 0.25 * σ * sin(l * y + m * z - σ * t)

gaus(y) = exp( -(y - gausW_center)^2 / (2 * gausW_width^2))

gauss_wave = [gaus(y) for y in yy, z in zz];

ywave = Fv_wave.(4500, zz, 9530) .* 2 .*1e6 .+ 4500;
zwave = zz[5:25:end]

dye_height = 20                               # height of initial dye
dye_smoothing = 10                            # tanh initial dye smoothing
dye_top = 300                                 # tanh initial dye going negative at
dye_centz = -250   
dye_centy = dye_centz/-Tanα
slope_end = Lz/Tanα

# along the slope
@inline dye(x, y, z, dH, sm, kt) = 0.5*( tanh( (dH-z) / sm ) + tanh(kt - z) )
# only dye in the fluid
@inline above_slope(y, z) = ifelse(z >= curvedslope(y), 1.0, 0.0)
# only as far as the slope goes
@inline slope_side(y) = ifelse(y <= slope_end, 1.0, 0.0)

# tanh IC
@inline c₀2(x, y, z) = slope_side(y) * above_slope(y, z) * dye(x, y, z, curvedslope(y)+dye_height, dye_smoothing, dye_top)

gaussforce_shape = [gaus(y) for y in yy, z in zz];
gaussforce_wave = [gaus(y)*  sin(l * y + m * z - σ * 9530) for y in yy, z in zz];

tanh_dye = [c₀2(x, y, z) for x in xx, y in 0:4:2000, z in zz];
tanh_dyeflat = [c₀2(x, slope_end, z) for x in 0:2:154, z in zz];
tanh_dyeflat2 = [c₀2(152, y, z) for y in 0:4:2000, z in zz];

c2 = log10.(clamp.(tanh_dye,7*1e-5,1));
c2_flat = log10.(clamp.(tanh_dyeflat,7*1e-5,1));
c2_flat2 = log10.(clamp.(tanh_dyeflat2,7*1e-5,1));

@inline mask2nd(X) = heaviside(X)* X^2 
@inline right_mask(y) = mask2nd((y-Ly+Sp_Region_right)/(Sp_Region_right))
right_maskshape = [right_mask(y) for y in yy, z in zz];

@inline left_mask(y) = mask2nd((Sp_Region_left-y)/Sp_Region_left)
@inline top_mask(z)= 0.5*(tanh((z+Sp_Center)/Sp_Region_z) + 1)
@inline corner_mask(y, z) = top_mask(z) * left_mask(y)
corner_maskshape = [corner_mask(y,z) for y in yy, z in zz];

f = Figure(resolution = (1700, 600),fontsize=26) #left, right, bottom, top
ga = f[1, 1] = GridLayout()

ax1 = Axis3(ga[1, 1] , azimuth = π/8, # rotation of plot
            elevation = 0.15, 
            aspect =(1,3,1), xtickwidth = 0,
            xticks = ([0,152], ["152", "0"]), xticklabelpad = -25, 
            zticks = ([-500, -250, 0]),
            xlabeloffset = 10, ylabeloffset = 20, zlabeloffset = 100,
            xlabel="x [m]", ylabel = "y [m]", zlabel = "z [m]",)
#Box(ga[1, 1], color = (:red, 0.2), strokewidth = 0)
limits!((0,152), (0,5500), (-500,0))

cb = Colorbar(ga[1, 2], limits = (-4,0), colormap = (:thermal), ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ),
    label = "Tracer Initialization",lowclip = :white, size = 25)

# right = 0 : right side of cb pulls protrusion content inward with an additional padding of 0.
cb.alignmode = Mixed(right = 0)
colsize!(ga, 2, Auto(0.01))
### adding plots 

surface!(ax1, -5:4:155, yy, land, shading = true, color=fill((:black, .7),100,100),)

# topography yz cut lower and uppbounds
lower = Point3f.(152, yy, -500);
upper = Point3f.(152, yy, curvedslope.(yy));
band!(ax1, lower, upper, color = (:black, 0.9), shading = true)

# arrows for wave propagation
# arrows directed at 0 units in the x,z direction and wave strength units in y  
xyz_vec = Point3f.(0, (Fv_wave.(4500, zwave, 9530) .* 2 .*1e6), 0)
wvalues = Fv_wave.(4500, zwave, 9530).* 2 .*1e6 .+ 4500
for x in xx[1:13:end]
    # position of tails at x=x, y=4500, z = each depth 
    xyz_tails = Point3f.(x, Fv_wave.(4500, zwave, 9530) .* 2 .*1e6 .+ 4500, zwave)
    arrows!(ax1, xyz_tails, xyz_vec, shading=false,
        arrowcolor = wvalues, linecolor = wvalues,
        linewidth = 7, arrowsize = Vec3f(10, 20, 80),
        align = :center, colormap = :balance)
end
# annotation
text!(ax1,Point.(152, 1200, -470), text = "α", align = (:right, :center), color = :white, 
    fontsize = 26, font = :bold)
text!(ax1,Point.(76, 4200, 70), text = "Forced Wave", align = (:left, :center), color = :black, 
    fontsize = 26, rotation = -π/60)
text!(ax1,Point.(76, 5400, 70), text = "Sponge", align = (:left, :center), color = :black, 
    fontsize = 26, rotation = -π/60)
text!(ax1, Point.(76, 750, -5), text = "Sponge", align = (:left, :center), color = :black, 
    fontsize = 26, rotation = -π/60)

# sponge regions
contourf!(ax1, yy, zz, right_maskshape; levels = (0.05:0.05:1), colormap = cgrad(:grayC, alpha=0.7), #:Greys_8, 
    transformation=(:yz, -5), linewidth = 5, alpha = 0.6)
contourf!(ax1, yy, zz, corner_maskshape; levels = (0.05:0.05:1), colormap = cgrad(:grayC, alpha=0.7), #:Greys_8, 
    transformation=(:yz, -5), linewidth = 5)
contourf!(ax1, yy, zz, gauss_wave; levels = (0.05:0.05:1), colormap = cgrad(:grayC, alpha=0.7), #:Greys_8, 
    transformation=(:yz, 0), linewidth = 3)

### adding plots to trcaer figures
contour!(ax1, xx, 0:4:2000, zz, c2, levels = (-4:0.05:0), colormap = :thermal)
contourf!(ax1, 0:2:154, zz, c2_flat; levels = (-4:0.05:0), colormap = :thermal,
    transformation=(:xz, slope_end+1))
contourf!(ax1, 0:4:2000, zz, c2_flat2; levels = (-4:0.05:0), colormap = :thermal,
    transformation=(:yz, 152))

resize_to_layout!(f)

save("Analysis/Plots/Poster_SimSetUp.png", f)

