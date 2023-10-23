using Measures
using Statistics
using Printf
using CurveFit
using JLD2
using CairoMakie

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

apath = "/glade/scratch/whitleyv/NewAdvection/Parameters/"

sn = "U250N100Lz100g100"

resS = 1.0

include("../parameters.jl")

pm = getproperty(SimParams(), Symbol(sn))

dzr = pm.dz * resS
dhr = pm.dh * resS

pm = merge(pm, (; dzr=dzr, dhr=dhr, Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2), nx = round(Int, pm.Lx/dhr),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ
Sp_extra = ifelse(z_start>0, 250.0, 0.0)
Ly = pm.Lyˢ+Sp_extra
ny = round(Int,Ly/pm.dhr)
slope_end = pm.Lzˢ/pm.Tanα

pm = merge(pm, (;Ly=Ly,ny=ny, slope_end=slope_end, Sp_extra=Sp_extra))

# if slope is in different spot than usual, need to move the curved part too!
const zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

apath = "Analysis/Plots/"
dpath = "Data/"
filescalename = dpath * "DeltavAll_mtest.jld2"
scale_file = jldopen(filescalename, "r+")

setnames = vcat(scale_file["setnames"])
sns = vcat(scale_file["sns"])

δ = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.45, 0.3, 0.4, 0.5, 0.55]./ pm.Ñ
δN = δ.^2 .* (pm.Ñ)^3  

thorpe_havg_tavg = scale_file["thorpe_havg_tavg"]
thorpe_hmax_tavg = scale_file["thorpe_hmax_tavg"]
Cheight_havg_tavg = scale_file["Cheight_havg_tavg"]
Cheight_hrms_trms = scale_file["Cheight_hrms_trms"]
eps_endAvg = scale_file["eps_endAvg"]
eps_endMax = scale_file["eps_endMaxAvg"]

LtN = thorpe_hmax_tavg.^2 .* (pm.Ñ)^3  

#### Compared to....

filescalename = dpath * "DeltavAllScale.jld2"
filesetnames =  "SetnamesList.jld2"
scale_file2 = jldopen(filescalename, "r+")

setnames = scale_file2["setnames"]
scale_file_sn = jldopen(filesetnames, "r+")
δs = scale_file_sn["δs"][1:22]

filescalename_ε = dpath * "DeltavDissip.jld2"
scale_file_ε = jldopen(filescalename_ε, "r+")

Ñ2 = zeros(length(setnames[1:22]))
for (m,setname) in enumerate(setnames[1:22])
    pm2 = getproperty(SimParams(), Symbol(setname))
    Ñ2[m] = pm2.Ñ
end

δN2 = δs.^2 .* Ñ2.^3

thorpe_hmax_tavg2 = scale_file2["thorpe_hmax_tavg"][1:22]
Cheight_havg_tavg2 =scale_file2["Cheight_havg_tavg"][1:22]
eps_endAvg2 = scale_file_ε["eps_endAvg"][1:22]

LtN2 = thorpe_hmax_tavg2.^2 .* Ñ2.^3

NoexD = 8:11

f = Figure(resolution = (1400, 1000), fontsize=26)
ga = f[1, 1] = GridLayout()

ax1 = Axis(ga[1, 1],  xlabel = "δ [m]", ylabel = "Tracer Thickness, Lₜᵣ [m]")
ax1.xticks = 0:40:160
ax1.yticks = 0:40:160
limits!(ax1, 0, 175, 0, 175)

ax2 = Axis(ga[2, 1],  xlabel = "δ [m]", ylabel = "Thorpe Scale, Lₜ [m]")
ax2.xticks = 0:40:160
ax2.yticks = 20:20:60
limits!(ax2, 0, 165, 0, 70)

ax3 = Axis(ga[1, 2],  xlabel = "δ²N³ [m²s⁻³]", ylabel = "Dissipation, ε [m²s⁻³]")
ax3.xticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
#ax3.yticks = (0:1e-7:3e-7, ["0", "1×10⁻⁷", "2×10⁻⁷", "3×10⁻⁷"])
limits!(ax3, -3e-5, 1.2*1e-3, 0, 4.5e-5)

ax4 = Axis(ga[2, 2],  xlabel = "Lₜ²N³ [m²s⁻³]", ylabel = "Dissipation, ε [m²s⁻³]")
ax4.xticks = (0:5*1e-4:1e-3, ["0", "5×10⁻⁴", "1×10⁻³"])
#ax4.yticks = (0:1e-7:3e-7, ["0", "1×10⁻⁷", "2×10⁻⁷", "3×10⁻⁷"])
limits!(ax4, -3e-5, 1.2*1e-3, 0, 4.5e-5)

m1 = scatter!(ax1, δs, Cheight_havg_tavg2, markersize = 25, marker = :utriangle, 
color =:dodgerblue2, strokewidth = 1, strokecolor = :black)
m2 = scatter!(ax1, δ, Cheight_havg_tavg, markersize = 25, marker = :circle, 
color =:firebrick2, strokewidth = 1, strokecolor = :black)
m3 = scatter!(ax1,  δ[NoexD], Cheight_havg_tavg[NoexD], markersize = 25, marker = :circle, 
color =:darkgreen, strokewidth = 1, strokecolor = :black)

scatter!(ax2, δs, thorpe_hmax_tavg2, markersize = 25, marker = :utriangle, 
color =:dodgerblue2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δ, thorpe_hmax_tavg, markersize = 25, marker = :circle, 
color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax2, δ[NoexD], thorpe_hmax_tavg[NoexD], markersize = 25, marker = :circle, 
color =:darkgreen, strokewidth = 1, strokecolor = :black)

scatter!(ax3, δN2, eps_endAvg2, markersize = 25, marker = :utriangle, 
color =:dodgerblue2, strokewidth = 1, strokecolor = :black)
scatter!(ax3, δN, eps_endAvg, markersize = 25, marker = :circle, 
color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax3, δN[NoexD], eps_endAvg[NoexD], markersize = 25, marker = :circle, 
color =:darkgreen, strokewidth = 1, strokecolor = :black)

scatter!(ax4, LtN2, eps_endAvg2, markersize = 25, marker = :utriangle, 
color =:dodgerblue2, strokewidth = 1, strokecolor = :black)
scatter!(ax4, LtN, eps_endAvg, markersize = 25, marker = :circle, 
color =:firebrick2, strokewidth = 1, strokecolor = :black)
scatter!(ax4, LtN[NoexD], eps_endAvg[NoexD], markersize = 25, marker = :circle, 
color =:darkgreen, strokewidth = 1, strokecolor = :black)

Legend(ga[1, 1, Top()], [m1, m2, m3], ["m = -π/Lz", "m = π/Lz", "No extra Lz"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 15), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :top, orientation = :horizontal)

savename = apath * "Delta_mtest"
save(savename * ".png", f)
