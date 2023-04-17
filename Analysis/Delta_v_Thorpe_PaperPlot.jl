using Plots
using Measures
using Statistics
using Printf
using Oceananigans
using CurveFit
using ArgParse
using JLD2

ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

apath = "/glade/scratch/whitleyv/NewAdvection/Parameters/"

sn = "U250Nfd250Lz100g100"

resS = 1.0

include("parameters.jl")

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
filescalename = dpath * "DeltavAllScale.jld2"

scale_file = jldopen(filescalename, "r+")

setnames = vcat(scale_file["setnames"],scale_file["setnames_oth"])
δ = vcat(scale_file["δ"],scale_file["δ_oth"])
thorpe_hmax_tavg = vcat(scale_file["thorpe_hmax_tavg"],scale_file["thorpe_hmax_tavg_oth"])

@info "Find Bestfit Lines"

(b_Ltmδ,m_Ltmδ) = linear_fit(δ, thorpe_hmax_tavg)


@info "Plotting!"

# varying N or U to change delat value in sim, varying N changes sigma and f respectively as well
VaryN = 1:11
VaryU = 12:22
# varying sigma to change criticality, holding N and f constant
Varysg = 23:26
# varying sigma but not changing criticaility, holding N and f const but changing topography
Varys = 27:30
# ones that are subcritical out of sigma varying set
Varysg_sub = 25:26

# simulations where N is varied to change delta
plot(δ, m_Ltmδ .* δ .+ b_Ltmδ, label =  "", color=:gray50, lw = 4)
scatter!(δ[Varysg], thorpe_hmax_tavg[Varysg], label ="", markersize = 8, marker=:pentagon, color =:green, markeralpha = 0.6,)
annotate!(δ[Varysg_sub], thorpe_hmax_tavg[Varysg_sub] .+ 2, text("Subcritical γ", 10))
# varying sigma without changing criticality, changing slope instead
scatter!(δ[Varys], thorpe_hmax_tavg[Varys], label ="", markersize = 8, marker=:pentagon, color =:green, markeralpha = 0.6,)
scatter!(δ[VaryN], thorpe_hmax_tavg[VaryN], label ="Vary N", markersize = 8, color =:green, xlabel = "δ [m]", ylabel = "Lₜ [m]",
tickfont = 15, guidefontsize = 20, titlefont = 22, title = "Maximum RMS Thorpe Scale (Lₜ) Scales Linearly with δ", 
bottom_margin=5mm, left_margin=10mm, right_margin = 10mm, legend = :topleft, marker=:c,
xlims = (0, 160), ylims = (0, 70), size = (1100, 700), yticks = 20:20:60, xticks = 0:40:160)
# simulations where U is varied to change delta
scatter!(δ[VaryU], thorpe_hmax_tavg[VaryU], label ="Vary V₀", markersize = 8, marker=:d, color = :green)
scatter!(-1 .* δ[Varysg], thorpe_hmax_tavg[Varysg], label ="Vary σ", markersize = 8, marker=:pentagon, color =:green, markeralpha = 0.6,)
# varying sigma to change the criticality 
plot!(δ, -1 .* m_Ltmδ .* δ .- b_Ltmδ, label =  @sprintf("%0.1f δ + %0.1f", m_Ltmδ, b_Ltmδ),
        color=:gray50, lw = 4)

savefig(apath * "Paper_Thorpe_v_Delta_All.png")
