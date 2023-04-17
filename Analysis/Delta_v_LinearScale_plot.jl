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
thorpe_havg_tavg = vcat(scale_file["thorpe_havg_tavg"],scale_file["thorpe_havg_tavg_oth"])
thorpe_hmax_tavg = vcat(scale_file["thorpe_hmax_tavg"],scale_file["thorpe_hmax_tavg_oth"])
Cheight_havg_tavg =vcat(scale_file["Cheight_havg_tavg"],scale_file["Cheight_havg_tavg_oth"])
Cheight_hrms_trms = vcat(scale_file["Cheight_hrms_trms"],scale_file["Cheight_hrms_trms_oth"])
Nheight_havg_tavg =vcat(scale_file["Nheight_havg_tavg"],scale_file["Nheight_havg_tavg_oth"])
Nheight_hrms_trms = vcat(scale_file["Nheight_hrms_trms"],scale_file["Nheight_hrms_trms_oth"])

@info "Find Bestfit Lines"

################################  THORPE SCALE v Delta
(b_Ltaδ,m_Ltaδ) = linear_fit(δ, thorpe_havg_tavg)
(b_Ltmδ,m_Ltmδ) = linear_fit(δ, thorpe_hmax_tavg)
################################ INTRUSIONS v Delta
(b_Chaδ,m_Chaδ) = linear_fit(δ, Cheight_havg_tavg)
(b_Chrδ,m_Chrδ) = linear_fit(δ, Cheight_hrms_trms)
################################ Isopycnals v Delta
(b_Nhaδ,m_Nhaδ) = linear_fit(δ, Nheight_havg_tavg)
(b_Nhrδ,m_Nhrδ) = linear_fit(δ, Nheight_hrms_trms)


@info "Plotting!"

VaryN = 1:11
VaryU = 12:22
Varysg = 23:26
Varys = 27:30
Varysg_sub = 25:26
################################ INTRUSIONS v Delta markerstrokecolor=:red
scatter(δ[VaryN], Cheight_havg_tavg[VaryN], label ="Mean, Vary N", markersize = 5, color =:purple, xlabel = "δ [m]", ylabel = "Lₕ [m]",
tickfont = 15, guidefontsize = 20, title = "Mean Dye Intrusion Height (Lₕ) vs δ", size = (1200,800), 
bottom_margin=5mm, left_margin=10mm, legend = :topleft, marker=:c)
scatter!(δ[VaryU], Cheight_havg_tavg[VaryU], label ="Mean, Vary U₀", markersize = 5, marker=:d, color =:purple)
scatter!(δ[Varysg], Cheight_havg_tavg[Varysg], label ="Mean, Vary σ, γ", markersize =7, marker=:star5, color =:purple)
scatter!(δ[Varys], Cheight_havg_tavg[Varys], label ="Mean, Vary σ, l", markersize = 7, marker=:rtriangle, color =:purple)
plot!(δ, m_Chaδ .* δ .+ b_Chaδ, label =  @sprintf("Mean Fit: %0.1f δ + %0.1f", m_Chaδ, b_Chaδ), 
    color=:gray35, lw = 4)
scatter!(δ[VaryN], Cheight_hrms_trms[VaryN], label ="RMS, Vary N", markersize = 5, marker=:c, color = :green)
scatter!(δ[VaryU], Cheight_hrms_trms[VaryU], label ="RMS, Vary U₀", markersize = 5, marker=:d, color = :green)
scatter!(δ[Varysg], Cheight_hrms_trms[Varysg], label ="RMS, Vary σ, γ", markersize = 7, marker=:star5, color =:green)
scatter!(δ[Varys], Cheight_hrms_trms[Varys], label ="RMS, Vary σ, l", markersize = 7, marker=:rtriangle, color =:green)
scatter!(δ[Varysg_sub], Cheight_havg_tavg[Varysg_sub], label ="Subcritical", markersize = 7, marker=:star5, color =:purple, markerstrokecolor=:red)
scatter!(δ[Varysg_sub], Cheight_hrms_trms[Varysg_sub], label ="", markersize = 7, marker=:star5, color =:green, markerstrokecolor=:red)
plot!(δ, m_Chrδ .* δ .+ b_Chrδ, label =  @sprintf("RMS Fit: %0.1f δ + %0.1f", m_Chrδ, b_Chrδ), 
    color=:gray50, lw = 4)
savefig(apath * "CHeightsvDelta_VaryAll.png")

################################  THORPE SCALE v Delta
scatter(δ[VaryN], thorpe_havg_tavg[VaryN], label ="Mean, Vary N", markersize = 5, color =:purple, xlabel = "δ [m]", ylabel = "Lₜ [m]",
tickfont = 15, guidefontsize = 20, title = "RMS Thorpe Scale (Lₜ) vs δ", size = (1100,700), 
bottom_margin=5mm, left_margin=10mm, legend = :topleft, marker=:c)
scatter!(δ[VaryU], thorpe_havg_tavg[VaryU], label ="Mean, Vary U₀", markersize = 5, marker=:d, color =:purple)
scatter!(δ[Varysg], thorpe_havg_tavg[Varysg], label ="Mean, Vary σ, γ", markersize =7, marker=:star5, color =:purple)
scatter!(δ[Varys], thorpe_havg_tavg[Varys], label ="Mean, Vary σ, l", markersize = 7, marker=:rtriangle, color =:purple)
plot!(δ, m_Ltaδ .* δ .+ b_Ltaδ, label =  @sprintf("Mean Fit: %0.1f δ + %0.1f", m_Ltaδ, b_Ltaδ), 
    color=:gray35, lw = 4)
scatter!(δ[VaryN], thorpe_hmax_tavg[VaryN], label ="Max, Vary N", markersize = 5, marker=:c, color = :green)
scatter!(δ[VaryU], thorpe_hmax_tavg[VaryU], label ="Max, Vary U₀", markersize = 5, marker=:d, color = :green)
scatter!(δ[Varysg], thorpe_hmax_tavg[Varysg], label ="Max, Vary σ, γ", markersize = 7, marker=:star5, color =:green)
scatter!(δ[Varys], thorpe_hmax_tavg[Varys], label ="Max, Vary σ, l", markersize = 7, marker=:rtriangle, color =:green)
scatter!(δ[Varysg_sub], thorpe_havg_tavg[Varysg_sub], label ="Subcritical", markersize = 7, marker=:star5, color =:purple, markerstrokecolor=:red)
scatter!(δ[Varysg_sub], thorpe_hmax_tavg[Varysg_sub], label ="", markersize = 7, marker=:star5, color =:green, markerstrokecolor=:red)
plot!(δ, m_Ltmδ .* δ .+ b_Ltmδ, label =  @sprintf("Max Fit: %0.1f δ + %0.1f", m_Ltmδ, b_Ltmδ),
        color=:gray50, lw = 4)
savefig(apath * "ThorpevDelta_VaryAll.png")

################################ Isopycnals v Delta
scatter(δ[VaryN], Nheight_havg_tavg[VaryN], label ="Mean, Vary N", markersize = 5, color =:purple, xlabel = "δ [m]", ylabel = "Lₙ₂ [m]",
tickfont = 15, guidefontsize = 20, title = "Mean N'²<0 Intrusion Height (Lₙ₂) vs δ", size = (1100,700), 
bottom_margin=5mm, left_margin=10mm, legend = :topleft, marker=:c)
scatter!(δ[VaryU], Nheight_havg_tavg[VaryU], label ="Mean, Vary U₀", markersize = 5, marker=:d, color =:purple)
scatter!(δ[Varysg], Nheight_havg_tavg[Varysg], label ="Mean, Vary σ, γ", markersize =7, marker=:star5, color =:purple)
scatter!(δ[Varys], Nheight_havg_tavg[Varys], label ="Mean, Vary σ, l", markersize = 7, marker=:rtriangle, color =:purple)
plot!(δ, m_Nhaδ .* δ .+ b_Nhaδ, label =  @sprintf("Mean Fit: %0.1f δ + %0.1f", m_Nhaδ, b_Nhaδ), 
    color=:gray35, lw = 4)
scatter!(δ[VaryN], Nheight_hrms_trms[VaryN], label ="RMS, Vary N", markersize = 5, marker=:c, color = :green)
scatter!(δ[VaryU], Nheight_hrms_trms[VaryU], label ="RMS, Vary U₀", markersize = 5, marker=:d, color = :green)
scatter!(δ[Varysg], Nheight_hrms_trms[Varysg], label ="RMS, Vary σ, γ", markersize = 7, marker=:star5, color =:green)
scatter!(δ[Varys], Nheight_hrms_trms[Varys], label ="RMS, Vary σ, l", markersize = 7, marker=:rtriangle, color =:green)
scatter!(δ[Varysg_sub], Nheight_havg_tavg[Varysg_sub], label ="Subcritical", markersize = 7, marker=:star5, color =:purple, markerstrokecolor=:red)
scatter!(δ[Varysg_sub], Nheight_hrms_trms[Varysg_sub], label ="", markersize = 7, marker=:star5, color =:green, markerstrokecolor=:red)
plot!(δ, m_Nhrδ .* δ .+ b_Nhrδ, label =  @sprintf("RMS Fit: %0.1f δ + %0.1f", m_Nhrδ, b_Nhrδ), 
    color=:gray50, lw = 4)
savefig(apath * "NHeightsvDelta_VaryAll.png")
