using Statistics
using Printf
using Oceananigans
using JLD2 
using CairoMakie

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"
apath = path_name * "Analysis/"

sn = "U250Nfd250Lz100g100"

include("parameters.jl")

pm = getproperty(SimParams(), Symbol(sn))

resS = 1.0
dzr = pm.dz * resS
dhr = pm.dh * resS

pm = merge(pm, (; dzr=dzr, dhr=dhr, Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2), nx = round(Int, pm.Lx/dhr),
                m = π/pm.Lz,
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

filesetnames = "SetList_mp.jld2"
scale_file = jldopen(filesetnames, "r+")
sns = scale_file["setnames"]
setnames = sns[1:22]
sfiles = scale_file["setfilenames"]

Lvals = length(3:2:9)
nut_avg = zeros(tlength, Lvals)

# choosing representative value
for (idx, m) in enumerate(3:2:9)
    @info "On $m.."
    setname = setnames[m]
    name_prefix = sfiles[m] * setname
    filepath = path_name * name_prefix * ".jld2"

    nu_timeseries = FieldTimeSeries(filepath, "νₑ");
    xn, yn, zn = nodes(nu_timeseries) #CCC

    pm2 = getproperty(SimParams(), Symbol(setname))

    uL = pm2.U₀ * pm2.Lzˢ

    tlength = length(nu_timeseries.times)
    #nut_std = zeros(tlength)
    #nut_count = zeros(tlength)

    threshold = 1e-4
    for i = 1:tlength
        nut = interior(nu_timeseries[i])
        nut_bool = (nut .> threshold)
        if sum(nut_bool) > 0
            nut_avg[i,idx] = mean(nut[nut_bool])
            #nut_std[i] = std(nut[nut_bool])
            #nut_count[i] = length(nut[nut_bool])
        end
    end

    #errorbars = 2 .* nut_std ./ sqrt.(nut_count)
    waves = nu_timeseries.times./ pm.Tσ
end

uL = pm.Lzˢ * .05 .* (3:2:9)
δ =  string.(round.(.05 .* (3:2:9) / pm.Ñ))

f = Figure(resolution = (800, 700), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1],  ylabel = rich("ν", subscript("t")), 
    xlabel = "Tσ" )
    # we want a log y axis with ticks at -8 through 0
    lines!(ax1, waves, ones(tlength) .* 1.05e-6, linewidth = 3, color = :black, label = "νₘ")
    for j = 1:4
        lines!(ax1, waves, nut_avg[:,j], label = δ[j], linewidth = 3)
    end
    #errorbars!(ax1, waves[3:end], nut_avg[3:end], errorbars[3:end], color = :red) 
    axislegend(ax1, label = vcat(["νₘ"], δ))
    savename = apath * "TurbulentViscosity"
save(savename * ".png", f) 

f = Figure(resolution = (800, 700), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1],  ylabel = rich("Re", subscript("t")), 
    xlabel = "Tσ", yscale = log10)
    # we want a log y axis with ticks at -8 through 0
    lines!(ax1, waves, ones(tlength) .* uL[2]/1.05e-6, linewidth = 3, color = :blue, label = "Reₘ, δ = 71")
    for j = 1:4
        lines!(ax1, waves[3:end], uL[j] ./ nut_avg[3:end,j], label = δ[j], linewidth = 3)
    end
    axislegend(ax1, label = vcat(["Reₘ"], δ))
    savename = apath * "ReynoldsNumber"
save(savename * ".png", f) 



