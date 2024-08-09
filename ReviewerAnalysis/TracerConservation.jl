using Statistics
using Printf
using JLD2
using CairoMakie
using Oceananigans

path_name =  "/glade/derecho/scratch/whitleyv/FROM_CHEYENNE/NewAdvection/Parameters/VaryU03C/wPE/"
ENV["GKSwstype"] = "nul" # if on remote HPC

apath = path_name * "Analysis/"
sn = "U250Nfd250Lz100g100"

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ

Sp_extra = ifelse(z_start>0, 250.0, 0.0)
Ly = pm.Lyˢ+Sp_extra
ny = round(Int,Ly/4)
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
sfiles = scale_file["setfilenames"]

setnames = sns[1:3:11]
setfiles = sfiles[1:3:11]

f = Figure(resolution = (800, 700), fontsize=26)
    ga = f[1, 1] = GridLayout()

    ax1 = Axis(ga[1, 1],  ylabel = "% Change in Total [%]", 
    xlabel = "t [hrs]", title = "Change in Total Tracer Concentration")

for (m, setname) in enumerate(setnames)
    
    name_prefix = setfiles[m] * setname
    filepath = path_name * name_prefix * ".jld2"
    @info "getting data from: " * setname

    c_timeseries = FieldTimeSeries(filepath, "Cs");
    xc, yc, zc = nodes(c_timeseries) #CCC

    ci = interior(c_timeseries)
    csum = sum(ci, dims = (1,2,3))[1,1,1,:]


    cerror = csum .- csum[1]
    percent_error = 100 .* (cerror ./ csum[1])

    pm2 = getproperty(SimParams(), Symbol(setname))

    δr = round(pm2.U₀ ./ pm2.Ñ)
    # we want a log y axis with ticks at -8 through 0
    lines!(ax1, c_timeseries.times ./ 3600, percent_error, linewidth = 3, label = rich("h", subscript("w"), " = $δr"))

end

    axislegend(ax1)

    savename = apath * "TracerConservation"

save(savename * ".png", f) 
