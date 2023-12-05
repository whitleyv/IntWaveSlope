using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

apath = path_name * "Analysis/"

include("parameters.jl")

setname = "U250N100Lz100g100"

# how many times were saved?
ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
m = -π/pm.Lz,
l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
Tf = 2*π/pm.f, 
Tσ = 2*π/pm.σ))

zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame
    
@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

name1_prefix = "IntWave_mp_noS_" * setname
name2_prefix = "IntWave_mp_noS_" * setname * "_budge"
filepath1 = path_name * name1_prefix * ".jld2"
filepath2 = path_name * name2_prefix * ".jld2"

@info "getting data from: " * setname

b_timeseries = FieldTimeSeries(filepath2,"b");
xb, yb, zb = nodes(b_timeseries) #CCC
Cg_timeseries = FieldTimeSeries(filepath2,"Cg");
∇κ∇Cg_timeseries = FieldTimeSeries(filepath2,"∇κ∇Cg");
∇κ∇B_timeseries = FieldTimeSeries(filepath2,"∇κ∇B");

ylength = 1000
CGi = interior(Cg_timeseries)[:,1:ylength,:,:];
Bi = interior(b_timeseries)[:,1:ylength,:,:];
∇κ∇Cgi = interior(∇κ∇Cg_timeseries)[:,1:ylength,:,:];
∇κ∇Bi = interior(∇κ∇B_timeseries)[:,1:ylength,:,:];

##################
# OFFLINE CALCULATION OF RHS Terms
##################

function GradFieldBoolean(xb, yb, zb, ylength, curvedslope)
    zlength = length(zb)
    xlength = length(xb)

    Ygrid = reshape(repeat(yb[3:ylength], (xlength-2)*(zlength-2)), ylength-2, zlength-2, xlength-2)
    SlopeGridY = curvedslope.(Ygrid);
    boolZY = (SlopeGridY .+ 4) .<= zb[3:end]';  # all the values greater than slope

    boolZYX = permutedims(boolZY, [3, 1, 2])

    return boolZYX
end

boolZYX = GradFieldBoolean(xb, yb, zb, ylength, curvedslope)

function GradField(xb, yb, zb, CGi, κi, boolZYX)
    Δx = 1/(xb[2] - xb[1])
    Δy = 1/(yb[2] - yb[1])
    Δz = 1/(zb[2] - zb[1])

    ∂xCg = (CGi[2:end, 2:end, 2:end, :] .- CGi[1:end-1, 2:end, 2:end, :]) .* Δx;
    ∂yCg = (CGi[2:end, 2:end, 2:end, :] .- CGi[2:end, 1:end-1, 2:end, :]) .* Δy;
    ∂zCg = (CGi[2:end, 2:end, 2:end, :] .- CGi[2:end, 2:end, 1:end-1, :]) .* Δz;

    κ∂xCg = κi[2:end, 2:end, 2:end, :] .* ∂xCg;
    κ∂yCg = κi[2:end, 2:end, 2:end, :] .* ∂yCg;
    κ∂zCg = κi[2:end, 2:end, 2:end, :] .* ∂zCg;

    ∂xκ∂xCg = (κ∂xCg[2:end, 2:end, 2:end, :] .- κ∂xCg[1:end-1, 2:end, 2:end, :]) .* Δx;
    ∂yκ∂yCg= (κ∂yCg[2:end, 2:end, 2:end, :] .- κ∂yCg[2:end, 1:end-1, 2:end, :]) .* Δy;
    ∂zκ∂zCg = (κ∂zCg[2:end, 2:end, 2:end, :] .- κ∂zCg[2:end, 2:end, 1:end-1, :]) .* Δz;

    ∇κ∇Cg = (∂xκ∂xCg .+ ∂yκ∂yCg .+ ∂zκ∂zCg) .* boolZYX;

    return ∇κ∇Cg
end

∇κ∇Cg = GradField(xb, yb, zb, CGi, κi, boolZYX)
∇κ∇B = GradField(xb, yb, zb, Bi, κi, boolZYX)

# if RHS terms offline, need to adjust the others to match the size and volume

# ∇ ⋅ κ∇b, ∇ ⋅ κ∇c
b∇κ∇Cg = ∇κ∇Cg .* Bi[3:end, 3:end, 3:end, :];
c∇κ∇B = ∇κ∇B .* CGi[3:end, 3:end, 3:end, :];
cb = (CGi .* Bi)[3:end, 3:end, 3:end, :] .* boolZYX;
csum = sum(CGi[3:end, 3:end, 3:end, :] .* boolZYX, dims=(1,2,3))[1,1,1,:];

##################
# ONLINE CALCULATION of RHS
##################

b∇κ∇Cg = Bi .* ∇κ∇Cgi;
c∇κ∇B = CGi .* ∇κ∇Bi;

cb = CGi .* Bi;

##################
# VOLUME INTEGRALS
##################

# Volume Integrating
∫b∇κ∇Cg = sum(b∇κ∇Cg, dims=(1,2,3))[1,1,1,:];
∫c∇κ∇B = sum(c∇κ∇B, dims=(1,2,3))[1,1,1,:];
∫cb = sum(cb, dims=(1,2,3))[1,1,1,:];

csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

# dividing by dye volume integration

∫b∇κ∇Cg_∫c = ∫b∇κ∇Cg ./ csum;
∫c∇κ∇B_∫c = ∫c∇κ∇B ./ csum;
∫cb_∫c = ∫cb ./ csum;
∫RHS_∫c = ∫b∇κ∇Cg_∫c .+ ∫cSGS_∫c;

##################
# TIME INTEGRALS of RHS (to match cb without derivs)
##################

tint_∫c∇κ∇B_∫c = zeros(161)
tint_∫b∇κ∇Cg_∫c = zeros(161)
tint_∫RHS_∫c = zeros(161)

for i = 1:161
    tint_∫b∇κ∇Cg_∫c[i] = sum(∫b∇κ∇Cg_∫c[1:i]) .* 600
    tint_∫RHS_∫c[i] = sum(∫RHS_∫c[1:i]) .* 600
    tint_∫c∇κ∇B_∫c[i] = sum(∫c∇κ∇B_∫c[2:i]) .* 600
end

filescalename = apath * "cwtdrhs_bothterms_" * setname *  ".jld2"
filescalename = apath * "cwtdrhs_bothterms_online_noS_" * setname *  ".jld2"

∂t_∫cb_∫c = (∫cb_∫c[2:end] .- ∫cb_∫c[1:end-1] )./ 600

jldsave(filescalename; 
∫c∇κ∇B_∫c, ∫b∇κ∇Cg_∫c, ∫cb_∫c, ∫RHS_∫c, ∂t_∫cb_∫c,
tint_∫b∇κ∇Cg_∫c, tint_∫RHS_∫c, tint_∫c∇κ∇B_∫c, 
)

##################
# PLOTTING
##################

wave_times = (0:600:(160*600))./pm.Tσˢ

# Plotting the terms in time integral form:

f = Figure(resolution = (1600, 800), fontsize=26)
    ga = f[1, 1] = GridLayout() 
    ax1 = Axis(ga[1, 1], ylabel = "∫cbdV / ∫cdV [ms⁻²]", xlabel = "Wave Periods [Tσ]")
    ax3 = Axis(ga[1, 2], ylabel = "Inidividual RHS Terms", xlabel = "Wave Periods [Tσ]")

    limits!(ax1, 0, 11, -1.5e-4, 1.5e-4)
    limits!(ax3, 0, 11, -1.5e-4, 1.5e-4)

    ax1.xticks = 0:2:10
    ax3.xticks = 0:2:10

    dbdtp = lines!(ax1, wave_times[2:end], ∫cb_∫c[2:end] .- ∫cb_∫c[2], color = :firebrick2, label = "∫cb dV/ ∫cdV", linewidth = 5)
        lines!(ax1, wave_times, tint_∫RHS_∫c .- tint_∫RHS_∫c[1], color = :dodgerblue2, label = "∫(∫c∇⋅κ∇b + b∇⋅κ∇cdV / ∫cdV) dt ", linewidth = 5)
        lines!(ax1, wave_times, 2 .* (tint_∫c∇κ∇B_∫c .- tint_∫c∇κ∇B_∫c[1]), color = :dodgerblue4, label = "∫(2∫c∇⋅κ∇b / ∫cdV) dt ", linewidth = 5)
    
    lines!(ax3, wave_times, tint_∫c∇κ∇B_∫c .- tint_∫c∇κ∇B_∫c[1], color = :dodgerblue1, linewidth = 5, label = "∫(∫c∇⋅κ∇bdV / ∫cdV)dt ")
        lines!(ax3, wave_times, tint_∫b∇κ∇Cg_∫c .- tint_∫b∇κ∇Cg_∫c[1], color = :dodgerblue4, linewidth = 5, label = "∫(∫b∇⋅κ∇cdV / ∫cdV)dt ")

    axislegend(ax1, position = :rt)
    axislegend(ax3, position = :rt)

savename = "cwtd_sumterm_intLHS_noS_" * setname 

save(apath * savename * ".png", f)

# Plotting the terms in time derivative form:

f = Figure(resolution = (1600, 800), fontsize=26)
    ga = f[1, 1] = GridLayout() 
    ax1 = Axis(ga[1, 1], ylabel = "∂t (∫cbdV / ∫cdV )[ms⁻²]", xlabel = "Wave Periods [Tσ]")
    ax3 = Axis(ga[1, 2], ylabel = "Inidividual RHS Terms", xlabel = "Wave Periods [Tσ]")

    limits!(ax1, 0, 11, -1.5e-8, 3.0e-8)
    limits!(ax3, 0, 11, -1.5e-8, 3.0e-8)

    ax1.xticks = 0:2:10
    ax3.xticks = 0:2:10

    dbdtp = lines!(ax1, wave_times[3:end], ∂t_∫cb_∫c[2:end], color = :firebrick2, label = "∂t [∫cb dV/ ∫cdV]", linewidth = 5)
        lines!(ax1, wave_times[2:end], ∫RHS_∫c[2:end], color = :dodgerblue2, label = "∫c∇⋅κ∇b + b∇⋅κ∇ dV / ∫cdV ", linewidth = 5)
        lines!(ax1, wave_times[2:end], 2 .* (∫c∇κ∇B_∫c[2:end] ), color = :dodgerblue4, label = "∫2c∇⋅κ∇b dV/ ∫cdV ", linewidth = 5)
    
    lines!(ax3, wave_times[2:end], ∫c∇κ∇B_∫c[2:end] , color = :dodgerblue1, linewidth = 5, label = "∫c∇⋅κ∇b dV / ∫cdV ")
        lines!(ax3, wave_times[2:end], ∫b∇κ∇Cg_∫c[2:end], color = :dodgerblue4, linewidth = 5, label = "∫b∇⋅κ∇c dV / ∫cdV ")

    axislegend(ax1, position = :rt)
    axislegend(ax3, position = :rt)

savename = "cwtd_sumterm_noS2_" * setname 

save(apath * savename * ".png", f)
