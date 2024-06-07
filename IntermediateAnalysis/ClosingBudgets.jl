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

function deriv(A, h)
    step = 1/h
    step2 = 1/(h*2)
    da = zeros(size(A));
    da[:,:,:,1] = (A[:,:,:,2] .- A[:,:,:,1]) .* step
    da[:,:,:,end] = (A[:,:,:,end] .- A[:,:,:,end-1]) .* step
    da[:,:,:,2:end-1] = (A[:,:,:,3:end] .- A[:,:,:,1:end-2]) .* step2
    return da
end

#########
#   CLOSING TRACER BUDGET
#########
#∂t(c) + u⋅∇c = ∇⋅κ∇c

Cg_timeseries = FieldTimeSeries(filepath1,"Cg");
∇κ∇Cg_timeseries = FieldTimeSeries(filepath2,"∇κ∇Cg");
u∇Cg_timeseries = FieldTimeSeries(filepath2,"u∇c");

CGi = interior(Cg_timeseries);
∇κ∇Cgi = interior(∇κ∇Cg_timeseries);
u∇Cgi = interior(u∇Cg_timeseries);

#∫ ∂C/∂t dV = 0
#∂t = Cg_timeseries.times[2:end] .-  Cg_timeseries.times[1:end-1];
#∂c = CGi[:,:,:,2:end] .- CGi[:,:,:,1:end-1];
#∫∂cdV = sum(∂c, dims = (1,2,3))[1,1,1,:]
#∫∂c∂tdV = ∫∂cdV ./ ∂t;

∂c∂t = deriv(CGi, 600.0);

∫cdV = sum(CGi, dims = (1,2,3))[1,1,1,:];

#∂∫cdV∂t = (∫cdV[2:end] .- ∫cdV[1:end-1]) ./ ∂t;
∫u∇CgdV = sum(u∇Cgi, dims = (1,2,3))[1,1,1,:]
∫∇κ∇CgdV = sum(∇κ∇Cgi, dims = (1,2,3))[1,1,1,:]
∫∂c∂tdV = sum(∂c∂t, dims = (1,2,3))[1,1,1,:]

### Plot both volume integrated LHS

##########
# CLOSING the BUOYANCY BUDGET
#########

b_timeseries = FieldTimeSeries(filepath1,"b");
xb, yb, zb = nodes(b_timeseries) #CCC

∇κ∇B_timeseries = FieldTimeSeries(filepath2,"∇κ∇B");
B_timeseries = FieldTimeSeries(filepath2,"b");
u∇B_timeseries = FieldTimeSeries(filepath2,"u∇b");

bi = interior(b_timeseries);
∇κ∇Bi = interior(∇κ∇B_timeseries);
u∇Bi = interior(u∇B_timeseries);
Bi = interior(B_timeseries);

const Sp_Region_right = 500                               # size of sponge region on RHS

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
@inline mask2nd(X) = heaviside(X)* X^2
@inline right_mask(y) = mask2nd((y-pm.Lyˢ+Sp_Region_right)/(Sp_Region_right))

xlength = length(xb)
ylength = length(yb)
zlength = length(zb)

# @inline b_sponge(x, y, z, b) =   0.001 * right_mask(y) * (pm.Ñ^2 * zb - b)

YZXgrid = reshape(repeat(0.001 .* right_mask.(yb), xlength*zlength), ylength, zlength, xlength)
ZXYgrid = permutedims(YZXgrid, [2,3,1]);
XYZgrid = permutedims(YZXgrid, [3,1,2]);
ZXYgrid_strat = pm.Ñ^2 .* zb .* ZXYgrid;

strat_sponge = sum(ZXYgrid_strat);
XYZgrid_b = XYZgrid .* Bi;
∫XYZgrid_bdV = sum(XYZgrid_b, dims = (1,2,3))[1,1,1,:]
∫_bSponge_dV = strat_sponge .- ∫XYZgrid_bdV;

# ∫∇ ⋅ κ ∇b dV 
∫∇κ∇BdV = sum(∇κ∇Bi, dims = (1,2,3))[1,1,1,:]

# ∂b/∂t
∂b∂t = deriv(bi, 600.0);

#∂b = bi[:,:,:,2:end] .- bi[:,:,:,1:end-1];
∫∂b∂tdV = sum(∂b∂t, dims = (1,2,3))[1,1,1,:];

∫u∇BdV = sum(u∇Bi, dims = (1,2,3))[1,1,1,:];
∫bdV = sum(bi, dims = (1,2,3));
#∂∫bdV∂t = (∫bdV[2:end] .- ∫bdV[1:end-1]) ./ ∂t;
RHS_b = ∫∇κ∇BdV .+ ∫_bSponge_dV

wave_times = b_timeseries.times./pm.Tσ

f = Figure(resolution = (1600, 2000), fontsize=26)
    ga = f[1, 1] = GridLayout() 
    ax1 = Axis(ga[1, 1], ylabel = "∫∂c∂tdV", xlabel = "Wave Periods [Tσ]")
    ax3 = Axis(ga[1, 2], ylabel = "∫∂b∂tdV", xlabel = "Wave Periods [Tσ]")
    ax2 = Axis(ga[2, 1], ylabel = "c RHS", xlabel = "Wave Periods [Tσ]")
    ax4 = Axis(ga[2, 2], ylabel = "b RHS", xlabel = "Wave Periods [Tσ]")
    ax5 = Axis(ga[3, 1], ylabel = "c: LHS - RHS", xlabel = "Wave Periods [Tσ]")
    ax6 = Axis(ga[3, 2], ylabel = "b: LHS - RHS", xlabel = "Wave Periods [Tσ]")

    #limits!(ax1, 0, 11, -1.5e-4, 1.5e-4)
    #limits!(ax3, 0, 11, -1.5e-4, 1.5e-4)

    #ax1.xticks = 0:2:10
    #ax3.xticks = 0:2:10

    pc = lines!(ax1, wave_times, ∫∂c∂tdV, color = :firebrick2, label = "∂c/∂t", linewidth = 5)

    pb = lines!(ax3, wave_times[2:end], (∫∂b∂tdV .+ ∫u∇BdV)[2:end], color = :firebrick2, label = "∂b/∂t + u⋅∇b", linewidth = 5)
    lines!(ax3, wave_times[2:end], RHS_b[2:end], color = :dodgerblue2, label = "∇⋅κ∇B + Sponge", linewidth = 5, linestyle = :dash)

    pc2 = lines!(ax2, wave_times, ∫∇κ∇CgdV, color = :dodgerblue2, label = "∇⋅κ∇c ", linewidth = 5)
    pc = lines!(ax2, wave_times, -1 .* ∫u∇CgdV, color = :firebrick4, label = "-u⋅∇c", linewidth = 5)

    #lines!(ax4, wave_times[2:end], (∫∂b∂tdV .- ∫_bSponge_dV)[2:end] .* 32, color = :dodgerblue2, label = "∂b/∂t - Sponge", linewidth = 5)
        lines!(ax4, wave_times, ∫∇κ∇BdV, color = :firebrick2, label = "∇⋅κ∇B ", linewidth = 5)
        lines!(ax4, wave_times, ∫u∇BdV, color = :firebrick4, label = "-u⋅∇b", linewidth = 5)

    pc3 = lines!(ax5, wave_times, ∫∂c∂tdV .+ ∫u∇CgdV .- ∫∇κ∇CgdV, color = :firebrick2, label = "e|LHS - RHS| = 0.19", linewidth = 5)
    
    pb4 = lines!(ax6, wave_times[2:end], (∫∂b∂tdV .- ∫_bSponge_dV)[2:end], color = :firebrick4, label = "e|∂b/∂t - Sponge| = 1.54", linewidth = 5)
    pb4 = lines!(ax6, wave_times[2:end], (∫u∇BdV .- ∫∇κ∇BdV)[2:end], color = :firebrick2, label = "e|u⋅∇b - ∇⋅κ∇c| = 9.0 × 10⁻⁵", linewidth = 5)
    pb4 = lines!(ax6, wave_times[2:end], (∫∂b∂tdV .- ∫_bSponge_dV .-  ∫∇κ∇BdV .+ ∫u∇BdV)[2:end], color = :dodgerblue2, label = "e|LHS - RHS| = 1.54", linewidth = 5)

    axislegend(ax1, position = :rb)
    axislegend(ax3, position = :rb)
    axislegend(ax2, position = :rb)
    axislegend(ax4, position = :rb)
    axislegend(ax5, position = :rb)
    axislegend(ax6, position = :rb)

savename = "cbClosures_" * setname 

save(apath * savename * ".png", f)

#########
#   CLOSING TRACER-BUOYANCY BUDGET
#########

Cg_timeseries = FieldTimeSeries(filepath1,"Cg");
b_timeseries = FieldTimeSeries(filepath1,"b");

Cg∇κ∇B_timeseries = FieldTimeSeries(filepath2,"Cg∇κ∇B");
B∇κ∇Cg_timeseries = FieldTimeSeries(filepath2,"B∇κ∇Cg");
u∇cb_timeseries = FieldTimeSeries(filepath2,"u∇cb");
B_timeseries = FieldTimeSeries(filepath2,"b");
CG_timeseries = FieldTimeSeries(filepath2,"Cg");

CGavgi = interior(CG_timeseries);
Bi = interior(B_timeseries);
CGi = interior(Cg_timeseries);
bi = interior(b_timeseries);
Cg∇κ∇Bi = interior(Cg∇κ∇B_timeseries);
B∇κ∇Cgi = interior(B∇κ∇Cg_timeseries);
u∇CgBi = interior(u∇cb_timeseries);

#∂t = Cg_timeseries.times[2:end] .-  Cg_timeseries.times[1:end-1];
cb = CGi .* bi;

∂cb∂t = deriv(cb, 600.0);
#∂cb = (cb[:,:,:,2:end] .- cb[:,:,:,1:end-1]); 

xb, yb, zb = nodes(b_timeseries) #CCC

const Sp_Region_right = 500                               # size of sponge region on RHS

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
@inline mask2nd(X) = heaviside(X)* X^2
@inline right_mask(y) = mask2nd((y-pm.Lyˢ+Sp_Region_right)/(Sp_Region_right))

xlength = length(xb)
ylength = length(yb)
zlength = length(zb)

# @inline b_sponge(x, y, z, b) =   0.001 * right_mask(y) * (pm.Ñ^2 * zb - b)

YZXgrid = reshape(repeat(0.001 .* right_mask.(yb), xlength*zlength), ylength, zlength, xlength)
ZXYgrid = permutedims(YZXgrid, [2,3,1]);
XYZgrid = permutedims(YZXgrid, [3,1,2]);
ZXYgrid_strat = pm.Ñ^2 .* zb .* ZXYgrid;
XYZgrid_strat = permutedims(ZXYgrid_strat, [2,3,1]);

XYZgrid_b = XYZgrid .* Bi;

bsponge = XYZgrid_strat .- XYZgrid_b;
cbsponge = CGavgi .* bsponge

∫_cbSponge_dV = sum(cbsponge, dims = (1,2,3))[1,1,1,:] .* 32

#∫∂cb∂tdV = (sum(∂cb∂t, dims = (1,2,3))[1,1,1,:]) ./ ∂t .* 32; 
∫∂cb∂tdV = sum(∂cb∂t, dims = (1,2,3))[1,1,1,:] .* 32; 
∫Cg∇κ∇BdV = sum(Cg∇κ∇Bi, dims = (1,2,3))[1,1,1,:] .* 32
∫B∇κ∇CgdV = sum(B∇κ∇Cgi, dims = (1,2,3))[1,1,1,:] .* 32
∫u∇CgBdV = sum(u∇CgBi, dims = (1,2,3))[1,1,1,:] .* 32

RHS_cb = ∫Cg∇κ∇BdV .+ ∫B∇κ∇CgdV .+ ∫_cbSponge_dV;
LHS_cb = ∫∂cb∂tdV .+ ∫u∇CgBdV

wave_times = b_timeseries.times./pm.Tσ

f = Figure(resolution = (2100, 1600), fontsize=26)
    ga = f[1, 1] = GridLayout() 
    ax1 = Axis(ga[1, 1], ylabel = "∫∂cb∂tdV", xlabel = "Wave Periods [Tσ]")
    ax3 = Axis(ga[1, 2], ylabel = "∫∂cb∂tdV", xlabel = "Wave Periods [Tσ]")
    ax2 = Axis(ga[2, 1], ylabel = "cb RHS Terms", xlabel = "Wave Periods [Tσ]")
    ax4 = Axis(ga[2, 2], ylabel = "cb Absolute Closure Error", xlabel = "Wave Periods [Tσ]")

    limits!(ax1, 0, 11, -1e-2, 1e-2)
    limits!(ax3, 0, 11, -1e-2, 1e-2)
    limits!(ax2, 0, 11, -5e-3, 5e-3)
    limits!(ax4, 0, 11, 0, 1e-2)

    #ax1.xticks = 0:2:10
    #ax3.xticks = 0:2:10
    pc = lines!(ax1, wave_times, ∫∂cb∂tdV, color = :firebrick2, label = "∂cb/∂t", linewidth = 5)
    lines!(ax1, wave_times, ∫u∇CgBdV, color = :firebrick4, label = "u⋅∇(cb)", linewidth = 5)

    pb = lines!(ax3, wave_times, LHS_cb, color = :firebrick2, label = "∂cb/∂t + u⋅∇(cb)", linewidth = 5)
    lines!(ax3, wave_times, RHS_cb, color = :dodgerblue4, label = "B∇⋅κ∇C + C∇⋅κ∇B + CSponge", linewidth = 5)

    pb2 =  lines!(ax2, wave_times, RHS_cb, color = :dodgerblue4, label = "B∇⋅κ∇C + C∇⋅κ∇B + CSponge", linewidth = 5,)
        lines!(ax2, wave_times, ∫Cg∇κ∇BdV, color = :firebrick2, label = "C∇⋅κ∇B", linewidth = 5)
        lines!(ax2,wave_times, ∫B∇κ∇CgdV, color = :green, label = "B∇⋅κ∇C ", linewidth = 5, linestyle = :dash)
    
    pb4 = lines!(ax4, wave_times, abs.(LHS_cb .- RHS_cb) , color = :firebrick2, label =  "|LHS- RHS|", linewidth = 5)
    pb4 = lines!(ax4, wave_times, abs.(∫∂cb∂tdV .- 2*∫Cg∇κ∇BdV) , color = :dodgerblue4, label =  "|∂cb/∂t- 2C∇⋅κ∇B|", linewidth = 5)

    axislegend(ax1, position = :rb)
    axislegend(ax3, position = :rb)
    axislegend(ax2, position = :rb)
    axislegend(ax4, position = :rt)

savename = "cb_mult_Closure_" * setname 

save(apath * savename * ".png", f)

