using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

apath = path_name * "Analysis/"

include("parameters.jl")

setname = "U300N100Lz100g100"

# how many times were saved?
zlength = 250
ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)
xlength = 38
tlength = 161 

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

name1_prefix = "IntWave_" * setname
name2_prefix = "IntWave_smdt_" * setname
name3_prefix = "IntWave_postspin_" * setname

filepath1 = path_name * name1_prefix * ".jld2"
filepath2 = path_name * name2_prefix * ".jld2"
filepath3 = path_name * "cgl_" * setname * ".jld2"
filepath4 = path_name * "cgl2_" * setname * ".jld2"

f3 = jldopen(filepath3)
CGli = f3["Cgli"];

f4 = jldopen(filepath4)
CGl2i = f4["CGl2i"];

@info "getting data from: " * setname

b_timeseries = FieldTimeSeries(filepath1,"b");
xb, yb, zb = nodes(b_timeseries) #CCC
# gaussian at slope but with increased time resolution
Cg_timeseries = FieldTimeSeries(filepath2,"Cg");
Cgr_timeseries = FieldTimeSeries(filepath1,"Cgr");

CGi = interior(Cg_timeseries)[:, 1:ylength, 1:zlength,1:2:end];
CGri = interior(Cgr_timeseries)[:, 1:ylength, 1:zlength,:];
Bi = interior(b_timeseries)[:, 1:ylength, 1:zlength,:];

spinlength = size(CGli)[4]
spinlength2 = size(CGl2i)[4]


#########################
#                   TRACER WEIGHTING Calculation
#########################

@info "Tracer Weighted Calculations.."

function calculate_tracer_weight(CGi, Bi)
    
    @info "Calculating Denominator of tracer wtd average..."
    # ∫ c dV
    csum = sum(CGi, dims=(1,2,3))[1,1,1,:];

    @info "Calculating first moment..."
    # ∫ cb dV
    cb = CGi .* Bi;
    cbsum = sum(cb, dims=(1,2,3))[1,1,1,:];

    # b̄ =  ∫ cb dV / ∫ c dV
    c_weighted_bavg = cbsum ./ csum;

    return c_weighted_bavg, csum
end

function calculate_cb_velocities(CGi, CGli, CGl2i, CGri, Bi, Δt)
    @info "Calculating Tracer Velocity..."
    ∂Cg_∂t = (CGi[:,:,:,2:end] .- CGi[:,:,:,1:end-1])./ Δt
    ∂Cgl_∂t = (CGli[:,:,:,2:end] .- CGli[:,:,:,1:end-1])./ Δt
    ∂Cgl2_∂t = (CGl2i[:,:,:,2:end] .- CGl2i[:,:,:,1:end-1])./ Δt
    ∂Cgr_∂t = (CGri[:,:,:,2:end] .- CGri[:,:,:,1:end-1])./ Δt

    @info "Calculating Buoynacy Velocity..."
    ∂b_∂t = (Bi[:,:,:,2:end] .- Bi[:,:,:,1:end-1])./ Δt

    return ∂Cg_∂t, ∂Cgl_∂t, ∂Cgl2_∂t, ∂Cgr_∂t, ∂b_∂t
end

function calculate_tracer_weight_velocity_terms(b̄_Cg, csum, ∂b_∂t, ∂C_∂t, CGi, Bi, Δt)
    
    @info "Calculating Tracer Wtd Velocity"
    ∂b̄_∂t =  (b̄_Cg[2:end] .- b̄_Cg[1:end-1])./ Δt

    @info "Calculating total two terms..."
    b∂C_∂t = ∂C_∂t .*Bi[:,:,:,2:end]
    C∂b_∂t = ∂b_∂t .*CGi[:,:,:,2:end]

    @info "Integrating..."
    ∫b∂C_∂t = sum(b∂C_∂t, dims = (1,2,3))[1,1,1,:]
    ∫C∂b_∂t = sum(C∂b_∂t, dims = (1,2,3))[1,1,1,:]

    ∫b∂C_∂t_csum = ∫b∂C_∂t ./ csum[2:end]
    ∫C∂b_∂t_csum = ∫C∂b_∂t ./ csum[2:end]

    return ∫b∂C_∂t_csum, ∫C∂b_∂t_csum, ∂b̄_∂t
end

Δt = 600.0

(b̄_Cg, Cg_sum) = calculate_tracer_weight(CGi, Bi)
(b̄_Cgr, Cgr_sum) = calculate_tracer_weight(CGri, Bi)
(b̄_Cgl, Cgl_sum) = calculate_tracer_weight(CGli, Bi[:,:,:,(tlength-spinlength+1):end])
(b̄_Cgl2, Cgl2_sum) = calculate_tracer_weight(CGl2i, Bi[:,:,:,(tlength-spinlength2+1):end])

(∂Cg_∂t, ∂Cgl_∂t, ∂Cgl2_∂t, ∂Cgr_∂t, ∂b_∂t) = calculate_cb_velocities(CGi, CGli, CGl2i, CGri, Bi, Δt)

(∫b∂Cg_∂t_Cg_sum, ∫Cg∂b_∂t_Cgsum, ∂b̄_∂t_Cg) = calculate_tracer_weight_velocity_terms(b̄_Cg, Cg_sum, ∂b_∂t, ∂Cg_∂t, CGi, Bi, 600.0)
(∫b∂Cgr_∂t_Cgr_sum, ∫Cgr∂b_∂t_Cgrsum, ∂b̄_∂t_Cgr) = calculate_tracer_weight_velocity_terms(b̄_Cgr, Cgr_sum, ∂b_∂t, ∂Cgr_∂t, CGri, Bi, 600.0)
(∫b∂Cgl_∂t_Cgl_sum, ∫Cgl∂b_∂t_Cglsum, ∂b̄_∂t_Cgl) = calculate_tracer_weight_velocity_terms(b̄_Cgl, Cgl_sum, ∂b_∂t[:,:,:,(tlength-spinlength+1):end], ∂Cgl_∂t, CGli, Bi[:,:,:,(tlength-spinlength+1):end], 600.0)
(∫b∂Cgl_∂t_Cgl2_sum, ∫Cgl∂b_∂t_Cgl2sum, ∂b̄_∂t_Cgl2) = calculate_tracer_weight_velocity_terms(b̄_Cgl2, Cgl2_sum, ∂b_∂t[:,:,:,(tlength-spinlength2+1):end], ∂Cgl2_∂t, CGl2i, Bi[:,:,:,(tlength-spinlength2+1):end], 600.0)

dt_times = b_timeseries.times[2:end]./pm.Tσ
dt_times_spin = b_timeseries.times[tlength-spinlength+2:end]./pm.Tσ
dt_times_spin2 = b_timeseries.times[tlength-spinlength2+2:end]./pm.Tσ

# plot statistics compared to delta values
filescalename = apath * "cWtd_Gauss_terms.jld2"

jldsave(filescalename;  
∫b∂Cg_∂t_Cg_sum, ∫Cg∂b_∂t_Cgsum, ∂b̄_∂t_Cg, 
∫b∂Cgr_∂t_Cgr_sum, ∫Cgr∂b_∂t_Cgrsum, ∂b̄_∂t_Cgr,
∫b∂Cgl_∂t_Cgl_sum, ∫Cgl∂b_∂t_Cglsum, ∂b̄_∂t_Cgl,
∫b∂Cgl_∂t_Cgl2_sum, ∫Cgl∂b_∂t_Cgl2sum, ∂b̄_∂t_Cgl2,
b̄_Cg, b̄_Cgr, b̄_Cgl, b̄_Cgl2,
dt_times, dt_times_spin, dt_times_spin2)