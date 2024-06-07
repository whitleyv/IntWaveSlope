using Statistics
using Printf
using Oceananigans
using Measures
using JLD2

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
apath = path_name * "Analysis/"
ENV["GKSwstype"] = "nul" # if on remote HPC

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

    return c_weighted_bavg
end

setnames = ["U150N100Lz100g100", "U250N100Lz100g100", "nExD_U300N100Lz100g100",]
Lvals = length(setnames)
∂ₜb̄ = zeros(Lvals, 161)

for (m, setname) in enumerate(setnames)
    name_mp_prefix = "IntWave_mp_" * setname
    filepath_mp = path_name * name_mp_prefix * ".jld2"

    Cg_timeseries = FieldTimeSeries(filepath_mp,"Cg");
    SGS_timeseries = FieldTimeSeries(filepath_mp, "SGS∇κ∇b");

    ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)

    Cg_mpi = interior(Cg_timeseries)[:,1:ylength,:,:];
    SGS_mpi = interior(SGS_timeseries)[:,1:ylength,:,:];

    @info "Tracer Weighted Calculations.."

    SGS_bar = calculate_tracer_weight(Cg_mpi, SGS_mpi)

    ∂ₜb̄[m, :] = 2 .* SGS_bar

end

filescalename = apath * "cwtdSGS.jld2"

jldsave(filescalename; 
    ∂ₜb̄, setnames,
)

dpath = "Data/"
apath = "Analysis/"
Lvals = 3

filescalename1 = dpath * "cwtdSGS.jld2"
scale_file1 = jldopen(filescalename1, "r+")

∂ₜb̄ = scale_file1["∂ₜb̄"]

∫∂ₜb̄ = zeros(Lvals, 161)
for i = 1:161
    ∫∂ₜb̄[:, i] = sum(∂ₜb̄[:,1:i], dims = 2) .* 600
end

∂ₜb̄_rWavg = zeros(Lvals, 147)
for i = 8:154
    ∂ₜb̄_rWavg[:, i-7] = mean(∂ₜb̄[:,i-7:i+7], dims = 2)
end

filescalename = dpath * "cwtdb1mom.jld2"
scale_file = jldopen(filescalename, "r+")

c_weighted_bavg_mp = hcat(hcat(scale_file["c_weighted_bavg_mp"][2,:], scale_file["c_weighted_bavg_mp"][4,:]),scale_file["c_weighted_bavg_mp"][7,:]) 


c_weighted_bavg_mp_rWavg = zeros(Lvals, 133)
for i = 15:147
    c_weighted_bavg_mp_rWavg[:, i-14] = mean(c_weighted_bavg_mp[i-14:i+14,:], dims = 1)
end

dt_c_weighted_bavg_mp = (c_weighted_bavg_mp[2:end,:] .- c_weighted_bavg_mp[1:end-1,:]) ./ 600.0
dt_c_weighted_bavg_mp_rWavgpost = (c_weighted_bavg_mp_rWavg[:, 2:end] .-  c_weighted_bavg_mp_rWavg[:, 1:end-1]) ./ 600.0
dt_c_weighted_bavg_mp_rWavg = zeros(Lvals, 132)
for i = 15:146
    dt_c_weighted_bavg_mp_rWavg[:, i-14] = mean(dt_c_weighted_bavg_mp[i-14:i+14,:], dims = 1)
end

using CairoMakie

include("parameters.jl")
sn = "U250N100Lz100g100"
pm = getproperty(SimParams(), Symbol(sn))

wave_times = (0:600:(160*600))./pm.Tσˢ

f = Figure(resolution = (2000, 1600), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "∂ₜb̄ = ⟨2∇⋅κ∇b⟩ [ms⁻³]")
ax2 = Axis(ga[2, 1], ylabel = "∂ₜb̄ = ⟨2∇⋅κ∇b⟩ [ms⁻³], Wavg", xlabel = "Wave Periods [Tσ]")
ax3 = Axis(ga[1, 2], ylabel = "b̄ Wavg [ms⁻²]")
ax4 = Axis(ga[2, 2], ylabel = "∂ₜb̄ Wavg [ms⁻³]", xlabel = "Wave Periods [Tσ]")

limits!(ax1, 0, 11, -2.8e-8, 2.8e-8)
limits!(ax2, 0, 11, -5.5e-9, 2e-9)
limits!(ax3, 0, 11, -2e-4, 1.2e-4)
limits!(ax4, 0, 11, -5.5e-9, 2e-9)

ax2.xticks = 0:2:10
ax4.xticks = 0:2:10

hidexdecorations!(ax1, grid = false)
hidexdecorations!(ax3, grid = false)

dbdtp = lines!(ax1, wave_times, ∂ₜb̄[1,:], color = :dodgerblue1, label = "42 m", linewidth = 5)
    lines!(ax1, wave_times, ∂ₜb̄[2,:], color = :dodgerblue3, label = "71 m", linewidth = 5)
    lines!(ax1, wave_times, ∂ₜb̄[3,:], color = :dodgerblue4, label = "85 m", linewidth = 5)
 
dbdtp = lines!(ax2, wave_times[8:154], ∂ₜb̄_rWavg[1,:], color = :dodgerblue1, label = "42 m", linewidth = 5)
    lines!(ax2, wave_times[8:154], ∂ₜb̄_rWavg[2,:], color = :dodgerblue3, label = "71 m", linewidth = 5)
    lines!(ax2, wave_times[8:154], ∂ₜb̄_rWavg[3,:], color = :dodgerblue4, label = "85 m", linewidth = 5)

lines!(ax3, wave_times[15:147], (c_weighted_bavg_mp_rWavg[1,:] .- c_weighted_bavg_mp_rWavg[1,1]), color = :dodgerblue1, linewidth = 5)
    lines!(ax3, wave_times[15:147], (c_weighted_bavg_mp_rWavg[2,:] .- c_weighted_bavg_mp_rWavg[2,1]), color = :dodgerblue3, linewidth = 5)
    lines!(ax3, wave_times[15:147], (c_weighted_bavg_mp_rWavg[3, :] .- c_weighted_bavg_mp_rWavg[3,1]), color = :dodgerblue4, linewidth = 5)

lines!(ax4, wave_times[15:146], dt_c_weighted_bavg_mp_rWavg[1,:], color = :dodgerblue1, linewidth = 5)
    lines!(ax4, wave_times[15:146], dt_c_weighted_bavg_mp_rWavg[2,:], color = :dodgerblue3, linewidth = 5)
    lines!(ax4, wave_times[15:146], dt_c_weighted_bavg_mp_rWavg[3,:], color = :dodgerblue4, linewidth = 5)

axislegend(ax1, position = :rt)

savename = "cwtd_sgs" 

save(apath * savename * ".png", f)


f = Figure(resolution = (1600, 800), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "Δb̄ = ∫⟨2∇⋅κ∇b⟩dt [ms⁻²]", xlabel = "Wave Periods [Tσ]")
ax3 = Axis(ga[1, 2], ylabel = "Δb̄ [ms⁻²]", xlabel = "Wave Periods [Tσ]")

limits!(ax1, 0, 11,  -2e-4, 1.2e-4)
limits!(ax3, 0, 11, -2e-4, 1.2e-4)

ax1.xticks = 0:2:10
ax3.xticks = 0:2:10

dbdtp = lines!(ax1, wave_times, ∫∂ₜb̄[1,:], color = :dodgerblue1, label = "42 m", linewidth = 5)
    lines!(ax1, wave_times, ∫∂ₜb̄[2,:], color = :dodgerblue3, label = "71 m", linewidth = 5)
    lines!(ax1, wave_times, ∫∂ₜb̄[3,:], color = :dodgerblue4, label = "85 m", linewidth = 5)
 
lines!(ax3, wave_times, (c_weighted_bavg_mp[:,1] .- c_weighted_bavg_mp[1,1]), color = :dodgerblue1, linewidth = 5)
    lines!(ax3, wave_times, (c_weighted_bavg_mp[:,2] .- c_weighted_bavg_mp[1,2]), color = :dodgerblue3, linewidth = 5)
    lines!(ax3, wave_times, (c_weighted_bavg_mp[:,3] .- c_weighted_bavg_mp[1,3]), color = :dodgerblue4, linewidth = 5)

axislegend(ax1, position = :rt)

savename = "cwtd_sgs_int" 

save(apath * savename * ".png", f)
