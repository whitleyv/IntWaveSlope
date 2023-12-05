using Statistics
using Printf
using JLD2
using CairoMakie

ENV["GKSwstype"] = "nul" # if on remote HPC

@info "Pulling Data..."

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/Analysis/"

sn1 = "U150N100Lz100g100"
sn2 = "U250N100Lz100g100"
sn3 = "U300N100Lz100g100"

timeindex = "beg"

filescalename1 = path_name * "BFluxTrip_" * timeindex * "_mp_" * sn1 * ".jld2"
filescalename2 = path_name * "BFluxTrip_" * timeindex * "_mp_" * sn2 * ".jld2"
filescalename3 = path_name * "BFluxTrip_" * timeindex * "_mp_" * sn3 * ".jld2"

SGSname2 = path_name * "BFluxTrip_SGS_mp_" * sn2 * ".jld2"

include("parameters.jl") 
sn = sn1
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                       nz = round(Int,pm.Lz/2),
                       m = -π/pm.Lz,
                       l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                       Tf = 2*π/pm.f, 
                       Tσ = 2*π/pm.σ))

scale_file1 = jldopen(filescalename1, "r+")
scale_file2 = jldopen(filescalename2, "r+")
scale_file3 = jldopen(filescalename3, "r+")
SGS_file2 = jldopen(SGSname2, "r+")                 

skeys = keys(scale_file1)

FluxDiv_turbWavg_beg1  = scale_file1[skeys[1]];
FluxDiv_phasedepWavg_beg1  = scale_file1[skeys[2]];
WaveAveragedVals_beg1  = scale_file1[skeys[5]];

yb = scale_file1[skeys[7]];
zb = scale_file1[skeys[8]];
phase_times = scale_file1[skeys[9]];

FluxDiv_turbWavg_beg2  = scale_file2[skeys[1]];
FluxDiv_phasedepWavg_beg2  = scale_file2[skeys[2]];
WaveAveragedVals_beg2  = scale_file2[skeys[5]];
SGS_Vals_beg2 = SGS_file2["WaveAveragedVals_" * timeindex]

FluxDiv_turbWavg_beg3  = scale_file3[skeys[1]];
FluxDiv_phasedepWavg_beg3  = scale_file3[skeys[2]];
WaveAveragedVals_beg3  = scale_file3[skeys[5]];

TurbVals_beg1 = scale_file1[skeys[4]];
TurbVals_beg2 = scale_file2[skeys[4]];
TurbVals_beg3 = scale_file3[skeys[4]];

PhaseVals_beg1 = scale_file1[skeys[3]];
PhaseVals_beg2 = scale_file2[skeys[3]];
PhaseVals_beg3 = scale_file3[skeys[3]];

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

land = curvedslope.(yb)
lin_land = linslope.(yb) # the actual z value the land hits for each value of y.
delta1 = 0.15/pm.Ñ
delta2 = 0.25/pm.Ñ
delta3 = 0.3/pm.Ñ

# c'w'
function fluxterms_dye(TurbVals_beg1, PhaseVals_beg1)
    @info "Multiplying fluxes"
    # (x,y,z,Ph,W)
    vc_turb = TurbVals_beg1.c_turb .* TurbVals_beg1.v_turb;
    wc_turb = TurbVals_beg1.c_turb .* TurbVals_beg1.w_turb;
    # (x,y,z,Ph)
    vc_phasedep = PhaseVals_beg1.c_Phdep .* PhaseVals_beg1.v_Phdep;
    wc_phasedep = PhaseVals_beg1.c_Phdep .* PhaseVals_beg1.w_Phdep;

    @info "Averaging fluxes..."
    # (x,y,z)
    vc_turbWavg = mean(vc_turb, dims = (4,5))[:,:,:,1,1];
    wc_turbWavg = mean(wc_turb, dims = (4,5))[:,:,:,1,1];
    # (x,y,z)
    vc_phasedepWavg = mean(vc_phasedep, dims = (4))[:,:,:,1];
    wc_phasedepWavg = mean(wc_phasedep, dims = (4))[:,:,:,1];

    dyeflux_values = (; vc_turbWavg, wc_turbWavg, vc_phasedepWavg, wc_phasedepWavg)

    return dyeflux_values
end

zlength = length(zb)
ylength = 880
xlength = 38

Ygrid = reshape(repeat(yb[2:ylength], xlength*(zlength-1)), ylength-1, zlength-1, xlength)
SlopeGridY = curvedslope.(Ygrid)
boolZY = (SlopeGridY .+ 4) .<= zb[2:end]';  # all the values greater than slope

boolZYX = permutedims(boolZY, [3, 1, 2])

Δy = yb[1]-yb[2]
Δz = zb[1]-zb[2]

function Flux_Divergence(Δy, Δz, boolZY, FluxVals)

    # (x,y,z)
    wc_turbWavg_dz = ((FluxVals.wc_turbWavg[:,2:end,1:end-1] .- FluxVals.wc_turbWavg[:,2:end,2:end])./Δz) .* boolZY;
    vc_turbWavg_dy = ((FluxVals.vc_turbWavg[:,1:end-1,2:end] .- FluxVals.vc_turbWavg[:,2:end,2:end])./Δy) .* boolZY;
    ∇c_turbWavg = wc_turbWavg_dz .+ vc_turbWavg_dy;
    
    FluxDiv_turbWavg = (; wc_turbWavg_dz, vc_turbWavg_dy, ∇c_turbWavg)

    wc_phasedepWavg_dz = ((FluxVals.wc_phasedepWavg[:,2:end,1:end-1] .- FluxVals.wc_phasedepWavg[:,2:end,2:end])./Δz) .* boolZY;
    vc_phasedepWavg_dy = ((FluxVals.vc_phasedepWavg[:,1:end-1,2:end] .- FluxVals.vc_phasedepWavg[:,2:end,2:end])./Δy) .* boolZY;
    ∇c_phasedepWavg = wc_phasedepWavg_dz .+ vc_phasedepWavg_dy;

    FluxDiv_phasedepWavg = (; wc_phasedepWavg_dz, vc_phasedepWavg_dy, ∇c_phasedepWavg)

    return FluxDiv_turbWavg, FluxDiv_phasedepWavg

end

dyeflux_values1 = fluxterms_dye(TurbVals_beg1, PhaseVals_beg1)
dyeflux_values2 = fluxterms_dye(TurbVals_beg2, PhaseVals_beg2)
dyeflux_values3 = fluxterms_dye(TurbVals_beg3, PhaseVals_beg3)

cFluxDiv_turbWavg1, cFluxDiv_phasedepWavg1 = Flux_Divergence(Δy, Δz, boolZYX, dyeflux_values1)
cFluxDiv_turbWavg2, cFluxDiv_phasedepWavg2 = Flux_Divergence(Δy, Δz, boolZYX, dyeflux_values2)
cFluxDiv_turbWavg3, cFluxDiv_phasedepWavg3 = Flux_Divergence(Δy, Δz, boolZYX, dyeflux_values3)

@info "Averaging in x..."
function xaverage(cFluxDiv_turbWavg1, cFluxDiv_phasedepWavg1, dyeflux_values1)
    ∇c_turbWavg_xavg = mean(cFluxDiv_turbWavg1.∇c_turbWavg, dims = 1)[1,:,:]
    ∇c_phasedepWavg_xavg = mean(cFluxDiv_phasedepWavg1.∇c_phasedepWavg, dims = 1)[1,:,:]

    vc_phasedepWavg_xavg = mean(dyeflux_values1.vc_phasedepWavg, dims = 1)[1,:,:]
    wc_phasedepWavg_xavg = mean(dyeflux_values1.wc_phasedepWavg, dims = 1)[1,:,:]

    vc_turbWavg_xavg = mean(dyeflux_values1.vc_turbWavg, dims = 1)[1,:,:]
    wc_turbWavg_xavg = mean(dyeflux_values1.wc_turbWavg, dims = 1)[1,:,:]

    xaveraged_fluxes = (; ∇c_turbWavg_xavg, ∇c_phasedepWavg_xavg, vc_phasedepWavg_xavg, wc_phasedepWavg_xavg, vc_turbWavg_xavg, wc_turbWavg_xavg)
    return xaveraged_fluxes

end

xaveraged_fluxes1 = xaverage(cFluxDiv_turbWavg1, cFluxDiv_phasedepWavg1, dyeflux_values1)
xaveraged_fluxes2 = xaverage(cFluxDiv_turbWavg2, cFluxDiv_phasedepWavg2, dyeflux_values2)
xaveraged_fluxes3 = xaverage(cFluxDiv_turbWavg3, cFluxDiv_phasedepWavg3, dyeflux_values3)

filescalename = path_name * "cFluxDeltaComp_" * timeindex * ".jld2"
jldsave(filescalename; 
xaveraged_fluxes1 = xaveraged_fluxes1,
xaveraged_fluxes2 = xaveraged_fluxes2,
xaveraged_fluxes3 = xaveraged_fluxes3,
)

∇_phasedepWavg_begxavg1 = mean(FluxDiv_phasedepWavg_beg1.∇_phasedepWavg, dims = 1)[1,:,:];
∇_turbWavg_begxavg1 = mean(FluxDiv_turbWavg_beg1.∇_turbWavg, dims = 1)[1,:,:];
c_Wavg_begxavg1 = mean(WaveAveragedVals_beg1.c_Wavg, dims = 1)[1,:,:];
b_Wavg_begxavg1 = mean(WaveAveragedVals_beg1.b_Wavg, dims = 1)[1,:,:];
SGS_Wavg_begxavg1 = mean(WaveAveragedVals_beg1.SGS_Wavg, dims = 1)[1,:,:];

∇_phasedepWavg_begxavg2 = mean(FluxDiv_phasedepWavg_beg2.∇_phasedepWavg, dims = 1)[1,:,:];
∇_turbWavg_begxavg2 = mean(FluxDiv_turbWavg_beg2.∇_turbWavg, dims = 1)[1,:,:];
c_Wavg_begxavg2 = mean(WaveAveragedVals_beg2.c_Wavg, dims = 1)[1,:,:];
b_Wavg_begxavg2 = mean(WaveAveragedVals_beg2.b_Wavg, dims = 1)[1,:,:];
SGS_Wavg_begxavg2 = mean(SGS_Vals_beg2.SGS_Wavg, dims = 1)[1,:,:];

∇_phasedepWavg_begxavg3 = mean(FluxDiv_phasedepWavg_beg3.∇_phasedepWavg, dims = 1)[1,:,:];
∇_turbWavg_begxavg3 = mean(FluxDiv_turbWavg_beg3.∇_turbWavg, dims = 1)[1,:,:];
c_Wavg_begxavg3 = mean(WaveAveragedVals_beg3.c_Wavg, dims = 1)[1,:,:];
b_Wavg_begxavg3 = mean(WaveAveragedVals_beg3.b_Wavg, dims = 1)[1,:,:];
SGS_Wavg_begxavg3 = mean(WaveAveragedVals_beg3.SGS_Wavg, dims = 1)[1,:,:];

filescalename = path_name * "FluxDeltaComp_" * timeindex * ".jld2"
    
jldsave(filescalename; 
∇_phasedepWavg_xavg = (; ∇_phasedepWavg_begxavg1, ∇_phasedepWavg_begxavg2, ∇_phasedepWavg_begxavg3),
∇_turbWavg_xavg = (; ∇_turbWavg_begxavg1, ∇_turbWavg_begxavg2, ∇_turbWavg_begxavg3),
c_Wavg_xavg = (;c_Wavg_begxavg1, c_Wavg_begxavg2, c_Wavg_begxavg3),
SGS_Wavg_xavg = (; SGS_Wavg_begxavg1, SGS_Wavg_begxavg2, SGS_Wavg_begxavg3)
)
@info "Rotating Data to match slope..."

function normal_profile(zht, lin_land, zc, bi, start_index, ylength)

    b_profile = zeros(zht)

    for hdx = 1:zht
        hab = hdx*2
        # where we want the hovmoller data to come from
        land_pz = lin_land .+ hab # where we want the hovmoller data to come from

        # for each y value, find the first z value above that number
        indexLength = length(start_index:ylength)
        depth_indices = zeros(Int, indexLength);

        for (zdx, i) in enumerate(start_index:ylength)
            depth_indices[zdx] = sum(zc .< land_pz[i]) +1
        end

        # the line at this height above bottom
        b_slopelines = zeros(indexLength)

        # at each y value
        for (zdx, j) in enumerate(start_index:ylength)
            b_slopelines[zdx] = bi[j, depth_indices[zdx]]
        end

        b_profile[hdx] = mean(b_slopelines)
    end

    return b_profile

end

st_y = 140 # 107
en_y = 212 # 250
zL = 75

c_Wavg_prof_beg1 =  normal_profile(zL, lin_land, zb, c_Wavg_begxavg1, st_y, en_y)
c_Wavg_prof_beg2 =  normal_profile(zL, lin_land, zb, c_Wavg_begxavg2, st_y, en_y)
c_Wavg_prof_beg3 =  normal_profile(zL, lin_land, zb, c_Wavg_begxavg3, st_y, en_y)

b_Wavg_prof_beg1 =  normal_profile(zL, lin_land, zb, b_Wavg_begxavg1, st_y, en_y)
b_Wavg_prof_beg2 =  normal_profile(zL, lin_land, zb, b_Wavg_begxavg2, st_y, en_y)
b_Wavg_prof_beg3 =  normal_profile(zL, lin_land, zb, b_Wavg_begxavg3, st_y, en_y)

SGS_Wavg_prof_beg1 =  normal_profile(zL, lin_land, zb, SGS_Wavg_begxavg1, st_y, en_y)
SGS_Wavg_prof_beg2 =  normal_profile(zL, lin_land, zb, SGS_Wavg_begxavg2, st_y, en_y)
SGS_Wavg_prof_beg3 =  normal_profile(zL, lin_land, zb, SGS_Wavg_begxavg3, st_y, en_y)

∇_phasedep_Wavg_prof_beg1 =  normal_profile(zL-1, lin_land, zb[2:end], ∇_phasedepWavg_begxavg1, st_y, en_y)
∇_phasedep_Wavg_prof_beg2 =  normal_profile(zL-1, lin_land, zb[2:end], ∇_phasedepWavg_begxavg2, st_y, en_y)
∇_phasedep_Wavg_prof_beg3 =  normal_profile(zL-1, lin_land, zb[2:end], ∇_phasedepWavg_begxavg3, st_y, en_y)

∇_turb_Wavg_prof_beg1 =  normal_profile(zL-1, lin_land, zb[2:end], ∇_turbWavg_begxavg1, st_y, en_y)
∇_turb_Wavg_prof_beg2 =  normal_profile(zL-1, lin_land, zb[2:end], ∇_turbWavg_begxavg2, st_y, en_y)
∇_turb_Wavg_prof_beg3 =  normal_profile(zL-1, lin_land, zb[2:end], ∇_turbWavg_begxavg3, st_y, en_y)

    ###[-∇p ][ -∇t ]
    ###[ c  ][ sgs ]

HAB = 2:2:zL*2

########## CENTER OF MASS?
cz1 = sum(HAB .* c_Wavg_prof_beg1) / sum(c_Wavg_prof_beg1)
cz2 = sum(HAB .* c_Wavg_prof_beg2) / sum(c_Wavg_prof_beg2)
cz3 = sum(HAB .* c_Wavg_prof_beg3) / sum(c_Wavg_prof_beg3)

cind1 = findfirst(HAB .>= cz1)
cind2 = findfirst(HAB .>= cz2)
cind3 = findfirst(HAB .>= cz3)

f1 = Figure(resolution = (1900, 1100), fontsize=26)
    ga = f1[1, 1] = GridLayout()    

    axtf = Axis(ga[1, 1], ylabel = "z' [m]", xlabel = "- ∇⋅⟨u⃗'b'⟩") #v
    axfs = Axis(ga[2, 1], ylabel = "z' [m]", xlabel = "- ∇⋅⟨u⃗'b'⟩ - ∇⋅⟨u⃗̃b̃⟩") #vh
    axpf = Axis(ga[1, 2], ylabel = "z' [m]", xlabel = "- ∇⋅⟨u⃗̃b̃⟩")
    axc = Axis(ga[2, 2], ylabel = "z' [m]", xlabel = "⟨c⟩") #dy v'b'
    axsgs =  Axis(ga[1, 3], ylabel = "z' [m]", xlabel = "⟨∇⋅κ∇b⟩") #dy v'b'
    axdbdt =  Axis(ga[2, 3], ylabel = "z' [m]", xlabel = "- ∇⋅⟨u⃗'b'⟩ - ∇⋅⟨u⃗̃b̃⟩ + ⟨∇⋅κ∇b⟩") #dy v'b'

    axtf.yticks = [50, 100]
    axfs.yticks = [50, 100]

    axpf.xticks = -4e-7:2e-7:0
    axtf.xticks = -1e-7:5e-8:0
    axfs.xticks = -4e-7:2e-7:0
    axc.xticks = 0:0.01:0.05
    axsgs.xticks = -8e-9:4e-9:4.1e-9
    axdbdt.xticks = -4e-7:2e-7:0

    limits!(axpf, -4.3e-7, 1e-7, 0, 150)
    limits!(axtf, -1e-7, 4e-8, 0, 150)
    limits!(axfs, -4.3e-7, 1e-7, 0, 150)
    limits!(axc, 0, 0.06, 0, 150)
    limits!(axsgs, -1.2e-8, 6e-9, 0, 150)
    limits!(axdbdt, -4.3e-7, 1e-7, 0, 150)

    hideydecorations!(axpf)
    hideydecorations!(axc)
    hideydecorations!(axsgs)
    hideydecorations!(axdbdt)

    deltaline1 = -5e-7:1e-7:-2e-7
    deltaline1L = length(deltaline1)

    lines!(axpf, -1 .* ∇_phasedep_Wavg_prof_beg1, HAB[2:end], color=:dodgerblue1, linewidth = 3, label= "δ = 42 m")
        lines!(axpf, -1 .* ∇_phasedep_Wavg_prof_beg2, HAB[2:end], color=:dodgerblue3, linewidth = 3, label = "δ = 71")
        lines!(axpf, -1 .* ∇_phasedep_Wavg_prof_beg3, HAB[2:end], color=:dodgerblue4, linewidth = 3, label = "δ = 85")
        lines!(axpf, deltaline1, delta1.*ones(deltaline1L), color = :firebrick1, linewidth = 2, linestyle = :dash)
        lines!(axpf, deltaline1, delta2.*ones(deltaline1L), color = :firebrick3, linewidth = 2, linestyle = :dash)
        lines!(axpf, deltaline1, delta3.*ones(deltaline1L), color = :firebrick4, linewidth = 2, linestyle = :dash)
        scatter!(axpf, Point((-1 .* ∇_phasedep_Wavg_prof_beg1)[cind1-1], HAB[cind1]), markersize = 10, color = :firebrick1)
        scatter!(axpf, Point((-1 .* ∇_phasedep_Wavg_prof_beg2)[cind2-1], HAB[cind2]), markersize = 10, color = :firebrick3)
        scatter!(axpf, Point((-1 .* ∇_phasedep_Wavg_prof_beg3)[cind3-1], HAB[cind3]), markersize = 10, color = :firebrick4)

    deltaline2 = -1.1e-7:1e-8:-5e-8
    deltaline2L = length(deltaline2)
    
    lines!(axtf, -1 .* ∇_turb_Wavg_prof_beg1, HAB[2:end], color=:dodgerblue1, linewidth = 3)
        lines!(axtf, -1 .* ∇_turb_Wavg_prof_beg2, HAB[2:end], color=:dodgerblue3, linewidth = 3)
        lines!(axtf, -1 .* ∇_turb_Wavg_prof_beg3, HAB[2:end], color=:dodgerblue4, linewidth = 3)
        lines!(axtf, deltaline2, delta1.*ones(deltaline2L), color = :firebrick1, linewidth = 2, linestyle = :dash)
        lines!(axtf, deltaline2, delta2.*ones(deltaline2L), color = :firebrick3, linewidth = 2, linestyle = :dash)
        lines!(axtf, deltaline2, delta3.*ones(deltaline2L), color = :firebrick4, linewidth = 2, linestyle = :dash)
     scatter!(axtf, Point((-1 .* ∇_turb_Wavg_prof_beg1)[cind1-1], HAB[cind1]), markersize = 10, color = :firebrick1)
     scatter!(axtf, Point((-1 .* ∇_turb_Wavg_prof_beg2)[cind2-1], HAB[cind2]), markersize = 10, color = :firebrick3)
     scatter!(axtf, Point((-1 .* ∇_turb_Wavg_prof_beg3)[cind3-1], HAB[cind3]), markersize = 10, color = :firebrick4)
 
    lines!(axfs, -1 .* ∇_turb_Wavg_prof_beg1 .- ∇_phasedep_Wavg_prof_beg1, HAB[2:end], color=:dodgerblue1, linewidth = 3)
        lines!(axfs,  -1 .* ∇_turb_Wavg_prof_beg2 .- ∇_phasedep_Wavg_prof_beg2, HAB[2:end], color=:dodgerblue3, linewidth = 3)
        lines!(axfs,  -1 .* ∇_turb_Wavg_prof_beg3 .- ∇_phasedep_Wavg_prof_beg3, HAB[2:end], color=:dodgerblue4, linewidth = 3)
        lines!(axfs, deltaline1, delta1.*ones(deltaline1L), color = :firebrick1, linewidth = 2, linestyle = :dash)
        lines!(axfs, deltaline1, delta2.*ones(deltaline1L), color = :firebrick3, linewidth = 2, linestyle = :dash)
        lines!(axfs, deltaline1, delta3.*ones(deltaline1L), color = :firebrick4, linewidth = 2, linestyle = :dash)
        scatter!(axfs, Point((-1 .* ∇_turb_Wavg_prof_beg1 .- ∇_phasedep_Wavg_prof_beg1)[cind1-1], HAB[cind1]), markersize = 10, color = :firebrick1)
        scatter!(axfs, Point((-1 .* ∇_turb_Wavg_prof_beg2 .- ∇_phasedep_Wavg_prof_beg2)[cind2-1], HAB[cind2]), markersize = 10, color = :firebrick3)
        scatter!(axfs, Point((-1 .* ∇_turb_Wavg_prof_beg3 .- ∇_phasedep_Wavg_prof_beg3)[cind3-1], HAB[cind3]), markersize = 10, color = :firebrick4)

    deltaline3 = 0.03:0.01:0.07
    deltaline3L = length(deltaline3)
   
    lines!(axc, c_Wavg_prof_beg1, HAB, color=:dodgerblue1, linewidth = 3)
        lines!(axc, c_Wavg_prof_beg2, HAB, color=:dodgerblue3, linewidth = 3)
        lines!(axc, c_Wavg_prof_beg3, HAB, color=:dodgerblue4, linewidth = 3)
        lines!(axc, deltaline3, delta1.*ones(deltaline3L), color = :firebrick1, linewidth = 2, linestyle = :dash)
        lines!(axc, deltaline3, delta2.*ones(deltaline3L), color = :firebrick3, linewidth = 2, linestyle = :dash)
        lines!(axc, deltaline3, delta3.*ones(deltaline3L), color = :firebrick4, linewidth = 2, linestyle = :dash)
        scatter!(axc, Point(c_Wavg_prof_beg1[cind1], HAB[cind1]), markersize = 10, color = :firebrick1)
        scatter!(axc, Point(c_Wavg_prof_beg2[cind2], HAB[cind2]), markersize = 10, color = :firebrick3)
        scatter!(axc, Point(c_Wavg_prof_beg3[cind3], HAB[cind3]), markersize = 10, color = :firebrick4)

    deltaline4 = -1.2e-8:1e-9:-4e-9
    deltaline4L = length(deltaline4)
    
    lines!(axsgs, SGS_Wavg_prof_beg1, HAB, color=:dodgerblue1, linewidth = 3)
        lines!(axsgs, SGS_Wavg_prof_beg2, HAB, color=:dodgerblue3, linewidth = 3)
        lines!(axsgs, SGS_Wavg_prof_beg3, HAB, color=:dodgerblue4, linewidth = 3)
        lines!(axsgs, deltaline4, delta1.*ones(deltaline4L), color = :firebrick1, linewidth = 2, linestyle = :dash)
        lines!(axsgs, deltaline4, delta2.*ones(deltaline4L), color = :firebrick3, linewidth = 2, linestyle = :dash)
        lines!(axsgs, deltaline4, delta3.*ones(deltaline4L), color = :firebrick4, linewidth = 2, linestyle = :dash)
        scatter!(axsgs, Point(SGS_Wavg_prof_beg1[cind1], HAB[cind1]), markersize = 10, color = :firebrick1)
        scatter!(axsgs, Point(SGS_Wavg_prof_beg2[cind2], HAB[cind2]), markersize = 10, color = :firebrick3)
        scatter!(axsgs, Point(SGS_Wavg_prof_beg3[cind3], HAB[cind3]), markersize = 10, color = :firebrick4)

    lines!(axdbdt, -1 .* ∇_turb_Wavg_prof_beg1 .- ∇_phasedep_Wavg_prof_beg1 .+ SGS_Wavg_prof_beg1[2:end], HAB[2:end], color=:dodgerblue1, linewidth = 3)
        lines!(axdbdt,  -1 .* ∇_turb_Wavg_prof_beg2 .- ∇_phasedep_Wavg_prof_beg2 .+ SGS_Wavg_prof_beg2[2:end], HAB[2:end], color=:dodgerblue3, linewidth = 3)
        lines!(axdbdt,  -1 .* ∇_turb_Wavg_prof_beg3 .- ∇_phasedep_Wavg_prof_beg3 .+ SGS_Wavg_prof_beg3[2:end], HAB[2:end], color=:dodgerblue4, linewidth = 3)
        lines!(axdbdt, deltaline1, delta1.*ones(deltaline1L), color = :firebrick1, linewidth = 2, linestyle = :dash)
        lines!(axdbdt, deltaline1, delta2.*ones(deltaline1L), color = :firebrick3, linewidth = 2, linestyle = :dash)
        lines!(axdbdt, deltaline1, delta3.*ones(deltaline1L), color = :firebrick4, linewidth = 2, linestyle = :dash)
        scatter!(axdbdt, Point((-1 .* ∇_turb_Wavg_prof_beg1 .- ∇_phasedep_Wavg_prof_beg1 .+ SGS_Wavg_prof_beg1[2:end])[cind1-1], HAB[cind1]), markersize = 10, color = :firebrick1)
        scatter!(axdbdt, Point((-1 .* ∇_turb_Wavg_prof_beg2 .- ∇_phasedep_Wavg_prof_beg2 .+ SGS_Wavg_prof_beg2[2:end])[cind2-1], HAB[cind2]), markersize = 10, color = :firebrick3)
        scatter!(axdbdt, Point((-1 .* ∇_turb_Wavg_prof_beg3 .- ∇_phasedep_Wavg_prof_beg3 .+ SGS_Wavg_prof_beg3[2:end])[cind3-1], HAB[cind3]), markersize = 10, color = :firebrick4)

    colgap!(ga, 15)
    rowgap!(ga, 5)

    axislegend(axpf, position = :lt)

#    savename = "SlopeNormalProfiles_SGSFluxes_" * sn
savename = "SlopeNormalProfiles_DeltaComp100m_wSGS_" * timeindex 

save(path_name * savename * ".png", f1)
