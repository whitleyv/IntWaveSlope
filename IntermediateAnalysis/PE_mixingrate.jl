using Printf
using Oceananigans
using JLD2
 
sn = "U250N100Lz100g100"

include("parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = -π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

zlength = pm.nz
ylength = 880
xlength = 38
tlength = 161

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

setname = sn

path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/"
apath = path_name * "Analysis/"

filescalename = apath * "PE_mixingrate.jld2"
filescalename2 = apath * "PE_mixingrateavg.jld2"

global Us = .05:.05:.55
global setnames=[]
for u = 50:50:550
    global setnames = [setnames ; @sprintf("U%dN100Lz100g100", u)]
end

Lvals = length(setnames)

ϕd_610avg = zeros(Lvals)
ϕi_610avg = zeros(Lvals)
ϕe_610avg = zeros(Lvals)
δ = zeros(Lvals)
Ñn = zeros(Lvals)

for (m, setname) in enumerate(setnames)
        
    pm2 = getproperty(SimParams(), Symbol(setname))

    pm2 = merge(pm2, (; Tanθ = sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
    Tanα = pm2.γ * sqrt((pm2.σ^2 - pm2.f^2)/(pm2.Ñ^2-pm2.σ^2)),
    m = -π/pm2.Lz,
    l = sqrt(((π/pm2.Lz)^2 * (pm2.f^2 - pm2.σ^2)) / (pm2.σ^2 - pm2.Ñ^2)),
    Tf = 2*π/pm2.f, 
    Tσ = 2*π/pm2.σ))

    δ[m] = pm2.U₀/pm2.Ñ
    Ñn[m] = pm2.Ñ

    # need to recalculate parameters each time because N, and wave length will change!
    @info "Pulling data..."
        
    name_prefix = "IntWave_" * setname
    filepath = path_name * name_prefix * ".jld2"

    κ_timeseries = FieldTimeSeries(filepath, "κₑ");
    b_timeseries = FieldTimeSeries(filepath, "b");
    N2_timeseries = FieldTimeSeries(filepath, "N2");

    xb, yb, zb = nodes(b_timeseries) #CCC

    bi = interior(b_timeseries)[:,1:ylength+1,1:zlength,1:tlength]; #ccc
    κi = interior(κ_timeseries)[:,1:ylength+1,1:zlength,1:tlength]; #ccc
    N2i = interior(N2_timeseries)[:,1:ylength+1,1:zlength+1,1:tlength]; #ccf

    Δx = 1/(xb[2] - xb[1])
    Δy = 1/(yb[2] - yb[1])
    Δz = 1/(zb[2] - zb[1])

    @info "Calculating buoyancy gradient..."
    ∂xb = bi[3:end,2:end-1,:,:] .- bi[1:end-2,2:end-1,:,:]; #CCC
    ∂yb =  bi[2:end-1,3:end,:,:] .- bi[2:end-1,1:end-2,:,:]; #CCC

    ∂b∂x = ∂xb.*Δx.*0.5;
    ∂b∂y = ∂xb.*Δy.*0.5;
    ∂b∂z = 0.5 .* (N2i[2:end-1,2:end-1,1:end-1,:] .+ N2i[2:end-1,2:end-1,2:end,:]);
    ∇b² = ∂b∂x.^2 .+ ∂b∂y.^2 .+ ∂b∂z.^2

    xlengthcut = xlength-2

    @info "Sorting to background PE..."

    Ygrid = reshape(repeat(yb[2:ylength], zlength), ylength-1, zlength)
    YXgrid = reshape(repeat(yb[2:ylength], zlength*xlengthcut), ylength-1, zlength, xlengthcut)
    SlopeGridYX = curvedslope.(YXgrid)
    SlopeGridY = curvedslope.(Ygrid)
    boolZYX = (SlopeGridYX .+ 4) .<= zb[1:zlength]' # all the values greater than slope
    non_boolZY = (SlopeGridY .+ 4) .> zb[1:zlength]' # all the values in the slope

    first_allowed_y_index = sum(non_boolZY, dims=1) .+ 1
    first_z_index_in_domain = sum(first_allowed_y_index .> ylength) + 1

    boolZ_perm = permutedims(boolZYX, [2, 1,3])
    #non_boolZ_perm = permutedims(non_boolZ, [2,1,3])

    background_APEb = zeros(xlengthcut, ylength-1, zlength, tlength);

    for i=1:tlength
        permdimsb = permutedims(bi[2:end-1,2:end-1,:, i], [3,2,1])

        sortb = sort(Array(permdimsb[boolZ_perm]))

        running_index_total = 1
        # to reshape you need to take into account the slope:
        for k = first_z_index_in_domain:zlength
            first_yidx = first_allowed_y_index[k]
            yL = length(first_yidx:ylength-1)
            num_indices_in_plane = yL*xlengthcut
            new_total = running_index_total + num_indices_in_plane
            background_APEb[:,first_yidx:ylength-1,k,i] = reshape(sortb[running_index_total:new_total-1], xlengthcut,yL)
            running_index_total = new_total
        end

    end

    @info "Calculating vertical deriv of background PE..."

    # taking the z derivative cuts our domain by one grid point on either end
    sorted_∂b∂z = Δz.*(background_APEb[:,:,3:end,:] .- background_APEb[:,:,1:end-2,:])

    @info "Calculating PE dissipation terms..."

    ϕi_preInt = κi[2:end-1,2:end-1,2:end-1,:] .* ∂b∂z[:,:,2:end-1,:]
    ϕd_preInt = κi[2:end-1,2:end-1,2:end-1,:] .* ∇b²[:,:,2:end-1,:] ./  sorted_∂b∂z     
    boolXYZ = permutedims(boolZYX, [3, 1, 2])[:,:,2:end-1]
    # boolean to remove values of things near slope weher zero under buoynacy val.
    boolposN = sorted_∂b∂z.>=0

    ϕi = zeros(tlength)
    ϕd = zeros(tlength)

    for i = 1:tlength
        ϕi_preInti = ϕi_preInt[i]
        ϕd_preInti = ϕd_preInt[i]
        boolposNi = boolposN[:,:,:,i]

        boolposXYZi = boolXYZ .& boolposNi
        
        ϕi[i] = (4*4*2) .* sum(ϕi_preInti[boolposXYZi])
        ϕd[i] = (4*4*2) .* sum(ϕd_preInti[boolposXYZi])
    end 

    ϕe = ϕd .- ϕi

    include("WaveValues.jl")
    wave_info=get_wave_indices(b_timeseries, pm2, tlength)

    Tσ6_idx = wave_info.T_Tσs[7] # completion of 6 waves
    Tσ10_idx = wave_info.T_Tσs[11] # completion of 6 waves

    ϕe_610avg[m] = mean(ϕe[Tσ6_idx:Tσ10_idx])
    ϕd_610avg[m] = mean(ϕd[Tσ6_idx:Tσ10_idx])
    ϕi_610avg[m] = mean(ϕi[Tσ6_idx:Tσ10_idx])

    @info "Saving to file..."

    jldopen(filescalename, "a+") do file
        mygroup = JLD2.Group(file, setname)
        mygroup["ϕd"] = ϕd
        mygroup["ϕi"] = ϕi
        mygroup["ϕe"] = ϕe
        mygroup["E_b"] = background_APEb
        mygroup["N_b"] = sorted_∂b∂z
    end


    # rolling wave average?
    @info "Computing Rolling Wave Averages..."
    # averages based on true wave times not consistent across N varying sims

    Wl = wave_info.Wl
    wav_arr_length = Wl*wave_info.nTσ
    wav_Avgarr_length = wav_arr_length - 2*Wl

    # global value in time
    ϕe_Wavg = zeros(wav_Avgarr_length);
    ϕd_Wavg = zeros(wav_Avgarr_length);
    ϕi_Wavg = zeros(wav_Avgarr_length);

    for (ki, tk) in enumerate(wave_info.WavePeriods[Wl+1:end-Wl])
        # glo2Wlbal value in time
        Widxs = wave_info.WavePeriods[ki:ki+2*Wl] 
        ϕe_Wavg[ki] = mean(ϕe[Widxs])
        ϕd_Wavg[ki] = mean(ϕd[Widxs])
        ϕi_Wavg[ki] = mean(ϕi[Widxs])

    end

    wavelength = 1+1/Wl:1/Wl:wave_info.nTσ-1

    ep2 = plot(wavelength, ϕe_Wavg, lw = 5, color = :green,
            label=@sprintf("<ϕe>₆_₁₀= %0.3f", ϕe_610avg[m]),
            xlabel = "Tσ", ylabel="ϕe [m²s⁻³]", legend_font = font(20), 
            guidefontsize = 20, titlefont=20, tickfont = 14, bottom_margin=10.0mm, left_margin=10.0mm, right_margin=20.0mm,
            legend = :topleft, size = (1000,800), title = big_title)
        plot!(wavelength, ϕd_Wavg, lw = 5, color = :gray50, label="ϕd")
        plot!(wavelength, ϕi_Wavg, lw = 5, color = :gray30, label="ϕi")


    savefig(ep2, apath * "PE_rateWav_" * setname * ".png")


end

jldsave(filescalename2; setnames, 
ϕe_610avg, ϕd_610avg, ϕi_610avg,
            δ, Ñn)
