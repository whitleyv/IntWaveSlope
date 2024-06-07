using Statistics
using Printf
using Oceananigans
using JLD2


path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

setnames = ["U100N100Lz100g100", "U150N100Lz100g100", "U200N100Lz100g100", "U250N100Lz100g100", 
"U300N100Lz100g100", "U450N100Lz100g100"]

setnames2 = ["U300N100Lz100g100", "U400N100Lz100g100", "U500N100Lz100g100"]

Lvals = length(setnames)
c_weighted_bavg_mp = zeros(Lvals+3, 161)
c_weighted_bavg_mn = zeros(Lvals, 161)


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

@info "Usual Set..."
for (m, setname) in enumerate(setnames)

    name_mp_prefix = "IntWave_mp_" * setname
    name_reg_prefix = "IntWave_" * setname

    filepath_mp = path_name * name_mp_prefix * ".jld2"
    filepath_reg = path_name * name_reg_prefix * ".jld2"

    b_mp_timeseries = FieldTimeSeries(filepath_mp,"b");
    b_reg_timeseries = FieldTimeSeries(filepath_reg,"b");

    c_reg_timeseries = FieldTimeSeries(filepath_reg, "Cg");
    c_mp_timeseries = FieldTimeSeries(filepath_mp, "Cg");


    #############
    # Gaussian trcaer release
    #############

    ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)

    c_mp_Fi = interior(c_mp_timeseries)[:,1:ylength,:,:];
    c_reg_Fi = interior(c_reg_timeseries)[:,1:ylength,:,:];

    b_mp_Fi = interior(b_mp_timeseries)[:,1:ylength,:,:];
    b_reg_Fi = interior(b_reg_timeseries)[:,1:ylength,:,:];

    @info "Tracer Weighted Calculations.."

    b̄_Cg_mp = calculate_tracer_weight(c_mp_Fi, b_mp_Fi)
    b̄_Cg_reg = calculate_tracer_weight(c_reg_Fi, b_reg_Fi)

    c_weighted_bavg_mn[m, :] = b̄_Cg_reg
    c_weighted_bavg_mp[m,:] = b̄_Cg_mp

end

@info "Change in domain set..."

for (m,setname) in enumerate(setnames2)

    name_mp_prefix = "IntWave_mp_nExD_" * setname

    filepath_mp = path_name * name_mp_prefix * ".jld2"

    b_mp_timeseries = FieldTimeSeries(filepath_mp,"b");

    c_mp_timeseries = FieldTimeSeries(filepath_mp, "Cg");


    #############
    # Gaussian trcaer release
    #############

    ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)

    c_mp_Fi = interior(c_mp_timeseries)[:,1:ylength,:,:];

    b_mp_Fi = interior(b_mp_timeseries)[:,1:ylength,:,:];

    @info "Tracer Weighted Calculations.."

    b̄_Cg_mp = calculate_tracer_weight(c_mp_Fi, b_mp_Fi)

    c_weighted_bavg_mp[m+Lvals,:] = b̄_Cg_mp


end

filescalename = apath * "cwtdb1mom.jld2"

jldsave(filescalename; 
c_weighted_bavg_mp, c_weighted_bavg_mn,
setnames, setnames2,
)