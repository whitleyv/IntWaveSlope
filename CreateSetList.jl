using JLD2
using Parameters

u = 50:50:550
setL = length(u)
Usets = "U" .* string.(u) .* "N100Lz100g100"
Nsets = "U250Nfd" .* string.(u) .* "Lz100g100"

# N repeat U250 so just replace with other one
@inline fun(z) = z==250
index250 = findfirst(fun, u)
Nsets[index250] = Usets[index250]

setnames = vcat(Usets, Nsets)

common_prefix = "IntWave_mp_"

Allfiles = split(repeat(common_prefix* " ", setL*2), " ")
deleteat!(Allfiles, length(Allfiles))

# all the N sets have same file name right now
Nprefix =  "noS_ "
Nfiles = split(repeat(Nprefix, setL), " ")
# remove last since it's just a space not a file
deleteat!(Nfiles, length(Nfiles))

Allfiles[setL+1:end] .= Allfiles[setL+1:end] .* Nfiles

# almost all the big U sets have no extra domain names right now

LargeUprefix =  "nExD_ "
Ubigfiles = split(repeat(LargeUprefix, length(u[index250+1:setL])), " ")
deleteat!(Ubigfiles, length(Ubigfiles))

Allfiles[index250+1:setL] .= Allfiles[index250+1:setL] .* Ubigfiles

# special cases:
Allfiles[index250+2] = common_prefix * "noS_"
Allfiles[index250+4] = common_prefix * "noS_"
Allfiles[setL+index250] = Allfiles[index250]

include("parameters.jl")

δs = zeros(setL*2)
Ns = zeros(setL*2)
Us = zeros(setL*2)
for (i,set) in enumerate(setnames)
    pm = getproperty(SimParams(), Symbol(set))
    δs[i] = pm.U₀/pm.Ñ
    Ns[i] = pm.Ñ
    Us[i] = pm.U₀
end

setnames_varyσ = ["U150N100Lz100g70s125", "U150N100Lz100g30s250", 
"U450N100Lz100g70s125", "U450N100Lz100g30s250", 
"U150N100Lz100g100s125", "U150N100Lz100g100s250",
 "U450N100Lz100g100s125", "U450N100Lz100g100s250"]
δ_VS= zeros(length(setnames_varyσ))
Ns_VS= zeros(length(setnames_varyσ))
σs_VS= zeros(length(setnames_varyσ))
γs_VS= zeros(length(setnames_varyσ))

for (m,setname) in enumerate(setnames_varyσ)
    pm2 = getproperty(SimParams(), Symbol(setname))
    
    δ_VS[m] = pm2.U₀/pm2.Ñ
    Ns_VS[m] = pm2.Ñ
    σs_VS[m] = pm2.σ
    γs_VS[m] = pm2.γ
end

σfiles = split(repeat(Nprefix, length(setnames_varyσ)), " ")
setfilenames_varyσ = common_prefix .* σfiles

# remove last since it's just a space not a file
deleteat!(σfiles, length(σfiles))

savename = "SetList_mp.jld2"
jldsave(savename;
setnames,
setnames_varyσ,
setfilenames = Allfiles,
setfilenames_varyσ,
δs,
Ns,
Us,
γ_varyσ=γs_VS,
N_varyσ=Ns_VS,
σ_varyσ=σs_VS,
δ_varyσ=δ_VS,
)