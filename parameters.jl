using Parameters

@with_kw struct SimParams
    Lx = 152
    Lzˢ = 500
    Lyˢ = 5500

    # parameters
    Ñˢ = 3.5 * 10^(-3)              # buoyancy frequency
    fˢ = Ñˢ/10.7                    # inertial frequency
    σˢ = 2.2*fˢ                     # tidal frequency
    V₀ˢ = 0.1                        # tide amplitude

    # topography parameters
    Tanθˢ = sqrt((σˢ^2 - fˢ^2)/(Ñˢ^2-σˢ^2))        # slope of internal wave energy propogation
    γᶜ = 1.9                               # when bulk slope is supercritical
    Tanαˢ = γᶜ * Tanθˢ                       # bulk topographic slope
    δˢ = V₀ˢ/Ñˢ

    dh = 4
    dz = 2
    dtˢ = 600.0 # sampling interval

    res1 = (dz=dz, dh=dh) 
    const_params = merge(res1, (; dye_height=20, dye_smoothing=10, Sp_R=500, g_width = 500/3, 
                    ycenter = 895, slope_endˢ = Lzˢ/Tanαˢ,Ñˢ =Ñˢ,Tanαˢ=Tanαˢ,
                    Lx=Lx, Lzˢ=Lzˢ, Lyˢ=Lyˢ, Tfˢ = 2*π/fˢ, Tσˢ = 2*π/σˢ))


    U5N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*0.05,     Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,     γ=γᶜ, dt=dtˢ))
    U10N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*0.1,     Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,     γ=γᶜ, dt=dtˢ))
    U50N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*0.5,     Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,     γ=γᶜ, dt=dtˢ))
    U100N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ,        Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,     γ=γᶜ, dt=dtˢ))
    U150N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,     γ=γᶜ, dt=dtˢ))
    U200N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*2,      Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,     γ=γᶜ, dt=dtˢ))
    U250N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*2.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,     γ=γᶜ, dt=dtˢ))
    U300N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*3,      Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U350N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*3.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U400N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*4,      Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U450N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U500N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*5,      Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U550N100Lz100g100 = merge(const_params,(;U₀=V₀ˢ*5.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ, γ=γᶜ, dt=dtˢ))

    U450N100Lz130g100 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=1.3*Lzˢ, γ=γᶜ, dt=dtˢ))
    U500N100Lz130g100 = merge(const_params,(;U₀=V₀ˢ*5,      Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=1.3*Lzˢ, γ=γᶜ, dt=dtˢ))
    U550N100Lz130g100 = merge(const_params,(;U₀=V₀ˢ*5.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=1.3*Lzˢ, γ=γᶜ, dt=dtˢ))

    U250N100Lz100g100Ra = merge(const_params,(;U₀=V₀ˢ*2.5,  Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=1.3*Lzˢ, γ=γᶜ, dt=dtˢ))
    U150N100Lz100g100Ra = merge(const_params,(;U₀=V₀ˢ*1.5,  Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=1.3*Lzˢ, γ=γᶜ, dt=dtˢ))
    U350N100Lz100g100Ra = merge(const_params,(;U₀=V₀ˢ*3.5,  Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=1.3*Lzˢ, γ=γᶜ, dt=dtˢ))


    # varying N and by proxy varying f and sigma with fixed U0
    U250Nfd50Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ,  Ñ=2.5*V₀ˢ/(0.5*δˢ), f=2.5*V₀ˢ/(0.5*δˢ*10.7), 
            σ=5.5*V₀ˢ/(0.5*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.2*dtˢ))
    U250Nfd100Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(δˢ),     f=2.5*V₀ˢ/(1.0*δˢ*10.7), 
            σ=5.5*V₀ˢ/(1.0*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.4*dtˢ))
    U250Nfd150Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(1.5*δˢ), f=2.5*V₀ˢ/(1.5*δˢ*10.7), 
            σ=5.5*V₀ˢ/(1.5*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.6*dtˢ))
    U250Nfd200Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(2.0*δˢ), f=2.5*V₀ˢ/(2.0*δˢ*10.7), 
            σ=5.5*V₀ˢ/(2.0*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.8*dtˢ))
    U250Nfd250Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(2.5*δˢ), f=2.5*V₀ˢ/(2.5*δˢ*10.7), 
            σ=5.5*V₀ˢ/(2.5*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U250Nfd300Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(3.0*δˢ), f=2.5*V₀ˢ/(3.0*δˢ*10.7), 
            σ=5.5*V₀ˢ/(3.0*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U250Nfd350Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(3.5*δˢ), f=2.5*V₀ˢ/(3.5*δˢ*10.7), 
            σ=5.5*V₀ˢ/(3.5*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U250Nfd400Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(4.0*δˢ), f=2.5*V₀ˢ/(4.0*δˢ*10.7), 
            σ=5.5*V₀ˢ/(4.0*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=dtˢ))
    U250Nfd450Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(4.5*δˢ), f=2.5*V₀ˢ/(4.5*δˢ*10.7), 
            σ=5.5*V₀ˢ/(4.5*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=1.7*dtˢ))
    U250Nfd500Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(5.0*δˢ), f=2.5*V₀ˢ/(5.0*δˢ*10.7), 
            σ=5.5*V₀ˢ/(5.0*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=1.9*dtˢ))
    U250Nfd550Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ, Ñ=2.5*V₀ˢ/(5.5*δˢ), f=2.5*V₀ˢ/(5.5*δˢ*10.7), 
            σ=5.5*V₀ˢ/(5.5*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=2.1*dtˢ))
   
    # varying N a little more closely at the small values
    U250Nfd60Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ,  Ñ=2.5*V₀ˢ/(0.6*δˢ), f=2.5*V₀ˢ/(0.6*δˢ*10.7), 
            σ=5.5*V₀ˢ/(0.6*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.2*dtˢ))
    U250Nfd70Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ,  Ñ=2.5*V₀ˢ/(0.7*δˢ), f=2.5*V₀ˢ/(0.7*δˢ*10.7),
            σ=5.5*V₀ˢ/(0.7*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.3*dtˢ))
    U250Nfd85Lz100g100 = merge(const_params,(;U₀=2.5*V₀ˢ,  Ñ=2.5*V₀ˢ/(0.85*δˢ), f=2.5*V₀ˢ/(0.85*δˢ*10.7), 
            σ=5.5*V₀ˢ/(0.85*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.3*dtˢ))

    # adding larger stratification values 
    U500Nfd51Lz100g100 = merge(const_params,(;U₀=5.0*V₀ˢ,  Ñ=5.0*V₀ˢ/(0.51*δˢ), f=5.0*V₀ˢ/(0.51*δˢ*10.7), 
            σ=11.0*V₀ˢ/(0.51*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.1*dtˢ))
    U500Nfd58Lz100g100 = merge(const_params,(;U₀=5.0*V₀ˢ,  Ñ=5.0*V₀ˢ/(0.58*δˢ), f=5.0*V₀ˢ/(0.58*δˢ*10.7), 
            σ=11.0*V₀ˢ/(0.58*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.1*dtˢ))
    U500Nfd66Lz100g100 = merge(const_params,(;U₀=5.0*V₀ˢ,  Ñ=5.0*V₀ˢ/(0.66*δˢ), f=5.0*V₀ˢ/(0.66*δˢ*10.7), 
            σ=11.0*V₀ˢ/(0.66*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.1*dtˢ))
    U500Nfd78Lz100g100 = merge(const_params,(;U₀=5.0*V₀ˢ,  Ñ=5.0*V₀ˢ/(0.78*δˢ), f=5.0*V₀ˢ/(0.78*δˢ*10.7), 
            σ=11.0*V₀ˢ/(0.78*δˢ*10.7), Lz=Lzˢ, γ=γᶜ, dt=0.2*dtˢ))

    # varying sigma on its own to change criticality Lx = Lxˢ
    U150N100Lz100g70s125 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=2.8*fˢ, Lz=Lzˢ,     γ=1.404, dt=0.8*dtˢ))
    U150N100Lz100g30s250 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=5.5*fˢ, Lz=Lzˢ,     γ=0.603, dt=0.4*dtˢ))
    U150N100Lz100g30s80 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=0.8*fˢ, Lz=Lzˢ,     γ=0.603, dt=dtˢ))
    U450N100Lz100g70s80 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=0.8*fˢ, Lz=Lzˢ,     γ=1.404, dt=dtˢ))
    U450N100Lz100g70s125 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=2.8*fˢ, Lz=Lzˢ,     γ=1.404, dt=0.8*dtˢ))
    U450N100Lz100g30s250 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=5.5*fˢ, Lz=Lzˢ,     γ=0.603, dt=0.4*dtˢ))

    # varying sigma on its own no change to criticality (changes topo)
    U150N100Lz100g100s125 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=2.8*fˢ, Lz=Lzˢ,   γ=γᶜ, dt=0.8*dtˢ))
    U150N100Lz100g100s250 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=5.5*fˢ, Lz=Lzˢ,   γ=γᶜ, dt=0.4*dtˢ))
    U150N100Lz100g100s80 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=1.8*fˢ, Lz=Lzˢ,    γ=γᶜ, dt=dtˢ))
    U450N100Lz100g100s80 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=1.8*fˢ, Lz=Lzˢ,    γ=γᶜ, dt=dtˢ))
    U450N100Lz100g100s125 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=2.8*fˢ, Lz=Lzˢ,   γ=γᶜ, dt=0.8*dtˢ))
    U450N100Lz100g100s250 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=5.5*fˢ, Lz=Lzˢ,   γ=γᶜ, dt=0.4*dtˢ))

    # varying criticality without varying sigma
    U150N100Lz100g70s100 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,         γ=1.4, dt=dtˢ))
    U150N100Lz100g30s100 = merge(const_params,(;U₀=V₀ˢ*1.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=0.6*Lzˢ,     γ=0.6, dt=dtˢ))
    U450N100Lz100g70s100 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=Lzˢ,         γ=1.4, dt=dtˢ))
    U450N100Lz100g30s100 = merge(const_params,(;U₀=V₀ˢ*4.5,    Ñ=Ñˢ, f=fˢ, σ=σˢ, Lz=0.6*Lzˢ,     γ=0.6, dt=dtˢ))

end
 
