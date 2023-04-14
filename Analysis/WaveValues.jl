function get_wave_indices(timeseries, p, tlength)
    
    @info "Calculating Wave Indices..."
    # we are going to average over each waves
    nTσ = round(Int,timeseries.times[tlength]/p.Tσ)

    isnada(val) = val == nothing

    T_Tσs = ones(Int,nTσ)
    for i = 2:nTσ
        abovTs(t) = t>=((i-1)*p.Tσ)
        idx = findfirst(abovTs, timeseries.times)
        # if there is not an index at that value use the last index
        T_Tσs[i] = isnada(idx) ? tlength : idx
    end

    Wl =  round(Int, mean(T_Tσs[2:end]-T_Tσs[1:nTσ-1]))
    WavePeriods= Array{Int64}(undef, Wl, nTσ)
    WavePeriods[1,:] = T_Tσs
    mT_Tσs = ones(Int,nTσ)
    for i=2:Wl
        fact = timeseries.times[i]
        mT_Tσs[1] = i
        for j = 1:nTσ-1
            abovWs(t) = t>=((j*p.Tσ) + fact)
            wt = findfirst(abovWs, timeseries.times)
            mT_Tσs[j+1] = isnada(wt) ? T_Tσs[j] : wt
        end
        WavePeriods[i,:] = mT_Tσs
    end

    wave_vals = (; WavePeriods=WavePeriods, nTσ=nTσ, Wl=Wl, T_Tσs=T_Tσs)

    return wave_vals

end
