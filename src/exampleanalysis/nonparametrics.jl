function mt_wavedirection(signal::AbstractMatrix{T}, n, n_overlap; fs=1,onesided=true, nw::Real=4, ntapers::Int=ceil(Int, 2nw)-1) where {T}
    config = DSP.MTSpectrogramConfig{T}(size(signal,2), n, n_overlap; fs=fs,onesided=onesided,nw=nw,ntapers=ntaper)
    configcross = DSP.MTCrossSpectraConfig{T}(size(signal,1), n; fs = fs,onesided=onesided,nw=nw,ntapers=ntaper)
    X = DSP.allocate_output(config)
    tempoutput = DSP.allocate_output(configcross)
    epochs = [DSP.arraysplit(r,n,n_overlap) for r in eachrow(signal)]
    shortsignal = zeros(eltype(epochs[1][1]), size(signal,1), length(epochs[1][1]))
    @views for (time_index,(z,x,y)) in enumerate(zip(epochs...))
        shortsignal[1,:] .= z
        shortsignal[2,:] .= x
        shortsignal[3,:] .= y
        mt_cross_power_spectra!(tempoutput, shortsignal, configcross)
        X[:,time_index] .= angle.(imag.(tempoutput[2,1,:]) .+ 1im .* imag.(tempoutput[3,1,:]))
    end
    return DSP.Periodograms.Spectrogram(X, config.mt_config.freq, config.time)
end

function compute_spectrogram(waveseries, wavetime)
    continuouswave = reduce(vcat, w[:,1] for w in waveseries)
    spec = mt_spectrogram(continuouswave, 2304, 2304÷2, fs = 1, nw=4, ntapers=7) # fs has to be int for some reason
    spectrogram = spec.power./(2π * 1.28) ./2 # equivalent to * Δ/2π and div 2 because onesided
    spectime = wavetime[1]:(wavetime[2]-wavetime[1])/2:wavetime[end]
    Ω = collect((2π*1.28) .* spec.freq)
    return spectrogram, Ω, spectime
end

function compute_errorgram(waveseries, wavetime)
    continuouswave = [reduce(vcat, w[:,i] for w in waveseries) for i in 1:3]
    spec = [mt_spectrogram(c, 2304, 0, nw=4, ntapers=7) for c in continuouswave]
    spectrogram = [s.power./(2π * 1.28) ./2 for s in spec] # equivalent to * Δ/2π and div 2 because onesided (DSP only allows spectrogram with integer fs so rescale is necessary)
    Ω = collect((2π*1.28) .* spec[1].freq)
    errorgram = @. log(spectrogram[2]+spectrogram[3])+2*log(tanh(OceanWaveSpectralFitting.approx_dispersion(Ω,40)))-log(spectrogram[1])
    return errorgram, Ω, wavetime
end