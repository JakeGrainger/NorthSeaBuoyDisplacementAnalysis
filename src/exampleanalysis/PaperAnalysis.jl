module PaperAnalysis
    using OceanWaveSpectralFitting, Optim, DSP
    include("nonparametrics.jl")
    
    const WindSea = JS_BWG_HNE_DL{0,40} # no aliasing and 40m water depth

    "
        movingaverage(vector, half_window_size)
    
    Compute the moving average of a vector"
    function movingaverage(vector::Vector{T}, half_window_size) where {T}
        out = similar(vector)
        for i = 1:length(vector)
            boundarydistance = min(length(vector) - i, i - 1)
            window = min(half_window_size, boundarydistance)
            @views out[i] = sum(vector[i-window:i+window])/(2window+1)
        end
        return out
    end

    "
        initialisepeak(ts, Δ, ωₗ, ωᵤ)

    Compute the peak frequency and direction for a certain frequency band."
    function initialisepeak(ts::Matrix{Float64}, Δ, ωₗ, ωᵤ)
        P = OceanWaveSpectralFitting.Periodogram(ts,Δ)
        S = P.ordinates[end÷2+1:end]
        Ω = P.Ω[end÷2+1:end]
        peakind = argmax([(ωₗ < ω < ωᵤ) ? real(s[1,1]) : -Inf for (s,ω) ∈ zip(S, Ω)])
        meandir = angle(imag(S[peakind][2,1])+1im*imag(S[peakind][3,1]))
        return Ω[peakind], mod2pi.(meandir.+pi) .- pi
    end

    "Function to generate an initial guess when we have a wind-sea only model."
    function initialguess(::Type{<:WindSea},ts,Δ,cutoffs)
        peak,mdir = initialisepeak(ts,Δ,maximum(cutoffs),4)
        return [0.7, peak, 2, 5, mdir, 4, 2.7, 0.55, 0.26]
    end

    "Type to store the timeseries, initial parameters, cutoff and proposed model."
    struct InitialiseData
        ts::Matrix{Float64}
        Δ::Float64
        cutoff::Float64
        model::Type{<:TimeSeriesModel}
        θ₀::Vector{Float64}
        function InitialiseData(model::Type{<:TimeSeriesModel},ts::Matrix{Float64},Δ::Float64,cutoff)
            θ₀ = initialguess(model,ts,Δ,cutoff)
            new(ts,Δ,cutoff,model,θ₀)
        end
    end

    "Make constraints for a wind-sea only model."
    function makeconstraints(::Type{<:WindSea},cut)
        lxW = [0, cut, 1,  3,-Inf, 0,  0,   0,  0  ]
        uxW = [4, 4,   Inf,8, Inf, 2π, Inf ,Inf,Inf]
        return TwiceDifferentiableConstraints(lxW,uxW)
    end

    "Type to store the results of fitting one choice of cutoff."
    struct FittingResults
        θ_est::Vector{Float64}
        varθ::Matrix{Float64}
        converged::Bool
        function FittingResults(ini,iterations)
            obj = DebiasedWhittleLikelihood(ini.model, ini.ts, ini.Δ; lowerΩcutoff = ini.cutoff, upperΩcutoff = 3.8)
            constraints = makeconstraints(ini.model, ini.cutoff)
            res = optimize(TwiceDifferentiable(Optim.only_fgh!(obj), ini.θ₀), constraints, ini.θ₀, IPNewton(), Optim.Options(iterations = iterations))
            obj = DebiasedWhittleLikelihood(ini.model, ini.ts, ini.Δ; lowerΩcutoff = ini.cutoff, upperΩcutoff = 3.8)
            EH = zeros(length(ini.θ₀),length(ini.θ₀))
            obj(nothing,nothing,EH,res.minimizer)
            varθ = nothing
            try 
                varθ = inv(EH)
            catch
                varθ = fill(NaN, size(EH))
            end
            new(res.minimizer, varθ, Optim.converged(res))
        end
    end

    "Type to store fitted models and the initial values."
    struct FittedModel
        ini::InitialiseData
        fits::FittingResults
        function FittedModel(ini::InitialiseData;iterations=3000)
            fits = FittingResults(ini,iterations)
            new(ini,fits)
        end
    end
end
