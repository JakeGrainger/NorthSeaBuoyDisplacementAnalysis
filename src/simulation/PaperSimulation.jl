module PaperSimulation
    using OceanWaveSpectralFitting, Optim, ProgressMeter
    import OceanWaveSpectralFitting: WhittleLikelihoodInference.GaussianProcess
    const WindSea = JS_BWG_HNE{0} # simulation with no aliasing
    include("../freqdirspectra/FrequencyDirectionSpectra.jl")
    include("leastsquares.jl")
    include("methodofmoments.jl")
    include("whittle.jl")

    ##
    "Helper function to initialise peak"
    function initialisepeak(ts::Matrix{Float64}, Δ, ωₗ, ωᵤ)
        P = OceanWaveSpectralFitting.Periodogram(ts,Δ)
        S = P.ordinates[end÷2+1:end]
        Ω = P.Ω[end÷2+1:end]
        peakind = argmax([(ωₗ < ω < ωᵤ) ? real(s[1,1]) : -Inf for (s,ω) ∈ zip(S, Ω)])
        meandir = angle(imag(S[peakind][2,1])+1im*imag(S[peakind][3,1]))
        return Ω[peakind], mod2pi.(meandir.+pi) .- pi
    end
    ##
    "Simulation study struct"
    struct SimulationStudy
        LSmlm::Vector{Vector{Float64}}
        LSmem::Vector{Vector{Float64}}
        M::Vector{Vector{Float64}}
        DW::Vector{Vector{Float64}}
        DWuni::Vector{Vector{Float64}}
        LSmlmflag::Vector{Bool}
        LSmemflag::Vector{Bool}
        Mflag::Vector{Bool}
        DWflag::Vector{Bool}
        DWuniflag::Vector{Bool}
        θ_true::Vector{Float64}
        cutoff::Float64
        function SimulationStudy(
            nreplications::Int,
            θ::Vector{Float64};
            ωₗ::Float64 = 0.0,
            name = "Simulation ",
        )
            @assert length(θ) == 9
            Δ = 1 / 1.28
            X = GaussianProcess(WindSea(θ), 2304, Δ)
            LSmlm = Vector{Vector{Float64}}(undef, nreplications)
            LSmlmflag = Vector{Bool}(undef, nreplications)
            LSmem = similar(LSmlm)
            M = similar(LSmlm)
            DW = similar(LSmlm)
            DWuni = similar(LSmlm)
            LSmemflag = similar(LSmlmflag)
            Mflag = similar(LSmlmflag)
            DWflag = similar(LSmlmflag)
            DWuniflag = similar(LSmlmflag)
            @showprogress name for rep ∈ 1:nreplications
                ts = rand(X).ts
                peak, dir = initialisepeak(ts, Δ, ωₗ, π / Δ)
                θ_ini = [0.7, peak, 2, 4.5, dir, 4.5, 3, 0.6, 0.3]
                LSmlmflag[rep], LSmlm[rep] =
                    LSfit(FrequencyDirectionSpectra.FDJONSWAP, ts, Δ, θ_ini, fdest = FrequencyDirectionSpectra.MLM; ωₗ = ωₗ)
                LSmemflag[rep], LSmem[rep] =
                    LSfit(FrequencyDirectionSpectra.FDJONSWAP, ts, Δ, θ_ini, fdest = FrequencyDirectionSpectra.MEM; ωₗ = ωₗ)
                Mflag[rep], M[rep] = momentsfit(FrequencyDirectionSpectra.FDJONSWAP, ts, Δ, θ_ini; ωₗ = ωₗ)
                DWflag[rep], DW[rep] = DWfit(WindSea, ts, Δ, θ_ini; ωₗ = ωₗ)
                DWuniflag[rep], DWuni[rep] = DWfit(JONSWAP{0}, ts[:,1], Δ, θ_ini[1:4]; ωₗ = ωₗ)
            end
            new(
                LSmlm,
                LSmem,
                M,
                DW,
                DWuni,
                LSmlmflag,
                LSmemflag,
                Mflag,
                DWflag,
                DWuniflag,
                θ,
                ωₗ,
            )
        end
    end
    Base.show(io::IO, sim::SimulationStudy) = print(io,"Simulation study results.")

end