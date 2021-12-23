"Return the Fourier coefficients of the bimodal wrapped Gaussian."
function GaussCoef(θ)
    φ₁ = θ[1]
    φ₂ = θ[2]
    σ = θ[3]
    a1 = exp(-(σ^2) / 2) * (cos(φ₁) + cos(φ₂))
    b1 = exp(-(σ^2) / 2) * (sin(φ₁) + sin(φ₂))
    a2 = exp(-2(σ^2)) * (cos(2φ₁) + cos(2φ₂))
    b2 = exp(-2(σ^2)) * (sin(2φ₁) + sin(2φ₂))
    return a1, b1, a2, b2
end
"Return the Fourier coefficients from the matrix of cross spectra."
function freqdir_fouriercoef2(S)
    denom = sqrt(abs((S[2, 2] + S[3, 3]) * S[1, 1]))
    a1 = imag(S[2, 1]) / denom
    b1 = imag(S[3, 1]) / denom
    a2 = real((S[2, 2] - S[3, 3]) / (S[2, 2] + S[3, 3]))
    b2 = 2 * real((S[2, 3]) / (S[2, 2] + S[3, 3]))
    return a1, b1, a2, b2
end
"Function to compute the least squares objective for method of moments from Ewans 1998."
function momentLS(θ, S)
    a1h, b1h, a2h, b2h = freqdir_fouriercoef2(S)
    a1, b1, a2, b2 = GaussCoef(θ)
    return (a1 - a1h)^2 + (b1 - b1h)^2 + (a2 - a2h)^2 + (b2 - b2h)^2
end
"Function to optimise moments objective at a given frequency."
function getPars(θ₀, S)
    res = optimize(x -> momentLS(x, S), θ₀)
    if Optim.converged(res)
        return [
            min(res.minimizer[1], res.minimizer[2]),
            max(res.minimizer[1], res.minimizer[2]),
            res.minimizer[3],
        ]
    else
        return fill(NaN, size(θ₀))
    end
end
"Function to compute the peak seperation under a given model at a frequency."
function truepeakseperation(θ, ω, ωₚ)
    β = θ[2]
    ν = θ[3]
    return β / exp(ν * ((ω > ωₚ) * ((ωₚ / ω) - 1.0) + 1.0))
end

"Function to compute the means and angular width of the bimodal spreading at a frequency."
function truemeansigma(θ, ω, ωₚ)
    σₗ = θ[4]
    σᵣ = θ[5]
    φₘ = θ[1]
    φₛ = truepeakseperation(θ, ω, ωₚ) / 2
    return φₘ - φₛ, φₘ + φₛ, (σₗ - (σᵣ / 3) * (4 * (ωₚ / ω)^2 - (ωₚ / ω)^8))
end

"Function to compute the means and angular width of the bimodal spreading at a frequency returning a vector."
function truemeansigmavec(θ, ω, ωₚ)
    φₘ₁, φₘ₂, σ = truemeansigma(θ, ω, ωₚ)
    return [φₘ₁, φₘ₂, σ]
end

##
"Objective function struct for moments. Requires known peak frequency from some other fit."
struct MomentObjective{T}
    moments::Vector{Vector{Float64}}
    Ω::T
    ωₚ::Float64
    function MomentObjective(ts, Δ, θ₀, ωₚ; ωₗ = 0, ωᵤ = Inf)
        B = BartlettPeriodogram(ts, Δ)[end÷2+1:end]
        S = B[ωₗ.<=B.Ω.<=ωᵤ]
        momentest = [getPars(truemeansigmavec(θ₀, ω, ωₚ), s) for (ω, s) ∈ zip(S.Ω, S.ordinates)]
        goodestind = [!isnan(m[1]) for m ∈ momentest]
        new{typeof(S.Ω)}(momentest[goodestind], S.Ω[goodestind], ωₚ)
    end
end
function momentobjective_summand(m, ω, θ, ωₚ)
    φₘ₁, φₘ₂, σ = truemeansigma(θ, ω, ωₚ)
    return (min(φₘ₁, φₘ₂) - m[1])^2 + (max(φₘ₁, φₘ₂) - m[2])^2 + (σ - m[3])^2
end
function (f::MomentObjective)(θ)
    return sum(momentobjective_summand(m, ω, θ, f.ωₚ) for (m, ω) ∈ zip(f.moments, f.Ω))
end

function momentsfit(model, ts, Δ, θ₀; ωₗ = 0, ωᵤ = Inf)
    ## fit LS
    B = BartlettPeriodogram(ts[:, 1], Δ, 128)[end÷2+1:end]
    Bz = B[ωₗ.<=B.Ω.<=ωᵤ]
    res1 = optimize(θ -> marginalbounds(θ,ωₗ,ωᵤ) ? squarediff(model(θ), Bz) : Inf, θ₀[1:4])
    ## fit moments
    if Optim.converged(res1)
        momobj = MomentObjective(ts, Δ, θ₀[5:end], res1.minimizer[2]; ωₗ = ωₗ, ωᵤ = ωᵤ)
        res2 = optimize(θ -> spreadingbounds(θ) ? momobj(θ) : Inf, θ₀[5:end])
        if Optim.converged(res2)
            return (true, [res1.minimizer; res2.minimizer])
        end
    end
    return (false, θ₀)
end
