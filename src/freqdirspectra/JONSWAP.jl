struct FDJONSWAP <: FDSpectralModel
    θ::Vector{Float64}
end

parameter(J::FDJONSWAP) = J.θ

function sdf(J::Union{FDJONSWAP,FDSwellJS}, ω)
    α = J.θ[1]
    ωₚ = J.θ[2]
    γ = J.θ[3]
    r = J.θ[4]
    ω = abs(ω)
    s = 4.0
    σ1² = 0.0049 + 0.0032 * (ω > ωₚ)
    δ = exp(-1 / (2 * σ1²) * (ω / ωₚ - 1)^2)
    sdf = α * ω^(-r) * exp(-(r / s) * (ω / ωₚ)^(-s)) * γ^δ / 2
    sdf *= (!isnan(sdf))
    return sdf
end

function spreading(J::FDJONSWAP, ω, ϕ)
    ω = abs(ω)
    halfPeakSep = J.θ[6] / exp(J.θ[7] * ((ω > J.θ[2])*((J.θ[2] / ω) - 1.0) + 1.0)) / 2
    φₘ₁ = J.θ[5] + halfPeakSep
    φₘ₂ = J.θ[5] - halfPeakSep
    σ = J.θ[8] - (J.θ[9]/3)*(4*(J.θ[2]/ω)^2 - (J.θ[2]/ω)^8)
    return 0.5 * (WrappedGaussian(ϕ, φₘ₁, σ) + WrappedGaussian(ϕ, φₘ₂, σ))
end