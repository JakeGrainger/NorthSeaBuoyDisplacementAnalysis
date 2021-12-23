struct FDSwell <: FDSpectralModel
    θ::Vector{Float64}
end

parameter(F::FDSwell) = F.θ

function sdf(F::FDSwell, ω)
    θ = F.θ
    α = θ[1]
    ωₚ = θ[2]
    κ = θ[3]
    sigω = sign(ω)
    ω = abs(ω)
    logbit = log(ω) - log(ωₚ)
    h = logbit/κ - κ
    sdf = α * exp(-0.5 * (h)^2)
    sdf *= (!isnan(sdf))
end

spreading(F::FDSwell, ω, ϕ) = WrappedGaussian(ϕ, F.θ[4], F.θ[5])


####

struct FDSwellJS <: FDSpectralModel
    θ::Vector{Float64}
end

parameter(F::FDSwellJS) = F.θ

spreading(F::FDSwellJS, ω, ϕ) = WrappedGaussian(ϕ, F.θ[5], F.θ[6])