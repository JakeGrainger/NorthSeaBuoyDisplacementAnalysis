"Make constraints for a wind-sea only model."
function makeconstraints(::Type{<:WindSea}, cut)
    lxW = [0, cut, 1, 3, -Inf, 0, 0, 0, 0]
    uxW = [4, 4, Inf, 8, Inf, 2π, Inf, Inf, Inf]
    return TwiceDifferentiableConstraints(lxW, uxW)
end
"Make constraints for a wind-sea only univariate model."
function makeconstraints(::Type{<:JONSWAP}, cut)
    lxW = [0, cut, 1, 3]
    uxW = [4, 4, Inf, 8]
    return TwiceDifferentiableConstraints(lxW, uxW)
end
"Fit with debiased Whittle"
function DWfit(model, ts, Δ, θ₀; ωₗ = 0, ωᵤ = Inf, iterations = 200)
    obj = DebiasedWhittleLikelihood(model, ts, Δ; lowerΩcutoff = ωₗ, upperΩcutoff = ωᵤ)
    constraints = makeconstraints(model, ωₗ)
    res = optimize(
        TwiceDifferentiable(Optim.only_fgh!(obj), θ₀),
        constraints,
        θ₀,
        IPNewton(),
        Optim.Options(iterations = iterations),
    )
    if Optim.converged(res)
        return (true, res.minimizer)
    else
        return (false, θ₀)
    end
end