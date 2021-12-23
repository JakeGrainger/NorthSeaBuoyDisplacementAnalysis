function WrappedGaussian(x, μ, σ)
    p = 0
    K = 10
    for k = -K:K
        p += exp(-0.5 * ((x - μ - 2pi * k) / σ)^2)
    end
    return p / (sqrt(2π) * σ)
end

function ω2k(ω, h)
    g = 9.81
    α = ω^2*(h/g)
    out = α*(tanh(α))^(-0.5)/h
    return out .* !isnan(out)
end