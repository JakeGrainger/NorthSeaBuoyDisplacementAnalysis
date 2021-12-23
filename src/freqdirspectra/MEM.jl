struct MEM{T,S}
    Ω::T
    Φ::S
    fdspec::Matrix{Float64}
    D::Matrix{Float64}
    function MEM(S, Ω, Φ)
        D = zeros(length(Φ), length(Ω)) # transpose on return because should be frequency direction
        fdspec = similar(D)
        for (i, s) ∈ enumerate(S)
            a1, b1, a2, b2 = freqdir_fouriercoef(s)
            c1 = a1 + 1im * b1
            c2 = a2 + 1im * b2
            F₁ = (c1 - c2 * conj(c1)) / (1 - abs(c1)^2)
            F₂ = c2 - c1 * F₁

            @. D[:, i] = real(
                (1 / (2pi)) * (1 - F₁ * conj(c1) - F₂ * conj(c2)) /
                (abs(1 - F₁ * exp(-1im * Φ) - F₂ * exp(-2im * Φ))^2),
            )

            for j = 1:length(Φ)
                fdspec[j, i] = D[j, i] * real.(s[1, 1])
            end
        end
        new{typeof(Ω), typeof(Φ)}(Ω, Φ, fdspec', D')
    end  
end
function freqdir_fouriercoef(S)
    denom = sqrt(abs((S[2, 2] + S[3, 3]) * S[1, 1]))
    a1 = imag(S[2,1]) / denom
    b1 = imag(S[3,1]) / denom
    a2 = real((S[2, 2] - S[3, 3]) / (S[2, 2] + S[3, 3]))
    b2 = 2 * real((S[2, 3]) / (S[2, 2] + S[3, 3]))
    return a1, b1, a2, b2
end
