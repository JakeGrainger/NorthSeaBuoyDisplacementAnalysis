struct MLM{T,S}
    Ω::T
    Φ::S
    fdspec::Matrix{Float64}
    D::Matrix{Float64}
    function MLM(S, Ω, Φ)
        length(S) == length(Ω) || throw(ArgumentError("Ω and S are not the same size"))
        Dtr = zeros(length(Φ), length(Ω)) # will transpose at the end
        G = zeros(ComplexF64,1,3)
        for (i, (ω,s)) ∈ enumerate(zip(Ω,S))
            for (j,ϕ) ∈ enumerate(Φ)
                allocateG!(G,ω,ϕ)
                Dtr[j, i] = 1.0/real((G) * inv(s) * (G'))[1]
            end
            k = length(Φ)/(sum(Dtr[:, i]) * 2pi)
            Dtr[:, i] .*= k
        end
        Dout = permutedims(Dtr)
        fdspec = zeros(length(Ω), length(Φ))
        for j ∈ 1:length(Φ), (i,s) ∈ enumerate(S)
            fdspec[i, j] = Dout[i, j] * real(s[1, 1])
        end
        new{typeof(Ω), typeof(Φ)}(Ω, Φ, fdspec, Dout)
    end
end

function allocateG!(G,ω,ϕ)
    G[1] = 1.0 
    G[2] = -1im * cos(ϕ)
    G[3] = -1im * sin(ϕ)
    return nothing
end
