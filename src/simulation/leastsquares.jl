##
"Compute the square difference between a model spreading and an estimate."
function squarediff(model::FrequencyDirectionSpectra.FDSpectralModel, fdest::Union{FrequencyDirectionSpectra.MEM,FrequencyDirectionSpectra.MLM})
    return squarediff(model, fdest.D, fdest.Ω, fdest.Φ)
end
"Compute the square difference between a model spreading and an estimate seperated into ordinates and frequencies."
function squarediff(model::FrequencyDirectionSpectra.FDSpectralModel, D, Ω, Φ)
    return sum((FrequencyDirectionSpectra.spreading(model, Ω[i], Φ[j]) - D[i, j])^2 for i ∈ 1:size(D, 1), j ∈ 1:size(D, 2))
end
"Compute the square difference between a model sdf and an estimate."
function squarediff(model::FrequencyDirectionSpectra.FDSpectralModel,B::BartlettPeriodogram)
    return squarediff(model, B.Ω, B.ordinates)
end
"Compute the square difference between a model sdf and an estimate seperated into ordinates and frequencies."
function squarediff(model::FrequencyDirectionSpectra.FDSpectralModel, Ω, ordinates)
    return sum((FrequencyDirectionSpectra.sdf(model, ω) - I)^2 for (ω, I) ∈ zip(Ω, ordinates))
end

##
function LSfit(
    model,
    ts,
    Δ,
    θ₀;
    fdest = FrequencyDirectionSpectra.MLM,
    ωₗ = 0,
    ωᵤ = Inf,
    window = 128,
    ndirection = 100,
)
    B = BartlettPeriodogram(ts, Δ, window)[end÷2+1:end]
    Bused = B[ωₗ.<=B.Ω.<=ωᵤ]
    S = fdest(Bused.ordinates, Bused.Ω, range(-π, π, length = ndirection))
    Bz = BartlettPeriodogram(ts[:, 1], Δ, window)[end÷2+1:end][ωₗ.<=B.Ω.<=ωᵤ]
    res1 = optimize(θ -> marginalbounds(θ,ωₗ,ωᵤ) ? squarediff(model(θ), Bz) : Inf, θ₀[1:4])
    if Optim.converged(res1)
        res2 = optimize(θ -> spreadingbounds(θ) ? squarediff(model([res1.minimizer; θ]), S) : Inf, θ₀[5:end])
        if Optim.converged(res2)
            return (true, [res1.minimizer; res2.minimizer])
        end
    end
    return (false, θ₀)
end


##
function marginalbounds(θ,ωₗ,ωᵤ)
    return 4>θ[1]>0 && ωᵤ>θ[2]>ωₗ && 20>θ[3]>1 && 8>θ[4]>3
end
function spreadingbounds(θ)
    return 2π>θ[2]>0 && θ[3]>0 && θ[4]>0 && θ[5]>0 
end