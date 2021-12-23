# Scenario 1, peak of the storm
θ₁ = [0.7, 0.8, 3.3, 5, π / 2, 4, 2.7, 0.55, 0.26]
# Scenario 2, higher frequency constant spread
θ₂ = [0.7, 1.1, 3.3, 5, π / 2, 4, 2.7, 0.55, 0]
# Scenario 3, PM style
θ₃ = [0.7, 1.0, 1, 5, π / 2, 4, 2.7, 0.55, 0.26]
#
Θ = [θ₁, θ₂, θ₃]

cutoff = [0.6, 0.75, 0.7]

##
include("src/simulation/PaperSimulation.jl")
import Random
Random.seed!(1234)

## Check on small example
@time sim_small = [
    PaperSimulation.SimulationStudy(1, θᵢ; ωₗ = c, name = "Scenario $i ") for
    ((i, θᵢ), c) ∈ zip(enumerate(Θ), cutoff)
]

## run simulation study
@time sim = [
    PaperSimulation.SimulationStudy(1000, θᵢ; ωₗ = c, name = "Scenario $i ") for
    ((i, θᵢ), c) ∈ zip(enumerate(Θ), cutoff)
]
using JLD2
@save "dat/simulation.jld" sim