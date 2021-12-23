module PaperFigures
    using   OceanWaveSpectralFitting,
            RecipesBase,
            LaTeXStrings,
            DSP,
            Statistics,
            JLD2,
            Dates
    import OceanWaveSpectralFitting.WhittleLikelihoodInference: HermitianPlot
    include("../freqdirspectra/FrequencyDirectionSpectra.jl")
    include("../exampleanalysis/PaperAnalysis.jl")
    include("../simulation/PaperSimulation.jl")
    include("generic.jl")
    for i in 1:11
        include("Figure$(i).jl")
    end
end