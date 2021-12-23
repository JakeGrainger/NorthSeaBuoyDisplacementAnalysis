module FrequencyDirectionSpectra

    using LinearAlgebra, RecipesBase, Distributions

    include("typestructure.jl")
    include("common.jl")
    include("swell.jl")
    include("JONSWAP.jl")
    include("MLM.jl")
    include("MEM.jl")
    include("plotting.jl")

end
