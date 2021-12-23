struct Figure9
    error::Matrix{Float64}
    Ω::Vector{Float64}
    dates::Vector{Float64}
    function Figure9()
        @load "dat/anonomous_data.jld"
        wavetime = wavetime ./24 ./3600
        errorgram, Ω, spectime = PaperAnalysis.compute_errorgram(waveseries,wavetime)
        new(errorgram, Ω, spectime)
    end
end

@recipe function f(f::Figure9)
    seriestype := :heatmap
    seriescolor := :balance
    clims := (-12,12)
    xticks := (1:5,map(i->"day $i",1:5))
    f.dates, f.Ω, f.error
end