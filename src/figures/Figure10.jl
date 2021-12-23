struct Figure10
    lowcutoff::Vector{Float64}
    highcutoff::Float64
    cutdates::Vector{Float64}
    power::Matrix{Float64}
    Ω::Vector{Float64}
    specdates::Vector{Float64}
    function Figure10()
        @load "dat/anonomous_data.jld"
        @load "dat/wavefits.jld"
        wavetime_days = wavetime ./ 24 ./ 3600
        spectrogram, Ω, specdates = PaperAnalysis.compute_spectrogram(waveseries,wavetime_days)
        new([f.ini.cutoff for f in fits], 3.8, wavetime_days[fitstart:fitend], spectrogram, Ω, specdates)
    end
end

@recipe function f(f::Figure10)
    xlims := (f.specdates[1],f.specdates[end])
    ylims := (0,1.28π)
    xticks := (1:5,map(i->"day $i",1:5))
    @series begin
        seriestype := :heatmap
        f.specdates, f.Ω, 10 .* log10.(f.power)
    end
    @series begin
        linestyle := :dot
        seriescolor := :black
        label := false
        linewidth := 1.3
        f.cutdates, _-> f.highcutoff
    end
    @series begin
        linestyle := :dot
        seriescolor := :black
        label := false
        linewidth := 1.3
        f.cutdates, f.lowcutoff
    end
    @series begin
        linestyle := :dash
        seriescolor := :black
        label := false
        linewidth := 1.3
        seriestype := :vline
        f.cutdates[[1,end]]
    end
end