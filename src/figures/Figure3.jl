struct Figure3a
    Hs::Vector{Float64}
    dates::Vector{Float64}
    function Figure3a(waveseries,wavedates)
        @assert length(waveseries) == length(wavedates)
        new(4*std.(waveseries),wavedates)
    end
end

struct Figure3b
    speed::Vector{Float64}
    dates::Vector{Float64}
end

struct Figure3c
    power::Matrix{Float64}
    Ω::Vector{Float64}
    dates::Vector{Float64}
    function Figure3c(waveseries,wavedates)
        spectrogram, Ω, specdates = PaperAnalysis.compute_spectrogram(waveseries,wavedates)
        new(spectrogram, Ω, specdates)
    end
end

struct Figure3d
    direction::Matrix{Float64}
    Ω::Vector{Float64}
    dates::Vector{Float64}
    function Figure3d(waveseries,wavedates)
        wavedir = PaperAnalysis.mt_wavedirection(permutedims(reduce(vcat, wave for wave in waveseries)),2304,2304÷2)
        specdates = wavedates[1]:(wavedates[2]-wavedates[1])/2:wavedates[end]
        Ω = collect((2π*1.28).*wavedir.freq)
        dir = wavedir.power
        new(mod2pi.(dir),Ω,specdates)
    end
end

struct Figure3e
    direction::Vector{Float64}
    dates::Vector{Float64}
    function Figure3e(direction,dates)
        @assert length(direction) == length(dates)
        new(mod2pi_nan.(deg2rad.(direction)), dates)
    end
end


struct Figure3
    fig3a::Figure3a
    fig3b::Figure3b
    fig3c::Figure3c
    fig3d::Figure3d
    fig3e::Figure3e
    function Figure3()
        @load "dat/anonomous_data.jld"
        new(
            Figure3a(waveseries, wavetime./24 ./3600),
            Figure3b(windspeed, windtime./24 ./3600),
            Figure3c(waveseries, wavetime./24 ./3600),
            Figure3d(waveseries, wavetime./24 ./3600),
            Figure3e(winddirection, windtime./24 ./3600)
        )
    end
end

## plotting
@recipe function f(f::Figure3a)
    yguide := L"$H_s$ (m)"
    @series begin
        seriestype := :heatmap
        xlims := (f.dates[1],f.dates[end])
        xticks := (1:5,map(i->"day $i",1:5))
        f.dates[end] .+ [2/24,3/24], 2:3, (x,y) -> 1 
    end
    @series begin
        label := false
        seriescolor := 1
        xlims := (f.dates[1],f.dates[end])
        xticks := (1:5,map(i->"day $i",1:5))
        ylims := (1.3,5.2)
        f.dates, f.Hs
    end 
end
@recipe function f(f::Figure3b)
    yguide := "wind speed (m/s)"
    @series begin    
        seriestype := :heatmap
        xlims := (f.dates[1],f.dates[end])
        xticks := (1:5,map(i->"day $i",1:5))
        f.dates[end] .+ [2/24,3/24], 2:3, (x,y) -> 1 
    end
    @series begin
        label := false
        xlims := (f.dates[1],f.dates[end])
        xticks := (1:5,map(i->"day $i",1:5))
        seriescolor := 2
        f.dates, f.speed
    end
end
@recipe function f(f::Figure3c)
    yguide := "angular frequency"
    seriestype := :heatmap
    xticks := (1:5,map(i->"day $i",1:5))
    f.dates, f.Ω, 10.0 .* log10.(f.power)
end
@recipe function f(f::Figure3d)
    yguide := "angular frequency"
    seriestype := :heatmap
    seriescolor := :hsv
    xticks := (1:5,map(i->"day $i",1:5))
    f.dates, f.Ω, f.direction
end

@recipe function f(f::Figure3e)
    @series begin
        seriestype := :heatmap
        seriescolor := :hsv
        f.dates, range(0,2π,length=100),(x,y)->y
    end
    @series begin
        yguide := "wind direction"
        label := false
        seriescolor := :black
        xlims := (f.dates[1], f.dates[end])
        xticks := (1:5,map(i->"day $i",1:5))
        f.dates, f.direction
    end
end

@recipe function f(f::Figure3)
    layout := (fieldcount(typeof(f)),1)
    for (i,field) in enumerate(fieldnames(typeof(f)))
        @series begin
            subplot := i
            getfield(f,field)
        end
    end
end