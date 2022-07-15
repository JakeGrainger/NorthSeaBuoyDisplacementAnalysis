struct Figure1a
    θ₁::Vector{Float64}
    function Figure1a(x=[0.7,0.8,3.3,5,π,4,2.7,0.55,0.26])
        new(x)
    end
end
struct Figure1b
    θ₂::Vector{Float64}
    function Figure1b(x= vcat([0.7,0.8,3.3,5,π,4,2.7,0.55,0.26], [0.2,0.45,0.05,2π/3,0.2]))
        new(x)
    end
end
struct Figure1c
    θ₃::Vector{Float64}
    function Figure1c(x = vcat([0.7,1.2,3.3,5,-π/4,4,2.7,0.55,0.26], [0.1,0.5,0.05,π/3,0.1], [0.2,0.6,0.05,5π/3,0.1]))
        new(x)
    end
end
struct Figure1{S}
    fig1a::Figure1a
    fig1b::Figure1b
    fig1c::Figure1c
    col::S
    function Figure1(col=:speed)
        new{typeof(col)}(Figure1a(),Figure1b(),Figure1c(),col)
    end
end

@recipe function f(f::Figure1a)
    ylims --> (0.2,1.8)
    legend --> true
    colorbar --> :right
    clims --> (0,1.25)
    colorbar_ticks --> (0:0.25:1.25)
    levels --> 50
    FrequencyDirectionSpectra.FDJONSWAP(f.θ₁),true
end

@recipe function f(f::Figure1b)
    ylims --> (0.2,1.8)
    legend --> true
    colorbar --> :right
    clims --> (0,1.25)
    colorbar_ticks --> (0:0.25:1.25)
    levels --> 50
    FrequencyDirectionSpectra.FDJONSWAP(f.θ₂[1:9]) + FrequencyDirectionSpectra.FDSwell(f.θ₂[10:end]),true
end

@recipe function f(f::Figure1c)
    ylims --> (0.2,1.8)
    legend --> true
    colorbar --> :right
    clims --> (0,1.25)
    colorbar_ticks --> (0:0.25:1.25)
    levels --> 50
    FrequencyDirectionSpectra.FDJONSWAP(f.θ₃[1:9]) + FrequencyDirectionSpectra.FDSwell(f.θ₃[10:14]) + FrequencyDirectionSpectra.FDSwell(f.θ₃[15:end]),true
end

@recipe function f(f::Figure1)
    ylims --> (0.2,1.8)
    layout := (1,3)
    clims --> (0,1.25)
    levels --> 50
    seriescolor := f.col
    @series begin
        subplot := 1
        seriestype := :contour
        FrequencyDirectionSpectra.FDJONSWAP(f.fig1a.θ₁),true
    end
    @series begin
        subplot := 2
        seriestype := :contour
        FrequencyDirectionSpectra.FDJONSWAP(f.fig1b.θ₂[1:9]) + FrequencyDirectionSpectra.FDSwell(f.fig1b.θ₂[10:end]),true
    end
    @series begin
        subplot := 3
        seriestype := :contour
        FrequencyDirectionSpectra.FDJONSWAP(f.fig1c.θ₃[1:9]) + FrequencyDirectionSpectra.FDSwell(f.fig1c.θ₃[10:14]) + FrequencyDirectionSpectra.FDSwell(f.fig1c.θ₃[15:end]),true
    end
end

struct ColorFig1
    fig1::Figure1
end

@recipe function f(fc::ColorFig1)
    ylims --> (0.2,1.8)
    layout := 1
    clims --> (0,1.25)
    levels --> 50
    legend := true
    seriescolor := fc.fig1.col
    seriestype := :contour
    FrequencyDirectionSpectra.FDJONSWAP(fc.fig1.fig1a.θ₁),true
end