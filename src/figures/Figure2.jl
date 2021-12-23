struct Figure2a
    θ::Vector{Float64}
    ϕₘ::Vector{Float64}
    palette::Vector{Symbol}
    function Figure2a(θ=[0.7,1.0,3.3,5, 0, 4,2.7,0.55,0.26],ϕₘ=0:π/4:4π/4,palette=[:Greens_9, :Greys_9, :Reds_9, :Blues_9])
        new(θ,ϕₘ,palette)
    end
end

struct Figure2b
    θ::Vector{Float64}
    ϕₘ::Vector{Float64}
    palette::Vector{Symbol}
    function Figure2b(θ=[0.7,1.0,3.3,5, 0, 4,2.7,0.55,0.26],ϕₘ=0:π/4:4π/4,palette=[:Greens_9, :Greys_9, :Reds_9, :Blues_9])
        new(θ,ϕₘ,palette)
    end
end

struct Figure2
    figa::Figure2a
    figb::Figure2b
    function Figure2(θ=[0.7,1.0,3.3,5, 0, 4,2.7,0.55,0.26],ϕₘ=0:π/4:4π/4,palette=[:Greens_9, :Greys_9, :Reds_9, :Blues_9])
        new(Figure2a(θ,ϕₘ,palette),Figure2b(θ,ϕₘ,palette))
    end
end

@recipe function f(f::Figure2a)
    for (dir,pal) ∈ zip(f.ϕₘ,f.palette)
        @series begin
            seriescolor := pal
            ylims := 0.5,1.7
            FrequencyDirectionSpectra.FDJONSWAP([f.θ[1:4];dir;f.θ[6:end]]), true
        end
    end
end

@recipe function f(f::Figure2b)
    Ω = range(0,π,length = 500)
    for (dir,pal) ∈ zip(f.ϕₘ, f.palette)
        @series begin
            seriescolor := pal
            HermitianPlot(Ω,sdf.(JS_BWG_HNE{0}([f.θ[1:4];dir;f.θ[6:end]]),Ω))
        end
    end
    for (i,(lab,part)) in enumerate(zip(["zz", "zx", "zy", "xz", "xx", "xy", "yz", "yx", "yy"],["","Im","Im","Re","","Im","Re","Re",""]))
        @series begin
            subplot := i
            if mod(i,3) != 1
                yformatter := _->""
            end
            if i < 7
                xformatter := _->""
            end
            label := false
            series_annotations := Main.Plots.series_annotations([L"%$(part)\;f_{%$(lab)}(\omega)"], Main.Plots.font(10))
            [2.3], [0.21]
        end
    end
end

