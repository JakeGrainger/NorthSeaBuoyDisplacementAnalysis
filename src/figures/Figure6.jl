struct Figure6
    θ1::Vector{Float64}
    θ2::Vector{Float64}
    θ3::Vector{Float64}
    function Figure6(
            θ1 = [0.7, 0.8, 3.3, 5, π / 2, 4, 2.7, 0.55, 0.26],
            θ2 = [0.7, 1.1, 3.3, 5, π / 2, 4, 2.7, 0.55, 0],
            θ3 = [0.7, 1.0, 1, 5, π / 2, 4, 2.7, 0.55, 0.26]
        )
        new(θ1,θ2,θ3)
    end
end

@recipe function f(f::Figure6)
    layout := (1,3)
    for (i,θ) in enumerate((f.θ1,f.θ2,f.θ3))
        @series begin
            ylims --> (0.4,1.9)
            yticks --> (0.6:0.3:1.5)
            clims --> (0,1.25)
            levels --> 50
            seriescolor := :thermal
            subplot := i
            FrequencyDirectionSpectra.FDJONSWAP(θ),true
        end
    end
end

##
struct ColorFig6
    fig6::Figure6
end
@recipe function f(fc::ColorFig6)
    legend := true
    layout := 1
    ylims --> (0.4,1.9)
    yticks --> (0.6:0.3:1.5)
    clims --> (0,1.25)
    levels --> 50
    seriescolor := :thermal
    FrequencyDirectionSpectra.FDJONSWAP(fc.fig6.θ1),true
end

