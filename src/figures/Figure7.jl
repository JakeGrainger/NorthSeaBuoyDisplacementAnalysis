struct Figure7{T}
    sim::Vector{T}
    parameternames::Vector{String}
    truncation::Vector{Tuple{Float64,Float64}}
    truncate::Bool
    function Figure7(;truncate=false)
        @load "dat/simulation.jld"
        truncation = [
            (-0.7,1.3), 
            (-0.1,0.1),
            (-1,2.7),
            (-2,3)
        ]
        new{eltype(sim)}(sim,makeparameternames(),truncation,truncate)
    end
end

@recipe function f(f::Figure7)
    count = 1
    goodindex = [s.DWflag .& s.DWuniflag .& s.LSmlmflag .& s.LSmemflag .& s.Mflag for s in f.sim]
    for (i,g) in enumerate(goodindex)
        if sum(g) != length(g)
            @warn "From scenario $i, some results failed. They were
                DW:     $(length(g)-sum(f.sim[i].DWflag))
                DW uni: $(length(g)-sum(f.sim[i].DWuniflag))
                LS mlm: $(length(g)-sum(f.sim[i].LSmlmflag))
                LS mem: $(length(g)-sum(f.sim[i].LSmemflag))
                moment: $(length(g)-sum(f.sim[i].Mflag))
            As a result, only $(sum(g)) results are used in the figure for Scenario $i.
            "
        end
    end
    color_palette := Main.Plots.palette([:royalblue, :green, :darkorange])
    layout := (4, 3)
    x = ["LS" "DW uni" "DW"]
    for i ∈ 1:4, k ∈ 1:3
        y = [getindex.(getfield(f.sim[k],t),i)[goodindex[k]] for t in [:LSmlm, :DWuni, :DW]]
        # @series begin
        #     seriestype := :violin
        #     subplot := count
        #     label := false
        #     x, y
        # end
        @series begin
            seriestype := :boxplot
            subplot := count
            label := false
            x, y
        end
        @series begin
            seriestype := :hline
            subplot := count
            seriescolor := :red
            linestyle := :dash
            label := false
            if k == 1
                yguide := L"%$(f.parameternames[i])"
                yguide_position := :left
            end
            if count ∈ 1:3
                xguide := "Scenario $k"
                xguide_position := :top
            end
            if f.truncate
                ylims := f.sim[k].θ_true[i] .+ f.truncation[i]
            end
            f.sim[k].θ_true[i:i]
        end
        count += 1
    end
end
