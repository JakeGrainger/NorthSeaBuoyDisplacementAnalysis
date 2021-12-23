struct Figure4
    θ::Vector{Float64}
    function Figure4(θ=[0.7,0.7,3.3,5,pi,4,2.7,0.55,0.26])
        new(θ)
    end
end

@recipe function f(f::Figure4)
    Ω = range(0,pi,length=200)
    Φ = range(0,2pi,length=200)
    model = FrequencyDirectionSpectra.FDJONSWAP(f.θ)
    layout := (1,3)
    @series begin
        subplot := 1
        xguide := "angular frequency"
        yguide := "direction"
        seriestype := :contourf
        seriescolor := :thermal
        levels := 30
        Ω, Φ, (ω,ϕ) -> FrequencyDirectionSpectra.freqdirspectra(model,ω,ϕ)
    end
    @series begin
        subplot := 2
        label := false
        xguide := "angular frequency"
        yguide := "spectral density"
        range(0,pi,length=1000), ω -> FrequencyDirectionSpectra.sdf(model, ω)
    end
    @series begin
        subplot := 3
        xguide := "angular frequency"
        yguide := "direction"
        seriestype := :contourf
        seriescolor := :acton
        levels := 30
        Ω, Φ, (ω,ϕ) -> FrequencyDirectionSpectra.spreading(model,ω,ϕ)
    end
end
