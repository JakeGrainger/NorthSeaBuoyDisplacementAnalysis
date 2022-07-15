struct Figure4{T,S1,S2,S3,S4}
    θ::Vector{Float64}
    levels::T
    col1::S1
    col2::S2
    col3::S3
    linecol::S4
    function Figure4(θ=[0.7,0.7,3.3,5,pi,4,2.7,0.55,0.26];levels=20,col1=:thermal,col2=1,col3=:acton,linecol=:acton)
        new{typeof(levels),typeof(col1),typeof(col2),typeof(col3),typeof(linecol)}(θ,levels,col1,col2,col3,linecol)
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
        seriescolor := f.col1
        levels := f.levels
        Ω, Φ, (ω,ϕ) -> FrequencyDirectionSpectra.freqdirspectra(model,ω,ϕ)
    end
    @series begin
        subplot := 1
        xguide := "angular frequency"
        yguide := "direction"
        seriestype := :contour
        seriescolor := f.linecol
        linewidth := 0.1
        levels := f.levels
        Ω, Φ, (ω,ϕ) -> FrequencyDirectionSpectra.freqdirspectra(model,ω,ϕ)
    end
    @series begin
        subplot := 2
        label := false
        xguide := "angular frequency"
        yguide := "spectral density"
        serieswidth := 2
        seriescolor := f.col2
        range(0,pi,length=1000), ω -> FrequencyDirectionSpectra.sdf(model, ω)
    end
    @series begin
        subplot := 3
        xguide := "angular frequency"
        yguide := "direction"
        seriestype := :contourf
        seriescolor := f.col3
        levels := 20
        Ω, Φ, (ω,ϕ) -> FrequencyDirectionSpectra.spreading(model,ω,ϕ)
    end
    @series begin
        subplot := 3
        xguide := "angular frequency"
        yguide := "direction"
        seriestype := :contour
        seriescolor := f.linecol
        linewidth := 0.1
        levels := 20
        Ω, Φ, (ω,ϕ) -> FrequencyDirectionSpectra.spreading(model,ω,ϕ)
    end
end
