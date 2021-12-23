@recipe function f(p::FrequencyDirectionSpectra.FDSpectralModel,compass=false;Δ=1,func=freqdirspectra)
    projection --> :polar
    seriestype --> :contour
    legend --> false
    Ω = range(0, π/Δ, length = 100)
    Φ = range(-pi, pi, length = 100)
    if compass
        xticks := (0:π/2:3π/2, ["E", "N", "W", "S"])
        return Φ, Ω, (ϕ,ω) -> func(p,ω,-ϕ+pi/2)
    else
        xticks := (0:π/4:7π/4, ["0", "π/4", "π/2", "3π/4", "π", "-3π/4", "-π/2", "-π/4"])
        return Φ, Ω, (ϕ,ω) -> func(p,ω,ϕ)
    end
end