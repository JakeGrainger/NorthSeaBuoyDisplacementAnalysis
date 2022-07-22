struct Figure5
    colors::Vector{Symbol}
    styles::Vector{Symbol}
    names::Vector{String}
    θ_default::Vector{Float64}
    α::Vector{Float64}
    ωₚ::Vector{Float64}
    γ::Vector{Float64}
    r::Vector{Float64}
    ϕₘ::Vector{Float64}
    β::Vector{Float64}
    ν::Vector{Float64}
    σₗ::Vector{Float64}
    σᵣ::Vector{Float64}
    function Figure5(;
            colors = parametercolors(),
            styles = [:dot, :dash, :dashdot, :solid],
            names = ["\\alpha", "\\omega_p", "\\gamma", "r", "\\omega_p", "\\beta", "\\nu", "\\sigma_l", "\\sigma_r"],
            θ_default = [0.7,0.8,3.3,5,pi,4,2.7,0.55,0.26],
            α = collect(0.4:0.2:1.0),
            ωₚ = collect(0.8:0.1:1.1),
            γ = [1,2.5,4,5.5],
            r = collect(4:0.5:5.5),
            ϕₘ = collect(range(π/2,3π/2,length=4)),
            β = collect(3:6),
            ν = collect(1:4),
            σₗ = collect(0.45:0.05:0.6),
            σᵣ = collect(0.15:0.05:0.3)
        )
        @assert all(length(styles) == length(a) for a in (α, ωₚ, γ, r, ϕₘ, β, ν, σₗ, σᵣ))
        new(colors,styles,names,θ_default,α, ωₚ, γ, r, ϕₘ, β, ν, σₗ, σᵣ)
    end
end
ϕₛ(ω,θ) = θ[6]*exp(-θ[7] * (abs(ω) > θ[2] ? θ[2]/abs(ω) : 1))
σ(ω,θ) = θ[8] - θ[9]/3 * (4*(θ[2]/abs(ω))^2 - (θ[2]/abs(ω))^8)
function makefulltheta(θ_default,val,ind)
    θ = copy(θ_default)
    θ[ind] = val
    return θ
end

@recipe function f(f::Figure5)
    layout := (3,3)
    Ω = range(0,2,length=200)
    for (i,parameter) ∈ enumerate([:α,:ωₚ,:γ,:r]), (j,par) ∈ enumerate(getfield(f,parameter))
        @series begin
            subplot := i
            seriescolor := f.colors[i]
            linestyle := f.styles[j]
            label := L"%$(f.names[i]) = %$par"
            if i < 4
                ylims := (0,1.8)
                yguide := L"f(\omega)"
            elseif i == 4
                yguide := L"10\,\log_{10}\;f(\omega)"
                yaxis := :log10
                ylims := (0.001,1.2^10)
                return Ω[2:end], ω -> (real(sdf(JONSWAP{0}(makefulltheta(f.θ_default,par,i)[1:4]),ω))) ^ 10 # for 10log10
            end
            return Ω, ω -> real(sdf(JONSWAP{0}(makefulltheta(f.θ_default,par,i)[1:4]),ω))
        end
    end
    for (i,parameter) ∈ enumerate([:ϕₘ,:β,:ν]), (j,par) ∈ enumerate(getfield(f,parameter))
        @series begin
            subplot := i+4
            seriescolor := f.colors[i+4]
            linestyle := f.styles[j]
            label := L"%$(f.names[i+4]) = %$par"
            Ω, ω -> makefulltheta(f.θ_default, par, i+4)[5] + ϕₛ(ω,makefulltheta(f.θ_default, par, i+4))/2
        end
        @series begin
            subplot := i+4
            ylims := (pi-2.1,pi+2.1)
            seriescolor := f.colors[i+4]
            linestyle := f.styles[j]
            label := false
            yguide := L"$\phi_{m2}(\omega)$ and $\phi_{m1}(\omega)$"
            if i == 3
                xguide := "angular frequency"
            end
            Ω, ω -> makefulltheta(f.θ_default, par, i+4)[5] - ϕₛ(ω,makefulltheta(f.θ_default, par, i+4))/2
        end
    end
    for (i,parameter) ∈ enumerate([:σₗ,:σᵣ]), (j,par) ∈ enumerate(getfield(f,parameter))
        @series begin
            subplot := i+7
            ylims := (0.08,1)
            seriescolor := f.colors[i+7]
            linestyle := f.styles[j]
            label := L"%$(f.names[i+7]) = %$par"
            yguide := L"\sigma(\omega)"
            xguide := "angular frequency"
            Ω, ω -> σ(ω,makefulltheta(f.θ_default,par,i+7))
        end
    end
end
