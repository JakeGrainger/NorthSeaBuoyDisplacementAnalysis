function parametercolors()
    α_color = :red
    ωₚ_color = :navy
    γ_color = :orange
    r_color = :purple
    ϕₘ_color = :blue
    β_color = :green
    ν_color = :magenta
    σₗ_color = :maroon
    σᵣ_color = :teal
    return [α_color,
            ωₚ_color,
            γ_color,
            r_color,
            ϕₘ_color,
            β_color,
            ν_color,
            σₗ_color,
            σᵣ_color]
end
function makeparameternames()
    return [
        "\\alpha", 
        "\\omega_p",
        "\\gamma",
        "r",
        "\\phi_m",
        "\\beta",
        "\\nu",
        "\\sigma_l",
        "\\sigma_r"
    ]
end

mod2pi_nan(x) = isnan(x) ? x : mod2pi(x)