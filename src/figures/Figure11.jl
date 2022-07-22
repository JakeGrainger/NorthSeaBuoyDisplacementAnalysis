struct Figure11
    θ_est::Vector{Vector{Float64}}
    CI::Vector{Vector{Float64}}
    pardates::Vector{Float64}
    parnames::Vector{String}
    colors::Vector{Symbol}
    plotlayout::Vector{Int}

    Hs::Vector{Float64}
    Hsdates::Vector{Float64}

    power::Matrix{Float64}
    Ω::Vector{Float64}
    specdates::Vector{Float64}

    winddirection::Vector{Float64}
    winddates::Vector{Float64}

    function Figure11()
        @load "dat/anonomous_data.jld"
        @load "dat/wavefits.jld"
        wavetime_days = wavetime ./ 24 ./ 3600
        windtime_days = windtime ./ 24 ./ 3600
        spectrogram, Ω, specdates = PaperAnalysis.compute_spectrogram(waveseries,wavetime_days)

        θ_est = [[f.fits.θ_est[i] for f ∈ fits] for i ∈ 1:9]
        CI = [[1.96 * sqrt(f.fits.varθ[i,i]) for f ∈ fits] for i ∈ 1:9]

        all(f.fits.converged for f in fits) || error("Some fits failed to converge")

        plotlayout = [2,3,5,7,4,6,6,8,8]
        colors = parametercolors()
        parnames = ["\\alpha", "\\omega_p", "\\gamma", "r", "\\phi_m", "\\beta", "\\nu", "\\sigma_l", "\\sigma_r"]
        new(θ_est, CI, wavetime_days[fitstart:fitend], parnames, colors, plotlayout, 
            4std.(waveseries), wavetime_days,
            spectrogram, Ω, specdates,
            winddirection,windtime_days)
    end
end

@recipe function f(f::Figure11;plotascontinuous=true)
    layout := (4,2)
    colorbar := false
    @series begin
        subplot := 1
        label := L"H_s"
        ylims := (1.5,5.5)
        xticks := (2:4,map(i->"day $i",2:4))
        xlims := (2-10/24, 5)
        f.Hsdates, f.Hs
    end
    @series begin
        seriestype := :heatmap
        subplot := 3
        f.specdates, f.Ω, 10 .* log10.(f.power)
    end
    @series begin
        seriescolor := :black
        subplot := 4
        label := "wind direction"
        f.winddates, mod2pi_nan.(deg2rad.(f.winddirection))
    end
    for i ∈ 1:9
        @series begin
            subplot := f.plotlayout[i]
            if f.plotlayout[i] == 3
                ylims := (0,2)
            end
            seriescolor := f.colors[i]
            if plotascontinuous
                ribbon := f.CI[i]
            else
                yerror := f.CI[i]
                seriestype := :scatter
                markersize := 3
                markerstrokecolor := f.colors[i]
                markerstrokewidth := 1
            end
            label := L"%$(f.parnames[i])"
            xticks := (2:4,map(i->"day $i",2:4))
            xlims := (2-10/24, 5)
            if i == 5
                return f.pardates, mod2pi.(f.θ_est[i])
            else
                return f.pardates, f.θ_est[i]
            end
        end
    end
end

