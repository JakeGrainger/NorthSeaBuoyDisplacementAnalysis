## code for figure gen
include("src/simulation/PaperSimulation.jl")
include("src/exampleanalysis/PaperAnalysis.jl")
include("src/figures/PaperFigures.jl")
using Plots, StatsPlots, Dates, OceanWaveSpectralFitting
theme(:vibrant)
fontsize = 10
default(xtickfontsize = fontsize, ytickfontsize = fontsize, colorbar_tickfontsize = fontsize, legendfontsize = fontsize, guidefontsize = fontsize)
pyplot(dpi = 300)
foldername = "fig/"
isdir(foldername) || mkdir(foldername)
@info "Loaded packages"

## Figure 1 - example S(ω,ϕ)
pyplot(size = (1000,300))
plot(PaperFigures.Figure1(cgrad(:speed,rev=true)), colorbar = false)
savefig(foldername * "figure1_nocolorbar.png")
plot(PaperFigures.ColorFig1(PaperFigures.Figure1(cgrad(:speed,rev=true))))
savefig(foldername * "figure1_colorbar.png")
@info "Finished figure 1"

## Figure 2 - S(ω,ϕ) vs f(ω)
pyplot(size = (1000,500))
plot(plot(PaperFigures.Figure2a()), plot(PaperFigures.Figure2b()), layout = (1,2))
savefig(foldername * "figure2.png")
@info "Finished figure 2"

## Figure 3 - buoy data summary
pyplot(size = (1000,1000))
plot(PaperFigures.Figure3())
savefig(foldername * "figure3.png")
@info "Finished figure 3"

## Figure 4 - S(ω,ϕ), f(ω), D(ω,ϕ)
pyplot(size = (1000,250))
plot(PaperFigures.Figure4(levels=[0;10.0.^(-14:1:-1);0.2:0.1:2.1],
    col3=cgrad(:speed,rev=true),
    col2=:black,
    col1=:speed,
    linecol=palette([:black,:black],2)),
    size=(800,200))
savefig(foldername * "figure4.png")
@info "Finished figure 4"

## Figure 5 - Vary parameters
pyplot(size = (1000,600))
plot(PaperFigures.Figure5(),legend=false)
savefig(foldername * "figure5.png")
@info "Finished figure 5"

## Figure 6 - Scenario S(ω,ϕ)
pyplot(size = (1000,300))
plot(PaperFigures.Figure6(col=cgrad(:speed,rev=true)))
savefig(foldername * "figure6_nocolorbar.png")
plot(PaperFigures.ColorFig6(PaperFigures.Figure6(col=cgrad(:speed,rev=true))))
savefig(foldername * "figure6_colorbar.png")
@info "Finished figure 6"

## Figure 7 - Marginal parameter estimates (simulation)
pyplot(size = (1000,800))
plot(PaperFigures.Figure7(truncate=true))
savefig(foldername * "figure7.png")
plot(PaperFigures.Figure7(truncate=false))
savefig(foldername * "figure7untruncated.png")
@info "Finished figure 7"

## Figure 8 - Directional parameter estimates (simulation)
pyplot(size = (1000,1000))
plot(PaperFigures.Figure8(truncate=true))
savefig(foldername * "figure8.png")
plot(PaperFigures.Figure8(truncate=false))
savefig(foldername * "figure8untruncated.png")
@info "Finished figure 8"

## Figure 9 - Heatmap of ̂R(ω,t)
pyplot(size = (1000,250))
plot(PaperFigures.Figure9())
savefig(foldername * "figure9.png")
@info "Finished figure 9"

## Figure 10 - Spectrogram with cutoffs
pyplot(size = (1000,250))
plot(PaperFigures.Figure10())
savefig(foldername * "figure10.png")
@info "Finished figure 10"

## Figure 11 - Parameter estimates (buoy data)
pyplot(size = (1000,750))
plot(PaperFigures.Figure11(),plotascontinuous=true) # set to false for scatter plot with CI bars
savefig(foldername * "figure11.png")
@info "Finished figure 11"
