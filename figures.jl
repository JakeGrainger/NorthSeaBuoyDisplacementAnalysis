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

## Figure 1 - example S(ω,ϕ)
pyplot(size = (1000,300))
plot(PaperFigures.Figure1(), colorbar = false)
savefig(foldername * "figure1_nocolorbar.png")
plot(PaperFigures.ColorFig1(PaperFigures.Figure1()))
savefig(foldername * "figure1_colorbar.png")
println("Finished figure 1")

## Figure 2 - S(ω,ϕ) vs f(ω)
pyplot(size = (1000,500))
plot(plot(PaperFigures.Figure2a()), plot(PaperFigures.Figure2b()), layout = (1,2))
savefig(foldername * "figure2.png")
println("Finished figure 2")

## Figure 3 - buoy data summary
pyplot(size = (1000,1000))
plot(PaperFigures.Figure3())
savefig(foldername * "figure3.png")
println("Finished figure 3")

## Figure 4 - S(ω,ϕ), f(ω), D(ω,ϕ)
pyplot(size = (1000,250))
plot(PaperFigures.Figure4())
savefig(foldername * "figure4.png")
println("Finished figure 4")

## Figure 5 - Vary parameters
pyplot(size = (1000,600))
plot(PaperFigures.Figure5(),legend=false)
savefig(foldername * "figure5.png")
println("Finished figure 5")

## Figure 6 - Scenario S(ω,ϕ)
pyplot(size = (1000,300))
plot(PaperFigures.Figure6())
savefig(foldername * "figure6_nocolorbar.png")
plot(PaperFigures.ColorFig6(PaperFigures.Figure6()))
savefig(foldername * "figure6_colorbar.png")
println("Finished figure 6")

## Figure 7 - Marginal parameter estimates (simulation)
pyplot(size = (1000,800))
plot(PaperFigures.Figure7(truncate=true))
savefig(foldername * "figure7.png")
plot(PaperFigures.Figure7(truncate=false))
savefig(foldername * "figure7untruncated.png")
println("Finished figure 7")

## Figure 8 - Directional parameter estimates (simulation)
pyplot(size = (1000,1000))
plot(PaperFigures.Figure8(truncate=true))
savefig(foldername * "figure8.png")
plot(PaperFigures.Figure8(truncate=false))
savefig(foldername * "figure8untruncated.png")
println("Finished figure 8")

## Figure 9 - Heatmap of ̂R(ω,t)
pyplot(size = (1000,250))
plot(PaperFigures.Figure9())
savefig(foldername * "figure9.png")
println("Finished figure 9")

## Figure 10 - Spectrogram with cutoffs
pyplot(size = (1000,250))
plot(PaperFigures.Figure10())
savefig(foldername * "figure10.png")
println("Finished figure 10")

## Figure 11 - Parameter estimates (buoy data)
pyplot(size = (1000,750))
plot(PaperFigures.Figure11())
savefig(foldername * "figure11.png")
println("Finished figure 11")
