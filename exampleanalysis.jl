using JLD2, ProgressMeter, DSP, Statistics
@load "dat/anonomous_data.jld"
Δ = 1/1.28
fitstart = 85
fitend = length(waveseries)-1

##
include("src/exampleanalysis/PaperAnalysis.jl")

## compute spectrogram, errorgram and mean direction for threshold choice
errorgram, Ω, spectime = PaperAnalysis.compute_errorgram(waveseries, wavetime)
wavedir = mod2pi.(PaperAnalysis.mt_wavedirection(permutedims(reduce(vcat, waveseries)), 2304, 0).power) # uses same Ω as above

##
"Compute the lower cutoff threshold from the errorgram."
function errorgramthreshold(R, Ω, tol, window)
    first_good = findprev(j->abs(mean(R[i] for i in j-window:j+window))>tol, window+1:length(R)-window, length(R)÷2)
    return isnothing(first_good) ? error("Failed to find cutoff") : Ω[first_good + window]
end

"Compute the threshold to remove the swell. Hard coded for the specific swell in this example."
function wavedir_threshold(dir,Ω)
    thresh = findprev(dir.>5,length(dir)÷2)
    return isnothing(thresh) ? NaN : Ω[thresh]
end

## choose threholds
swellend = fitend-50 # final record which actually contains swell
errorthreshold = [errorgramthreshold(a,Ω, 2, 20) for a in eachcol(errorgram[:,fitstart:fitend])]
swellthreshold = [wavedir_threshold(w,Ω) for w in eachcol(wavedir[:,fitstart:swellend])]
threshold = max.(errorthreshold,[swellthreshold;zeros(length(errorthreshold)-length(swellthreshold))])

## initialise fits
ini = [PaperAnalysis.InitialiseData(PaperAnalysis.WindSea,waveseries[i],Δ,threshold[k]) for (k,i) ∈ enumerate(fitstart:fitend)]

## fit models
fits = @showprogress [PaperAnalysis.FittedModel(i) for i ∈ ini]
@save "dat/wavefits.jld" fits fitstart fitend