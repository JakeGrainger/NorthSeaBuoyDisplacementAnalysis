abstract type FDSpectralModel end

struct AdditiveFDModel{M1,M2} <: FDSpectralModel
    model1::M1
    model2::M2
end

freqdirspectra(FD::FDSpectralModel, ω, ϕ) = spreading(FD, ω, ϕ) * sdf(FD, ω)
freqdirspectra(aFD::AdditiveFDModel, ω, ϕ) = freqdirspectra(aFD.model1, ω, ϕ) + freqdirspectra(aFD.model2, ω, ϕ)

Base.:+(FD1::FDSpectralModel,FD2::FDSpectralModel) = AdditiveFDModel(FD1,FD2)

