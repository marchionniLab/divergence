% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{computeMultivariateBinaryMatrix}
\alias{computeMultivariateBinaryMatrix}
\title{Compute the binary matrix with digitized divergence coding}
\usage{
computeMultivariateBinaryMatrix(seMat, Baseline)
}
\arguments{
\item{seMat}{SummarizedExperiment with assay to be digitized, in [0, 1], with each column corresponding to a sample and each 
row corresponding to a feature; usually in quantile form.}

\item{Baseline}{A Baseline object; this corresponds to the output of findMultivariateGammaWithSupport() or 
computeMultivariateSupport()}
}
\value{
A matrix with the same columns as Mat, with rows being the multivariate features, containing the binary form data.
}
\description{
Function for obtaining the binary form for a matrix for multivariate divergence of data given a baseline range
}
\examples{
baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
seMat.base = SummarizedExperiment(assays=list(data=baseMat))
assays(seMat.base)$quantile = computeQuantileMatrix(seMat.base)
baseline = computeMultivariateSupport(seMat=seMat.base, FeatureSets=msigdb_Hallmarks)
dataMat = breastTCGA_Mat[, breastTCGA_Group != "NORMAL"]
seMat = SummarizedExperiment(assays=list(data=dataMat))
assays(seMat)$quantile = computeQuantileMatrix(seMat)
Mat.div = computeMultivariateBinaryMatrix(seMat=seMat, Baseline=baseline)

}
\keyword{binary}
\keyword{digitization}
