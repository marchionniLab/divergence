% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{computeUnivariateDigitization}
\alias{computeUnivariateDigitization}
\title{Perform ternary digitization}
\usage{
computeUnivariateDigitization(seMat, seMat.base, computeQuantiles = TRUE,
  gamma = c(1:9/100, 1:9/10), beta = 0.95, alpha = 0.01,
  parallel = TRUE, verbose = TRUE, findGamma = TRUE, Groups = NULL,
  classes = NULL)
}
\arguments{
\item{seMat}{SummarizedExperiment with assay to be digitized, in [0, 1], with each column corresponding to a sample and each 
row corresponding to a feature; usually in quantile form.}

\item{seMat.base}{SummarizedExperiment with baseline assay in [0, 1], with each column corresponding to a sample and each 
row corresponding to a feature}

\item{computeQuantiles}{Logical; apply quantile transformation to both data and baseline matrices (TRUE or FALSE; defaults to TRUE).}

\item{gamma}{Range of gamma values to search through. 
By default gamma = {0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9}.}

\item{beta}{Parameter for eliminating outliers (0 < beta <= 1). By default beta=0.95.}

\item{alpha}{Expected proportion of divergent features per sample to be estimated. The optimal gamma providing
this level of divergence in the baseline data will be searched for.}

\item{parallel}{Logical indicating whether to compute features parallelly with mclapply on
Unix based systems (defaults to TRUE, switched to FALSE if parallel package is not available).}

\item{verbose}{Logical indicating whether to print status related messages during computation (defaults
to TRUE).}

\item{findGamma}{Logical indicating whether to search for optimal gamma values through the given gamma values (defaults to 
TRUE). If FALSE, the first value given in gamma will be used.}

\item{Groups}{Factor indicating class association of samples (optional).}

\item{classes}{Vector of class labels (optional).}
}
\value{
A list with elements:
			Mat.div: divergence coding of data matrix in ternary (-1, 0, 1) form, of same dimensions at seMat
		baseMat.div: divergence coding of base matrix in ternary (-1, 0, 1) form, of same dimensions at seMat.base
			div: data frame with the number of divergent features in each sample, including upper and lower divergence
		features.div: data frame with the divergent probability of each feature; divergence probability for each 
				phenotype in included as well if 'Groups' and 'classes' inputs were provided.
		Baseline: a list containing a "Ranges" data frame with the baseline interval  for each feature, and a "Support" 
				binary matrix of the same dimensions as Mat indicating whether each sample was a support or a feature or not
				(1=support, 0=not in the support), gamma: selected gamma value, alpha: the expected number of divergent features per 
					sample computed over the baseline data matrix, optimal: logical indicaing whether the selected gamma value provided 
					the necessary alpha requirement, alpha_space: a data frame with alpha values for each gamma searched
}
\description{
Function for obtaining the digitized form, along with other relevant statistics and measures 
given a data matrix and a baseline matrix
}
\examples{
baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
dataMat = breastTCGA_Mat[, breastTCGA_Group != "NORMAL"]
seMat.base = SummarizedExperiment(assays=list(data=baseMat))
seMat = SummarizedExperiment(assays=list(data=dataMat))
div = computeUnivariateDigitization(
  seMat = seMat,
  seMat.base = seMat.base,
 parallel = TRUE
)
assays(seMat)$div = div$Mat.div

}
\keyword{digitize}
\keyword{ternary}
