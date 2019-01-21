# divergence

### Wikum Dinalankara, Qian Ke, and Luigi Marchionni | Johns Hopkins University

This package provides functionality for performing divergence analysis. For more information on the divergence framework, see [Dinalankara et al, "Digitizing omics profiles by divergence from a baseline", PNAS 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5939095/).

***

Recently we introduced divergence analysis, a framework for analysing high dimensional omics data. The principle behind this approach is to encode deviation from a baseline profile as a simple digitized coding. Given a matrix of high dimensional omics data (RNA-Seqexpression, microarray expression, CpG level methylation, etc) where one group of samples are designated as the baseline cohort, the other sample profiles can be converted to a binary or ternary string indicating whether each omics feature of the sample is outside the baseline range or not. This can be performed at the individual feature level (in which
case the result will be a string of 0, 1, and -1 indicating no divergence, deviation by being above the baseline range, and deviation by being below the baseline range, respectively), or it can be perfomed at the multivariate feature level. Multivariate features are sets of single features (e.g. a functional gene set) and if divergence analysis is performed at this level, the result will be binary (0 indicating no deviation, 1 indicating deviation). The method applies no special normalization procedures other than a sample specific quantile transformation which is applied to each sample before baseline computation and divergence coding. 

The package contains a small data set (a subset of the TCGA breast cancer dataset), which is used in the vignette to demonstrate the workflow.
