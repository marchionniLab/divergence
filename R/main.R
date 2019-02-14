#' Compute quantile transformations
#'
#' Function for computing the quantile transformation for one or more samples
#' supplied as columns of a matrix.
#'
#' @param Mat Matrix of data, with each column corresponding to a sample and each row corresponding to a feature.
#'
#' @return A matrix of the same dimensions as Mat with the quantile data.
#'
#' @export
#'
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' baseMat.q = computeQuantileMatrix(baseMat)
#'

computeQuantileMatrix <- function(Mat){
	getQuantileMat(Mat=Mat)
}

# ==================================================================================================
#' Estimate the baseline support
#'
#' Function for computing the basline support for univariate features given gamma
#' and beta parameters.
#'
#' @param Mat Matrix of data in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param gamma Parameter for selecting radius around each support point (0 < gamma < 1).
#' By default gamma = 0.1.
#' @param beta Parameter for eliminating outliers (0 < beta <= 1). By default beta=0.95.
#' @param parallel Logical indicating whether to compute features parallelly with mclapply on
#' Unix based systems (defaults to TRUE, switched to FALSE if parallel package is not available).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements "Ranges": data frame with the baseline interval 
#' for each feature, "Support": binary matrix of the same dimensions as Mat indicating whether each sample was a support 
#' for a feature or not (1=support, 0=not in the support), "gamma": gamma value, and "alpha": the expected number of divergent 
#' features per sample estimated over the samples.
#' 
#' @keywords baseline, support
#' @export
#' 
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' baseMat.q = computeQuantileMatrix(baseMat)
#' baseline = computeUnivariateSupport(Mat=baseMat.q)
#'


computeUnivariateSupport <- function(Mat, gamma=0.1, beta=0.95, parallel=TRUE, verbose=TRUE){

	parallel = check_parallel(parallel)

	computeRanges(Mat=Mat, gamma=gamma, beta=beta, parallel=parallel, verbose=verbose)
}

# ==================================================================================================
#' Search for optimal gamma and associated support
#'
#' Function for searching through a range of gamma values for finding the smallest gamma that 
#' provides expected proportion of divergent features per sample less than or equal to alpha.
#'
#' @param Mat Matrix of data in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param gamma Range of gamma values to search through. 
#' By default gamma = \{0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9\}.
#' @param beta Parameter for eliminating outliers (0 < beta <= 1). By default beta=0.95.
#' @param alpha Expected proportion of divergent features per sample to be estimated
#' over the samples in Mat. By default alpha = 0.01; i.e. search for the smallest gamma that provides
#' 1\% or less number of divergent features per sample. 
#' @param parallel Logical indicating whether to compute features parallelly with mclapply on
#' Unix based systems (defaults to TRUE, switched to FALSE if parallel package is not available).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements "Ranges": data frame with the baseline interval 
#' for each feature, "Support": binary matrix of the same dimensions as Mat indicating whether each sample was a support 
#' for a feature or not (1=support, 0=not in the support), "gamma": gamma value, and "alpha": the expected number of divergent 
#' features per sample estimated over the samples, "optimal": logical indicaing whether the 
#' selected gamma value provided the necessary alpha requirement, and "alpha_space": a data frame with alpha values 
#' for each gamma searched.
#' 
#' @keywords gamma
#' @export
#' 
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' baseMat.q = computeQuantileMatrix(baseMat)
#' baseline = findUnivariateGammaWithSupport(Mat=baseMat.q)
#'

findUnivariateGammaWithSupport <- function(Mat, gamma=c(1:9/100, 1:9/10), beta=0.95, alpha=0.01, parallel=TRUE, verbose=TRUE){

	parallel = check_parallel(parallel)

  findGamma(Mat=Mat, gamma=gamma, beta=beta, alpha=alpha, parallel=parallel, verbose=verbose)

}

# ==================================================================================================
#' 
#' Compute the ternary matrix with digitized divergence coding
#'
#' Function for obtaining the ternary form for a matrix of data given a baseline range
#'
#' @param Mat Matrix of data to be digitized, in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param Baseline A list with a data frame element "Ranges" containing the baseline range of each features; 
#' this corresponds to the output of findUnivariateGammaWithSupport() or computeUnivariateSupport()
#'
#' @return A matrix with the same dimensions as Mat containing the ternary form data.
#'
#' @keywords ternary digitization
#' @export
#'
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' baseMat.q = computeQuantileMatrix(baseMat)
#' baseline = computeUnivariateSupport(Mat=baseMat.q)
#' dataMat = breastTCGA_Mat[, breastTCGA_Group != "NORMAL"]
#' dataMat.q = computeQuantileMatrix(dataMat)
#' divMat = computeUnivariateTernaryMatrix(Mat=dataMat.q, Baseline=baseline)
#'

computeUnivariateTernaryMatrix <- function(Mat, Baseline){
			computeTernary(Mat=Mat, Baseline=Baseline)
}

# ==================================================================================================
#' 
#' Perform ternary digitization
#'
#' Function for obtaining the digitized form, along with other relevant statistics and measures 
#' given a data matrix and a baseline matrix 
#'
#' @param Mat Matrix of data to be digitized, in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param baseMat Matrix of baseline data in [0, 1], (usually in quantiles), with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param computeQuantiles Apply quantile transformation to both data and baseline matrices (TRUE or FALSE; defaults to TRUE).
#' @param gamma Range of gamma values to search through. 
#' By default gamma = {0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9}.
#' @param beta Parameter for eliminating outliers (0 < beta <= 1). By default beta=0.95.
#' @param alpha Expected proportion of divergent features per sample to be estimated. The optimal gamma providing
#' this level of divergence in the baseline data will be searched for.
#' @param parallel Logical indicating whether to compute features parallelly with mclapply on
#' Unix based systems (defaults to TRUE, switched to FALSE if parallel package is not available).
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#' @param findGamma Logical indicating whether to search for optimal gamma values through the given gamma values (defaults to 
#' TRUE). If FALSE, the first value given in gamma will be used.
#' @param Groups Factor indicating class association of samples (optional).
#' @param classes Vector of class labels (optional).
#'
#' @return A list with elements:
#' 			Mat.div: divergence coding of data matrix in ternary (-1, 0, 1) form, of same dimensions at Mat
#'			baseMat.div: divergence coding of base matrix in ternary (-1, 0, 1) form, of same dimensions at Mat
#' 			div: data frame with the number of divergent features in each sample, including upper and lower divergence
#'			features.div: data frame with the divergent probability of each feature; divergence probability for each 
#'					phenotype in included as well if 'Groups' and 'classes' inputs were provided.
#'			Baseline: a list containing a "Ranges" data frame with the baseline interval  for each feature, and a "Support" 
#'					binary matrix of the same dimensions as Mat indicating whether each sample was a support or a feature or not
#'					(1=support, 0=not in the support), gamma: selected gamma value, alpha: the expected number of divergent features per 
#' 					sample computed over the baseline data matrix, optimal: logical indicaing whether the selected gamma value provided 
#' 					the necessary alpha requirement, alpha_space: a data frame with alpha values for each gamma searched
#' 
#' @keywords digitize ternary
#' @export
#'
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' dataMat = breastTCGA_Mat[, breastTCGA_Group != "NORMAL"]
#' div = computeUnivariateDigitization(
#'   Mat = dataMat,
#'   baseMat = baseMat,
#'	 parallel = TRUE
#')
#'

computeUnivariateDigitization <- function(Mat, baseMat, 
															computeQuantiles=TRUE,
                              gamma=c(1:9/100, 1:9/10),
                              beta=0.95, 
                              alpha=0.01,
                              parallel=TRUE,
                              verbose=TRUE,
                              findGamma=TRUE, 
                              Groups=NULL,                             
                              classes=NULL
){

	parallel = check_parallel(parallel)

	computeTernaryDigitization(Mat=Mat, baseMat=baseMat, computeQuantiles=computeQuantiles,
		gamma=gamma, beta=beta, alpha=alpha, parallel=parallel, verbose=verbose, 
		findGamma=findGamma, Groups=Groups, classes=classes)

}













# ==================================================================================================
#'
#' Compute chi-squared test
#'
#' Given a binary or ternary data matrix with class associations of samples, computes chi-squared tests
#' for each feature between given groups
#'
#' @param Mat Matrix of digitized binary or ternary data with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param Groups Factor indicating class association of samples
#' @param classes Vector of class labels; the test will be applied between the classes given.
#'
#' @return A data frame with columns 'statistic' and 'pval'.
#'
#' @keywords chi-squared
#' @export
#'
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' dataMat = breastTCGA_Mat[, breastTCGA_Group != "NORMAL"]
#' div = computeUnivariateDigitization(
#'   Mat = dataMat,
#'   baseMat = baseMat,
#' 	 parallel = TRUE
#')
#' sel = which(colnames(breastTCGA_Mat) %in% colnames(dataMat))
#' div.chi = computeChiSquaredTest(Mat=div$Mat.div, 
#'                                 Groups=breastTCGA_ER[sel],
#'                                 classes=c("Positive", "Negative"))
#'
#'

computeChiSquaredTest <- function(Mat, Groups, classes){

	chiSquaredTest(Mat=Mat, Groups=Groups, classes=classes)

}














# ==================================================================================================
#' Estimate the baseline support
#'
#' Function for computing the basline support for multivariate features given gamma
#' and beta parameters.
#'
#' @param Mat Matrix of data in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param FeatureSets The multivariate features in list or matrix form. In list form, each list element
#' should be a vector of individual features; in matrix form, it should be a binary matrix with rownames
#' being individual features and column names being the names of the feature sets.
#' @param gamma Parameter for selecting radius around each support point (0 < gamma < 1).
#' By default gamma = 0.1.
#' @param beta Parameter for eliminating outliers (0 < beta <= 1). By default beta=0.95.
#' @param distance Type of distance to be calculated between points. Any type of distance that can be passed on
#' to the dist function can be used (default 'euclidean'). 
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements:
#'		Support: a matrix indicating which samples were included in the support.
#'		Baseline_list: a list where each element is the baseline of a multivariate feature. 
#'		featureMat: the multivariate features in matrix form.
#'		alpha: the expected number of divergent multivariate features per sample.
#'		gamma: the gamma parameter used for baseline computation.
#'		distance: the type of distance used for baselien computation.
#' 
#' @keywords baseline, support
#' @export
#' 
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' baseMat.q = computeQuantileMatrix(baseMat)
#' baseline = computeMultivariateSupport(Mat=baseMat.q, FeatureSets=msigdb_Hallmarks)
#'

computeMultivariateSupport <- function(Mat, FeatureSets, gamma=0.1, beta=0.95, distance="euclidean", verbose=TRUE){


	computeFeatureSetSupport(Mat=Mat, FeatureSets=FeatureSets, gamma=gamma, beta=beta, distance=distance, verbose=verbose)

}

# ==================================================================================================
#' Find optimal gamma and corresponding support for list of feature sets
#'
#' Function for searching through a range of gamma values for finding the smallest gamma and support that 
#' provides expected proportion of divergent features per sample less than or equal to alpha.
#'
#' @param Mat Matrix of data in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param FeatureSets The multivariate features in list or matrix form. In list form, each list element
#' should be a vector of individual features; in matrix form, it should be a binary matrix with rownames
#' being individual features and column names being the names of the feature sets.
#' @param gamma Range of gamma values to search through. 
#' By default gamma = \{0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9\}.
#' @param beta Parameter for eliminating outliers (0 < beta <= 1). By default beta=0.95.
#' @param alpha Expected proportion of divergent features per sample to be estimated
#' over the samples in Mat. By default alpha = 0.01; i.e. search for the smallest gamma that provides
#' 1\% or less number of divergent features per sample. 
#' @param distance Type of distance to be calculated between points. Any type of distance that can be passed on
#' to the dist function can be used (default 'euclidean'). 
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#'
#' @return A list with elements:
#'		Support: a matrix indicating which samples were included in the support.
#'		Baseline: a list where each element is the baseline of a multivariate feature. 
#'		featureMat: the multivariate features in matrix form.
#'		alpha: the expected number of divergent multivariate features per sample.
#'		gamma: the gamma parameter selected.
#'		distance: the type of distance used for baselien computation.
#' 		optimal: TRUE or FALSE indicating whether the alpha criteria was met
#' 		alpha_space: the alpha values correspinding to the gamma values searched through
#'
#' @keywords gamma
#' @export
#' 
#' @examples
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' baseMat.q = computeQuantileMatrix(baseMat)
#' baseline = findMultivariateGammaWithSupport(Mat=baseMat.q, FeatureSets=msigdb_Hallmarks)
#'

findMultivariateGammaWithSupport <-function(Mat, FeatureSets, gamma=1:9/10, beta=0.95, alpha=0.01, distance="euclidean", verbose=TRUE){


	findFeatureSetGammaAndSupport(Mat=Mat, FeatureSets=FeatureSets, gamma=gamma, beta=beta, alpha=alpha, distance=distance, verbose=verbose)

}

# ==================================================================================================
#' 
#' Compute the binary matrix with digitized divergence coding
#'
#' Function for obtaining the binary form for a matrix for multivariate divergence of data given a baseline range
#'
#' @param Mat Matrix of data to be digitized, in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param Baseline A Baseline object; this corresponds to the output of findMultivariateGammaWithSupport() or 
#' computeMultivariateSupport()
#'
#' @return A matrix with the same columns as Mat, with rows being the multivariate features, containing the binary form data.
#'
#' @keywords binary digitization
#' @export
#'
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' baseMat.q = computeQuantileMatrix(baseMat)
#' baseline = computeMultivariateSupport(Mat=baseMat.q, FeatureSets=msigdb_Hallmarks)
#' dataMat = breastTCGA_Mat[, breastTCGA_Group != "NORMAL"]
#' dataMat.q = computeQuantileMatrix(dataMat)
#' divMat = computeMultivariateBinaryMatrix(Mat=dataMat.q, Baseline=baseline)
#'

computeMultivariateBinaryMatrix <- function(Mat, Baseline){

	computeFeatureSetBinaryMatrix(Mat=Mat, Baseline=Baseline)

}

# ==================================================================================================
#' 
#' Perform binary digitization
#'
#' Function for obtaining the digitized form, along with other relevant statistics and measures 
#' given a data matrix and a baseline matrix with multivariate features of interest
#'
#' @param Mat Matrix of data to be digitized, in [0, 1], with each column corresponding to a sample and each 
#' row corresponding to a feature; usually in quantile form.
#' @param baseMat Matrix of baseline data in [0, 1], (usually in quantiles), with each column corresponding to a sample and each 
#' row corresponding to a feature
#' @param FeatureSets The multivariate features in list or matrix form. In list form, each list element
#' should be a vector of individual features; in matrix form, it should be a binary matrix with rownames
#' being individual features and column names being the names of the feature sets.
#' @param computeQuantiles Apply quantile transformation to both data and baseline matrices (TRUE or FALSE; defaults to TRUE).
#' @param gamma Range of gamma values to search through. 
#' By default gamma = {0.01, 0.02, ... 0.09, 0.1, 0.2, ..., 0.9}.
#' @param beta Parameter for eliminating outliers (0 < beta <= 1). By default beta=0.95.
#' @param alpha Expected proportion of divergent features per sample to be estimated. The optimal gamma providing
#' this level of divergence in the baseline data will be searched for.
#' @param distance Type of distance to be calculated between points. Any type of distance that can be passed on
#' to the dist function can be used (default 'euclidean'). 
#' @param verbose Logical indicating whether to print status related messages during computation (defaults
#'  to TRUE).
#' @param findGamma Logical indicating whether to search for optimal gamma values through the given gamma values (defaults to 
#' TRUE). If FALSE, the first value given in gamma will be used.
#' @param Groups Factor indicating class association of samples
#' @param classes Vector of class labels
#'
#' @return A list with elements:
#' 			Mat.div: divergence coding of data matrix in ternary (-1, 0, 1) form, of same dimensions at Mat
#'			baseMat.div: divergence coding of base matrix in binary form, of same column names at Mat, rows being multivariate
#'					features.
#' 			div: data frame with the number of divergent features in each sample
#'			features.div: data frame with the divergent probability of each feature; divergence probability for each 
#'					phenotype in included as well if 'Groups' and 'classes' inputs were provided.
#'			Baseline: a list containing a "Ranges" data frame with the baseline interval  for each feature, and a "Support" 
#'					binary matrix of the same dimensions as Mat indicating whether each sample was a support or a feature or not
#'					(1=support, 0=not in the support), 
#'			gamma: selected gamma value
#'			alpha: the expected number of divergent features per sample computed over the baseline data matrix
#' @keywords digitize ternary
#' @export
#'
#' @examples
#' baseMat = breastTCGA_Mat[, breastTCGA_Group == "NORMAL"]
#' dataMat = breastTCGA_Mat[, breastTCGA_Group != "NORMAL"]
#' div = computeMultivariateDigitization(
#'   Mat = dataMat,
#'   baseMat = baseMat,
#'   FeatureSets = msigdb_Hallmarks
#' )
#'

computeMultivariateDigitization <- function(Mat, baseMat, FeatureSets,
                              computeQuantiles=TRUE,
                              gamma=c(1:9/100, 1:9/10),
                              beta=0.95, 
                              alpha=0.01,
                              distance="euclidean",
                              verbose=TRUE,
                              findGamma=TRUE, 
                              Groups=NULL,                             
                              classes=NULL){


	computeFeatureSetDigitization(
		Mat=Mat, 
		baseMat=baseMat, 
		FeatureSets=FeatureSets,
		computeQuantiles=computeQuantiles,
		gamma=gamma, beta=beta, alpha=alpha,
		distance=distance, verbose=verbose,
		findGamma=findGamma, Groups=Groups, classes=classes
	)

	
}





