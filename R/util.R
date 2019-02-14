
# Utility functions

feature_list_to_matrix = function(feature_list){

  # rownames will be features, colnames will be set names
  mat = t(sapply(unique(unlist(feature_list)), function(x) sapply(feature_list, function(y) x %in% y ) )) * 1

  if(length(feature_list) == 1){
    if(nrow(mat) == 1)
      mat = t(mat)
      rownames(mat) = feature_list[[1]]
      colnames(mat) = names(feature_list)
  }

  mat

}

check_gamma_beta = function(gamma, beta){

	if(any(gamma <= 0) | any(gamma >= 1) )
		stop("gamma must be in (0, 1)")

	if(beta <= 0 | beta > 1)
		stop("beta must be in (0, 1]")

}

check_gammas_beta = function(gammas, beta){

	if(! is.numeric(gammas)){
		stop("Invalid gamma values provided")
	}else{
		check_gamma_beta(gammas, beta)
	}

}

check_alpha = function(alpha){

	if(! (alpha > 0 & alpha < 1))
		stop("alpha must be in (0, 1)")

}

check_feature_set = function(FeatureSets){

  if(is.list(FeatureSets)){

    if(is.null(names(FeatureSets))){
        stop("Feature set names not provided")
    }

    featureMat = feature_list_to_matrix(FeatureSets)

  }else if(is.matrix(FeatureSets)){

    # check if matrix is binary
    if(! all(FeatureSets == 0 | FeatureSets == 1)){
      stop("Invalid feature sets")
    }

    featureMat = FeatureSets
  }else{
    stop("Invalid feature sets")
  }
  featureMat

}

check_parallel = function(parallel){

  if(! is.logical(parallel)){

    warning("Ivalid parallel parameter; setting parallel = FALSE")
    parallel = FALSE

  }else if(parallel && identical(.Platform$OS.type, "windows")){

  	# user has set parallel = TRUE but is on windows
    warning("Parallel support not available on Windows; setting parallel = FALSE")
    parallel = FALSE


  }

  parallel

}

matrix_to_SE = function(mat){

  message("Convering input matrix to SummarizedExperiment")

  SummarizedExperiment(assays=list(data=mat))

}




