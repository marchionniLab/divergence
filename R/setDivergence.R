

# compute support for a single feature set
computeSingleFeaturesetSupport = function(Mat, features, gamma, beta, distance, j=NULL){

  check_gamma_beta(gamma, beta)

  features = intersect(rownames(Mat), features)

  if(length(features) == 0)
      stop("No features available in given data")

  Mat = Mat[features, , drop=FALSE]

  if(is.null(j))
    j = max(floor(gamma * ncol(Mat)), 1)

  # get distance to j'th nearest neighbor from each point in x

  Mat_Dist = as.matrix(dist(t(Mat), method = distance, diag = FALSE))

  xjn = apply(Mat_Dist, 1, function(y){
    sort(y, decreasing=FALSE)[j+1]
  })

  ### sel  indicate whether a sample is used to build the support
  sel = which(xjn <= quantile(xjn, beta))

  Radius = xjn[sel]
  Centermatrix = Mat[, sel]
  Mat_Dist_m2c = Mat_Dist[, sel]

  ### compare with radius

  v = apply(Mat_Dist_m2c,1, function(y){
    as.numeric(all(y > Radius))
  })

  prob = sum(v==1)/length(v)

  list(Radius=Radius, Centers=Centermatrix, features=features, alpha=prob, distance=distance)

}

computeSingleFeaturesetBinaryVector = function(Mat, Baseline){

  if(! all(Baseline$features %in% rownames(Mat)))
    stop("Some features in the support not available in the given data")

  Mat = Mat[Baseline$features, , drop=FALSE]

  Centers = Baseline$Centers
  Radius = Baseline$Radius
  distance = Baseline$distance

  mn_m <- dim(Mat)
  mn_c <- dim(Centers)

  Mat_Matrix <- cbind(Centers, Mat)

  Mat_Dist <- as.matrix(dist(t(Mat_Matrix), method=distance, diag=FALSE))

  ### get distance from Mat to Centers
  Mat_Dist_m2c <-  Mat_Dist[c((mn_c[2]+1):(mn_c[2]+mn_m[2])),c(1:mn_c[2])]

  ### compare with radius
  v <- apply(Mat_Dist_m2c,1, function(y){
    as.numeric(all(y > Radius))
  })

  v

}

computeFeatureSetSupport = function(seMat, FeatureSets, gamma=0.1, beta=0.95, distance="euclidean", verbose=TRUE){

  check_gamma_beta(gamma, beta)

  featureMat = check_feature_set(FeatureSets)

  Mat = get_mat_from_SE(seMat)

  Baseline_list = lapply(1:ncol(featureMat), function(j){
      features = rownames(featureMat)[which(featureMat[, j] == 1)]
      computeSingleFeaturesetSupport(Mat=Mat, features=features, gamma=gamma, beta=beta, distance=distance, j=NULL)
  })
  names(Baseline_list) = colnames(featureMat)

  Support = 1 * t(sapply(Baseline_list, function(B) computeSingleFeaturesetBinaryVector(Mat=Mat, Baseline=B)) == 0)

  alpha = mean(sapply(Baseline_list, function(x) x$alpha))

  list(Support=Support, Baseline_list=Baseline_list, featureMat=featureMat, alpha=alpha, gamma=gamma, distance=distance)

}

findFeatureSetGammaAndSupport = function(seMat, FeatureSets, gamma=1:9/10, beta=0.95, alpha=0.01, distance="euclidean", verbose=TRUE){

  check_gammas_beta(gamma, beta)
  check_alpha(alpha)

  if(verbose)
    message(sprintf("Searching optimal support for alpha threshold=%g\n", alpha))

  fgs = list()
  optimal = FALSE

  Mat = get_mat_from_SE(seMat)

  for(i in seq_along(gamma)){

      S = computeFeatureSetSupport(
        seMat=seMat, 
        FeatureSets=FeatureSets, 
        gamma=gamma[i], 
        beta=beta, 
        distance=distance, 
        verbose=verbose)

      if(verbose){
        message(sprintf("\t[gamma=%g, beta=%g, alpha=%g]\n", gamma[i], beta, S$alpha))
      }

      fgs[[i]] = S

      if(S$alpha <= alpha){
        optimal = TRUE
        break
      }

  }

  alpha_space = data.frame(gamma=gamma, alpha=NA)
  alpha_space$alpha[seq_along(fgs)] = sapply(fgs, function(x) x$alpha)

  S$gamma = gamma[i]
  S$optimal = optimal
  S$alpha_space = alpha_space

  S

}

computeFeatureSetBinaryMatrix = function(seMat, Baseline){

  Mat = get_mat_from_SE(seMat)

  binMat = t(sapply(Baseline$Baseline_list, function(B){
      computeSingleFeaturesetBinaryVector(Mat=Mat, Baseline=B)
  }))

  binMat

}

computeFeatureSetDigitization = function(seMat, seMat.base, FeatureSets,
                              computeQuantiles=TRUE,
                              gamma=c(1:9/100, 1:9/10),
                              beta=0.95, 
                              alpha=0.01,
                              distance="euclidean",
                              verbose=TRUE,
                              findGamma=TRUE, 
                              Groups=NULL,                             
                              classes=NULL){


  stopifnot(rownames(seMat) == rownames(seMat.base))

  if(! is.null(Groups)){

    stopifnot(length(Groups) == ncol(seMat))

    if(! is.factor(Groups))
      Groups = factor(Groups)

    if(is.null(classes))
       classes = levels(Groups)
 
  }

  check_alpha(alpha)
  check_gammas_beta(gamma, beta)

  featureMat = check_feature_set(FeatureSets)

  if(computeQuantiles){

    if(verbose)
      message(sprintf("Computing quantiles..\n"))

    assays(seMat.base)$quantile = getQuantileMat(seMat.base)
    assays(seMat)$quantile = getQuantileMat(seMat)

  }

  if(findGamma){

    B = findFeatureSetGammaAndSupport(
      seMat=seMat.base, 
      FeatureSets=FeatureSets, 
      gamma=gamma, beta=beta, alpha=alpha, distance=distance, verbose=verbose
    )

  }else{
    
    if(verbose)
      message(sprintf("Using gamma=%g\n", gamma[1]))
    
    B = computeFeatureSetSupport(
      seMat=seMat.base, 
      FeatureSets=FeatureSets, 
      gamma=gamma[1], beta=beta, distance=distance, verbose=verbose
    )

  }

  DMat = computeFeatureSetBinaryMatrix(seMat=seMat, Baseline=B)
  baseDMat = computeFeatureSetBinaryMatrix(seMat=seMat.base, Baseline=B) 

  D = rowMeans(DMat)
  N = colSums(DMat)
    
  df = data.frame(feature=rownames(DMat), prob.div=D)

  if(! is.null(Groups)){

    classDiv = sapply(classes, function(x) rowMeans(DMat[, which(Groups == x)]))

    df = data.frame(df, classDiv)
    colnames(df) = c("feature", "prob.div", paste("prob.div.", classes, sep=""))

  }

  list(
    Mat.div=DMat,
    baseMat.div = baseDMat,
    div = data.frame(sample=colnames(DMat), count.div=N),
    features.div = df,
    Baseline = B,
    gamma = B$gamma,
    alpha = B$alpha
  )
    
}






