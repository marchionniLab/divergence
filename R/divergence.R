
###
### divergence computation functions
###

### ====================================================
### quantiles
### ====================================================

quantileTransform = function(x){
  
  #if(rank_samples)
  #  x = rank(x, ties.method="average")
  
  # start with rank 0, i.e. all values with min(x) will have rank 0
  # the percentile transformation depends on the non-zero ranks
  x = x - min(x)
  
  y = x
  # separate out the positive entries and apply transformation
  z = rank(x[x > 0], ties.method="min")
  y[x > 0] = (z - min(z))/sum(x > 0)
  y

}


### apply quantiles transformation to a matrix
###     rank_samples: TRUE if ranking should be applied before converting to percentiles
getQuantileMat = function(Mat){
  
  newMat = apply(Mat, 2, function(x) quantileTransform(x))
  
  rownames(newMat) = rownames(Mat)
  colnames(newMat) = colnames(Mat)
  
  newMat
}

### ====================================================
### get expected proportion of divergent feature per sample
### ====================================================

getAlpha = function(Mat, Baseline){
  
  mean(colSums(abs(computeTernary(Mat=Mat, Baseline=Baseline)))/nrow(Mat))
  
}

### ====================================================
### computing ranges
### ====================================================

computeSingleRange = function(x, gamma, beta, j=NULL){
  

  if(is.null(j))
    j = max(floor(gamma * length(x)), 1)
  
  # get distance to j'th nearest neighbor from each point in x
  xjn = sapply(1:length(x), function(i){
    sort(abs(x[-i]-x[i]), partial=j)[j]
  })

  sel = which(xjn <= quantile(xjn, beta))
  
  x = x[sel]
  xjn = xjn[sel]
  
  r = c(
    max(x[which.min(x)] - xjn[which.min(x)], 0),
    min(x[which.max(x)] + xjn[which.max(x)], 1)
  )
  
  list(range=r, support=sel)
  
}

getRangeList = function(Mat, gamma=0.1, beta=0.95, par=TRUE, nmax=200, mmax=1000, verbose=TRUE){

	L = NULL

  # nmax = cols
  # mmax = rows

  # for testing:
  # nmax = 50

	if(par){

    # check_parallel should have already checked if parallel package is available by this point
		require(parallel)

		if(nrow(Mat) > mmax && ncol(Mat) > nmax){

			# for large matrices, break-up into chunks and process
			# each setion separately with mclapply

			L = getRangesBySections(Mat=Mat, gamma=gamma, beta=beta, par=par, nmax=nmax, mmax=mmax, verbose=verbose)

		}else{

			# process in parallel by each row

      tryCatch({
        L = mclapply(1:nrow(Mat), function(i) computeSingleRange(Mat[i, ], gamma=gamma, beta=beta, j=NULL),
          mc.cores=detectCores()-1
        )
      }, error = function(e){warning(e)})

    	if(sum(sapply(L, is.null)) > 0 || length(L) < nrow(Mat)){

    		# did not retrurn all features; re-run without parallel
        if(verbose)
      		cat("Not all features returned; re-running without parallelization\n")

    		L = getRangeList(Mat=Mat, gamma=gamma, beta=beta, par=FALSE, nmax=nmax, mmax=mmax)

    	}

		}

	}else{

		# run without any parallelization

		L = lapply(1:nrow(Mat), function(i){
      	tempx = NULL
      	tryCatch({
        	tempx = computeSingleRange(Mat[i, ], gamma=gamma, beta=beta, j=NULL)
      	}, error = function(e){warning(e)})
      	tempx
    	})

	}

	names(L) = rownames(Mat)

	L

}

getRangesBySections = function(Mat, gamma=0.1, beta=0.95, par=TRUE, nmax=200, mmax=1000, verbose=TRUE){

  dirnum = paste(sapply(1:5, function(i) sample(0:9, 1)), collapse="")
  dirname = sprintf("ranges_%s", dirnum)

  if(dir.exists(dirname)){
    dirnum = paste(sapply(1:5, function(i) sample(0:9, 1)), collapse="")
    dirname = sprintf("ranges_%s", dirnum)
  }

  if(verbose)
    cat(sprintf("Creating directory %s for saving files\n", dirname))

  if(dir.exists(dirname)){
    cat(sprintf("Directory %s already exits. Contents will be overwritten.\n", dirname))
    tryCatch({
      system(sprintf("rm %s/*.rda", dirname))
    }, error = function(e){print(e)})
  }else{
    dir.create(dirname)
  }

  # apply to each feature in the baseline data matrix

  nc = detectCores()-1
  #nc = floor(detectCores()/2)
  
  n = nrow(Mat)
  k = ceiling(n/mmax)
  
  slots = lapply(1:k, function(i){ 
    x = 1:mmax + ((i-1) * mmax) 
    x[x <= n]
  })
  
  #print(sapply(slots, function(x) c(x[1:2], x[(length(x)-1):length(x)], length(x)) ))
  
  for(j in 1:length(slots)){
    
    if(verbose)
    	cat(sprintf("Processing features %d to %d\n", min(slots[[j]]), max(slots[[j]]) ))
    
    partialRanges = getRangeList(Mat[slots[[j]], ], gamma=gamma, beta=beta, par=par, nmax=nmax, mmax=mmax)
    
    if(verbose)
      cat("Saving..\n")
    save(partialRanges, file=sprintf("%s/%d.rda", dirname, j))
    
    rm(partialRanges)
    gc()
  
  }
  
  # assemble
  if(verbose) 
    cat("Assembling..\n")

  L = list()
  for(j in 1:length(slots)){
    
    ll = load(sprintf("%s/%d.rda", dirname, j))
    
    L[[j]] = partialRanges
    
    rm(list=ll)
    gc()
    
  }

  rL = Reduce(c, L)
  if(! all(names(rL) == rownames(Mat))){
      warning("Warning: rownames after assembly not in correct order; will try to re-order them..\n")
      tryCatch({
        rL = rL[rownames(Mat)]
      }, error = function(e){print(e)})
  }

  if(verbose)
    cat("Deleting temporary files..")

  unlink(dirname, recursive=TRUE)

  if(verbose)
    cat("done.\n")

  rL
}

computeRanges = function(Mat, gamma=0.1, beta=0.95, parallel=TRUE, verbose=TRUE){
  
  if(! is.matrix(Mat)){
    stop("Input data is not in matrix form")
  }

  check_gamma_beta(gamma, beta)

  if(verbose){
    cat(sprintf("Computing ranges from %d reference samples for %d features\n", 
                ncol(Mat), nrow(Mat)))
    cat(sprintf("[beta=%g, gamma=%g]\n", beta, gamma))
  }

  # apply to each feature in the baseline data matrix
  L = getRangeList(Mat=Mat, gamma=gamma, beta=beta, par=parallel, verbose=verbose)
  
  if(! all(sapply(L, function(x) "range" %in% names(x)))){
    stop("Not all feature level supports were computed")
  }

  R = data.frame(t(sapply(L, function(x) x$range)))
  colnames(R) = c("baseline.low", "baseline.high")
  
  S = t(sapply(L, function(x) 1*(1:ncol(Mat) %in% x$support) ))
  colnames(S) = colnames(Mat)
  
  Baseline = list(Ranges=R, Support=S)

  alpha=getAlpha(Mat=Mat, Baseline=Baseline)
  if(verbose)
    cat(sprintf("[Expected proportion of divergent features per sample=%g]\n", alpha))
  
  
  list(Baseline=Baseline,
       gamma=gamma,
       alpha=alpha)
  
}

### ====================================================
### search for gamma
### ====================================================

findGamma = function(Mat,
	gamma=c(1:9/100, 1:9/10),
  beta=0.95, 
  alpha=0.01,
  parallel=TRUE,
  verbose=TRUE){
 
  if(! is.matrix(Mat)){
    stop("Input data is not in matrix form")
  }

  check_gammas_beta(gamma, beta)

 if(verbose)
    cat(sprintf("Searching optimal gamma for alpha=%g\n", alpha))
   
  optimal_gamma = -1
  
  gamma = sort(gamma)
  names(gamma) = paste("g", 1:length(gamma), sep="")
  
  g = c()
  e = rep(NA, length(gamma))
  names(e) = names(gamma)
    
  RangesList = list()
  for(i in 1:length(gamma)){
    L = computeRanges(Mat=Mat, gamma=gamma[i], beta=beta, parallel=parallel, verbose=verbose)
    
    RangesList[[i]] = L
    e[i] = L$alpha
    
    g = c(g, gamma[i])
    if(e[i] <= alpha)
      break
  }
  names(RangesList) = names(gamma)[1:length(RangesList)]
  names(g) = names(gamma)[1:length(g)]

  optimal = FALSE
  if(length(which(e <= alpha)) < 1){
    sel = which.max(g)
  }else{
    sel = which(g == min(g[which(e <= alpha)]))
  }
  optimal_gamma = g[sel]
  R_star = RangesList[[ sel ]]
  e_star = e[sel]
  
  if(e_star <= alpha)
    optimal = TRUE
  
  temp_e = e
  names(temp_e) = g
  #print(temp_e)

  if(verbose)
    cat(sprintf("Search results for alpha=%g: gamma=%g, expectation=%g, optimal=%s\n", alpha, optimal_gamma, e_star, optimal))
  
  list(Baseline=R_star$Baseline,
       gamma=optimal_gamma,
       alpha=e_star,
       optimal=optimal,
       alpha_space=data.frame(gamma=gamma, alpha=e)
       )
  
}

### ====================================================
### compute ternary form
### ====================================================

computeTernary = function(Mat, Baseline){
  
  R = Baseline$Ranges

  lower = "baseline.low"
  upper = "baseline.high"

  if(nrow(Mat) != nrow(R))
    cat(sprintf("WARNING [%d, %d]\n", nrow(Mat), nrow(R)))
  
  DMat = ((Mat < R[, lower]) * (-1)) + ((Mat > R[, upper]) * 1)
  rownames(DMat) = rownames(Mat)
  colnames(DMat) = colnames(Mat)
  DMat
  
}

### ====================================================
### compute divergences
### ====================================================

computeTernaryDigitization = function(Mat, baseMat, 
                              computeQuantiles=TRUE,
                              gamma=c(1:9/100, 1:9/10),
                              beta=0.95, 
                              alpha=0.01,
                              parallel=TRUE,
                              verbose=TRUE,
                              findGamma=TRUE, 
                              Groups=NULL,                             
                              classes=NULL){

  if( ! all(rownames(Mat) == rownames(baseMat)) ){
    stop("Feature names different in data and baseline matrices")
  }

  if(! is.null(Groups)){

    stopifnot(length(Groups) == ncol(Mat))

    if(! is.factor(Groups))
      Groups = factor(Groups)

    if(is.null(classes))
       classes = levels(Groups)
 
  }

  check_alpha(alpha)

  check_gammas_beta(gamma, beta)

  if(computeQuantiles){
    if(verbose)
      cat(sprintf("Computing quantiles..\n"))
    baseMat = getQuantileMat(baseMat)
    Mat = getQuantileMat(Mat)
  }

  if(findGamma){
    L = findGamma(Mat=baseMat, gamma=gamma, beta=beta, alpha=alpha, parallel=parallel, verbose=verbose)
  }
  else{
    if(verbose)
      cat(sprintf("Using gamma=%g\n", gamma[1]))
    L = findGamma(Mat=baseMat, gamma=gamma[1], beta=beta, alpha=alpha, parallel=parallel, verbose=FALSE)
  }


  DMat_ternary = computeTernary(Mat=Mat, Baseline=L$Baseline)
  
  baseMat_ternary = computeTernary(Mat=baseMat, Baseline=L$Baseline)

  DMat = abs(DMat_ternary)
  
  D = rowMeans(DMat)
  N = colSums(DMat)
  
  Npos = colSums(DMat_ternary > 0)
  Nneg = colSums(DMat_ternary < 0)
  
  df = data.frame(feature=rownames(DMat), prob.div=D)

  if(! is.null(Groups)){

    classDiv = sapply(classes, function(x) rowMeans(DMat[, which(Groups == x)]))

    df = data.frame(df, classDiv)
    colnames(df) = c("feature", "prob.div", paste("prob.div.", classes, sep=""))

  }

  list(Mat.div=DMat_ternary,
      baseMat.div = baseMat_ternary,
      div = data.frame(sample=colnames(Mat), count.div=N, count.div.upper=Npos, count.div.lower=Nneg),
      features.div = df,
      Baseline = L$Baseline,
      gamma = L$gamma,
      alpha = L$alpha,
      optimal = L$optimal,
      alpha_space = L$alpha_space)
    
}

