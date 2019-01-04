
### =====================================================
### compute chi-squared test
### =====================================================

chiSquaredTest = function(Mat, Groups, classes){
  
  if(! all(classes %in% Groups))
    stop("Not all given classes are available in the groups variable")
  
  ch = apply(Mat, 1, function(x){
    
    if(length(unique(x)) > 1){
      u = sapply(classes, function(y) sapply(unique(x), function(z) sum(x[Groups == y] == z)))
      chisq.test(u)
    }else{
      list(statistic=NA, p.value=NA)
    }

  })
  
  ch2 = data.frame(t(sapply(ch, function(x) c(x$statistic, x$p.value))))
  colnames(ch2) = c("statistic", "pval")
  ch2[order(ch2$pval, na.last=TRUE), ]
  
}
