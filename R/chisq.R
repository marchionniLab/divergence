
### =====================================================
### compute chi-squared test
### =====================================================

rowChiSq = function(obs, exp) {
    if(length(unique(obs)) > 1) {
        chisq.test(xtabs(~obs + exp))
    } else {
      list(statistic=NA, p.value=NA)
    }
}

chiSquaredTest = function(Mat, Groups, classes){
    
    if(! all(classes %in% Groups))
    		stop("Not all given classes are available in the groups variable")

    ch = apply(Mat, 1, rowChiSq, Groups)
    
    ch2 = data.frame(t(sapply(ch, function(x) c(x$statistic, x$p.value))))
  	colnames(ch2) = c("statistic", "pval")
  	ch2[order(ch2$pval, na.last=TRUE), ]

}
