xvalidate<-function (X, y, method=sdda, fold = NULL, trace = FALSE,...)
{
  if(is.null(colnames(X))) colnames(X)<-1:ncol(X)  # added patch
  
  
  n <- nrow(X)
  C <- rep(0, n)
  if (is.null(fold)) {
    grp <- 1:n
  }
  else {
    grp <- rep(1:fold, n/fold + 1)[1:n]
    grp[order(as.numeric(y) + runif(n))] <- grp
  }
  for (i in sort(unique(grp))) {
    if (trace)
      cat("X-validate group", i, " of ", max(grp), "\n")
    tt <- (grp == i)
    tmp <- method(X[!tt, ], y[!tt],...)
   
    C[tt] <- predict(tmp, X[tt, , drop = FALSE],...)
  }
  return(C)
}

