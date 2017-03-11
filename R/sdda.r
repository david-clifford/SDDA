
dlda <- function(X,y,priors) {
  n <- nrow(X)
  p <- ncol(X)

  ## minor changes to allow y to be a factor
  if(is.factor(y)) {
    G <- nlevels(y)
    yLevels <- levels(y)
  } else {
    G <- max(y)
    yLevels <- 1:G
  }

  inty <- y
  if(is.factor(y)) inty <- unclass(y)

  pp <- .C("dldawrap",as.integer(n),as.integer(p),as.integer(G),
           as.double(X),as.integer(inty),
           means=matrix(0.0,G,p),
           vars=double(p),
           counts=integer(G),
           DUP=FALSE,
           package="SDDA")

  if(missing(priors))   {
    priors <- rep(1.0,G)
    priors[which.max(pp$counts)] <- 1.0001
    priors <- priors/sum(priors)
  }

  pp <- c(pp[c("means","vars","counts")],list(priors=priors))

  rownames(pp$means) <- yLevels
  colnames(pp$means) <- colnames(X)
  names(pp$vars) <- colnames(X)
  names(pp$counts) <- yLevels

  class(pp) <- "dlda"
  return(pp)
}

dqda <- function(X,y,priors) {
  n <- nrow(X)
  p <- ncol(X)

  ## minor changes to allow y to be a factor
  if(is.factor(y)) {
    G <- nlevels(y)
    yLevels <- levels(y)
  } else {
    G <- max(y)
    yLevels <- 1:G
  }

  inty <- y
  if(is.factor(y)) inty <- unclass(y)

  pp <- .C("dqdawrap",as.integer(n),as.integer(p),as.integer(G),
           as.double(X),as.integer(inty),
           means=matrix(0.0,G,p),
           vars=matrix(0.0,G,p),
           counts=integer(G),
           DUP=FALSE,
           package="SDDA")

  if(missing(priors))   {
    priors <- rep(1.0,G)
    priors[which.max(pp$counts)] <- 1.0001
    priors <- priors/sum(priors)
  }

  pp <- c(pp[c("means","vars","counts")],list(priors=priors))

  rownames(pp$means) <- yLevels
  colnames(pp$means) <- colnames(X)
  colnames(pp$vars) <- colnames(X)
  rownames(pp$vars) <- yLevels
  names(pp$counts) <- yLevels

  class(pp) <- "dqda"
  return(pp)
}


predict.dlda <- function(object, newdata,
                         type=c("class","prob"),...)
{
  type <- match.arg(type)
  G <- length(object$counts)
  n <- nrow(newdata)
  p <- ncol(object$means)

  pp <- .C("preddldawrap",as.integer(p),as.integer(G),
           object$means, object$vars, object$counts,
           object$priors,
           as.integer(n), newdata, R=matrix(0.0,n,G),
           DUP=FALSE,
           package="SDDA")

  yLevels <- rownames(object$means)
  probs <- pp$R
  rownames(probs) <- rownames(newdata)
  colnames(probs) <- yLevels
  if (type=="prob") return(probs)
  else {
    ## return the level of y with the largest probability
    ret <- apply(probs,1,which.max)
    ret <- yLevels[ret]
    return(ret)
  }
}

predict.dqda <- function(object, newdata,
                         type=c("class","prob"),...)
{
  type <- match.arg(type)
  G <- length(object$counts)
  n <- nrow(newdata)
  p <- ncol(object$means)

  pp <- .C("preddldawrap",as.integer(p),as.integer(G),
           object$means, object$vars, object$counts,
           object$priors,
           as.integer(n), newdata, R=matrix(0.0,n,G),
           DUP=FALSE,
           package="SDDA")

  yLevels <- rownames(object$means)
  probs <- pp$R
  rownames(probs) <- rownames(newdata)
  colnames(probs) <- yLevels
  if (type=="prob") return(probs)
  else {
    ## return the level of y with the largest probability
    ret <- apply(probs,1,which.max)
    ret <- yLevels[ret]
    return(ret)
  }
}

sdda <- function(X,y,priors,start=rep(FALSE,ncol(X)),never=rep(FALSE,ncol(X)),method="lda",...) {

  if(method!="lda" & method!="qda") {
    print("Please choose method='lda' or method='qda'\n")
    stop()
  }

  if(method=="lda") {
    if(missing(priors)) tt <- dlda(X,y)
    else tt <- dlda(X,y,priors)
  }

  if(method=="qda") {
    if(missing(priors)) tt <- dqda(X,y)
    else tt <- dqda(X,y,priors)
  }

  return(sdda.dda(tt,X,y,start,never,...))
}

sdda.dda <- function(obj, X, y, start=rep(FALSE,p), never=rep(FALSE,p),
                 useprob=TRUE, usecache=TRUE, usexval=TRUE, jmax= -20) {
  p <- ncol(X)
  n <- nrow(X)
  G <- length(obj$counts)

  if( inherits(obj,"dlda")) lda <- 1
  else lda <- 0

  pflg <- ifelse(useprob,1,0)
  cflg <- ifelse(usecache,1,0)
  xflg <- ifelse(usexval,1,0)

  if(any(start)) never[start] <- FALSE

  jj <- max(1,abs(jmax))

  ## minor changes to allow y to be a factor
  inty <- y
  if(is.factor(y)) inty <- unclass(y)

  pp <- .C("sddawrap", as.integer(n), as.integer(p), as.integer(G),
           X, as.integer(inty),
           obj$means, obj$vars, obj$counts, obj$priors,
           start, never, S=logical(p),
           as.integer(lda), as.integer(pflg), as.integer(cflg),
           as.integer(xflg),as.integer(jmax),
           ecrit = double(jj), pcrit=double(jj),
           DUP=FALSE,
           package="SDDA")

  obj$S <- pp$S
  names(obj$S) <- colnames(X)
  obj$ecrit <- pp$ecrit
  obj$pcrit <- pp$pcrit

  class(obj) <- c("sdda",class(obj))
  return(obj)
}

predict.sddanull <- function(object,newdata,
                             type=c("class","prob"),...)
{

  type <- match.arg(type)
  G <- length(object$counts)
  priors <- object$priors

  n <- nrow(newdata)
  R <- matrix(priors,nrow=n,ncol=G,byrow=TRUE)
  rownames(R) <- rownames(newdata)
  colnames(R) <- rownames(object$means)

  if (type=="prob") return(R)
  else return(apply(R,1,which.max))
}


predict.sdda <- function(object, newdata, ...)
{
  S <- object$S
  if(all(!S)) return(predict.sddanull(object,newdata,...))

  if(inherits(object,"dlda")) VARS <- object$vars[S]
  else VARS <-   object$vars[,S,drop=F]

  objT <- list(means=object$means[,S,drop=F], vars= VARS,
               counts=object$counts, priors=object$priors)

  class(objT) <- class(object)[-1]

  predict(objT, newdata[,S,drop=F],...)


}

plotdiag <- function(obj) {
  k <- min(sum(obj$S)+1, length(obj$ecrit))
  op <- par(mfrow=c(2,1))
  plot(1:k, obj$pcrit[1:k])
  plot(1:k, obj$ecrit[1:k])
  par(op)
  invisible()
}

plot.sdda <- function(obj, X, y, ...) {
  S <- obj$S
  j <- sum(S)
  if(j==0) {
    warning("No genes selected - cannot plot")
    return(invisible())
  }
  if(j==1) {
    plot(y, X[,S])
  } else {
    require(MASS)
    ll <- lda(X[,S],unclass(y))
    xx <- X[,S]%*%ll$scaling
    pairs(xx,lab=unclass(y),col=unclass(y)+1)
  }
  invisible()
}

which.genes <- function(obj) return(which(obj$S))

summary.sdda <- function(object, ...) print.sdda(object,...)

print.sdda <- function(x,plot=FALSE,...) {
  cat("SDDA using",class(x)[-1],".\n")

  cat("n =",nrow(x$means),"samples and p =",ncol(x$means),"variables.\n")

  cat("Group levels are:",levels(as.factor(rownames(x$means))),".\n")

  cat(length(which.genes(x)),"variables are chosen in total.\n")

  cat("Variables are",names(which.genes(x)),"\n")

  if(plot) {
    cat("See plot diagnostics\n")
    plotdiag(x)
  }
}
