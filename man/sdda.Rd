\name{sdda}
\alias{sdda}

\title{Stepwise Diagonal Discriminant Analysis}

\description{ Uses a forwards stepwise strategy to add variables one at
  a time to a diagonal discriminant rule.}

\usage{sdda(X, y, priors, start = rep(FALSE, ncol(X)), never = rep(FALSE, ncol(X)),method="lda",...)  }

\arguments{
  \item{X}{Training data matrix - rows are observations, columns are
    variables.}
  \item{y}{A factor of true class labels, or a numeric vector with
    values 1, 2, 3, ... G where G is the number of classes.}
  \item{priors}{Prior probabilities for the different classes, if left
    unspecified these default to equal probability to belong to each group}
  \item{start}{Logical vector indicating which variables (if any) should
    start in the rule.}
  \item{never}{Logical vector indicating which variables (if any) should
    never be in the rule.}
   \item{method}{Choose linear or quadratic discriminant analysis, "lda"
     or "qda"}
   \item{...}{Other terms as appropriate}
 }
 
\details{ Uses a forwards stepwise strategy to add variables one at a
  time to a diagonal discriminant rule. At each step the variable that
  most decreases the (leave-one-out) cross-validated error rate is
  added. Adding stops when the cross-validated error rate cannot be
  decreased.  }

\value{
  \item{means}{The means for each variable within each group}
  \item{vars}{Under DLDA (default) this is a vector of the variances of
    each variable, under DQDA this is a matrix of variances for each
    variable within each group}
  \item{counts}{The number of samples in each group}
  \item{priors}{The prior probability of being in each group}
  \item{S}{Vector of flags of the variables in the rule, including
    \code{start}}
  \item{ecrit}{X-validated error rate criterion. Used for diagnostics, see \code{\link{plotdiag}}}
  \item{pcrit}{X-validated likelihood criterion. Used for diagnostics, see \code{\link{plotdiag}}}
}

\author{Glenn Stone}

\examples{
  data(iris)
  Y <- iris[,5]
  noise <- matrix(rnorm(50*nrow(iris)),nrow=nrow(iris))
  colnames(noise) = paste("noise",1:ncol(noise),sep="")
  X <- cbind(as.matrix(iris[,1:4]),noise)
  dim(X)

  s1 <- sdda(X,Y)

  which.genes(s1)
  plotdiag(s1)
  summary(s1)
  
  ## Probabilities for each class
  predict(s1,newdata=X,type="prob")
  rs <- predict(s1,X)
  table(Y,rs)

  ## Leave-one-out crossvalidation
  xv <- xvalidate(X,Y,trace=TRUE)
  table(Y,xv)
}

\seealso{\code{\link{xvalidate}},\code{\link{predict.sdda}}, \code{\link{plotdiag}}}

\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
