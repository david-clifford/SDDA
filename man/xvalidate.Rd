\name{xvalidate}
\alias{xvalidate}

\title{ Produces cross validated fitted values for sdda.}

\description{Ten-fold cross-validation, or leave one out cross
  validation (default).}

\usage{xvalidate(X, y, method=sdda, fold = NULL, trace = FALSE, ...)}

\arguments{
  \item{X}{real matrix with n rows (samples) and p columns (variables or
genes)}
  \item{y}{vector of class labels}
  \item{method}{sdda is the only option in this version of the library}
  \item{fold}{an integer specifying the number of folds to use in the
    cross validation. Default value NULL gives leave one out cross
    validation.}
  \item{trace}{a logical variable with value TRUE if information about
    cross validation progress to be output and value FALSE if no output
    required.}
  \item{\dots}{Other terms as appropriate} }

\details{ The data is divided into fold groups. The fitted value for
  each group is obtained by using a model built from the data for the
  remaining groups.  For sdda the usual cross validated class labels are
  produced.  }

\value{ An n by 1 vector of cross validated fitted values.}

\author{ Glenn Stone }

\note{  Missing values in X and y are not allowed}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
