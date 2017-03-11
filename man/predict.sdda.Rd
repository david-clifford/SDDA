\name{predict.sdda}
\alias{predict.sdda}

\title{Predict functionality for stepwise diagonal discriminant analysis}

\description{Produces predicted classes or class probabilities for an
  \code{sdda} rule.}

\usage{ predict.sdda(object, newdata, ...)}

\arguments{
  \item{object}{A fitted \code{sdda} rule object.}
  \item{newdata}{A dataset from which we wish to make
    predictions. Unlike many predict methods this is required, as sdda
    does not store the original data due to potential memory problems.}
  \item{\dots}{parameters passed to the based \code{dlda} or \code{dqda}
    predict method. Specify \code{type} to be "class" or "prob" to
    return the predicted class or the predicted probabilities for each
    class, by default the predicted class is returned}
}

\value{The predicted class (by default), or the probabilites of being
from each class} 

\author{Glenn Stone}

\keyword{models}
\keyword{multivariate}
