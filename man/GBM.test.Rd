\name{GBM.test}
\alias{GBM.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test GBM predictor
}
\description{
This function tests a regression model for a given \code{X.test} feature matrix, \code{Y.test} response vector, and working parameters.
}
\usage{
GBM.test(model, X.test, Y.test, M.test)
}
\arguments{
  \item{model}{
Model returned by \code{\link{GBM.train}} function.
}
  \item{X.test}{
Input N-by-p feature matrix of unseen samples. Columns correspond to features, rows correspond to samples.
}
  \item{Y.test}{
Input N-element response vector of unseen samples.
}
  \item{M.test}{
Number of extensions of boosting model to take when predicting response. Must be not greater than \code{M.train} used when training boosting model.
}
}
\value{
Result of regression
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

\seealso{
\code{\link{GBM.train}}
}
