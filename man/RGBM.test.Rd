\name{RGBM.test}
\alias{RGBM.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test rgbm predictor
}
\description{
This function tests a regression model for a given \code{X.test} feature matrix, \code{Y.test} response vector, and working parameters.
}
\usage{
RGBM.test(model, X.test, Y.test, M.test)
}
\arguments{
  \item{model}{
Model returned by \code{\link{RGBM.train}} function.
}
  \item{X.test}{
Input S-by-P feature matrix of unseen samples. Columns correspond to features, rows correspond to samples.
}
  \item{Y.test}{
Input S-element response vector of unseen samples.
}
  \item{M.test}{
Number of extensions of boosting model to take when predicting response. Must be not greater than \code{M.train} used when training boosting model.
}
}
\value{
Result of regression
}
\author{
Raghvendra Mall <raghvendra5688@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
