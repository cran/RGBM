\name{GBM.train}
\alias{GBM.train}
\title{
Train GBM predictor
}
\description{
This function trains a regression model for a given \code{X.train} feature matrix, \code{Y.train} response vector, and working parameters. A model returned by this function can be used to predict response for unseen data with \code{\link{GBM.test}} function.
}
\usage{
GBM.train(X.train, Y.train, s_f = 0.3, s_s = 1, lf =1, M.train = 5000, nu = 0.001)
}
\arguments{
  \item{X.train}{
Input N-by-p feature matrix of training samples. Columns correspond to features, rows correspond to samples.
}
  \item{Y.train}{
Input N-element response vector of training samples.
}
  \item{s_f}{
Sampling rate of features, 0<s_f<=1. Fraction of columns from X.train, which will be sampled without replacement to calculate each extesion in boosting model. By default it's 0.3.
}
  \item{s_s}{
Sampling rate of samples, 0<s_s<=1. Fraction of rows from X.train, which will be sampled with replacement to calculate each extension in boosting model. By default it's 1.
}
\item{lf}{
Loss function: 1-> Least Squares and 2 -> Least Absolute Deviation
}
  \item{M.train}{
Number of extensions in boosting model, e.g. number of iterations of the main loop of RGBM algorithm. By default it's 5000.
}
  \item{nu}{
Shrinkage factor, learning rate, 0<nu<=1. Each extension to boosting model will be multiplied by the learning rate. By default it's 0.001.
}
}
\value{
Regression model is a structure containing all the information needed to predict response for unseen data
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
\seealso{
\code{\link{GBM.test}}
}

