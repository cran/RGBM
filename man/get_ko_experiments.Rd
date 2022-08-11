\name{get_ko_experiments}
\alias{get_ko_experiments}
\title{
Get indices of experiments where knockout or knockdown happened
}
\description{
This function provides the indices of all those samples (out of N) where it is known apriori that a gene was either knocked-out or was knocked-down. This information is useful for the \code{\link{null_model_refinement_step}} which utilizes the \code{\link{z_score_effect}} technique (with the help of this information).
}
\usage{
get_ko_experiments(K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{
  N-by-p initial perturbation matrix. It directly corresponds to E matrix, e.g. if K[i,j] is equal to 1, it means that gene j was knocked-out in experiment i. Single gene knock-out experiments are rows of K with only one value 1. Colnames of K is set to be the set of all genes.  By default it's a matrix of zeros of the same size as E, e.g. unknown initial perturbation state of genes.
}
}
\value{
Return a vector containing the indices of all the samples where a gene was knocked-out/down.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
\seealso{
\code{\link{null_model_refinement_step}}, \code{\link{z_score_effect}}
}
