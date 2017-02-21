\name{second_GBM_step}
\alias{second_GBM_step}
\title{
Re-iterate through the core GBM model building with optimal set of Tfs for each target gene 
}
\description{
This function re-performs the core GBM model building (only one time) using the optimal set of transcription factors obtained from \code{\link{select_ideal_k}} followed by \code{\link{get_colids}} for individual target gene to return a regularized GRN. 
}
\usage{
second_GBM_step(E, K, df_colids, tfs, targets, Ntfs, Ntargets, lf,  M,  nu, s_f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
  N-by-p expression matrix. Columns correspond to genes, rows correspond to experiments. E is expected to be already normalized using standard methods, for example RMA. Colnames of E is the set of all genes. 
}
  \item{K}{
  N-by-p initial perturbation matrix. It directly corresponds to E matrix, e.g. if K[i,j] is equal to 1, it means that gene j was knocked-out in experiment i. Single gene knock-out experiments are rows of K with only one value 1. Colnames of K is set to be the set of all genes. By default it's a matrix of zeros of the same size as E, e.g. unknown initial perturbation state of genes.
}
  \item{df_colids}{
  A matrix made up of column vectors where each column vector represents the optimal set of active Tfs which regulate each target gene and obtained from \code{\link{get_colids}}. Some column vectors are just made up of zeros indicating that corresponding target genes are isolated and not regulated by any Tf
}
  \item{tfs}{
  List of names of transcription factors.
}
  \item{targets}{
  List of names of target genes.
}
  \item{Ntfs}{
  Total number of transcription factors used in the experiment.
}
  \item{Ntargets}{
  Total number of target genes used in the experiment
}
  \item{lf}{
  Loss Function: 1 -> Least Squares and 2 -> Least Absolute Deviation  
}
  \item{M}{
Number of extensions in boosting model, e.g. number of iterations of the main loop of RGBM algorithm. By default it's 5000.
}
  \item{nu}{
Shrinkage factor, learning rate, 0<nu<=1. Each extension to boosting model will be multiplied by the learning rate. By default it's 0.001.
}
  \item{s_f}{
Sampling rate of transcription factors, 0<s_f<=1. Fraction of transcription factors from E, as indicated by \code{tfs} vector, which will be sampled without replacement to calculate each extesion in boosting model. By default it's 0.3.
}
}
\value{
Returns a regularized GRN of the form Ntfs-by-Ntargets
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

\seealso{
\code{\link{first_GBM_step}}
}
