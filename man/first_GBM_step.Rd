\name{first_GBM_step}
\alias{first_GBM_step}
\title{
Perform either LS-Boost or LAD-Boost (\code{\link{GBM}}) on expression matrix E followed by the \code{\link{null_model_refinement_step}}
}
\description{
This function utilizes the core gradient boosting machine model (\code{\link{GBM}}) followed by the refinement step to generate the first adjacency matrix A of size p x p using the list of Tfs and the set of target genes. Several such adjacency matrices (A) are obtained based on the number of iterations to be performed. All these adjacency matrices are averaged to reduce the noise in the inferred intermediate GRN.
}
\usage{
first_GBM_step(E, K, tfs, targets, Ntfs, Ntargets, lf, M, nu,s_f, no_iterations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
  N-by-p expression matrix. Columns correspond to genes, rows correspond to experiments. E is expected to be already normalized using standard methods, for example RMA. Colnames of E is the set of all genes. 
}
  \item{K}{
  N-by-p initial perturbation matrix. It directly corresponds to E matrix, e.g. if K[i,j] is equal to 1, it means that gene j was knocked-out in experiment i. Single gene knock-out experiments are rows of K with only one value 1. Colnames of K is set to be the set of all genes.  By default it's a matrix of zeros of the same size as E, e.g. unknown initial perturbation state of genes.

}
  \item{tfs}{
  List of names of transcription factors. In case of presence of prior mechanistic network it is a subset of all the p genes whereas in absence of such a mechanistic network it is a list of names of all the p genes. 
}
  \item{targets}{
  List of names of target genes. In case of presence of prior mechanistic network it is a subset of all the p genes whereas in absence of such a mechanistic network it is a list of names of all the p genes. 
}
  \item{Ntfs}{
  Total number of transcription factors used in the experiment.
}
  \item{Ntargets}{
  Total number of target genes used in the experiment.
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
  \item{no_iterations}{
  Number of iterations to perform equivalent to building that many core LS-Boost/LAD-Boost models and then averaging them to have smooth edge-weights in the inferred intermediate GRN. 
}
}
\value{
Intermediate Gene Regulatory Network in form of a Ntfs-by-Ntargets adjacency matrix.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}


\seealso{
\code{\link{second_GBM_step}}
}
