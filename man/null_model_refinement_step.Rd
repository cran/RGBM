\name{null_model_refinement_step}
\alias{null_model_refinement_step}
\title{
Perform the null model refinement step
}
\description{
We used this function for refining the edge-weights in an inferred GRN (A) by utilizing matrix (S2) obtained from null-mutant zscore effect (\code{\link{z_score_effect}}) as shown in \cite{Slawek J, Arodz T} i.e. A = A x S2.
}
\usage{
null_model_refinement_step(E, A, K,tfs, targets, Ntfs, Ntargets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
  N-by-p expression matrix. Columns correspond to genes, rows correspond to experiments. E is expected to be already normalized using standard methods, for example RMA. Colnames of E is the set of all genes. 
}
  \item{A}{
  Intermediate GRN network in the form of a p-by-p adjacency matrix.
}
  \item{K}{
  N-by-p initial perturbation matrix. It directly corresponds to E matrix, e.g. if K[i,j] is equal to 1, it means that gene j was knocked-out in experiment i. Single gene knock-out experiments are rows of K with only one value 1. Colnames of K is set to be the set of all genes.  By default it's a matrix of zeros of the same size as E, e.g. unknown initial perturbation state of genes.
}
  \item{tfs}{
  List of names of transcription factors
}
  \item{targets}{
  List of names of target genes
}

  \item{Ntfs}{
  Number of transcription factors used while building the GBM (\code{\link{GBM}}) model.
}
  \item{Ntargets}{
  Number of targets used while building the GBM (\code{\link{GBM}}) model.
}
}
\value{
  Returns a refined adjacency matrix A in the form of a Ntfs-by-Ntargets matrix.
}
\references{
Slawek J, Arodz T. ENNET: inferring large gene regulatory networks from expression data using gradient boosting. BMC systems biology. 2013 Oct 22;7(1):1.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}


\seealso{
\code{\link{z_score_effect}}
}
