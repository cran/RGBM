\name{z_score_effect}
\alias{z_score_effect}
\title{
Generates a matrix S2 of size Ntfs x Ntargets using the null-mutant zscore algorithm \cite{Prill, Robert J., et al}
}
\description{
This function generates a matrix of the form Ntfs-by-Ntargets using the steps proposed in null-mutant zscore method and acts as a refinement step for the inferred GRN where this matrix is multiplied element by element with the inferred adjacency matrix A. However, this step is only effective in presence of additional source of information like knockout, knockdown or which genes are intially perturbed in time-series expression data.
}
\usage{
z_score_effect(E, K, tfs, targets, Ntfs, Ntargets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
  N-by-p expression matrix. Columns correspond to genes, rows correspond to experiments. E is expected to be already normalized using standard methods, for example RMA. Colnames of E is the set of all genes. 
}
  \item{K}{
N-by-p initial perturbation matrix. It directly corresponds to E matrix, e.g. if K[i,j] is equal to 1, it means that gene j was knocked-out in experiment i. Single gene knock-out experiments are rows of K with only one value 1. Colnames of K is set to be the set of all genes. By default it's a matrix of zeros of the same size as E, e.g. unknown initial perturbation state of genes.
}
  \item{tfs}{
  List of names of transcription factors
}
  \item{targets}{
  List of names of target genes
}

  \item{Ntfs}{
  Total number of transcription factors used in the experiment.
}
  \item{Ntargets}{
  Total number of target genes used in the experiment.
}
}
\value{
Returns an S2 matrix of form Ntfs-by-Ntargets. In absence of any additional knockout/knockdown/perturbation information the S2 matrix is a matrix of ones.
}
\references{
Prill, Robert J., et al. "Towards a rigorous assessment of systems biology models: the DREAM3 challenges." PloS one 5.2 (2010): e9202.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{null_model_refinement_step}}
}
