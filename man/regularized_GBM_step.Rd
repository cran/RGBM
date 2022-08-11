\name{regularized_GBM_step}
\alias{regularized_GBM_step}
\title{
Perform the regularized GBM modelling once the initial GRN is inferred
}
\description{
This function undertakes all the proposed steps for regularizing the list of transcription factors for individual target gene followed by re-iterating through the core GBM model and the refinement step to produce the final reverse engineered GRN.
}
\usage{
regularized_GBM_step(E, A_prev, K, tfs, targets, Ntfs, Ntargets, lf, M, nu, s_f,
                     experimentid, outputpath, sample_type, mink=0,real=0)
}
\arguments{
  \item{E}{
  N-by-p expression matrix. Columns correspond to genes, rows correspond to experiments. E is expected to be already normalized using standard methods, for example RMA. Colnames of E is the set of all genes. 
}
  \item{A_prev}{
  An intermediate inferred GRN obtained from \code{\link{first_GBM_step}}
}
  \item{K}{
  N-by-p initial perturbation matrix. It directly corresponds to E matrix, e.g. if K[i,j] is equal to 1, it means that gene j was knocked-out in experiment i. Single gene knock-out experiments are rows of K with only one value 1. Colnames of K is set to be the set of all genes. By default it's a matrix of zeros of the same size as E, e.g. unknown initial perturbation state of genes.
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
  \item{experimentid}{
  The id of the experiment being conducted. It takes natural numbers like 1,2,3 etc. By default it's 1.
}
  \item{outputpath}{
  Location where the Adjacency_Matrix and Images folder will be created.
}
  \item{sample_type}{
  String arguement representing a label for the experiment i.e. in case of DREAM3 challenge sample_type="DREAM3".
}
  \item{mink}{
   User specified threshold i.e. the minimum number of Tfs to be considered while optimizing the L-curve criterion. By default it's 0. 
}
  \item{real}{
    Numeric value 0 or 1 corresponding to simulated or real experiment respectively.
  }
}
\value{
Returns the final inferred GRN in form of Ntfs-by-Ntargets matrix
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}


\seealso{
\code{\link{first_GBM_step}}
}
