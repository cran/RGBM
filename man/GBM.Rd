\name{GBM}
\alias{GBM}
\title{
Calculate Gene Regulatory Network from Expression data using either LS-TreeBoost or LAD-TreeBoost
}
\description{
This function calculates a Ntfs-by-Ntargets adjacency matrix A from N-by-p expression matrix E. E is expected to be given as input. E is assumed to have p columns corresponding to all the genes, Ntfs represents the number of transcription factors and Ntargets represents the number of target genes and N rows corresponding to different experiments. Additionally, GBM function takes matrix of initial perturbations of genes K of the same size as E, and other parameters including which loss function to use (LS = 1, LAD = 2). As a result, GBM returns a squared matrix A of edge confidences of size Ntfs-by-Ntargets. A subset of known transcription factors can be defined as a subset of all p genes.
}
\usage{
GBM(E = matrix(rnorm(100), 10, 10), K = matrix(0, nrow(E), ncol(E)), 
     tfs = paste0("G",c(1:10)), targets = paste0("G",c(1:10)), 
     s_s = 1, s_f = 0.3, lf = 1, 
     M = 5000,nu = 0.001, scale = TRUE,center = TRUE, optimization.stage = 2)
}
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
  \item{s_s}{
Sampling rate of experiments, 0<s_s<=1. Fraction of rows of E, which will be sampled with replacement to calculate each extension in boosting model. By default it's 1.
}
  \item{s_f}{
Sampling rate of transcription factors, 0<s_f<=1. Fraction of transcription factors from E, as indicated by \code{tfs} vector, which will be sampled without replacement to calculate each extesion in boosting model. By default it's 0.3.
}
  \item{lf}{
Loss function: 1 -> Least Squares, 2 -> Least Absolute deviation
}
  \item{M}{
Number of extensions in boosting model, e.g. number of iterations of the main loop of RGBM algorithm. By default it's 5000.
}
  \item{nu}{
Shrinkage factor, learning rate, 0<nu<=1. Each extension to boosting model will be multiplied by the learning rate. By default it's 0.001.
}
  \item{scale}{
Logical flag indicating if each column of E should be scaled to be unit standard deviation. By default it's TRUE.
}
  \item{center}{
Logical flag indicating if each column of E should be scaled to be zero mean. By default it's TRUE.
}
  \item{optimization.stage}{
Numerical flag indicating if re-evaluation of edge confidences should be applied after calculating initial V, optimization.stage=\{0,1,2\}. If optimization.stage=0, no re-evaluation will be applied. If optimization.stage=1, variance-based optimization will be applied. If optimization.stage=2, variance-based and z-score based optimizations will be applied.
}
}
\value{
  \item{A}{
Gene Regulatory Network in form of a Ntfs-by-Ntargets adjacency matrix.
}
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{GBM.train}}, \code{\link{GBM.test}}, \code{\link{v2l}}
}
\examples{
# load RGBM library
library("RGBM")
# this step is optional, it helps speed up calculations, run in parallel on 2 processors
library(doParallel)
cl <- makeCluster(2)
# run network inference on a 100-by-100 dummy expression data.
V = GBM()
stopCluster(cl)
}
