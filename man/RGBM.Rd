\name{RGBM}
\alias{RGBM}
\title{
Regularized Gradient Boosting Machine for inferring GRN
}
\description{
This function performs the proposed regularized gradient boosting machines for reverse engineering GRN. It allows the user to provide prior information in the form of a mechanistic network g_M and after generation of an initially inferred GRN using the core GBM model undergoes a pruning step. Here we detect and remove isolated nodes using the \code{\link{select_ideal_k}} function along with identification of the optimal set of transcription factors for each target gene. We then re-iterate through the GBM followed by the refinement step to generate the final re-constructed GRN.
}
\usage{
RGBM(E = matrix(rnorm(100), 10, 10), K = matrix(0, nrow(E), ncol(E)), 
     g_M = matrix(1, 10, 10), tfs = paste0("G", c(1:10)), 
     targets = paste0("G", c(1:10)), lf = 1, M = 5000, nu = 0.001, s_f = 0.3,
     no_iterations = 2, mink = 0, experimentid = 1, outputpath= "DEFAULT",
     sample_type = "Exp1_", real = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
  N-by-p expression matrix. Columns correspond to genes, rows correspond to experiments. E is expected to be already normalized using standard methods, for example RMA. Colnames of E is the set of all p genes and Ntfs represents the number of transcription factors and Ntargets represents the number of target genes.
}
  \item{K}{
  N-by-p initial perturbation matrix. It directly corresponds to E matrix, e.g. if K[i,j] is equal to 1, it means that gene j was knocked-out in experiment i. Single gene knock-out experiments are rows of K with only one value 1. Colnames of K is set to be the set of all genes.  By default it's a matrix of zeros of the same size as E, e.g. unknown initial perturbation state of genes.
}
  \item{g_M}{
  Initial mechanistic network in the form of an adajcency matrix (Ntf-by-Ntargets). Here each column is a binary vector where only those elements are 1 when the corresponding transcription factor has a connection with that target gene. Colnames of g_M should be same as names of targets and Rownames of g_M should be same as names of Tfs. By default it's a matrix of ones of size Ntfs x Ntargets.
}
  \item{tfs}{
  List of names of transcription factors
}
  \item{targets}{
  List of names of target genes
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
  Number of times initial GRN to be constructed and then averaged to generate smooth edge weights for the initial GRN as shown in \code{\link{first_GBM_step}}
}
  \item{mink}{
  specified threshold i.e. the minimum number of Tfs to be considered while optimizing the L-curve criterion. By default it's 0. 
}
  \item{experimentid}{
  The id of the experiment being conducted. It takes natural numbers like 1,2,3 etc. By default it's 1.
}
  \item{outputpath}{
  Location where intermediate Adjacency_Matrix and Images folder will be created.  By default 
  it's a temp directory (e.g. /tmp/Rtmp...)
}
  \item{sample_type}{
    String arguement representing a label for the experiment i.e. in case of DREAM3 challenge sample_type="DREAM3".
}
  \item{real}{
    Numeric value 0 or 1 corresponding to simulated or real experiment respectively.
  }
}
\value{
Returns the final inferred GRN of form Ntfs-by-Ntargets adjacency matrix.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}


\seealso{
\code{\link{select_ideal_k}}, \code{\link{first_GBM_step}}
}
\examples{
# load RGBM library
library("RGBM")
# this step is optional, it helps speed up calculations, run in parallel on 2 processors
library(doParallel)
cl <- makeCluster(2)
# run network inference on a 100-by-100 dummy expression data.
A = RGBM()
stopCluster(cl)

}