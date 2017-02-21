\name{get_colids}
\alias{get_colids}
\title{
Get the indices of recitifed list of Tfs for individual target gene 
}
\description{
This function is used to identify the recitified list of transcription factors for individual target genes after analysing the variable importance scores (where non-essential Tfs are pruned). These list of Tfs are usually different for individual target genes. Hence we maintain this in the form an adjacency matrix where the rownames correspond to all the Tfs and colnames correspond to all the target genes. Each column is a binary vector where all the values corresponding to the rectified Tfs active for that target are 1 while rest of the values are zeros.
}
\usage{
get_colids(A, ideal_k, tfs, targets, Ntfs, Ntargets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
  Adjacency Matrix A obtained after the GBM and refinement step.
}
  \item{ideal_k}{
  A vector containing the optimal value of k (no of active TFs) for each target gene obtained from \code{\link{select_ideal_k}}.
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
  Total number of target genes used in the experiment.
}
}
\value{
The function returns an adjacency matrix where the rownames correspond to all the Tfs and colnames correspond to all the target genes. Each column is a binary vector where all the values corresponding to the rectified Tfs active for that target are 1 while rest of the values are zeros.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}


\seealso{
\code{\link{get_tf_indices}}
}

