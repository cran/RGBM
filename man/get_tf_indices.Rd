\name{get_tf_indices}
\alias{get_tf_indices}
\title{
Get the indices of all the TFs from the data
}
\description{
This function provides the indices of all the transcription factors which are present in the expression matrix. In case of DREAM Challenges it will return the indices as 1,...,p for all the p genes in the data as the transcription factors are not known beforehand.
}
\usage{
get_tf_indices(E, tfs, Ntfs)
}
\arguments{
  \item{E}{
  E is the expression matrix of size N x p where N is number of examples and p is the number of genes. Here the column names of expression matrix is the list of all the genes present in the E matrix. Colnames of E is the set of all genes. 
}
  \item{tfs}{
  List of names of transcription factors.
}
  \item{Ntfs}{
  Total number of transcription factors used in the experiment.
}
}
\value{
Returns the indices of all the transcription factors present in E matrix.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

\seealso{
\code{\link{get_colids}}
}

