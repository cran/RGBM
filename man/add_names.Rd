\name{add_names}
\alias{add_names}
\title{
Add row and column names to the adjacency matrix A
}
\description{
Here we add the names of the transcription factors (Tfs) as rownames and names of the target genes as column names to the adjacency matrix A.
}
\usage{
add_names(A, tfs, targets)
}
\arguments{
  \item{A}{
  Adjacency matrix A obtained as a result of GBM procedure.
}
  \item{tfs}{
  List of names of transcription factors.
}
  \item{targets}{
  List of names of target genes.
}
}
\details{
In case of DREAM Challenge datasets list of transcription factors is same as list of target genes and are referred as G1, ..., G100.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
