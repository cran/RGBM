\name{normalize_matrix_colwise}
\alias{normalize_matrix_colwise}
\title{
Column normalize the obtained adjacency matrix
}
\description{
We perform a column normalization on an adjacency matrix A equivalent to inferred GRN 
}
\usage{
normalize_matrix_colwise(A,Ntargets)
}
\arguments{
  \item{A}{
  Inferred GRN in the form of Ntfs-by-Ntargets matrix
}
  \item{Ntargets}{
  Total number of target genes used in the experiment
}
}
\value{
Column Normalized GRN of size Ntfs-by-Ntargets
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
