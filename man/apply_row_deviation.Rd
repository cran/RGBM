\name{apply_row_deviation}
\alias{apply_row_deviation}
\title{
Apply row-wise deviation on the inferred GRN
}
\description{
This function performs a row-wise standard deviation of network A to generate an S1 matrix which is then used to modify the weights in network A}
\usage{
apply_row_deviation(A,Ntfs,Ntargets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
  Inferred GRN in the form of Ntfs-by-Ntargets matrix
}
   \item{Ntfs}{
  Total number of transcription factors used in the experiment.
}
  \item{Ntargets}{
  Total number of target genes used in the experiment
}
}
\value{
Refined adjacency matrix A in the form of Ntfs-by-Ntargets matrix
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
