\name{regulate_regulon_size}
\alias{regulate_regulon_size}
\title{
Regulate the size of the regulon for each TF
}
\description{
We control the size of the regulon for each TF by using a heuristic to remove the edges whose weights are small
}
\usage{
regulate_regulon_size(A)
}
\arguments{
  \item{A}{
  Inferred GRN in the form of Ntfs-by-Ntargets matrix
}
}
\value{
Refined adjacency matrix A in the form of Ntfs-by-Ntargets matrix
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
