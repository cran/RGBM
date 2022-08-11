\name{transform_importance_to_weights}
\alias{transform_importance_to_weights}
\title{
Log transforms the edge-weights in the inferred GRN
}
\description{
This function performs an inverse absolute log-transformation of the non-zero edge weights in the final inferred GRN (A) to make the edge-weights more comprehensible and understandable.
}
\usage{
transform_importance_to_weights(A)
}
%- maybe also 'usage' for other objects documented here.
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

