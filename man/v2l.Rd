\name{v2l}
\alias{v2l}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert adjacency matrix to a list of edges
}
\description{
This function converts adjacency matrix A to a sorted list of edges, e.g. a list in which edges are sorted by decreasing confidence.
}
\usage{
v2l(A, max = 1e+05, check.names = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Input adjacency matrix.
}
  \item{max}{
Maximal length of the resulting list. This number may be lower than the number of all the edges from adjacency matrix. Then only top \code{max} edges will be returned.
}
  \item{check.names}{
  Checks name of the gene ids
  }
}
\value{
A data frame of sorted edges: (1) list of sources (2) list of destinations (3) list of confidences. Elements in all the lists correspond to each other.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}


