\name{get_filepaths}
\alias{get_filepaths}
\title{
Generate filepaths to maintain adjacency matrices and images
}
\description{
This function generates a set of filepaths which are used to keep the adjacency matrix A obtained after the \code{\link{first_GBM_step}} + \code{\link{null_model_refinement_step}}. It also generates a path where an image of the variable importance curves for several target genes can be kept.
}
\usage{
get_filepaths(A_prev, experimentid, outputpath, sample_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A_prev}{
  Adjacency matrix A obtained after \code{\link{first_GBM_step}} + \code{\link{null_model_refinement_step}}.
}
  \item{experimentid}{
  The id of the experiment being conducted. It takes natural numbers like 1,2,3 etc. By default it's 1.
}
  \item{outputpath}{
  Location where the Adjacency_Matrix and Images folder will be created.
}
  \item{sample_type}{
  String arguement representing a label for the experiment i.e. in case of DREAM3 challenge sample_type="DREAM3".
}
}
\value{
Returns a data frame where the first element in the data frame is the location where the Adjacency_Matrix folder is located in the filesystem, second element represents the location where the Images folder is located in the filesystem, third element represents the path to the file where the Adjacency_Matrix will be written.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
