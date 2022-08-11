\name{select_ideal_k}
\alias{select_ideal_k}
\title{
Identifies the optimal value of k i.e. top k Tfs for each target gene
}
\description{
This function detects the optimal number of transcription factors which are regulating each target gene. This number is different for different target genes. It utilizes a heuristic to also detect the isolated targets which are not regulated by any transcription factor. To the detect the optimal number of Tfs for each target gene, it uses a notion similar to that used for optimization of the L-curve criterion for Tikonov regularization by evaluating the variable importance curve for each target gene.
}
\usage{
select_ideal_k(experimentid, mink, filepath, imagepath, adjacency_matrix_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experimentid}{
  The id of the experiment being conducted. It takes natural numbers like 1,2,3 etc. By default it's 1.
}
  \item{mink}{
  User specified threshold i.e. the minimum number of Tfs to be considered while optimizing the L-curve criterion. By default it's 0. 
}
  \item{filepath}{
  Path where some intermediate files will be written and provided by the function \code{\link{get_filepaths}}.
}
  \item{imagepath}{
    Path where an image of the variable importance curves for first 16 target genes will be written and provided by the function \code{\link{get_filepaths}}.

}
  \item{adjacency_matrix_path}{
    Path where an intermediate adjacency matrix will be written and provided by the function \code{\link{get_filepaths}}.
}
}
\value{
Returns a vector where each element represents the optimal number of transcription factors for each target gene.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}


