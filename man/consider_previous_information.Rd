\name{consider_previous_information}
\alias{consider_previous_information}
\title{
Remember the intermediate inferred GRN while generating the final inferred GRN
}
\description{
This function combines the adjacency matrix A_prev obtained as a result of \code{\link{first_GBM_step}} with the adjacency matrix A obtained as a result of \code{\link{second_GBM_step}}. All the edges in the matrix A which have non-zero weights are given machine precision weights initially. We then perform a harmonic mean for each element of A_prev and A to obtain a regularized adjacency matrix (A_final). As a result of this procedure transcriptional regulations which were strong and present in both A_prev and A end up getting highest weights in A_final. We finally remove all edges whose weights are less than machine precision from A_final. 
}
\usage{
consider_previous_information(A, A_prev,real)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
  Inferred GRN from the \code{\link{second_GBM_step}}
}
  \item{A_prev}{
  Inferred GRN from the \code{\link{first_GBM_step}}
}
  \item{real}{
  Numeric value 0 or 1 corresponding to simulated or real experiment respectively.
}
}
\value{
Returns an adjacency matrix A_final of the form Ntfs-by-Ntargets
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{first_GBM_step}}, \code{\link{second_GBM_step}}
}

\examples{
## The function is currently defined as
function (A, A_prev) 
{
  #Utilize Past Information also to not remove true positives
  A_prev[A_prev==0] <- .Machine$double.eps;
  A_prev <- transform_importance_to_weights(A_prev);
  A[A==0] <- .Machine$double.eps;
  epsilon <- 1/log(1/.Machine$double.eps);
  A <- transform_importance_to_weights(A);
  A_final <- 2*A*A_prev/(A+A_prev);
  A_final <- A_final - epislon;
  A_final[A_final<0] <- 0.0;
  return(A_final);
}
}
