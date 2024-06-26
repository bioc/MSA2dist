% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aa2selfscore.R
\name{aa2selfscore}
\alias{aa2selfscore}
\title{aa2selfscore}
\usage{
aa2selfscore(aa, scorematrix = "BLOSUM62")
}
\arguments{
\item{aa}{\code{AAStringSet} [mandatory]}

\item{scorematrix}{score matrix to use [default: BLOSUM62]}
}
\value{
\code{data.frame}
}
\description{
This function return the selfscore from an \code{AAStringSet}.
}
\examples{
data(woodmouse, package="ape")
#cds2aa(dnabin2dnastring(woodmouse), shorten=TRUE,
#genetic.code=Biostrings::getGeneticCode("2"))
woodmouse |> dnabin2dnastring() |> cds2aa(shorten=TRUE,
genetic.code=Biostrings::getGeneticCode("2")) |> aa2selfscore()
}
\seealso{
\code{\link[Biostrings]{XStringSet-class}},
\code{\link[pwalign]{substitution_matrices}}
}
\author{
Kristian K Ullrich
}
