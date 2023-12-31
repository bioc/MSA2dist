% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/region.R
\name{region}
\alias{region}
\title{region}
\usage{
region(seq)
}
\arguments{
\item{seq}{\code{DNAStringSet} or \code{AAStringSet} [mandatory]}
}
\value{
region \code{IRanges} object from \code{metadata}
}
\description{
This function shows the region slot from a
\code{DNAStringSet} or an \code{AAStringSet}
\code{metadata} information.\cr
}
\examples{
## load example sequence data
data(iupac, package="MSA2dist")
iupac.aa <- iupac |> cds2aa(shorten = TRUE)
## create region
region1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
## add region
iupac.aa <- iupac.aa |> addregion2string(region=region1)
iupac.aa |> region()
}
\seealso{
\code{\link[MSA2dist]{addpop2string}}
}
\author{
Kristian K Ullrich
}
