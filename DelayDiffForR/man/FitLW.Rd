\name{FitLW}
\alias{FitLW}
\title{
A function to fit a length-weight relationship
}
\description{
Fitting is done by minimizing the chisquare sum([(obs-pred)/sd]^2)
}
\usage{FitLW(length, weight, sd)}
\examples{
data(lw)
FitLW(lw[,1], lw[,2], lw[,3])
}
