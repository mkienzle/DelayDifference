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

# customize synthetic dataset
length<-runif(1000,20,80);
a<-0.75; b<-2 ;
weight<-a*length^b;lengthbin<-floor(length);
data<-cbind(as.numeric(names(tapply(weight,lengthbin+0.5,mean))),tapply(weight,lengthbin+0.5,mean),sqrt(tapply(weight,lengthbin+0.5,var)));
FitLW(data[,1], data[,2], data[,3])

}
