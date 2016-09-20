pkgname <- "DelayDiffForR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('DelayDiffForR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("FitLW")
### * FitLW

flush(stderr()); flush(stdout())

### Name: FitLW
### Title: A function to fit a length-weight relationship
### Aliases: FitLW

### ** Examples

data(lw)
FitLW(lw[,1], lw[,2], lw[,3])

# customize synthetic dataset
length<-runif(1000,20,80);
a<-0.75; b<-2 ;
weight<-a*length^b;lengthbin<-floor(length);
data<-cbind(as.numeric(names(tapply(weight,lengthbin+0.5,mean))),tapply(weight,lengthbin+0.5,mean),sqrt(tapply(weight,lengthbin+0.5,var)));
FitLW(data[,1], data[,2], data[,3])




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
