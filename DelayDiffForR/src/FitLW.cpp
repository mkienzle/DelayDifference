// CREATED  14 Sep 2016
// MODIFIED 19 Sep 2016

// Example of passing a vector
// From https://www.r-bloggers.com/calling-c-from-r-using-rcpp

// Compile using
// PKG_CXXFLAGS=`Rscript -e 'Rcpp:::CxxFlags()'` PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` R CMD SHLIB FitLW.cpp

// g++ -I/usr/local/lib/R/include -DNDEBUG  -I/usr/local/include   -I/home/mkienzle/R/x86_64-unknown-linux-gnu-library/3.2/Rcpp/include -fpic  -g -O2  -c FitLW.cpp -o FitLW.o
// g++ -shared -L/usr/local/lib -o FitLW.so FitLW.o


// in R
// library(Rcpp)
// .libPaths(c( .libPaths(), "/home/mkienzle/mystuff/Software/root-6.06.06/lib/", "/usr/local/lib64/"))
// dyn.load('~/mystuff/Software/root-6.06.06/lib/libMinuit2.so')
// dyn.load('FitLW.so')
// my.data <- read.table("lw.dat")
// .Call('FitLW', my.data$V1, my.data$V2, my.data$V3)

# include <Rcpp.h>
# include <cstdlib>
# include <iostream>

#include "lib_facilities.h"
#include "ObjectiveFCN.h"

//using namespace std;
using namespace ROOT::Minuit2;

RcppExport SEXP FitLW(SEXP x, SEXP y, SEXP z){
  int i,n;
  Rcpp::NumericVector length(x);
  Rcpp::NumericVector weight(y);
  Rcpp::NumericVector sd(z);
  n=length.size();

  int status=0;
  
  // Declare the objective function (also called FCN function in documentation)  
  ObjectiveFCN fFCN(length, weight, sd);

  // create Minuit parameters with names
  MnUserParameters upar;
  upar.Add("a", 0.001, 0.001);
  upar.Add("b", 3, 1);

  // create MIGRAD minimizer with MnStrategy 0 (strategy to calculate first and second derivative
  // with fewer function calls -- less precise result)
  MnMigrad migrad(fFCN, upar, 0);

  // Minimize
  FunctionMinimum min = migrad();

  // output
  std::cout<<"minimum: "<< min << std::endl;

  return(Rcpp::wrap(status));
  
     }



