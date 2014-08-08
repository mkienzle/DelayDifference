// @(#)root/minuit2:$Id: GaussFcn.cxx 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "LogLikLinearReg.h"
#include <stdio.h>
#include <math.h>

#include <cassert>

namespace ROOT {

   namespace Minuit2 {


     // log-likelihood of a linear model with iid Gaussian errors that depends on 3 parameters
     // intercept and slope of a line and the standard error of a Gaussian random variable with mean = 0
     double LogLikLinearReg::operator()(const std::vector<double>& par) const {
  
       assert(par.size() == 3);

       double ss = 0., loglik = 0., sigma = 0.;

       for(unsigned int n = 0; n < fYvalues.size(); n++) 
	 {
	   // Calculate the squared difference between the linear model of X-values and the observed Y-values
	   ss += pow(line( &fXvalues[n], &par[0])- fYvalues[n],2);
	 }

       // Std. dev. of the error around the line (sigma)
       sigma = par[2];

       // Quantity to maximize, the log-likelihood of a Gaussian distributed variable with mean given by a linear function
       loglik = fYvalues.size() * log(sqrt(2 * 3.1416) * sigma) + 1/(2 * pow(sigma,2)) * ss;
  return loglik;
}


  }  // namespace Minuit2

}  // namespace ROOT
