// CREATED  10 Sep 2013
// MODIFIED 16 Sep 2016

// PURPOSE log-likelihood function to fit a delay difference model to a time-series of catch
//         assuming the predicted square root of catches is distributed as a Gaussian with mean square root of observed catched and 
//         standard deviation (sigma) to be estimated


// MODIFIED FROM @(#)root/minuit2:$Id: DemoGaussSim.cxx 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  
// Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT 

# include "ObjectiveFCN.h"
# include <stdio.h>
# include <math.h>
# include <iostream>

# include <cassert>

       // Global variables
       //extern double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
       //extern double M;

namespace ROOT {

   namespace Minuit2 {


     // log-likelihood 
     double ObjectiveFCN::operator()(const std::vector<double>& par) const {
  
       double CHISQ=0.;

       // Calculate the difference between model and dataCalculate the sum of square of differnce between *** square root *** of estimated and observed catch
       for(unsigned int i = 0; i < fLength.size(); i++) 
	 {
	   CHISQ+= pow((fWeight[i] - par[0] * pow(fLength[i], par[1])) / fsd[i], 2.0);
	   //std::cout << "The value of the chisq is" << CHISQ << std::endl;
	 }

       return CHISQ;
}

  }  // namespace Minuit2

}  // namespace ROOT
