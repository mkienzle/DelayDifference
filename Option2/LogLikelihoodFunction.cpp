// CREATED  10 Sep 2013
// MODIFIED 26 May 2014

// PURPOSE log-likelihood function to fit a delay difference model to a time-series of catch
//         assuming the predicted square root of catches is distributed as a Gaussian with mean square root of observed catched and 
//         standard deviation (sigma) to be estimated


// MODIFIED FROM @(#)root/minuit2:$Id: DemoGaussSim.cxx 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  
// Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT 

# include "LogLikelihoodFunction.h"
# include "DelayDifference.h"
# include <stdio.h>
# include <math.h>
# include <iostream>

# include <cassert>

       // Global variables
       extern double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
       extern double M;

namespace ROOT {

   namespace Minuit2 {


     // log-likelihood 
     double LogLikelihoodFunction::operator()(const std::vector<double>& par) const {
  
       //       assert(par.size() == 5);
       double ss = 0., loglik = 0., sigma = 0.;
       
       // Scale parameter
       double Targeted_catchability = par[0] * CatchabilityScalingFactor;
       double Nontargeted_catchability = par[1] * CatchabilityScalingFactor;
       sigma = par[2]; // The Std. dev. of the Gaussian distribution of predicted sqrt catch around observed sqrt catch (sigma)

       // Calculate biomasses according to the delay difference model
       std::vector<double> Biomass(fCatchValues.size(), 0.0), FishingMortality(fCatchValues.size(), 0.0), EstimatedCatch(fCatchValues.size(), 0.0);

       DelayDifference(fTargetedEffortValues, fNontargetedEffortValues, Biomass, par);

       // Calculate estimated catch using (a) biomasses from the delay difference and (b) fishing mortalities
       for(unsigned int counter=0; counter < fCatchValues.size(); counter++)
      {

	// Estimated catch
	FishingMortality[counter] = Targeted_catchability * fTargetedEffortValues[counter] + Nontargeted_catchability * fNontargetedEffortValues[counter];

	EstimatedCatch[counter] = FishingMortality[counter] / (M + FishingMortality[counter]) * \
	  Biomass[counter] * (1 - exp(- (M + FishingMortality[counter])));
}

       // Calculate the sum of square of differnce between *** square root *** of estimated and observed catch
       for(unsigned int n = 0; n < fCatchValues.size(); n++) 
	 {
	   	   ss+=pow( sqrt(EstimatedCatch[n]) - sqrt(fCatchValues[n]), 2);
	 }

       //// objective function: negative log-likelihood of a Gaussian distributed variable
       loglik = fCatchValues.size() * log(sqrt(2 * M_PI) * sigma) + 1/(2 * pow(sigma,2)) * ss;
       //std::cout << "Sigma is  " << sigma << "\n";
       
       // NOTE: the residual sum of square (sigma) is (by some scientists) calculated rather than estimated
       //loglik = fCatchValues.size() * log(sqrt(2 * M_PI) * ss / (double) (fCatchValues.size() -2 )) + 1/(2 * pow(sigma,2)) * ss;
       //std::cout << "Residual sum of square is  " << ss / (double) ( fCatchValues.size() - 2) << "\n";

       return loglik;
}

  }  // namespace Minuit2

}  // namespace ROOT
