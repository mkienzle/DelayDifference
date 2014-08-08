// CREATED  10 Sep 2013
// MODIFIED  7 Nov 2013

// PURPOSE log-likelihood objective function to fit a weekly delay difference model to a time-series of catch
//         assuming the predicted square root of catches is distributed as a Gaussian with mean square root of observed catched and 
//         standard deviation (sigma) to be estimated


// MODIFIED FROM @(#)root/minuit2:$Id: DemoGaussSim.cxx 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  
// Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT 

# include "LogLikelihood_2par_Function.h"
# include "2par_WeeklyDelayDifference.h"
# include <stdio.h>
# include <math.h>

# include "FixPar_2parVersion.h"

# include <cassert>

namespace ROOT {

   namespace Minuit2 {


     // log-likelihood 
     double LogLikelihoodFunction::operator()(const std::vector<double>& par) const {
  
       assert(par.size() == 2);

       double ss = 0., loglik = 0., catchability, sigma = 0.;

       // Scale parameter
       catchability = par[0] * CatchabilityScalingFactor;
       sigma = par[1]; // The Std. dev. of the Gaussian distribution of predicted catch around observed catch (sigma)

       // Calculate biomasses according to the delay difference model
       std::vector<double> Biomass(fCatchValues.size(), 0.0), FishingMortality(fCatchValues.size(), 0.0), EstimatedCatch(fCatchValues.size(), 0.0);
       //double Biomass[fCatchValues.size()];
       //double FishingMortality[fCatchValues.size()];
       //double EstimatedCatch[fCatchValues.size()];

       for(unsigned int counter=0; counter < fCatchValues.size(); counter++)
      {
	FishingMortality.at(counter) = catchability * fEffortValues[counter];
	Biomass.at(counter) = WeeklyDD_2par(counter, fCatchValues.size(), fEffortValues, par);
	EstimatedCatch.at(counter) = FishingMortality[counter] / (M + FishingMortality[counter]) * \
	  Biomass[counter] * (1 - exp(- (M + FishingMortality[counter])));
}

       // Calculate the sum of square of differnce between *** square root *** of estimated and observed catch
       for(unsigned int n = 0; n < fCatchValues.size(); n++) 
	 {
	   ss+=pow( sqrt(EstimatedCatch[n]) - sqrt(fCatchValues[n]), 2);
	 }

       // negative log-likelihood of a Gaussian distributed variable (objective function)
       loglik = fCatchValues.size() * log(sqrt(2 * M_PI) * sigma) + 1/(2 * pow(sigma,2)) * ss;

       return loglik;
}

  }  // namespace Minuit2

}  // namespace ROOT
