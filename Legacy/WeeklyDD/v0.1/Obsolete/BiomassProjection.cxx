// CREATED   10 Sept 2013
// MODIFIED  19 Nov 2013

// PURPOSE calculate biomass using a delay difference model with weekly timesteps

# include <math.h>
# include "FixPar.h"
# include <vector>
# include <iostream>
# include <cassert>

std::vector<double> vonMisesRecDist(double a, double b);

// STATUS working, converging to the simulated parameters

// ARGUMENTS 1. the time step to evaluate the biomass at
//           2. the total number of time steps considered in the analysis
//           3. a vector of targeted effort
//           4. a vector of non-targeted effort 
//           5. a vector of 11 parameters estimated by MINUIT, in order: targeted catchability; non-targeted catchability; sigma; biomass at step 1 and 2; von mises mean and variance, recruitment in year 1..4)

int BiomassProjection(const long unsigned int max_timestep, const double &TotTargetedEffort, const std::vector<double> &NontargetedEffort, const std::vector<double> &par, std::vector<double> &Biomass, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability){

  // for(unsigned int my_count = 0; my_count < par.size(); ++my_count)
  //   std::cout << par[my_count] << "\t";
  // std::cout << "\n";

  // Check for wrong input
  //assert(TargetedEffort.size() == max_timestep);
  assert(FishingPattern.size() == NWPY);

  // Define necessary containers
  //std::vector< double > Biomass(max_timestep, 0.0), Rec(max_timestep, 0.0), survival(max_timestep, 0.0);
  std::vector< double > SSB(max_timestep, 0.0), Rec(max_timestep, 0.0), survival(max_timestep, 0.0);
  std::vector< double > FishingMortality(max_timestep, 0.0), PredictedCatch(max_timestep, 0.0);

  std::vector< double > TargetedEffort(max_timestep, 0.0);

  // Remember that the optimizer work on parameters of the same order of magnitude
  double Targeted_catchability = par[0] * CatchabilityScalingFactor;
  double Nontargeted_catchability = par[1] * CatchabilityScalingFactor;

  // Initialise the vector of biomass
  Biomass[0] = par[3] *  BiomassScalingFactor;
  Biomass[1] = par[4] *  BiomassScalingFactor;

  // von mises parameters
  double vm_mean = par[5]; // make sure it varies between -M_PI and +M_PI
  double vm_sigma = par[6]; // should be positive

  // Initialise the vector of Effort
  for(unsigned int i=0; i < max_timestep; i += 1){
    TargetedEffort.at(i) = FishingPattern.at(i%NWPY) * TotTargetedEffort * Availability.at(i%NWPY);
    std::cout << "Targeted effort is " << TargetedEffort.at(i) << "\n";
 }
  // Calculate survival given the vector of effort
  for(unsigned int i=0; i < max_timestep; i += 1) {
    FishingMortality[i] = Targeted_catchability * TargetedEffort[i] + Nontargeted_catchability * NontargetedEffort[i];
    survival[i] = exp(-(M + FishingMortality[i]));
  }


  //  for(unsigned int i=12; i < max_timestep; i += 52) {
  //    int year_number = i / 52;
  //    //Rec.at(i) = 1.0e9;
  //    //std::cout << " rec par" << 4 + year_number << " is " << par[4 + year_number] << "\n";
  //    Rec.at(i) = par[4 + year_number] * RecruitmentScalingFactor;
  //    //if( i == 12) Rec.at(i) = par.at(4) * RecruitmentScalingFactor;
  //    //else Rec.at(i) = 1.0e9;
  //  }

  // Calculate the proportion of recruitment in each week 
  std::vector<double> RecDist(52, 0.0);
  RecDist = vonMisesRecDist(vm_mean, vm_sigma);
  
  // Create the vector of recruitment for the first year
    for(unsigned int counter = 0; counter < NWPY; counter++){
      Rec[counter] = par[7 + counter / NWPY] * RecruitmentScalingFactor * RecDist[counter % NWPY];
    }
  // Initialise recruitment in the first year
  //for(unsigned int k = 0; k < NWPY; ++k)
  //  Rec.at(k) = par[6] * RecruitmentScalingFactor * RecDist.at(k);


  // Example how to print the vector of parameters
  //std::cout << "This is par[0] " << par[0] << "\n";
  //std::cout << "This is par[1] " << par[1] << "\n";
  //std::cout << "Number of parameters " << par.size() << "\n";

  // Recursive calculation of biomass
  for(unsigned int counter = 2; counter < max_timestep; counter++){

    // compute recruitment for the year ahead
    if(counter % 52 == 0){
      double TotalSSB = 0.0;

      // Calculate the spawning stock biomass during the past year
      for(unsigned int i = 0; i < NWPY; ++i){
	SSB.at(counter - i - 1) = PropMature.at(NWPY - i - 1) * Biomass.at(counter - i - 1);
	TotalSSB += SSB.at(counter - i - 1);
}
      // Use a Beverton and Holt function to calculate the associated recruitment
      double PredRec = SrPar[0] * TotalSSB / (SrPar[1] + TotalSSB);

      // Distribute this recruitment over the coming year
      for(unsigned int k = 0; k < NWPY; ++k)
	Rec.at(counter + k) = PredRec * RecDist.at(k);
    }
    // calculate biomass
    Biomass[counter] = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter]; 

    PredictedCatch[counter] = FishingMortality[counter] / (M + FishingMortality[counter]) * Biomass[counter] * ( 1 - survival[counter]);

  }

  for(unsigned int counter = 2; counter < max_timestep; counter++)
    std::cout << counter << "\t" << Rec[counter] << "\t" << Biomass[counter] << "\t" << SSB.at(counter) << "\t" << PredictedCatch[counter] << "\n";

  // Calculate the total catch in the last year
  double TotalCatch = 0.0;
  for(unsigned int k = 0; k < NWPY; ++k)
    TotalCatch += PredictedCatch[max_timestep - 1 - k];


  // Return biomass >= 0
  //if(Biomass[time_index] < 0) 
  return TotalCatch;
  //else
  //  return Biomass[time_index];

}
