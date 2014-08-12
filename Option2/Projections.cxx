// CREATED   10 Sep 2013
// MODIFIED  29 Apr 2014

// PURPOSE calculate biomass using a delay difference model with weekly timesteps

#include <math.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <random>
#include <fstream>
#include <sstream>

std::vector<double> vonMisesRecDist(double a, double b);

// STATUS working

// ARGUMENTS 1. the total number of time steps considered in the analysis
//           2. the total quantity of targeted effort applied to the fishery over 1 year
//           3. a vector of non targeted effort (set to zero for the moment)
//           4. a vector of parameter of the delay difference estimated by fitting the model to data
//           5. a vector density of maturity (its sum = 1) based on the proportion of biomass mature in each week from biological sampling. 
//           6. 3 parameters of the linearized form of the Ricker stock-recruitment relationship
//           7. a vector giving the proportion of effort in each week (it sum = 1)
//           8. a vector giving the fraction ( 1 <= x <= 0) of stock available to fishing in each week

int Projections(const long unsigned int max_timestep, const double &TotTargetedEffort, const std::vector<double> &NontargetedEffort, const std::vector<double> &par, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability){

  // Global variables
  extern int NIPY;
  extern double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
  extern double rho, wk, wk_1, M;

  // Normal random generator to simulate recruitment deviations from the mean
  std::default_random_engine generator;
  generator.seed(5); // NOTE that the sequence of random number is the same between all simulations which effectively allows to compare just the effect of varying effort by opposition of using generator.seed(time(NULL)); which would change the time series of recruitment at each level of effort

  std::normal_distribution<double> distribution(0.0, SrPar[2]);

  // Check for wrong input
  assert(NontargetedEffort.size() == max_timestep);
  assert(FishingPattern.size() == NIPY);
  
  // Define necessary containers
  std::vector< double > Biomass(max_timestep, 0.0);
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

  // Construct the vector of Effort which is constant from year to year
  for(unsigned int i=0; i < max_timestep; i += 1){
    TargetedEffort[i] = FishingPattern[i%NIPY] * TotTargetedEffort * Availability[i%NIPY];
 }

  // Calculate survival given the vector of fishing effort and fixed natural mortality
  for(unsigned int i=0; i < max_timestep; i += 1) {
    FishingMortality[i] = Targeted_catchability * TargetedEffort[i] + Nontargeted_catchability * NontargetedEffort[i];
    survival[i] = exp(-(M + FishingMortality[i]));
  }

  // Calculate the proportion of recruitment in each week (kept constant between years)
  std::vector<double> RecDist(NIPY, 0.0);
  RecDist = vonMisesRecDist(vm_mean, vm_sigma);
  
  // Create the vector of recruitment for the first year
    for(unsigned int counter = 0; counter < NIPY; counter++){
      Rec[counter] = par[7 + counter / NIPY] * RecruitmentScalingFactor * RecDist[counter % NIPY];
    }

  // Recursive calculation of biomass
  for(unsigned int counter = 2; counter < max_timestep; counter++){

    // compute recruitment for the year ahead
    if(counter % NIPY == 0){
      double TotalSSB = 0.0;

      // Calculate the spawning stock biomass during the past year
      for(unsigned int i = 0; i < NIPY; ++i){
	SSB.at(counter - i - 1) = PropMature.at(NIPY - i - 1) * Biomass.at(counter - i - 1);
	TotalSSB += SSB.at(counter - i - 1);
}
      // Use the linearized Ricker function to simulate recruitment
      double RandNumber = distribution(generator);

      // Random recruitment
      double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB + RandNumber) * TotalSSB;
      // Deterministic recruitment
      //double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB) * TotalSSB;
      
      //if(PredRec < 0.0) throw;

      // Distribute this recruitment over the coming year
      for(unsigned int k = 0; k < NIPY; ++k)
	Rec[counter + k] = PredRec * RecDist[k];
    }
    // calculate biomass
    Biomass[counter] = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter]; 

    PredictedCatch[counter] = FishingMortality[counter] / (M + FishingMortality[counter]) * Biomass[counter] * ( 1 - survival[counter]);

  }

  // Output into files
  std::ofstream EFQ;
  EFQ.open ("Results/Simulation/Projections.txt");
  EFQ << "timestep,biomass,catch,effort,recruitment,SSB\n";    

  // Recursive calculation of biomass
  for(unsigned int counter = 0; counter < max_timestep; counter++){
    EFQ << counter << "," << Biomass[counter] << "," << PredictedCatch[counter] << "," << TargetedEffort[counter] << "," << Rec[counter] << "," << SSB[counter] << "\n";
  }
  EFQ.close();

  // Calculate the total catch in the last year
  double TotalCatch = 0.0;
  for(unsigned int k = 0; k < NIPY; ++k)
    TotalCatch += PredictedCatch[max_timestep - 1 - k];

  return TotalCatch;

}
