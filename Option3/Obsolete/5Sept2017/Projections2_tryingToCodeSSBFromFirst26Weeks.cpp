// CREATED   10 Sep 2013
// MODIFIED  25 Aug 2017

// NOTE using Beverton and Hold stock-recruitment relationship instead of Ricker since 14 Feb. 2017

// PURPOSE calculate biomass using a delay difference model with weekly timesteps

#include <math.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <random>
#include <fstream>
#include <sstream>

#include "../UsefulFunctions.h"

//std::vector<double> vonMisesRecDist(double a, double b);

// STATUS working

// ARGUMENTS 1. the total number of time steps considered in the analysis
//           2. the total quantity of targeted effort applied to the fishery over 1 year
//           3. a vector of non targeted effort (set to zero for the moment)
//           4. a vector of parameter estimated by fitting the delay difference model
//           5. a vector density of maturity (its sum = 1) based on the proportion of biomass mature in each week from biological sampling. 
//           6. 3 parameters of the linearized form of the Ricker stock-recruitment relationship
//           7. a vector giving the proportion of effort in each week (it sum = 1)
//           8. a vector giving the fraction ( 1 <= x <= 0) of stock available to fishing in each week

int Projections2(const long unsigned int max_timestep, const double &TotTargetedEffort, const std::vector<double> &NontargetedEffort, const std::vector<Parameter> &par, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability){

  // Global variables
  extern int NIPY;

  // Read parameters
  double rho = GetParameterValueAccordingToSymbol(par, "rho");
  double wk = GetParameterValueAccordingToSymbol(par, "wk");
  double wk_1 = GetParameterValueAccordingToSymbol(par, "wk_1");
  double CatchabilityScalingFactor = GetParameterValueAccordingToSymbol(par, "CatchabilityScalingFactor");
  double BiomassScalingFactor = GetParameterValueAccordingToSymbol(par, "BiomassScalingFactor");
  double RecruitmentScalingFactor = GetParameterValueAccordingToSymbol(par, "RecruitmentScalingFactor");
  //int NIPY = (int) GetParameterValueAccordingToSymbol(par, "NIPY");

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
  std::vector< double > SSB2(max_timestep, 0.0);
  std::vector< double > FishingMortality(max_timestep, 0.0), PredictedCatch(max_timestep, 0.0);

  std::vector< double > TargetedEffort(max_timestep, 0.0);

  // Remember that the optimizer work on parameters of the same order of magnitude
  double M = GetParameterValueAccordingToSymbol(par, "M");
  double Targeted_catchability = GetParameterValueAccordingToSymbol(par, "q1") * GetParameterValueAccordingToSymbol(par, "CatchabilityScalingFactor");
  double Nontargeted_catchability = GetParameterValueAccordingToSymbol(par, "q2") * GetParameterValueAccordingToSymbol(par, "CatchabilityScalingFactor");

  // Initialise the vector of biomass
  Biomass[0] = GetParameterValueAccordingToSymbol(par, "B1") * GetParameterValueAccordingToSymbol(par, "BiomassScalingFactor");
  Biomass[1] = GetParameterValueAccordingToSymbol(par, "B2") * GetParameterValueAccordingToSymbol(par, "BiomassScalingFactor");

  // von mises parameters
  double vm_mean = GetParameterValueAccordingToSymbol(par, "vm_mean");
  double vm_sigma = GetParameterValueAccordingToSymbol(par, "vm_sigma"); // should be positive

  // Construct the vector of Effort which is constant from year to year
  for(unsigned int i=0; i < max_timestep; i += 1){
    TargetedEffort[i] = FishingPattern[i%NIPY] * TotTargetedEffort;
 }

  // Calculate survival given the vector of fishing effort and fixed natural mortality
  for(unsigned int i=0; i < max_timestep; i += 1) {
    FishingMortality[i] = (Targeted_catchability * TargetedEffort[i] + Nontargeted_catchability * NontargetedEffort[i]) * Availability[i%NIPY];
    survival[i] = exp(-(M + FishingMortality[i]));
  }

  // Calculate the proportion of recruitment in each week (kept constant between years)
  std::vector<double> RecDist(NIPY, 0.0);
  RecDist = vonMisesRecDist(vm_mean, vm_sigma,NIPY);
  
  // Create the vector of recruitment for the first year (using estimate of magnitude of recruitment from the first year)
    for(unsigned int counter = 0; counter < NIPY; counter++){
      Rec[counter] = GetParameterValueAccordingToShortName(par, "Recruit year 1") * GetParameterValueAccordingToSymbol(par, "RecruitmentScalingFactor") * RecDist[counter % NIPY];
    }

  // Recursive calculation of biomass
  for(unsigned int counter = 2; counter < max_timestep; counter++){

//     // compute recruitment for the year ahead
//         if(counter % NIPY == 0){
// 	  double TotalSSB = 0.0;

//       // Calculate the spawning stock biomass during the past year
//       for(unsigned int i = 0; i < NIPY; ++i){
//       SSB.at(counter - i - 1) = PropMature.at(NIPY - i - 1) * Biomass.at(counter - i - 1);
//       TotalSSB += SSB.at(counter - i - 1);
// }
//       // Use the linearized Ricker function to simulate recruitment
//       double RandNumber = distribution(generator);

//       // Random recruitment
//       //double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB + RandNumber) * TotalSSB;
//       double PredRec = exp(SrPar[0] - log(SrPar[1] + TotalSSB) + RandNumber) * TotalSSB;

//       // Deterministic recruitment
//       //double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB) * TotalSSB;
      
//       //if(PredRec < 0.0) throw;

//       // Distribute this recruitment over the coming year
//       for(unsigned int k = 0; k < NIPY; ++k)
// 	Rec[counter + k] = PredRec * RecDist[k];
//     }

    // This was changed to calculate SSB from the first 26 weeks of each biological years
    
    // compute recruitment for the year ahead
    if((counter+26) % NIPY == 0){
    std::cout << "The counter is equal to " << counter << std::endl;
    //    if(counter % NIPY == 0){
    double TotalSSB2 = 0.0;

    // Calculate the spawning stock biomass during the past year
    //  for(unsigned int i = 0; i < NIPY; ++i){
    for(unsigned int i = 0; i < 26; ++i){
    //  SSB.at(counter - i - 1) = PropMature.at(NIPY - i - 1) * Biomass.at(counter - i - 1);
      std::cout << "Biomass at ( " << counter - i - 1 << ") is = " << Biomass.at(counter - i - 1) << std::endl;
     SSB2.at(counter - i - 1) = PropMature.at(i) * Biomass.at(counter - i - 1);
    //   TotalSSB += SSB.at(counter - i - 1);
     TotalSSB2 += SSB2.at(counter - i - 1);
}
    std::cout << "At counter=" << counter << ", SSB2=" << TotalSSB2 << std::endl;
      // Use the linearized Ricker function to simulate recruitment
      double RandNumber2 = distribution(generator);

      // Random recruitment
      //double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB + RandNumber) * TotalSSB;
      //double PredRec = exp(SrPar[0] - log(SrPar[1] + TotalSSB2) + RandNumber2) * TotalSSB2;
      double PredRec = exp(SrPar[0] - log(SrPar[1] + TotalSSB2)) * TotalSSB2;
      std::cout << "Predicted recruitment = " << PredRec << std::endl;
      
      // Deterministic recruitment
      //double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB) * TotalSSB;
      
      if(PredRec < 0.0) throw;

      // Distribute this recruitment over the present year
      for(unsigned int k = 0; k < NIPY; ++k){
      Rec[counter + k - 26] = PredRec * RecDist[k];
      std::cout << counter + k - 26 << " ";
      }
    }
	
    // calculate biomass
    Biomass[counter] = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter]; 

    if(Biomass[counter] < 0) Biomass[counter] = 0;

    std::cout << "Biomass at ( " << counter << ") is = " << Biomass.at(counter) << std::endl;
      
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
