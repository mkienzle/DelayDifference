// CREATED   10 Sept 2013
// MODIFIED   8 Nov 2013

// PURPOSE calculate biomass using a delay difference model with weekly timesteps

# include <math.h>
# include "FixPar_2parVersion.h"
# include <vector>
# include <iostream>
# include <cassert>

// ARGUMENTS 1. the time step to evaluate the biomass at
//           2. the total number of time steps considered in the analysis
//           3. a vector of effort 
//           4. a vector of 2 parameters estimated by MINUIT (catchability and sigma)

double WeeklyDD_2par(const unsigned int time_index, const long unsigned int max_timestep, const std::vector<double> &Effort, const std::vector<double> &par){

  // Check for wrong input
  assert(Effort.size() == max_timestep);
  assert(par.size() == 2);

  // Define necessary containers
  std::vector< double > Biomass(max_timestep, 0.0), Rec(max_timestep, 0.0), survival(max_timestep, 0.0);

  // Remember that the optimizer work on parameters of the same order of magnitude
  double catchability = par[0] * CatchabilityScalingFactor;

  // Calculate survival given the vector of effort
  for(unsigned int i=0; i < max_timestep; i += 1) {
    survival.at(i) = exp(-(M + catchability * Effort[i]));
  }

  // Fill recruitment vector according to timing of
  for(unsigned int i=12; i < max_timestep; i += 52) Rec.at(i) = 1.0e9;

  // Initialise the vector of biomass
  Biomass.at(0) = B1;
  Biomass.at(1) = B2;

  // Example how to print the vector of parameters
  //std::cout << "This is par[0] " << par[0] << "\n";
  //std::cout << "This is par[1] " << par[1] << "\n";
  //std::cout << "Number of parameters " << par.size() << "\n";

  // Recursive calculation of biomass
  for(unsigned int counter = 2; counter < max_timestep; counter++){

    Biomass.at(counter) = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter]; 
  }

  // Return biomass >= 0
  if(Biomass[time_index] < 0) 
    return 0;
  else
    return Biomass[time_index];

}
