// CREATED   10 Sept 2013
// MODIFIED   3 Feb 2014

// PURPOSE calculate biomass using a delay difference model with monthly timesteps

# include <math.h>
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

void MonthlyDD(const std::vector<double> &TargetedEffort, const std::vector<double> &NontargetedEffort, std::vector<double> &Biomass, const std::vector<double> &par){

  // Global variables
  extern int NMPY;
  extern double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
  extern double rho, wk, wk_1, M;

  //std::cout << "In MonthlyDD function, rho " << rho << std::endl;
  //std::cout << "In MonthlyDD function, wk " << wk << std::endl;
  //std::cout << "In MonthlyDD function, wk_1 " << wk_1 << std::endl;
  //std::cout << "In MonthlyDD function, M " << M << std::endl;

  //  for(unsigned int i=0; i < par.size(); ++i){
  //  std::cout << "In MonthlyDD function, Parameter " << i+1 << " is " << par.at(i) << std::endl;
  //}

  assert(TargetedEffort.size() == Biomass.size());

  unsigned int max_timestep = TargetedEffort.size();

  // Check for wrong input
  //assert(TargetedEffort.size() == max_timestep);
  //assert(par.size() == 11);

  // Define necessary containers
  std::vector< double > Rec(max_timestep, 0.0), survival(max_timestep, 0.0);

  // Remember that the optimizer work on parameters of the same order of magnitude
  double Targeted_catchability = par[0] * CatchabilityScalingFactor;
  double Nontargeted_catchability = par[1] * CatchabilityScalingFactor;

  // Initialise the vector of biomass
  Biomass[0] = par[3] *  BiomassScalingFactor;
  Biomass[1] = par[4] *  BiomassScalingFactor;

  // von mises parameters
  double vm_mean = par[5]; // make sure it varies between -M_PI and +M_PI
  double vm_sigma = par[6]; // should be positive

  // Calculate survival given the vector of effort
  for(unsigned int i=0; i < max_timestep; i += 1) {
    survival[i] = exp(-(M + Targeted_catchability * TargetedEffort[i] + Nontargeted_catchability * NontargetedEffort[i]));
    //std::cout << "In MonthlyDelaydifference, survival at step " << i+1 << " is " << survival[i] << std::endl;
 }

  // Calculate the proportion of recruitment in each month
  std::vector<double> RecDist(NMPY, 0.0);
  RecDist[6] = 1.0;
  //RecDist = vonMisesRecDist(vm_mean, vm_sigma);

  
  // Create the vector of recruitment
    for(unsigned int counter = 0; counter < max_timestep; counter++){
      Rec[counter] = par[7 + counter / NMPY] * RecruitmentScalingFactor * RecDist[counter % NMPY];
      //std::cout << "Recruitment at time " << counter << " is " << Rec[counter] << std::endl;
    }

  // Example how to print the vector of parameters
  //std::cout << "This is par[0] " << par[0] << "\n";
  //std::cout << "This is par[1] " << par[1] << "\n";
  //std::cout << "Number of parameters " << par.size() << "\n";

  // Recursive calculation of biomass
  for(unsigned int counter = 2; counter < max_timestep; counter++){

    // calculate biomass
    Biomass[counter] = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter]; 

    // For some reason, it seems that omiting - rho * survival[counter-1] * wk_1 *  Rec[counter-1] fits the simulated data better
    //Biomass[counter] = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] + wk * Rec[counter]; 

    if(Biomass[counter] < 0){Biomass[counter] = 0;} 
 
    //std::cout << "Biomass in step " << counter << " is " << Biomass[counter] << std::endl;
  }

  // Return biomass >= 0
  //if(Biomass[time_index] < 0) 
  //  return 0;
  //else
  //  return Biomass[time_index];

}
