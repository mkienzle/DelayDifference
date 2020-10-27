// CREATED   10 Sept 2013
// MODIFIED  25 Oct 2020

// Copyright (c) 2013, 2014 Queensland Government, Department of Agriculture, Forestry, and Fisheries
// Programmed by Marco Kienzle
// This code is distributed under the GNU GPL v.3 license (see at the bottom for more info)

// PURPOSE calculate biomass using a delay difference model with weekly timesteps

# include <math.h>
# include <vector>
# include <iostream>
# include <cassert>

std::vector<double> vonMisesRecDist(double a, double b);

// STATUS working, converging to the simulated parameters

// ARGUMENTS 1. a vector of targeted effort (CONSTANT)
//           2. a vector of non-targeted effort (CONSTANT)
//           3. a vector of biomass (TO BE FILLED)
//           4. a vector of 12 parameters estimated by MINUIT, in order: natural mortality, targeted catchability; non-targeted catchability; sigma; biomass at step 1 and 2; von mises mean and variance, recruitment in year 1..4)

void WeeklyDD(const std::vector<double> &TargetedEffort, const std::vector<double> &NontargetedEffort,  std::vector<double> &Biomass, const std::vector<double> &par){

  // Global variables
  extern int NIPY, EndOfSpawningWeekNumber;
  extern double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
  extern double rho, wk, wk_1;//, M;

  assert(TargetedEffort.size() == Biomass.size());

  unsigned int max_timestep = TargetedEffort.size();

  // Define necessary containers
  std::vector< double > Rec(max_timestep, 0.0), survival(max_timestep, 0.0);

  // Remember that the optimizer work on parameters of the same order of magnitude
  double M = par[0];
  double Targeted_catchability = par[1] * CatchabilityScalingFactor;
  double Nontargeted_catchability = par[2] * CatchabilityScalingFactor;

  // Initialise the vector of biomass
  Biomass[0] = par[4] *  BiomassScalingFactor;
  Biomass[1] = par[5] *  BiomassScalingFactor;

  // von mises parameters
  double vm_mean = par[6]; // make sure it varies between -M_PI and +M_PI
  double vm_sigma = par[7]; // should be positive

  // Calculate survival given the vector of effort
  for(unsigned int i=0; i < max_timestep; i += 1) {
    survival[i] = exp(-(M + Targeted_catchability * TargetedEffort[i] + Nontargeted_catchability * NontargetedEffort[i]));
  }

  // Calculate the proportion of recruitment in each week 
  std::vector<double> RecDist(NIPY, 0.0);
  RecDist = vonMisesRecDist(vm_mean, vm_sigma);
  
  // Create the vector of recruitment
    for(unsigned int counter = 0; counter < max_timestep; counter++){
      Rec[counter] = par[8 + counter / NIPY] * RecruitmentScalingFactor * RecDist[counter % NIPY];
    }

  // Example how to print the vector of parameters
  //std::cout << "This is par[0] " << par[0] << "\n";
  //std::cout << "This is par[1] " << par[1] << "\n";
  //std::cout << "Number of parameters " << par.size() << "\n";

  // Recursive calculation of biomass
  for(unsigned int counter = 2; counter < max_timestep; counter++){

    // set biomass to zero every year in week 43 (1st of August) to
    // represent all squid die after spawning (senescence), spawning season assumed to be June--July
    // 52 allowed for convergence
    if(EndOfSpawningWeekNumber != 0) {
      //      if( (counter % EndOfSpawningWeekNumber) == 0) {Biomass[counter-2] = 0; Biomass[counter-1] = 0;}
      // Remember that vector indices start at 0 in C++
      if( (counter % (EndOfSpawningWeekNumber - 1)) == 0) {

	////////////////////////////////////////////////////////////////
	// Calculate biomass assuming senescence (previous biomass died)
	////////////////////////////////////////////////////////////////
	
	// The biomass in the week after the end of spawning is only given by the recruitment
	Biomass[counter] = wk * Rec[counter]; } else {

	// The biomass 2 weeks after the end of spawning is given by recruitment + biomass in the previous week (i.e. the week after spawning)
	if(counter % (EndOfSpawningWeekNumber - 1) == 1) {
	Biomass[counter] = survival[counter - 1] * (wk - rho * wk_1) * Biomass[counter - 1] / wk + rho * survival[counter-1] * Biomass[counter-1] + wk * Rec[counter];}
      
      else {

	      
      Biomass[counter] = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter]; }
	    
      }}
     else {
    // calculate biomass
       Biomass[counter] = survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter]; }

    if(Biomass[counter] < 0){Biomass[counter] = 0;}
    
  } // End of the for loop
    
  }
  

///// Licensing agreement
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
