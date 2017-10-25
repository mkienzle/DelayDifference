// CREATED   10 Sep 2013
// MODIFIED  25 August 2017

// NOTE using Beverton and Holt stock-recruitment relationship instead of Ricker since 14 Feb. 2017

// PURPOSE calculate biomass using a delay difference model with weekly timesteps

#include <math.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <random>
#include <fstream>
#include <sstream>
#include <chrono>

#include "../UsefulFunctions.h"

//std::vector<double> vonMisesRecDist(double a, double b);

// STATUS working

// ARGUMENTS 1. the total number of time steps considered in the analysis

//           2. a vector of non targeted effort (set to zero for the moment)
//           3. a vector of parameter estimated by fitting the delay difference model
//           4. a vector density of maturity (its sum = 1) based on the proportion of biomass mature in each week from biological sampling. 
//           5. 3 parameters of the linearized form of the Ricker stock-recruitment relationship
//           6. a vector giving the proportion of effort in each week (it sum = 1)
//           7. a vector giving the fraction ( 1 <= x <= 0) of stock available to fishing in each week

//int Projections2(const long unsigned int max_timestep, const double &TotTargetedEffort, const std::vector<double> &NontargetedEffort, const std::vector<Parameter> &par, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability){
  
int Projections2_forMDPWithEffortLimits(const long unsigned int max_timestep, const double &LowerTargetedEffort, const double &UpperTargetEffort, const std::vector<double> &NontargetedEffort, const std::vector<Parameter> &par, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability){

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
  std::default_random_engine generator(std::random_device{}()); // random seed

  //std::default_random_engine generator;
  //generator.seed(5); // NOTE that the sequence of random number is the same between all simulations which effectively allows to compare just the effect of varying effort by opposition to random seed above
  
  // Random generator for a Normal distribution with mean 0 and standard deviation read from file
  std::normal_distribution<double> distribution(0.0, SrPar[2]);

  // Uniform integer generator
  //unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  
  //std::default_random_engine generator2(seed1);
  //std::default_random_engine generator2(std::random_device{}());
  //std::uniform_int_distribution<int> distribution2(1,50000);
  
  // Check for wrong input
  assert(NontargetedEffort.size() == max_timestep);
  assert(FishingPattern.size() == NIPY);
  
  // Define necessary containers
  std::vector< double > Biomass(max_timestep, 0.0);
  std::vector< double > SSB(max_timestep, 0.0), Rec(max_timestep, 0.0), survival(max_timestep, 0.0);
  std::vector< double > FishingMortality(max_timestep, 0.0), PredictedCatch(max_timestep, 0.0);

  std::vector< int > TotYearlyTargetedEffort(max_timestep / NIPY, 0);
  std::vector< double > TargetedEffort(max_timestep, 0.0);

  std::vector< double > YearlyRec(max_timestep, 0.0);
  std::vector< double > YearlySSB(max_timestep, 0.0);
  
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

  double TotalSSB = 0.0;
  
  // Draw TotYearlyTargetedEffort using a uniform distribution draw between lower and upper boundary
  // for as many years as there are in the simulation
  std::default_random_engine UnifDistGenerator(std::random_device{}());
  std::uniform_real_distribution<double> Unifdistribution(LowerTargetedEffort, UpperTargetEffort);

  for(int i = 0; i < max_timestep/NIPY; i +=1){
    TotYearlyTargetedEffort[i] = Unifdistribution(UnifDistGenerator);
    //std::cout << "Random effort value is: " << TotYearlyTargetedEffort[i] << std::endl;
  }

  // Construct the vector of Effort which is constant from year to year
  //for(unsigned int i=0; i < max_timestep; i += 1){
  //  TargetedEffort[i] = FishingPattern[i%NIPY] * TotTargetedEffort;
  //}

  // Create a randomly generated vector of total annual targeted effort
  //for(unsigned int i=0; i < max_timestep / NIPY; i +=1){
  //  double RandNumber2 = distribution2(generator2);

  // if( (rand() % 51) < 1 ){ // ensure that 1 out of 51 yearly effort is zero
  //    TotTargetedEffort[i] = 0; }
  //  else{ 
  //    TotTargetedEffort[i] = RandNumber2;}
  //  //std::cout << "Random effort is " << TotTargetedEffort[i] << std::endl;
  //    }

  // Construct the vector of Effort in each week given total effort varies each year
  for(unsigned int i=0; i < max_timestep; i += 1){
    TargetedEffort[i] = FishingPattern[i%NIPY] * TotYearlyTargetedEffort[i/NIPY];
    //std::cout << "Random weekly effort is " << TargetedEffort[i] << std::endl;
  }
  
  // Calculate survival given the vector of fishing effort and fixed natural mortality
  for(unsigned int i=0; i < max_timestep; i += 1) {
    FishingMortality[i] = (Targeted_catchability * TargetedEffort[i] + Nontargeted_catchability * NontargetedEffort[i]) * Availability[i%NIPY];
    //std::cout << "Fishing mortality is " << FishingMortality[i] << std::endl;
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

    // compute recruitment using the biomass of adults in the first 26 weeks of each biological years
    //if( counter % 52 > 0 && counter % 26 == 0){
    //  std::cout << "Is counter a multiple of 26 ? " << counter << std::endl;


      //for(unsigned int i = 0; i < 25; ++i) TotalSSB += SSB[counter - i];
    //}

    // When you are at the beginning of the year
    if( (counter == 2) || (counter % NIPY == 0) ){
      //std::cout << "counter is equal to " << counter << ", I am in the beginning of year loop" << std::endl;
      TotalSSB = 0.0;
      
      // Create a vector to hold reproducers 
      std::vector< double > Reproducer(NIPY + 2, 0.0);
      std::vector< double > TmpSSB(NIPY, 0.0);

      // The biomass of reproducer is initialised with biomass
      Reproducer[0] = Biomass[counter - 2];
      //std::cout << "Biomass[" << counter -2 << "] = " << Biomass[counter-2] << std::endl;
      Reproducer[1] = Biomass[counter - 1];
      //std::cout << "Biomass[" << counter -1 << "] = " << Biomass[counter-1] << std::endl;
      
      // Survivor and growth of reproducer is iterated over the next year
      for(int i = 2; i < NIPY + 2; i++){
	
	Reproducer[i] = survival[counter-1+i-2] * Reproducer[i - 1] + rho * survival[counter-1 + i - 2] * Reproducer[i -1] - rho * survival[counter-1 + i - 2] *survival[counter -2 + i - 2] * Reproducer[i - 2];
	//std::cout << "The biomass of reproducer at index i=" << i << " is " << Reproducer[i] << std::endl;
	TmpSSB[i-2] = Reproducer[i] * PropMature.at(i-2);
	SSB[counter + i - 2] = TmpSSB[i-2];
	//std::cout << "The Tmp SSB is " << TmpSSB[i-2] << std::endl;
	TotalSSB += TmpSSB[i-2];
      }
      
      //std::cout << "This year, the spawning stock biomass (SSB) is:" << TotalSSB << std::endl;
      // save TotalSSB
      //std::cout << "Stored at the index " << counter/52 <<  std::endl;
      YearlySSB[counter/52] = TotalSSB;
      
      // Stock-recruitment relationship
      double RandNumber = distribution(generator);
      //std::cout << "Random N(0,0.4) realisation is " << RandNumber << " and log(10) = " << log(10) << std::endl;

      // Random recruitment
      //double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB + RandNumber) * TotalSSB; // Ricker SRR
      double PredRec = exp(SrPar[0] - log(SrPar[1] + TotalSSB) + RandNumber) * TotalSSB; // Beverton-Holt SRR
      YearlyRec[counter / 52] = PredRec;
      
      // Deterministic recruitment
      //double PredRec = exp(SrPar[0] + SrPar[1] * TotalSSB) * TotalSSB; // Ricker 
      //double PredRec = exp(SrPar[0] - log(SrPar[1] + TotalSSB)) * TotalSSB; // Beverton-Holt
      //std::cout << "Recruitment is " << PredRec << std::endl;
     
      //if(PredRec < 0.0) throw;

      // Distribute this recruitment over the coming year
      for(unsigned int k = 0; k < NIPY; k++)
	Rec[counter + k] = PredRec * RecDist[k];
      }
    // calculate biomass
    Biomass[counter] = std::max(survival[counter-1] * Biomass[counter-1] + rho * survival[counter-1] * Biomass[counter-1] - rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - rho * survival[counter-1] * wk_1 *  Rec[counter-1] + wk * Rec[counter],0.0); 

    PredictedCatch[counter] = FishingMortality[counter] / (M + FishingMortality[counter]) * Biomass[counter] * ( 1 - survival[counter]);

    //std::cout << "Biomass at counter " << counter << "is " << Biomass[counter] << " and prop mature is " << PropMature.at(counter % NIPY) << std::endl;
  }

  // Output into files
  time_t currentTime = time(0);
  tm* currentDate = localtime(&currentTime);
  char filename[256] = {0};

  strftime(filename, 256, "Results/Simulation/Projections_forMDPWithEffortLimits.txt", currentDate);
  
  std::ofstream EFQ;
  EFQ.open (filename);
  EFQ << "timestep,biomass,catch,effort,recruitment,SSB\n";    

  std::ofstream EFQ2;
  char filename2[256] = {0};
  
  strftime(filename2, 256, "Results/Simulation/Projections_forMDP_byYearWithEffortLimits.txt", currentDate);

  EFQ2.open (filename2);
  EFQ2 << "year,biomass,Catch,Effort,Recruitment,SSB,YearRec,YearSSB,B0,B1\n";    


  double YearBiomass = 0.0;
  double YearCatch = 0.0;
  double YearTargetEffort = 0.0;
  double YearRec = 0.0;
  double YearSSB = 0.0;

  // Go through each vectors
  for(unsigned int counter = 0; counter < max_timestep; counter++){
    EFQ << counter << "," << Biomass[counter] << "," << PredictedCatch[counter] << "," << TargetedEffort[counter] << "," << Rec[counter] << "," << SSB[counter] << "\n";


    YearBiomass += Biomass[counter];
    YearCatch += PredictedCatch[counter];
    YearTargetEffort += TargetedEffort[counter];
    YearRec += Rec[counter];
    YearSSB += SSB[counter];

    // At the end of each year
    if( (counter+1) % 52 == 0){
      //std::cout << "Counter = " << counter << std::endl;
      
      EFQ2 << (counter+1)/52 << "," << YearBiomass << "," << YearCatch << "," << YearTargetEffort << "," << YearRec << "," << YearSSB << "," << YearlyRec[(counter+1)/52 - 1] << "," << YearlySSB[(counter+1)/52 - 1] << "," << Biomass[counter - 51] << "," << Biomass[counter - 50] << "\n";
      YearBiomass = 0.0;
      YearCatch = 0.0;
      YearTargetEffort = 0.0;
      YearRec = 0.0;
      YearSSB = 0.0;
    }
    
  }
    EFQ.close();
    EFQ2.close();
    
  // Calculate the total catch in the last year
  double TotalCatch = 0.0;
  for(unsigned int k = 0; k < NIPY; ++k)
    TotalCatch += PredictedCatch[max_timestep - 1 - k];

  return TotalCatch;

}