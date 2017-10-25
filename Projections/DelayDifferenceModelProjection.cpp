// file DelayDifferenceModelProjection.cpp contains definitions of functions applicable to the DelayDifferenceModelProjection class

// CREATED  25 Aug 2017
// MODIFIED  1 Sep 2017
// AUTHOR Marco.Kienzle@gmail.com

#include "DelayDifferenceModel.h"
#include "DelayDifferenceModelProjection.h"
#include "../lib_facilities2.h"
#include "../UsefulFunctions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>

// Constructor
DelayDifferenceModelProjection::DelayDifferenceModelProjection(std::string filename)
{

  // Read information about input files descriptions and paths
  //std::cout << "I am going to read projection data from " << filename << std::endl;
  std::vector<ProjectionInputFiles> ProjInputFiles = ReadProjectionInputFileDescription(filename);

  // Get the parameters of the delay difference model from file
  //std::cout << "Parameters of the DD model to be read from " << GetProjectionFilePathAccordingToShortName(ProjInputFiles, "DDpar") << std::endl;  
  DDmodel.Parameterize(GetProjectionFilePathAccordingToShortName(ProjInputFiles, "DDpar"));

  // Get the stock recruitment relationship parameters from file
  //std::cout << "Warning: assuming Beverton and Holt stock recruitment relationship " << std::endl;
  SRR.Name = "Beverton and Holt";
  SRR.par = ReadParameterFromSingleColumnFile(GetProjectionFilePathAccordingToShortName(ProjInputFiles, "SRRpar"));
  assert(SRR.par.size() == 3);
  //std::cout << "The size of parameter vector for SRR is " << SRR.par.size() << std::endl;
  
  // Get the Weekly effort pattern
  WeeklyEffortPattern = ReadParameterFromSingleColumnFile(GetProjectionFilePathAccordingToShortName(ProjInputFiles, "EffortPat"));
  //std::cout << "The size of the weekly effort pattern " << WeeklyEffortPattern.size() << std::endl;
  assert(WeeklyEffortPattern.size() == DDmodel.NIPY);
  
  // Get the Weekly availability
  WeeklyAvailability = ReadParameterFromSingleColumnFile(GetProjectionFilePathAccordingToShortName(ProjInputFiles, "Avail"));
  //std::cout << "The size of the weekly availability " << WeeklyAvailability.size() << std::endl;
  assert(WeeklyAvailability.size() == DDmodel.NIPY);

  // Get the Weekly percentage of mature adults
  WeeklyPercMature = ReadParameterFromSingleColumnFile(GetProjectionFilePathAccordingToShortName(ProjInputFiles, "PercSpawners"));
  //std::cout << "The size of the weekly percentage of mature adults " << WeeklyPercMature.size() << std::endl;
  assert(WeeklyPercMature.size() == DDmodel.NIPY);

  // Calculate the weekly distribution of recruitment according to the von Mises distribution
  RecDist = vonMisesRecDist(DDmodel.vm_mean, DDmodel.vm_sd,DDmodel.NIPY);
  
  // Initialise the biomass vector using parameter B1 and B2
  Biomass.push_back(DDmodel.B1);
  Biomass.push_back(DDmodel.B2);
  
  // Create a vector of parameter using the file
  //std::vector<Parameter> ModelInputParameters = ReadOutputParameterDescription(filename);

}

// Set the period over which to project dynamics of the stock
void DelayDifferenceModelProjection::SetNbOfYear(const int Years)
{
  //std::cout << "I will set the number of years over which to project to " << Years << std::endl;
  NbOfYearOfProjection = Years;
}

// Set effort to be used in the projection
void DelayDifferenceModelProjection::DefineEffort(std::string Type, const int value)
{
  //std::cout << "I will set effort to " << Type << " using " << value << " boat days" << std::endl;
  if( Type.compare("constant") == 0 ){
    for(unsigned int i=0; i < NbOfYearOfProjection * DDmodel.NIPY; i++)
      Effort.push_back(value * WeeklyEffortPattern[ i % 52 ]);
  }
  
}

// Draw effort randomly between LowerBound and UpperBound
void DelayDifferenceModelProjection::DefineEffort(std::string Type, const int LowerBound, const int UpperBound)
{

  //std::cout << "I will set effort to " << Type << " randomly choosing between " << LowerBound << " and " << UpperBound << std::endl;

  // using a uniform distribution
  if( Type.compare("runif") == 0 ){

    // Uniform integer generator
    std::default_random_engine generator2(std::random_device{}());
    std::uniform_int_distribution<int> distribution2(LowerBound, UpperBound);

    // Store total effort in each year into a vector
    std::vector<double> YearlyEffort;
    for(unsigned int i=0; i < NbOfYearOfProjection; i++){
      YearlyEffort.push_back( (double) distribution2(generator2));
      //std::cout << "Random effort is " << YearlyEffort[i] << std::endl;
    }

    // Distribute yearly effort to each week using the weekly effort pattern
    for(unsigned int i=0; i < NbOfYearOfProjection * DDmodel.NIPY; i++)
      Effort.push_back(YearlyEffort[i/52] * WeeklyEffortPattern[ i % 52 ]);
  } // End of runif

    // automatic insert a certain proportiof zero effort automatically
  if( Type.compare("auto") == 0 ){

    // Uniform integer generator
    std::default_random_engine generator2(std::random_device{}());
    std::uniform_int_distribution<int> distribution2(LowerBound, UpperBound);

    // Store total effort in each year into a vector
    std::vector<double> YearlyEffort;
    for(unsigned int i=0; i < NbOfYearOfProjection; i++){

      // ensure that 1 out of 31 yearly effort is zero
      if( (rand() % 31) < 1 ){ 
	YearlyEffort.push_back( 0.0);}
      else{ 
	YearlyEffort.push_back( (double) distribution2(generator2));}
      //std::cout << "Random effort is " << YearlyEffort[i] << std::endl;
    }
  
  // Distribute yearly effort to each week using the weekly effort pattern
  for(unsigned int i=0; i < NbOfYearOfProjection * DDmodel.NIPY; i++)
    Effort.push_back(YearlyEffort[i/52] * WeeklyEffortPattern[ i % 52 ]);
  } // End of auto

}

// Draw effort randomly between LowerBound and UpperBound
void DelayDifferenceModelProjection::Project()
{
  //std::cout << "Performing the projection" << std::endl;

  // resize vectors to required size
  Biomass.resize(NbOfYearOfProjection * DDmodel.NIPY);

  SSB.resize(NbOfYearOfProjection * DDmodel.NIPY);
  YearlySSB.resize(NbOfYearOfProjection);

  Rec.resize(NbOfYearOfProjection * DDmodel.NIPY);
  YearlyRec.resize(NbOfYearOfProjection);

  survival.resize(NbOfYearOfProjection * DDmodel.NIPY);
  FishingMortality.resize(NbOfYearOfProjection * DDmodel.NIPY);
  PredictedCatch.resize(NbOfYearOfProjection * DDmodel.NIPY);

  // Define some variables
  double TotalSSB = 0.0;
 
  /////// Gaussian random generator for recruitment
  // Normal random generator to simulate recruitment deviations from the mean
  std::default_random_engine generator(std::random_device{}()); // random seed

  //std::default_random_engine generator;
  //generator.seed(5); // NOTE that the sequence of random number is the same between all simulations which effectively allows to compare just the effect of varying effort by opposition to random seed above
  
  // Random generator for a Normal distribution with mean 0 and standard deviation read from file
  std::normal_distribution<double> distribution(0.0, SRR.par[2]);

  // Calculate survival given the vector of fishing effort and fixed natural mortality
  //std::cout << "The size of the fishing mortality vector is " << FishingMortality.size() << std::endl;
  
  for(unsigned int i=0; i < FishingMortality.size(); i += 1) {
    FishingMortality[i] = DDmodel.Targeted_q * Effort[i] * WeeklyAvailability[i % DDmodel.NIPY];
    //std::cout << "Fishing mortality is " << FishingMortality[i] << std::endl;
    survival[i] = exp(-(DDmodel.M + FishingMortality[i]));
  }

  // Create the vector of recruitment for the first year (using estimate of magnitude of recruitment from the first year)
    for(unsigned int counter = 0; counter < DDmodel.NIPY; counter++){
      Rec[counter] = DDmodel.EstimatedYearlyRecruitment[0] * RecDist[counter % DDmodel.NIPY];
    }

      // Recursive calculation of biomass
    for(unsigned int counter = 2; counter < Biomass.size(); counter++){

      // When you are at the beginning of the year
      if( (counter == 2) || (counter % DDmodel.NIPY == 0) ){
      //std::cout << "counter is equal to " << counter << ", I am in the beginning of year loop" << std::endl;
      TotalSSB = 0.0;
      
      // Create a vector to hold reproducers 
      std::vector< double > Reproducer(DDmodel.NIPY + 2, 0.0);
      std::vector< double > TmpSSB(DDmodel.NIPY, 0.0);

      // The biomass of reproducer is initialised with biomass
      Reproducer[0] = Biomass[counter - 2];
      //std::cout << "Biomass[" << counter -2 << "] = " << Biomass[counter-2] << std::endl;
      Reproducer[1] = Biomass[counter - 1];
      //std::cout << "Biomass[" << counter -1 << "] = " << Biomass[counter-1] << std::endl;
      
      // Survivor and growth of reproducer is iterated over the next year
      for(int i = 2; i < DDmodel.NIPY + 2; i++){
        
        Reproducer[i] = survival[counter-1+i-2] * Reproducer[i - 1] + DDmodel.rho * survival[counter-1 + i - 2] * Reproducer[i -1] - DDmodel.rho * survival[counter-1 + i - 2] * survival[counter -2 + i - 2] * Reproducer[i - 2];
        //std::cout << "The biomass of reproducer at index i=" << i << " is " << Reproducer[i] << std::endl;
        TmpSSB[i-2] = Reproducer[i] * WeeklyPercMature.at(i-2);
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
      //double PredRec = exp(SRR.par[0] + SRR.par[1] * TotalSSB + RandNumber) * TotalSSB; // Ricker SRR
      double PredRec = exp(SRR.par[0] - log(SRR.par[1] + TotalSSB) + RandNumber) * TotalSSB; // Beverton-Holt SRR
      YearlyRec[counter / 52] = PredRec;
      
      // Deterministic recruitment
      //double PredRec = exp(SRR.par[0] + SRR.par[1] * TotalSSB) * TotalSSB; // Ricker 
      //double PredRec = exp(SRR.par[0] - log(SRR.par[1] + TotalSSB)) * TotalSSB; // Beverton-Holt
      //std::cout << "Recruitment is " << PredRec << std::endl;
     
      //if(PredRec < 0.0) throw;
       //if(PredRec < 0.0) throw;

      // Distribute this recruitment over the coming year
      for(unsigned int k = 0; k < DDmodel.NIPY; k++)
        Rec[counter + k] = PredRec * RecDist[k];
      } // End of the loop over when you are at the beginning of the year
      
     
    // calculate biomass
    Biomass[counter] = std::max(survival[counter-1] * Biomass[counter-1] + DDmodel.rho * survival[counter-1] * Biomass[counter-1] - DDmodel.rho * survival[counter-1] *survival[counter -2] * Biomass[counter-2] - DDmodel.rho * survival[counter-1] * DDmodel.wk_1 *  Rec[counter-1] + DDmodel.wk * Rec[counter],0.0); 

    PredictedCatch[counter] = FishingMortality[counter] / (DDmodel.M + FishingMortality[counter]) * Biomass[counter] * ( 1 - survival[counter]);

    //std::cout << "Biomass at counter " << counter << "is " << Biomass[counter] << " and prop mature is " << WeeklyPercMature.at(counter % DDmodel.NIPY) << std::endl;

    } // End of the loop recursively calculating Biomass
    
}

// Print fisheries quantities to "Screen", to a generic "File" name or to "FileWithDateTag"
void DelayDifferenceModelProjection::Print(std::string Option)
{
  //std::cout << "Print the projection's results" << std::endl;

  // Print to screen
  if ( Option.compare("Screen") == 0 ) { 

    // Print column names
    std::cout << "year,biomass,Catch,Effort,Recruitment,SSB,YearRec,YearSSB,B0,B1\n";
    
    // Initialize the yearly variable
    double YearBiomass = 0.0;
    double YearCatch = 0.0;
    double YearTargetEffort = 0.0;
    double YearRec = 0.0;
    double YearSSB = 0.0;

    // Go through each vectors
    for(unsigned int counter = 0; counter < Biomass.size(); counter++){

      // Sum over a year
      YearBiomass += Biomass[counter];
      YearCatch += PredictedCatch[counter];
      YearTargetEffort += Effort[counter];
      YearRec += Rec[counter];
      YearSSB += SSB[counter];

      // At the end of each year
      if( (counter+1) % 52 == 0){
      //std::cout << "Counter = " << counter << std::endl;
	std::cout << (counter+1)/52 << "," << YearBiomass << "," << YearCatch << "," << YearTargetEffort << "," << YearRec << "," << YearSSB << "," << YearlyRec[(counter+1)/52 - 1] << "," << YearlySSB[(counter+1)/52 - 1] << "," << Biomass[counter - 51] << "," << Biomass[counter - 50] << "\n";

      // reset
      YearBiomass = 0.0;
      YearCatch = 0.0;
      YearTargetEffort = 0.0;
      YearRec = 0.0;
      YearSSB = 0.0;
      }
    }
  } // End of option to print to screen
  
  
  // Print to file
  if ( Option.compare("FileWithDateTag") == 0 ) { 
    time_t currentTime = time(0);
    tm* currentDate = localtime(&currentTime);
    char filename2[256] = {0};

    strftime(filename2, 256, "Results/Simulation/Projections_forMDP_byYear-%Y-%m-%d_at_%H-%M-%S.csv", currentDate);
    std::ofstream EFQ2;
    EFQ2.open(filename2);

    // Header of file
    EFQ2 << "year,biomass,Catch,Effort,Recruitment,SSB,YearRec,YearSSB,B0,B1\n";
    
    // Initialize the yearly variable
    double YearBiomass = 0.0;
    double YearCatch = 0.0;
    double YearTargetEffort = 0.0;
    double YearRec = 0.0;
    double YearSSB = 0.0;

    // Go through each vectors
    for(unsigned int counter = 0; counter < Biomass.size(); counter++){

      // Sum over a year
      YearBiomass += Biomass[counter];
      YearCatch += PredictedCatch[counter];
      YearTargetEffort += Effort[counter];
      YearRec += Rec[counter];
      YearSSB += SSB[counter];

      // At the end of each year
      if( (counter+1) % 52 == 0){
      //std::cout << "Counter = " << counter << std::endl;
	 EFQ2 << (counter+1)/52 << "," << YearBiomass << "," << YearCatch << "," << YearTargetEffort << "," << YearRec << "," << YearSSB << "," << YearlyRec[(counter+1)/52 - 1] << "," << YearlySSB[(counter+1)/52 - 1] << "," << Biomass[counter - 51] << "," << Biomass[counter - 50] << "\n";
	
      // reset
      YearBiomass = 0.0;
      YearCatch = 0.0;
      YearTargetEffort = 0.0;
      YearRec = 0.0;
      YearSSB = 0.0;
      }
    }
  } // End of option to print to file with date tag

  // Print to file
  if ( Option.compare("File") == 0 ) { 
    char filename2[256] = "Results/Simulation/Projections_forMDP_byYear.csv";

    std::ofstream EFQ2;
    EFQ2.open(filename2);

    // Header of file
    EFQ2 << "year,biomass,Catch,Effort,Recruitment,SSB,YearRec,YearSSB,B0,B1\n";
    
    // Initialize the yearly variable
    double YearBiomass = 0.0;
    double YearCatch = 0.0;
    double YearTargetEffort = 0.0;
    double YearRec = 0.0;
    double YearSSB = 0.0;

    // Go through each vectors
    for(unsigned int counter = 0; counter < Biomass.size(); counter++){

      // Sum over a year
      YearBiomass += Biomass[counter];
      YearCatch += PredictedCatch[counter];
      YearTargetEffort += Effort[counter];
      YearRec += Rec[counter];
      YearSSB += SSB[counter];

      // At the end of each year
      if( (counter+1) % 52 == 0){
      //std::cout << "Counter = " << counter << std::endl;
	 EFQ2 << (counter+1)/52 << "," << YearBiomass << "," << YearCatch << "," << YearTargetEffort << "," << YearRec << "," << YearSSB << "," << YearlyRec[(counter+1)/52 - 1] << "," << YearlySSB[(counter+1)/52 - 1] << "," << Biomass[counter - 51] << "," << Biomass[counter - 50] << "\n";
	
      // reset
      YearBiomass = 0.0;
      YearCatch = 0.0;
      YearTargetEffort = 0.0;
      YearRec = 0.0;
      YearSSB = 0.0;
      }
    }
  } // End of option to print to file
  
} // End of Print function

