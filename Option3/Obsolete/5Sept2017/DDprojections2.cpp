// CREATED  19 November 2013
// MODIFIED 14 February 2017

// version 0.1

// NOTE using Beverton and Hold stock-recruitment relationship instead of Ricker since 14 Feb. 2017

// USAGE DD_Option3Projections Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedBevertonAndHoltCoef Data/AverageEffortPattern.txt Data/Availability.txt 100

// USAGE INTO A LOOP (see DD_Option3Projections.sh) DD_Option3Projections.sh Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedBevertonAndHoltCoef Data/AverageEffortPattern.txt Data/Availability.txt 


// COMPILE make DD_Option3Projections

// PURPOSE project biomass long into the future using a delay difference model according to
//         1. a set of parameters for the delay difference model 
//         2. a set of proportion of biomass spawning in each interval per year
//         3. estimates of the stock-recruitment relationship
//         4. proportion of total effort (yearly) allocated to each interval during the year
//         5. pattern of prawn availability
//         6. total yearly effort 

#include <fstream>
#include <sstream>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <cstdlib>
#include <cassert>

#include "../lib_facilities.h"
#include "../UsefulFunctions.h"

using namespace std;

  // Global variables
  int NIPY;

int Projections2(const long unsigned int max_timestep, const int &TargetedEffort, const std::vector<double> &NontargetedEffort, const std::vector<Parameter> &par, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability);

 
int main(int argc, char *argv[]){

  //  std::cout << std::endl;
  //  std::cout << "*** BEWARE this program does not implement SSB calculations from first 26 weeks ***" << std::endl;
  //  std::cout << "*** IS INCONSISTENT WITH HOW THE Beverton and Holt Stock-Recruitment was calculated ***" << std::endl;
  //  std::cout << std::endl;
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Read the values of fixed parameters from file in local directory into global variables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Read into a vector of class Parameter from CSV file containing values for the delay difference parameters
std::string arg1(argv[1]);
std::vector<Parameter> ModelInputParameters = ReadOutputParameterDescription(arg1);

double rho = GetParameterValueAccordingToSymbol(ModelInputParameters, "rho");
double wk = GetParameterValueAccordingToSymbol(ModelInputParameters, "wk");
double wk_1 = GetParameterValueAccordingToSymbol(ModelInputParameters, "wk_1");
double CatchabilityScalingFactor = GetParameterValueAccordingToSymbol(ModelInputParameters, "CatchabilityScalingFactor");
double BiomassScalingFactor = GetParameterValueAccordingToSymbol(ModelInputParameters, "BiomassScalingFactor");
double RecruitmentScalingFactor = GetParameterValueAccordingToSymbol(ModelInputParameters, "RecruitmentScalingFactor");
NIPY = (int) GetParameterValueAccordingToSymbol(ModelInputParameters, "NIPY");

  const int NbYear = 50;
  double a;

  // Read the fraction of biomass sexually mature at any time step 
  std::ifstream MaturityFile(argv[2]);
  std::vector< double > PropMature;

while(MaturityFile >> a)
  {
PropMature.push_back(a);
}

assert(PropMature.size() == NIPY);
MaturityFile.close();

    //cout << "Read in " << PropMature.size() << " proportion of sexually mature biomass\n";

    // Read single column file containing the stock recruitment parameters
    std::ifstream SRFile(argv[3]);
    std::vector< double > SrPar;

    while(SRFile >> a){
      SrPar.push_back(a);}
    assert(SrPar.size() == 3);
    SRFile.close();

    //cout << "Read in " << SrPar.size() << " stock recruitment parameters\n";

    // Read single column file the fishing pattern
    std::ifstream FishingPatternFile(argv[4]);
    std::vector< double > FishingPattern;

    while(FishingPatternFile >> a){
      FishingPattern.push_back(a);}
    assert(FishingPattern.size() == NIPY);
    FishingPatternFile.close();

    //cout << "Read in " << FishingPattern.size() << " fishing pattern values\n";

    // Read single column file the availability
    std::ifstream AvailabilityFile(argv[5]);
    std::vector< double > Availability;

    while(AvailabilityFile >> a){
      Availability.push_back(a);}
    assert(Availability.size() == NIPY);
    AvailabilityFile.close();

    //cout << "Read in " << Availability.size() << " availability values\n";

    // Total targeted effort in a year is given as an input
    int TotTargetedEffort = stoi(argv[6],NULL);
    //std::cout << "Total effort per year is  " << TotTargetedEffort << std::endl;

    // Non-targeted effort is assumed to be null
    std::vector< double >   NontargetedEffort(NbYear * NIPY, 0.0);//, proj_biomass(NbYear * NIPY, 0.0);

 
    // Perform projections and write results to file
    Projections2(NbYear*NIPY, TotTargetedEffort, NontargetedEffort, ModelInputParameters, PropMature, SrPar, FishingPattern, Availability);
    
  return 0;

}