// CREATED  19 November 2013
// MODIFIED 18 August 2017

// version 0.1

// NOTE using Beverton and Hold stock-recruitment relationship instead of Ricker since 14 Feb. 2017

// USAGE DD_Option3Projections_forMDP Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedBevertonAndHoltCoef Data/AverageEffortPattern.txt Data/AverageAvailability.txt 1500000 1000 2000 218665 205087

// RESULTS in files:
//                1. Results/Simulation/Projections_byYear_forMDP.txt
//                2. Results/Simulation/Projections_forMDP.txt

// COMPILE make DD_Option3Projections_forMDP

// PURPOSE project biomass long into the future using a delay difference model according to
//         1. a set of parameters for the delay difference model 
//         2. a set of proportion of biomass spawning in each interval per year
//         3. estimates of the stock-recruitment relationship
//         4. proportion of total effort (yearly) allocated to each interval during the year
//         5. pattern of prawn availability
//         6. Recruitment in the first year
//         7. lower boundary of the range of effort to draw a value using an uniform distribution
//         8. upper boundary of the range of effort
//         9. Biomass at time step 0
//        10. Biomass at time step 1
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

int Projections2_forMDP(const long unsigned int max_timestep, const double &LowerTargetedEffort, const double &UpperTargetEffort, const std::vector<double> &NontargetedEffort, const std::vector<Parameter> &par, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability, const double &B0, const double &B1, const double &RecruitFirstYear);

 
int main(int argc, char *argv[]){

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

 const int NbYear = 2;
 const double Recruitment1stYear = strtod(argv[6], NULL);
 //std::cout << "The 1st year recruitment is:" << Recruitment1stYear << std::endl;
 
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

//std::cout << "Read in " << PropMature.size() << " proportion of sexually mature biomass\n";

    // Read single column file containing the stock recruitment parameters
    std::ifstream SRFile(argv[3]);
    std::vector< double > SrPar;

    while(SRFile >> a){
      SrPar.push_back(a);}
    assert(SrPar.size() == 3);
    SRFile.close();

    //std::cout << "Read in " << SrPar.size() << " stock recruitment parameters\n";

    // Read single column file the fishing pattern
    std::ifstream FishingPatternFile(argv[4]);
    std::vector< double > FishingPattern;

    while(FishingPatternFile >> a){
      FishingPattern.push_back(a);}
    assert(FishingPattern.size() == NIPY);
    FishingPatternFile.close();

    //std::cout << "Read in " << FishingPattern.size() << " fishing pattern values\n";

    // Read single column file the availability
    std::ifstream AvailabilityFile(argv[5]);
    std::vector< double > Availability;

    while(AvailabilityFile >> a){
      Availability.push_back(a);}
    assert(Availability.size() == NIPY);
    AvailabilityFile.close();

    //std::cout << "Read in " << Availability.size() << " availability values\n";

    // Total targeted effort in a year is given as an input
    //double n = strtod(argv[6],NULL);
    //printf("Total effort per year is: %f\n",n);

    // Total targeted effort in a year is given as an input
    double Lowern = strtod(argv[7],NULL);
    double Uppern = strtod(argv[8],NULL);
    //std::cout << "Total effort per year range between: " << Lowern << " and " << Uppern << std::endl;

    // Read initial biomasses passed as arguments
    double B0 = strtod(argv[9],NULL);
    double B1 = strtod(argv[10],NULL);
    
    // Non-targeted effort is assumed to be null
    std::vector< double >   NontargetedEffort(NbYear * NIPY, 0.0);//, proj_biomass(NbYear * NIPY, 0.0);

 
    // Perform projections and write results to file
    Projections2_forMDP(NbYear*NIPY, Lowern, Uppern, NontargetedEffort, ModelInputParameters, PropMature, SrPar, FishingPattern, Availability, B0, B1, Recruitment1stYear);
    
  return 0;

}
