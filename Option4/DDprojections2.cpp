// CREATED  19 November 2013
// MODIFIED 29 July 2014

// version 0.1 option 4

// USAGE (see DD_Option4Projections.sh) DD_Option4Projections FixedWeeklyParameters.txt /tmp/ParEstimates.txt Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/AverageAvailability.txt 100

// COMPILE make DD_v0.3Projections

// PURPOSE project biomass long into the future using a delay difference model according to
//         1. a set of fixed parameters also used in the delay difference model
//         2. a set of delay difference parameters estimated during the fitting process
//         3. a set of proportion of biomass spawning in each interval per year
//         4. estimates of the stock-recruitment relationship
//         5. proportion of total effort (yearly) allocated to each interval during the year
//         6. pattern of prawn availability
//         7. total yearly effort 

#include <fstream>
#include <sstream>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
//#include "FixPar.h"
#include <cstdlib>
#include <cassert>

#include "../lib_facilities.h"
#include "../UsefulFunctions.h"

using namespace std;

// Global variables
double rho, wk, wk_1;
double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
int NIPY;

int Projections2(const long unsigned int max_timestep, const double &TargetedEffort, const std::vector<double> &NontargetedEffort, const std::vector<double> &par, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability);

int main(int argc, char *argv[]){

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Read the values of fixed parameters from file in local directory into global variables
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //string FixParfilename = "FixParameters.txt";
  std::string arg1(argv[1]);
  //display_file(arg1);

  rho = fill_from_file(arg1, "Brody growth coefficient");
  wk = fill_from_file(arg1, "Estimated weight at recruitment");
  wk_1 = fill_from_file(arg1, "Parameter defining weight one timestep before recruitment");
  CatchabilityScalingFactor = fill_from_file(arg1, "Catchability scaling factor");
  BiomassScalingFactor = fill_from_file(arg1, "Biomass scaling factor");
  RecruitmentScalingFactor = fill_from_file(arg1, "Recruitment scaling factor");
  NIPY = (int) fill_from_file(arg1, "Number of intervals in a year");

  const int NbYear = 150;

    // Read single column file containing parameters 
    std::ifstream ParFile(argv[2]);
    std::vector< double > Par;
    double a;

    while(ParFile >> a){
      Par.push_back(a);}
    ParFile.close();
    
    cout << "Read in " << Par.size() << " estimated parameters from catch and effort data\n";

    // Read single column file containing the proportion of biomass sexually mature in each week
    std::ifstream MaturityFile(argv[3]);
    std::vector< double > PropMature;

    while(MaturityFile >> a){
      PropMature.push_back(a);}
    
    assert(PropMature.size() == NIPY);
    MaturityFile.close();

    cout << "Read in " << PropMature.size() << " proportion of sexually mature biomass\n";

    // Read single column file containing the stock recruitment parameters
    std::ifstream SRFile(argv[4]);
    std::vector< double > SrPar;

    while(SRFile >> a){
      SrPar.push_back(a);}
    assert(SrPar.size() == 3);
    SRFile.close();

    cout << "Read in " << SrPar.size() << " stock recruitment parameters\n";

    // Read single column file the fishing pattern
    std::ifstream FishingPatternFile(argv[5]);
    std::vector< double > FishingPattern;

    while(FishingPatternFile >> a){
      FishingPattern.push_back(a);}
    assert(FishingPattern.size() == NIPY);
    FishingPatternFile.close();

    cout << "Read in " << FishingPattern.size() << " fishing pattern values\n";

    // Read single column file the availability
    std::ifstream AvailabilityFile(argv[6]);
    std::vector< double > Availability;

    while(AvailabilityFile >> a){
      Availability.push_back(a);}
    assert(Availability.size() == NIPY);
    AvailabilityFile.close();

    cout << "Read in " << Availability.size() << " availability values\n";

    // Total targeted effort in a year is given as an input
    double n = strtod(argv[7],NULL);
    printf("Total effort per year is: %f\n",n);

    // Non-targeted effort is assumed to be null
    std::vector< double >   NontargetedEffort(NbYear * NIPY, 0.0);//, proj_biomass(NbYear * NIPY, 0.0);

 
    // Perform projections and write results to file
    Projections2(NbYear*NIPY, n, NontargetedEffort, Par, PropMature, SrPar, FishingPattern, Availability);
    
  return 0;

}
