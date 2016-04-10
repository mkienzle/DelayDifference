// CREATED  19 Nov 2013
// MODIFIED  6 Apr 2016

// COMPILE make all or make DD_Option2Projections

// VERSION 0.1

// USAGE DD_Option2Projections Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt 100

// COMPILE make DD_Option2Projections

// PURPOSE project biomass long into the future using a weekly delay difference model according to
//         1. a set of fixed and estimated parameters from fitting the delay difference model to data
//         2. a set of proportion of biomass spawning in each interval per year
//         3. estimates of the stock-recruitment relationship
//         4. proportion of total effort (yearly) allocated to each interval during the year
//         5. pattern of prawn availability
//         6. total yearly effort (in boat-days) 

// REFERENCE Kienzle et al. Fisheries Research 155 (2014) p. 138-148
//           Hilborn and Walters (1992) Quantitative Fisheries Stock Assessment

// OBSOLETE COMPILE g++ -std=c++0x -g -o WeeklyDDprojections -I/usr/include /usr/lib/prob.o Projections.cxx vonMisesRecDist.cxx WeeklyDDprojections.cpp

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
double rho, wk, wk_1, M;
double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
int NIPY;

int Projections(const long unsigned int max_timestep, const double &TargetedEffort, const std::vector<double> &NontargetedEffort, const std::vector<Parameter> &ParameterVector, std::vector<double> &PropMature, std::vector<double> &SrPar, std::vector<double> &FishingPattern, std::vector<double> &Availability);

int main(int argc, char *argv[]){

  std::cout << "Beginning of projection program\n";
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Read the values of fixed parameters from file in local directory into global variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //string FixParfilename = "FixParameters.txt";
  std::string arg1(argv[1]);
  //display_file(arg1);

  // Read file into a vector of class Parameter

  std::vector<Parameter> DelayDifferenceParameters = ReadOutputParameterDescription(arg1);
  std::cout << "After reading the parameter values from file\n";

  // Assign fixed model's parameter set according to the information provided by the user in the csv file passed as 2nd argument

  rho = GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "rho");
  wk = GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "wk");
  wk_1 = GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "wk_1");
  M = GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "M");
  CatchabilityScalingFactor = GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "CatchabilityScalingFactor");
  BiomassScalingFactor = GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "BiomassScalingFactor");
  RecruitmentScalingFactor = GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "RecruitmentScalingFactor");
  NIPY = (int) GetParameterValueAccordingToSymbol(DelayDifferenceParameters, "NIPY");

  const int NbYear = 150;

    // Read single column file containing Delay Difference parameters estimated by fitting the model to data 
    //std::ifstream ParFile(argv[2]);
    //std::vector< double > Par;
    //double a;

    //while(ParFile >> a){
    //  Par.push_back(a);}
    //ParFile.close();
    
    //cout << "Read in " << Par.size() << " estimated parameters from catch and effort data\n";

    // Read single column file containing the proportion of biomass sexually mature in each week
    std::ifstream MaturityFile(argv[2]);
    std::vector< double > PropMature;
    double a;

    while(MaturityFile >> a){
      PropMature.push_back(a);}
    assert(PropMature.size() == NIPY);
    MaturityFile.close();

    cout << "Read in " << PropMature.size() << " proportion of sexually mature biomass\n";

    // Read single column file containing the stock recruitment parameters estimated by linear regression on transformed data see Hilborn and Walters (1992)
    std::ifstream SRFile(argv[3]);
    std::vector< double > SrPar;

    while(SRFile >> a){
      SrPar.push_back(a);}
    assert(SrPar.size() == 3);
    SRFile.close();

    cout << "Read in " << SrPar.size() << " stock recruitment parameters\n";

    // Read single column file the fishing pattern
    std::ifstream FishingPatternFile(argv[4]);
    std::vector< double > FishingPattern;

    while(FishingPatternFile >> a){
      FishingPattern.push_back(a);}
    assert(FishingPattern.size() == NIPY);
    FishingPatternFile.close();

    cout << "Read in " << FishingPattern.size() << " fishing pattern values\n";

    // Read single column file containing availability
    std::ifstream AvailabilityFile(argv[5]);
    std::vector< double > Availability;

    while(AvailabilityFile >> a){
      Availability.push_back(a);}
    assert(Availability.size() == NIPY);
    AvailabilityFile.close();

    cout << "Read in " << Availability.size() << " availability values\n";

    // Total targeted effort in a year is given as an input
    double n = strtod(argv[6],NULL);
    printf("Total effort per year is: %f\n",n);


    // Non-targeted effort is assumed to be null
    std::vector< double >   NontargetedEffort(NbYear * NIPY, 0.0);//, proj_biomass(NbYear * NIPY, 0.0);

     // Perform projections and write results to file
    Projections(NbYear*NIPY, n, NontargetedEffort, DelayDifferenceParameters, PropMature, SrPar, FishingPattern, Availability);
    
  return 0;

}
