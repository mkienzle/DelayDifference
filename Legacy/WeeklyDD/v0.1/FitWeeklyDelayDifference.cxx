// CREATED  10 Sep 2013
// MODIFIED 23 Apr 2014

// VERSION 0.1

// Copyright (c) 2013, 2014 Queensland Government, Department of Agriculture, Forestry, and Fisheries
// Programmed by Marco Kienzle
// This code is distributed under the GNU GPL v.3 license (see at the bottom for more info)

// PURPOSE fit a delay difference model (with weekly time-steps) to a series of catch and effort data
//         by maximum likelihood

// USAGE ./FitWeeklyDelayDifference Data/SimData4.txt

// ARGUMENT the path to a file containing formatted data into 7 columns format containing: (1) timestep (2) label for the type of year used (3) numeric value of year (4) numeric value for week (5) catch (6) targeted and (7) un-targeted effort data from file

// DATA input data were simulated using the script SimulatePopDynamic.R in ~/mystuff/Work/DEEDI/Moreton Bay Prawn Trawl fishery/Analysis/Scripts/DelayDifferenceModel/1989-2010/Weekly/Test the code 

// COMPILE make

// OBSOLETE g++ -std=c++0x -g -o FitWeeklyDelayDifference -I/usr/include /usr/lib/prob.o vonMisesRecDist.cxx FitWeeklyDelayDifference.cxx LogLikelihoodFunction.cxx WeeklyDelayDifference.cxx `root-config --glibs --cflags` -lMinuit2

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "Riostream.h"

# include "WeeklyDelayDifference.h"
# include "LogLikelihoodFunction.h"
# include "FixPar.h"

using namespace ROOT::Minuit2;

int main( int argc, char *argv[]) {

  // Starting message
  std::cout << "Down under delay difference model version 0.1 (2013-12-19)\n";
  std::cout << "Copyright (C) 2013 Queensland Government, Department of Agriculture, Forestry and Fisheries\n";
  std::cout << "This code is distributed under the GNU GPL v.3 license (http://www.gnu.org/licenses/)\n";

  // We assume argv[1] is a filename to read the data from
  std::ifstream ifs(argv[1]);

  std::vector< double > timestep, Year, Week, FishCatch, Targeted_Effort, Nontargeted_Effort;
  std::vector<string> YearType;
  double a, c, d, e, f, g;
  string b;
  long unsigned counter=0;

  while( ifs >> a >> b >> c >> d >> e >> f >> g){
    timestep.push_back( a ); YearType.push_back (b); 
    Year.push_back(c); Week.push_back(d); 
    FishCatch.push_back(e); Targeted_Effort.push_back(f); Nontargeted_Effort.push_back(g);
    //printf(" Reading %.2f %.2f %.2f from file \n", a, b, c, d);
    counter++;
 }
  ifs.close();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Standard minimization using MIGRAD
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Declare the objective function (also called FCN function in documentation)  
  LogLikelihoodFunction fFCN(FishCatch, Targeted_Effort, Nontargeted_Effort);

  // create Minuit parameters with names
  MnUserParameters upar;
  upar.Add("Targeted q", 4, 2);
  upar.Add("Nontargeted q", 2, 2);
  upar.Add("sigma", 10, 0.1);
  upar.Add("Biomass1", 1, 2);
  upar.Add("Biomass2", 1, 2);
  upar.Add("vm_mean", 0, 1);
  upar.Add("vm_sigma", 5, 2);

  std::string RecVarName="Recruit year";
  for(unsigned int i = 1; i <= counter / 52; ++i){
    upar.Add(RecVarName + " " + std::to_string(i), 0.2 * (i+1), .5);
  }
  // Assert parameters domain
  upar.SetLimits("Targeted q", 0, 1e2);
  upar.SetLimits("Nontargeted q", 0, 1e2);
  upar.SetLimits("sigma", 0, 1e2);
  upar.SetLimits("Biomass1", 0, 1e3);
  upar.SetLimits("Biomass2", 0, 1e3);
  upar.SetLimits("vm_mean", -M_PI, M_PI);
  upar.SetLimits("vm_sigma", 1e-3, 80);

  for(unsigned int i = 1; i <= counter/52; ++i){
    upar.SetLimits(RecVarName + " " + std::to_string(i), 1e-2, 1e2);
  }

  cout << "The number of variable is " << upar.Params().size() << "\n";

  // create MIGRAD minimizer with MnStrategy 0 (strategy to calculate first and second derivative with fewer function calls -- less precise result)
  MnMigrad migrad(fFCN, upar, 0);

  // Fix a parameter
  migrad.Fix("Targeted q");
  migrad.Fix("Nontargeted q");
  migrad.Fix("sigma");

  // Minimize
  FunctionMinimum min = migrad();

  // output
  std::cout<<"minimum: "<< min << std::endl;
  
  migrad.Release("Targeted q");
  migrad.Release("Nontargeted q");
  migrad.Release("sigma");
  
  FunctionMinimum min2 = migrad();
  std::cout<<"minimum2: "<< min2 << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output results to file
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ofstream FitFile;
  FitFile.open ("Results/FitOutcome.txt");
  FitFile << "# Minimum negative log likelihood\n" << min2 << "\n";
  FitFile.close();

  ofstream FitParamFile;
  FitParamFile.open ("Results/ParameterEstimates.txt");
  for(unsigned int i = 0; i < upar.Params().size(); ++i)
    FitParamFile << upar.Name(i) << "," << min2.UserState().Value(i) << "," << min2.UserState().Error(i) << "\n";
  FitParamFile.close();

  // Output estimates of fisheries quantities: fishing mortality, catch and biomass 

  ofstream EFQ;
  EFQ.open ("Results/EstimatedFisheriesQuantities.txt");
  EFQ << "timestep,YearType,Year,Week,EstimatedFishingMort,EstimatedBiomass,EstimatedCatches\n";    
  
  std::vector<double> EstimatedCatches(counter, 0.0), EstimatedBiomass(counter, 0.0), EstimatedFishingMort(counter, 0.0);
  std::vector<double> Residuals(counter, 0.0);

  // Get parameter estimates into a vector
  std::vector<double> EstimatedPar;
  for(unsigned int i = 0; i < upar.Params().size(); i++) EstimatedPar.push_back(min.UserState().Value(i));


  // Calculate biomass
  WeeklyDD(Targeted_Effort, Nontargeted_Effort, EstimatedBiomass, EstimatedPar);

  // Calculate estimated catch
  for(unsigned int iii=0; iii < counter; iii++)
    {
      EstimatedFishingMort.at(iii) = min.UserState().Value("Targeted q") * CatchabilityScalingFactor * Targeted_Effort[iii] + min.UserState().Value("Nontargeted q") * CatchabilityScalingFactor * Nontargeted_Effort[iii];
      //EstimatedBiomass.at(iii) = WeeklyDD(iii, counter, Targeted_Effort, Nontargeted_Effort, EstimatedPar);
      EstimatedCatches.at(iii) = EstimatedFishingMort[iii] / (M + EstimatedFishingMort[iii]) * EstimatedBiomass[iii] * \
	(1 - exp(- (M + EstimatedFishingMort[iii])));
      Residuals.at(iii) = (sqrt(EstimatedCatches[iii]) - sqrt(FishCatch[iii])) / min.UserState().Value("sigma");

      EFQ << timestep[iii] << "," << YearType[iii] << "," << Year[iii] << "," << Week[iii] << "," << EstimatedFishingMort[iii] << "," << EstimatedBiomass[iii] << "," << EstimatedCatches[iii] << "\n";  
    }
  EFQ.close();

  return 0;
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
