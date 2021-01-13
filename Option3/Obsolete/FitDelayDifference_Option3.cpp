// CREATED  10 Sep 2013
// MODIFIED  6 Apr 2016

// VERSION 0.1

// Copyright (c) 2013-17 Queensland Government, Department of Agriculture and Fisheries
// Programmed by Marco Kienzle
// This code is distributed under the GNU GPL v.3 license (see at the bottom for more info)

// PURPOSE fit a delay difference model (with weekly time-steps) to a series of catch and effort data
//         by maximum likelihood

// USAGE DelayDifference_Option3 Data/SimData4.txt InputParameterDescription.csv

// ARGUMENT(1) the path to a file containing formatted data into 7 columns format containing: (1) timestep (2) label for the type of year used (3) numeric value of year (4) numeric value for week (5) catch (6) targeted and (7) un-targeted effort data from file
// ARGUMENT(2) this program reads input parameters from a CSV file passed as the second argument. This file contains the following information: LongName,ShortName,Symbol,Type,Value,LowerLimit,UpperLimit,Unit

// DATA an example of input data format is shown SimulatePopDynamic.R Test the code 

// COMPILE make

// CHANGES 2014-05-27 vonMisesRecDist was modified so that its input "a" is chnaged into ("a" %modulo% M_PI) so "a" can vary between -Inf and +Inf while the mean of von Mises is [-M_PI; M_PI]

#include "../lib_facilities.h"
#include "DelayDifference.h"
#include "LogLikelihoodFunction.h"
#include "../UsefulFunctions.h"

// Send error message if problem reading input datafile
class ReadingError { };

using namespace ROOT::Minuit2;

// Global variables
double rho, wk, wk_1; //, M;
double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
int NIPY;

int main( int argc, char *argv[]) {

  // Print a starting message
  banner("0.1, option 3", "2016-04-06");

  // We assume argv[1] is a filename to open
  std::ifstream ifs(argv[1]);

  std::vector< double > timestep, Year, Week, FishCatch, Targeted_Effort, Nontargeted_Effort;
  std::vector<string> YearType;
  double a, c, d, e, f, g;
  string b;
  long unsigned counter=0;

  try{   // Assess if data file is read entirely
  while( ifs >> a >> b >> c >> d >> e >> f >> g){
    timestep.push_back( a ); YearType.push_back (b); 
    Year.push_back(c); Week.push_back(d); 
    FishCatch.push_back(e); Targeted_Effort.push_back(f); Nontargeted_Effort.push_back(g);
    //printf(" Reading %.2f %.2f %.2f from file \n", a, b, c, d);
    counter++;
 }
  // When finishing reading, if EOF not reached throw an error
  if(!ifs.eof()){ throw ReadingError();}
  ifs.close();
  }catch(ReadingError){ std::cerr << "Error: could not read the entire file\n"; return 1;}


  std::cout << "In main, the size of FishCatch is " << FishCatch.size() << "\n";

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Read the values of fixed parameters from file in local directory into global variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read file into a vector of class Parameter
  std::string arg2(argv[2]);
  std::vector<Parameter> ModelInputParameters = ReadParameterDescription(arg2);

  // Assign fixed model's parameter set according to the information provided by the user in the csv file passed as 2nd argument

  rho = GetParameterValueAccordingToSymbol(ModelInputParameters, "rho");
  wk = GetParameterValueAccordingToSymbol(ModelInputParameters, "wk");
  wk_1 = GetParameterValueAccordingToSymbol(ModelInputParameters, "wk_1");
  CatchabilityScalingFactor = GetParameterValueAccordingToSymbol(ModelInputParameters, "CatchabilityScalingFactor");
  BiomassScalingFactor = GetParameterValueAccordingToSymbol(ModelInputParameters, "BiomassScalingFactor");
  RecruitmentScalingFactor = GetParameterValueAccordingToSymbol(ModelInputParameters, "RecruitmentScalingFactor");
  NIPY = (int) GetParameterValueAccordingToSymbol(ModelInputParameters, "NIPY");

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Standard minimization using MIGRAD
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Declare the objective function (also called FCN function in documentation)  
  LogLikelihoodFunction fFCN(FishCatch, Targeted_Effort, Nontargeted_Effort);

  // create Minuit parameters with names and expected uncertainty
  MnUserParameters upar;
  upar.Add("Natural mortality", GetParameterValueAccordingToShortName(ModelInputParameters, "Natural mortality"), 0.001);
  upar.Add("Targeted q", GetParameterValueAccordingToShortName(ModelInputParameters, "Targeted q"), 0.5);
  upar.Add("Nontargeted q", GetParameterValueAccordingToShortName(ModelInputParameters, "Nontargeted q"), 0.5);
  upar.Add("sigma", GetParameterValueAccordingToShortName(ModelInputParameters, "sigma"), 1);
  upar.Add("Biomass1", GetParameterValueAccordingToShortName(ModelInputParameters, "Biomass1"), 0.5);
  upar.Add("Biomass2", GetParameterValueAccordingToShortName(ModelInputParameters, "Biomass2"), 0.5);
  upar.Add("vm_mean", GetParameterValueAccordingToShortName(ModelInputParameters, "vm_mean"), 0.1);
  upar.Add("vm_sigma", GetParameterValueAccordingToShortName(ModelInputParameters, "vm_sigma"), 1);

  std::string RecVarName="Recruit year";
  for(unsigned int i = 1; i <= counter / NIPY; ++i){
    //upar.Add(RecVarName + " " + std::to_string(i), 0.2 * (i+1), .5);
    upar.Add(RecVarName + " " + std::to_string(i), GetParameterValueAccordingToShortName(ModelInputParameters, RecVarName + " " + std::to_string(i)), .5);
  }
  // Assert parameters domain
  upar.SetLimits("Natural mortality", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "Natural mortality"), 
GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "Natural mortality"));
  upar.SetLimits("Targeted q", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "Targeted q"),
GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "Targeted q"));
  upar.SetLimits("Nontargeted q", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "Nontargeted q"),
 GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "Nontargeted q"));
  upar.SetLimits("sigma", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "sigma"),
GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "sigma"));
  upar.SetLimits("Biomass1", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "Biomass1"), 
GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "Biomass1"));
  upar.SetLimits("Biomass2", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "Biomass2"), 
GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "Biomass2"));
  upar.SetLimits("vm_mean", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "vm_mean"),
GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "vm_mean"));
  upar.SetLimits("vm_sigma", GetParameterLowerLimitAccordingToShortName(ModelInputParameters, "vm_sigma"), 
GetParameterUpperLimitAccordingToShortName(ModelInputParameters, "vm_sigma"));

  for(unsigned int i = 1; i <= counter/NIPY; ++i){
    //upar.SetLimits(RecVarName + " " + std::to_string(i), 1e-2, 1e2);
        upar.SetLimits(RecVarName + " " + std::to_string(i), 
		   GetParameterLowerLimitAccordingToShortName(ModelInputParameters, RecVarName + " " + std::to_string(i)),
		   GetParameterUpperLimitAccordingToShortName(ModelInputParameters, RecVarName + " " + std::to_string(i)));	
  }

  cout << "The number of variable estimated by Minuit will be " << upar.Params().size() << "\n";

  // create MIGRAD minimizer with MnStrategy 0 (strategy to calculate first and second derivative with fewer function calls -- less precise result)
  MnMigrad migrad(fFCN, upar, 0);

  // Fix a parameter
  migrad.Fix("Natural mortality");
  migrad.Fix("Targeted q");
  migrad.Fix("Nontargeted q");
  migrad.Fix("sigma");

  // Minimize
  FunctionMinimum min = migrad();

  // output
  std::cout<<"minimum: "<< min << std::endl;

  migrad.Release("Natural mortality");
  migrad.Release("Targeted q");
  migrad.Release("Nontargeted q");
  migrad.Release("sigma");
  
  FunctionMinimum min2 = migrad();
  std::cout<<"minimum2: "<< min2 << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output results to file
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<Parameter> ModelOutputParameters=ModelInputParameters;

  ofstream FitFile;
  FitFile.open ("Results/FitOutcome.txt");
  FitFile << "# Minimum negative log likelihood\n" << min2 << "\n";
  FitFile.close();

  ofstream FitParamFile;
  FitParamFile.open ("Results/ParameterEstimates.txt");
  for(unsigned int i = 0; i < upar.Params().size(); ++i){

    SetParameterValueAccordingToShortName(ModelOutputParameters, upar.Name(i), min2.UserState().Value(i));
    SetParameterUncertaintyAccordingToShortName(ModelOutputParameters, upar.Name(i),min2.UserState().Error(i));

    FitParamFile << upar.Name(i) << "," << min2.UserState().Value(i) << "," << min2.UserState().Error(i) << "\n";
  }
  FitParamFile.close();

  // Write parameters to csv file
  WriteParameterToCSV("Results/DelayDifferenceModelParameters.csv", ModelOutputParameters);

  // Output estimates of fisheries quantities: fishing mortality, catch and biomass 

  ofstream EFQ;
  EFQ.open ("Results/EstimatedFisheriesQuantities.txt");
  EFQ << "timestep,YearType,Year,Week,EstimatedFishingMort,EstimatedBiomass,EstimatedCatches\n";    
  
  std::vector<double> EstimatedCatches(counter, 0.0), EstimatedBiomass(counter, 0.0), EstimatedFishingMort(counter, 0.0);
  std::vector<double> Residuals(counter, 0.0);

  // Get parameter estimates into a vector
  std::vector<double> EstimatedPar;
  for(unsigned int i = 0; i < upar.Params().size(); i++) EstimatedPar.push_back(min2.UserState().Value(i));

  // Calculate estimated catch
  WeeklyDD(Targeted_Effort, Nontargeted_Effort, EstimatedBiomass, EstimatedPar);

  for(unsigned int iii=0; iii < counter; iii++)
    {
      EstimatedFishingMort.at(iii) = min2.UserState().Value("Targeted q") * CatchabilityScalingFactor * Targeted_Effort[iii] + min2.UserState().Value("Nontargeted q") * CatchabilityScalingFactor * Nontargeted_Effort[iii];

      EstimatedCatches.at(iii) = EstimatedFishingMort[iii] / (min2.UserState().Value("Natural mortality") + EstimatedFishingMort[iii]) * EstimatedBiomass[iii] * \
	(1 - exp( -(min2.UserState().Value("Natural mortality") + EstimatedFishingMort[iii] ) ));

      Residuals.at(iii) = (sqrt(EstimatedCatches[iii]) - sqrt(FishCatch[iii])) / min2.UserState().Value("sigma");

      EFQ << timestep[iii] << "," << YearType[iii] << "," << Year[iii] << "," << Week[iii] << "," << EstimatedFishingMort[iii] << "," << EstimatedBiomass[iii] << "," << EstimatedCatches[iii] << "\n";  
    }
  EFQ.close();


  // Define the value for the program to return on exit consistent with Linux standard i.e. 0=success 1=failure
  // if final state of minimizer is invalid, return failure
  // if any parameter at one of its boundaries, return failure
  bool failed=false, AtUpperBoundary = false, AtLowerBoundary;
  // check each parameters if at upper/lower boundaries

  // for all parameters
  for(unsigned int i=0; i < upar.Params().size(); i++){

    // An exception for non-targeted catchability which parameter estimation is often switched-off by using zeroes for non-targeted effort resulting in huge error estimates associated with this parameter
    if(!(i == 2 && (min2.UserState().Value(i) - 2.0) < 1e-5 )){

      AtLowerBoundary = fabs(min2.UserState().Value(i) - min2.UserState().Parameter(i).LowerLimit()) < min2.UserState().Error(i); 
      AtUpperBoundary = fabs(min2.UserState().Value(i) - min2.UserState().Parameter(i).UpperLimit()) < min2.UserState().Error(i);
	if( AtLowerBoundary || AtUpperBoundary )
	  { 
	    failed=true;
	  }
    }

  }
  // check if minimizer status is valid
  if(! min2.IsValid()){failed=true;}

  // The convention of linux software is 0=success and 1=failure
  //std::cout << "What is the value returned by the program ? " << failed << std::endl;
  return failed;
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