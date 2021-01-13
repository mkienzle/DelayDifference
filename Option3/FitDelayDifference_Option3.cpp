// CREATED  10 Sep 2013
// MODIFIED 25 Oct 2020

// VERSION 0.3

// Copyright (c) 2013, 2014 Queensland Government, Department of Agriculture, Forestry, and Fisheries
// Programmed by Marco Kienzle
// This code is distributed under the GNU GPL v.3 license (see at the bottom for more info)

// PURPOSE fit a delay difference model (with weekly time-steps) to a series of catch and effort data
//         by maximum likelihood

// USAGE    DelayDifference_Option3 Data/SimData4.txt FixParameters.txt
// OPTIONAL DelayDifference_Option3 Data/SimData4.txt FixParameters.txt --ZeroBiomassAtWeek 42

// ARGUMENTS
// First argument: path to a file containing formatted data into 7 columns format containing: (1) timestep (2) label for the type of year used (3) numeric value of year (4) numeric value for week (5) catch (6) targeted and (7) un-targeted effort data from file
// Second argument: path to a file containing values of FixParameters
// Third (optional argument): a number of week in the fishing year at which the biomass will be set to zero (to represent senescence, total mortality after spawning)

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
double rho, wk, wk_1, M;
double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
int NIPY;
int EndOfSpawningWeekNumber = 0;

int main( int argc, char *argv[]) {

  // Print a starting message
  banner("0.1, option 3", "2020-10-25");

  //if ( argc != 2) {// argc should be 2 for correct execution
  //  cout << "usage: " << argv[0] << " <filename>\n";
  //  else {
    // We assume argv[1] is a filename to open
    std::ifstream ifs(argv[1]);

  // std::ifstream ifs( "Data/SimData4.txt" );
  //}
  // works well on real data
  //std::ifstream ifs( "Data/TigerWeeklyData1989-2010");


  // Grab the ZeroBiomassAtWeek number if passed as an argument
    std::string WeekNumber_strg;

    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--ZeroBiomassAtWeek") {
	  if ((i + 1) < argc) { // Make sure we aren't at the end of argv!
                WeekNumber_strg = argv[i+1]; // Increment 'i' so we don't get the argument as the next argv[i].
		EndOfSpawningWeekNumber = stoi(WeekNumber_strg);

            }             
        }
    }
    
    
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Read the values of fixed parameters from file in local directory into global variables
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::string arg2(argv[2]);
  //string FixParfilename = "FixParameters.txt";
  //display_file(FixParfilename);

  rho = fill_from_file(arg2, "Brody growth coefficient");
  wk = fill_from_file(arg2, "Estimated weight at recruitment");
  wk_1 = fill_from_file(arg2, "Parameter defining weight one timestep before recruitment");
  M = fill_from_file(arg2, "Natural mortality"); // Natural mortality is estimated in this version of the delay difference
  CatchabilityScalingFactor = fill_from_file(arg2, "Catchability scaling factor");
  BiomassScalingFactor = fill_from_file(arg2, "Biomass scaling factor");
  RecruitmentScalingFactor = fill_from_file(arg2, "Recruitment scaling factor");
  NIPY = (int) fill_from_file(arg2, "Number of intervals in a year");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Standard minimization using MIGRAD
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Declare the objective function (also called FCN function in documentation)  
  LogLikelihoodFunction fFCN(FishCatch, Targeted_Effort, Nontargeted_Effort);

  // create Minuit parameters with names and expected uncertainty
  MnUserParameters upar;
  //upar.Add("Natural mortality", 2.34 / (double) NIPY , 0.001);
  upar.Add("Natural mortality", M, 0.001);
  upar.Add("Targeted q", 1, 0.5);
  upar.Add("Nontargeted q", 2, 0.5);
  upar.Add("sigma", 10, 1);
  upar.Add("Biomass1", 1, 0.5);
  upar.Add("Biomass2", 1, 0.5);
  upar.Add("vm_mean", -0.6, 0.1);
  upar.Add("vm_sigma", 1, 1);

  std::string RecVarName="Recruit year";
  for(unsigned int i = 1; i <= counter / NIPY; ++i){
    upar.Add(RecVarName + " " + std::to_string(i), 0.2 * (i+1), .5);
  }
  // Assert parameters domain
  upar.SetLimits("Natural mortality", 0, 0.3);
  upar.SetLimits("Targeted q", 0, 1e2);
  upar.SetLimits("Nontargeted q", 0, 1e2);
  upar.SetLimits("sigma", 0, 1e2);
  upar.SetLimits("Biomass1", 0, 1e3);
  upar.SetLimits("Biomass2", 0, 1e3);
  upar.SetLimits("vm_mean", -M_PI, M_PI);
  upar.SetLimits("vm_sigma", 1e-3, 80);

  for(unsigned int i = 1; i <= counter/NIPY; ++i){
    upar.SetLimits(RecVarName + " " + std::to_string(i), 1e-4, 1e4);
  }

  cout << "The number of variable is " << upar.Params().size() << "\n";

  // create MIGRAD minimizer with MnStrategy 0 (strategy to calculate first and second derivative with fewer function calls -- less precise result)
  MnMigrad migrad(fFCN, upar, 0);

  // Fix a parameter
  migrad.Fix("Natural mortality");
  migrad.Fix("Targeted q");
  migrad.Fix("Nontargeted q");
  migrad.Fix("sigma");
  

  // Minimize
  FunctionMinimum min = migrad(50000);

  // output
  std::cout<<"minimum: "<< min << std::endl;

  migrad.Release("Natural mortality");
  migrad.Release("Targeted q");
  migrad.Release("Nontargeted q");
  migrad.Release("sigma");
  
  FunctionMinimum min2 = migrad(50000);
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
