// CREATED   1 April 2016
// MODIFIED  5 April 2016
// COMPILE g++ -std=c++0x -g -o test -I/usr/include /usr/lib/prob.o ../UsefulFunctions.cpp test.cpp

// USAGE ./test InputParameterDescription.csv

#include "../lib_facilities.h"
#include "../UsefulFunctions.h"

// Global variables
double rho, wk, wk_1, M;
double CatchabilityScalingFactor, BiomassScalingFactor,RecruitmentScalingFactor;
int NIPY;

int main( int argc, char *argv[]) 
{

  // We assume argv[1] is a filename to read parameters from 
    std::string arg1(argv[1]);
    //cout << "Input file name is " << arg1 << "\n";
  //display_file(arg1);

  std::vector<Parameter> ModelInputParameters = ReadParameterDescription(arg1);

  
  std::cout << "The starting value of vm_mean is " << GetParameterValueAccordingToShortName(ModelInputParameters, "vm_mean") << " bound from " << GetParameterLowerLimitAccordingToSymbol(ModelInputParameters, "vm_mean") << " to " << GetParameterUpperLimitAccordingToSymbol(ModelInputParameters, "vm_mean") << "\n"; 


  std::vector<Parameter> Par2 = ReadOutputParameterDescription(argv[2]);

  //display_file(arg1);

  return 0;
}
