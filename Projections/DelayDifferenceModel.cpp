// file DelayDifferenceModel.cpp contains definitions of functions applicable to the DelayDifferenceModel class

// CREATED  25 August 2017
// MODIFIED 25 August 2017
// AUTHOR Marco.Kienzle@gmail.com

#include "DelayDifferenceModel.h"
#include "../lib_facilities2.h"
#include "../UsefulFunctions.h"
//#include "../UsefulClasses.h"

// Initialise an instance of this object by reading parameters values from file
// and check
void DelayDifferenceModel::Parameterize(std::string filename)
{
  //std::cout << "I am going to read data from " << filename << std::endl;

  // Create a vector of parameter using the file
  std::vector<Parameter> ModelInputParameters = ReadOutputParameterDescription(filename);

  // Temporal coverage of the model
  NIPY = GetParameterValueAccordingToSymbol(ModelInputParameters, "NIPY");

  // Biological parameters
  B1 = GetParameterValueAccordingToSymbol(ModelInputParameters, "B1") * GetParameterValueAccordingToSymbol(ModelInputParameters, "BiomassScalingFactor");
  B2 = GetParameterValueAccordingToSymbol(ModelInputParameters, "B2") * GetParameterValueAccordingToSymbol(ModelInputParameters, "BiomassScalingFactor");

  vm_mean = GetParameterValueAccordingToSymbol(ModelInputParameters, "vm_mean");
  vm_sd = GetParameterValueAccordingToSymbol(ModelInputParameters, "vm_sigma");
  
  rho = GetParameterValueAccordingToSymbol(ModelInputParameters, "rho");

  wk = GetParameterValueAccordingToSymbol(ModelInputParameters, "wk");
  wk_1 = GetParameterValueAccordingToSymbol(ModelInputParameters, "wk_1");

  M = GetParameterValueAccordingToSymbol(ModelInputParameters, "M");

  EstimatedYearlyRecruitment.push_back(GetParameterValueAccordingToShortName(ModelInputParameters, "Recruit year 1") * GetParameterValueAccordingToSymbol(ModelInputParameters, "RecruitmentScalingFactor"));
    
  // Fisheries parameters  
  Targeted_q = GetParameterValueAccordingToSymbol(ModelInputParameters, "q1") * GetParameterValueAccordingToSymbol(ModelInputParameters, "CatchabilityScalingFactor");
  NonTargeted_q = GetParameterValueAccordingToSymbol(ModelInputParameters, "q2") * GetParameterValueAccordingToSymbol(ModelInputParameters, "CatchabilityScalingFactor");

  // Check for error
  if (B1 < 0) error("DelayDifferenceModel object constructor: B1 is negative ");
  if (B2 < 0) error("DelayDifferenceModel object constructor: B2 is negative ");

  if (vm_sd < 0) error("DelayDifferenceModel object constructor: von Mises standard deviation is negative ");
  if (rho < 0 || rho > 1) error("DelayDifferenceModel object constructor: Brody growth coef. outside [0;1] range ");
  if (wk < 0) error("DelayDifferenceModel object constructor: weight at recruitment age is negative ");
  if (wk_1 < 0) error("DelayDifferenceModel object constructor: weight at 1 timestep before recruitment age is negative ");

  if (M < 0) error("DelayDifferenceModel object constructor: natural mortality is negative ");
  if (Targeted_q < 0) error("DelayDifferenceModel object constructor: targeted catchability ");
  if (NonTargeted_q < 0) error("DelayDifferenceModel object constructor: non targeted catchability ");  
}

// Print object to screen
void DelayDifferenceModel::Print(std::string Option)
{
  if( Option.compare("Screen") == 0 ){

    std::cout << "# Number of timesteps in a year in the model" << std::endl << NIPY << std::endl;
    std::cout << "# Brody growth coefficient" << std::endl << rho << std::endl;
    std::cout << "# Biomass at timestep 1" << std::endl << B1 << std::endl;
    std::cout << "# Biomass at timestep 2" << std::endl << B2 << std::endl;
    std::cout << "# Natural mortality" << std::endl << M << std::endl;

    std::cout << "# Weight at recruitment age" << std::endl << wk << std::endl;
    std::cout << "# Weight 1 timestep before recruitment age" << std::endl << wk_1 << std::endl;

    std::cout << "# von Mises mean" << std::endl << vm_mean << std::endl;
    std::cout << "# von Mises standard deviation " << std::endl << vm_sd << std::endl;
    
    std::cout << "# Targeted catchability" << std::endl << Targeted_q << std::endl;
    std::cout << "# Non Targeted catchability" << std::endl << NonTargeted_q << std::endl;
  } else {
      // Open file for reading
    ofstream file(Option.c_str());    // open file for reading
    if (!file) error("can't open input file ",Option);

    file << "# Brody growth coefficient" << std::endl << rho << std::endl;
    file << "# Biomass at timestep 1" << std::endl << B1 << std::endl;
    file << "# Biomass at timestep 2" << std::endl << B2 << std::endl;
    file << "# Natural mortality" << std::endl << M << std::endl;

    file << "# Weight at recruitment" << std::endl << wk << std::endl;
    file << "# Weight 1 timestep before recruitment" << std::endl << wk_1 << std::endl;

    file << "# Targeted catchability" << std::endl << Targeted_q << std::endl;
    file << "# Non Targeted catchability" << std::endl << NonTargeted_q << std::endl;

    
    file.close();
  }
  
}


