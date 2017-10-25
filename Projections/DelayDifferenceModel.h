// file DelayDifferenceModel.h

// CREATED  25 August 2017
// MODIFIED 31 August 2017
// AUTHOR Marco.Kienzle@gmail.com

#include "../lib_facilities2.h"

class DelayDifferenceModel {

 public:

  // Initialise an instance of this object by reading parameters values from file
  void Parameterize(std::string filename);

  // Print object to screen (if Option is "Screen") otherwise to file (Option gives the path to the file)
  void Print(std::string Option);

  //// Temporal coverage of the model
  // Number of time steps per year
  unsigned int NIPY;

   //// Biological parameters
  
  // Brody growth parameter
  double rho;

  // Weight at recruitment
  double wk, wk_1;
 
  // Biomasses at timestep 1 and 2
  double B1, B2;

  // von Mises distribution of recruitment
  double vm_mean, vm_sd;

  // Estimated yearly recruitment
  std::vector<double> EstimatedYearlyRecruitment;
  
  //// Fishing fleet parameters
  
  // Targeted and non-targeted catchability
  double Targeted_q, NonTargeted_q;  

  // Natural mortality
  double M;

 private:


};

