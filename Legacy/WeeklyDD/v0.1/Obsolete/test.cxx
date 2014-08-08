// a test to see how to access to minuit parameter vector

#include <iostream>

//#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
//#include "Minuit2/MnMigrad.h"
//#include "Minuit2/MnMinos.h"
//#include "Minuit2/MnContours.h"
//#include "Minuit2/MnPlot.h"
//#include "Minuit2/MinosError.h"
//#include "Minuit2/ContoursError.h"

using namespace ROOT::Minuit2;

void access_parameter(std::vector<double> &par){

  std::cout << "first par is " << par[0] << "\n";
  std::cout << "second par is " << par[1] << "\n";
  
}
