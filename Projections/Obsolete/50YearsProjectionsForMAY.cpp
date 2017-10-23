// CREATED  25 August 2017
// MODIFIED  5 Sept 2017

// PURPOSE perform projections of a delay difference model over 50 years varying fishing effort each year
//         to create data to calculate a transition matrix for a Markov Decision Process 

// USAGE 50YearsProjectionsForMAY 1500

// BACKGROUND a DelayDifferenceModel class defines a delay-difference model object. This object contains all parameters estimated
//            from fitting a delay-difference model to logbook data.
//            the DelayDifferenceModelProjection is made of a DelayDifferenceModel object + other type required to perform a projection
//            of the delay-difference model forward in time

#include "DelayDifferenceModel.h"
#include "DelayDifferenceModelProjection.h"

#include "../UsefulFunctions.h"

int main(int argc, char *argv[]){

  //// Creating a Delay Difference Model Projection
  DelayDifferenceModelProjection MBTigerPrawnProjection("Data/ConfigFileForProjection1.csv");

  // Set number of years to project the stock
  MBTigerPrawnProjection.SetNbOfYear(50);
  
  // Generate constant effort
  MBTigerPrawnProjection.DefineEffort("constant", stoi(argv[1]));
  
  // Project the stock 
  MBTigerPrawnProjection.Project();
  
  // Printing the parameter of the delay diffence
  MBTigerPrawnProjection.Print("File");

  return 0;
};
