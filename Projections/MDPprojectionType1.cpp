// CREATED  25 August 2017
// MODIFIED  1 Sept 2017

// PURPOSE perform projections of a delay difference model over 1e6 years varying fishing effort each year
//         to create data used later to calculate a transition matrix for a Markov Decision Process 

// BACKGROUND a DelayDifferenceModel class defines a delay-difference model object. This object contains all parameters estimated
//            from fitting a delay-difference model to logbook data.
//            the DelayDifferenceModelProjection is made of a DelayDifferenceModel object + other type required to perform a projection
//            of the delay-difference model forward in time

#include "DelayDifferenceModel.h"
#include "DelayDifferenceModelProjection.h"

#include "../UsefulFunctions.h"

int main(){

  //// Creating a Delay Difference Model Projection
  DelayDifferenceModelProjection MBTigerPrawnProjection("Data/ConfigFileForProjection1.csv");

  // Set number of years to project the stock
  MBTigerPrawnProjection.SetNbOfYear(1e6);
  
  // Generate random effort using the auto option which generate automatically
  // generate about 3% of draw in the bins spanning 1-30.000 with width 1.000, including 3% of zeros (algorithm written with Martin Peron)
  MBTigerPrawnProjection.DefineEffort("auto", 1, 30000);
  
  // Project the stock 
  MBTigerPrawnProjection.Project();
  
  // Printing the parameter of the delay diffence
  MBTigerPrawnProjection.Print("FileWithDateTag");

  return 0;
};