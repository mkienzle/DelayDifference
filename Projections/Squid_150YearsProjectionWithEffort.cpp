// CREATED  25 Aug 2017
// MODIFIED 18 Nov 2020

// PURPOSE perform projections for NZ squid
// BASED ON new programs to perform projections using classes and methods

// BACKGROUND a DelayDifferenceModel class defines a delay-difference model object. This object contains all parameter estimated
//            from fitting a delay-difference model to logbook data.
//            the DelayDifferenceModelProjection is made of a DelayDifferenceModel object + other type required to perform a projection
//            of the delay-difference model forward in time

#include "DelayDifferenceModel.h"
#include "DelayDifferenceModelProjection.h"

#include "../UsefulFunctions.h"

int main(int argc, char *argv[]){

  // Creating a Delay Difference model object and printing to screen
  //DelayDifferenceModel MBTigerPrawn;
  //MBTigerPrawn.Initialize("Data/1990-2010Model10DelayDifferenceModelParameters.csv");
  //MBTigerPrawn.Print("Screen"); //MBTigerPrawn.Print("/tmp/file.txt");

  //// Creating a Delay Difference Model Projection
  DelayDifferenceModelProjection MBTigerPrawnProjection("Data/ConfigFileForProjection.csv");
  //std::cout << "I am here " << std::endl;
  
  //std::cout << "The proportion of effort in the first week is " << MBTigerPrawnProjection.WeeklyEffortPattern[0] << std::endl;
  
  MBTigerPrawnProjection.SetNbOfYear(150);
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Various way of generating a vector of effort

  // to calculate MAY, constant effort over long period of time was used
  MBTigerPrawnProjection.DefineEffort("constant", stoi(argv[1]));

  // J. Filar wanted to draw effort from an interval 
  //MBTigerPrawnProjection.DefineEffort("runif", 0, 0);

  // An algorithm to simulate a certain number of zero was written with Martin Peron
  //MBTigerPrawnProjection.DefineEffort("auto", 1, 30000);
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  MBTigerPrawnProjection.Project();
  
  // Printing the parameter of the delay diffence
  MBTigerPrawnProjection.Print("Screen");
  //MBTigerPrawnProjection.Print("File");
  //MBTigerPrawnProjection.Print("FileWithDateTag");

  // MBTigerPrawnProjection.DelayDifferenceParameters.Print("Screen");

  
  //std::vector<ProjectionInputFiles> ProjectionsDataFile = ReadProjectionInputFileDescription("Data/ConfigFileForProjection1.csv");

  //  for(int i = 0; i < ProjectionInputFiles.size(); i++)
  //  std::cout << "The number of projection input file are " << ProjectionsDataFile.size() << std::endl;

  return 0;
};
