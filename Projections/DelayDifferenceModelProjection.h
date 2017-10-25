// file DelayDifferenceModelProjection.h

// CREATED  25 August 2017
// MODIFIED  1 Sept. 2017
// AUTHOR Marco.Kienzle@gmail.com

#include "../UsefulClasses.h"
#include "../lib_facilities2.h"

// This DelayDifferenceModelProjection class inherit from a DelayDifferenceModel which contains parameters estimated from fitting a delay-difference model to logbook data. Variable required for projecting a delay-difference (such as the number of year to project forward, the stock-recruitment relationship to close the loop, ...) are combined with the delay difference model into this class.


class DelayDifferenceModelProjection {

 public:

  // Initialise an instance of this object by reading information from file
  DelayDifferenceModelProjection(std::string filename);

  // A Delay Difference Model
  DelayDifferenceModel DDmodel;

  // Set the period over which to project dynamics of the stock
  void SetNbOfYear(const int Years);

  //// Function to set effort to be used in the projection
  // Set a single constant effort, Type = "constant"
  void DefineEffort(std::string Type, const int value);
  // Set randomly chose value between LowerBound and UpperBound using various method : Type = "runif" or "auto" 
  void DefineEffort(std::string Type, const int LowerBound, const int UpperBound);

  // Calculate future biomass, catch, SSB, ... in the fishery according to the delay-difference model
  void Project();
  
  // Print object to screen (if Option is "Screen") otherwise to file (Option gives the path to the file)
  void Print(std::string Option);

 private:

  // Temporal coverage of the projections
  unsigned int NbOfYearOfProjection;

  //// Biological variables
  // Stock-recruitment relationship
  Stock_Recruitment_Relationship SRR;

  // Weekly availability
  std::vector<double> WeeklyAvailability;

  // Weekly percentage of mature adults
  std::vector<double> WeeklyPercMature;

  // Weekly distribution of recruitment according to von Mises
  std::vector<double> RecDist;
  
  //// Fishing variable
  // Vector effort
  std::vector<double> Effort;
   
  // Pattern of fishing effort (Weekly)
  std::vector<double> WeeklyEffortPattern;

  // Vector of biomass
  std::vector<double> Biomass;

  // Vector of Spawning Stock Biomass
  std::vector<double> SSB;
  std::vector<double> YearlySSB;
  
  // Vector of recruitment
  std::vector<double> Rec;
  std::vector<double> YearlyRec;
  
  // Vector of survival
  std::vector<double> survival;
  
  // Vector of fishing mortalities
  std::vector<double> FishingMortality;
  
  // Vector of catches
  std::vector<double> PredictedCatch;
  
};

