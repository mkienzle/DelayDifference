// CREATED  29 April 2014
// MODIFIED 30 Aug 2017

// AUTHOR Marco.Kienzle@gmail.com

// PURPOSE some useful classes

// header guards to avoid class redefinition error
#ifndef USEFUL_CLASSES_H
#define USEFUL_CLASSES_H

// A class to hold parameters used in the model
class Parameter
{
   public:
     std::string LongName = "NA";
     std::string ShortName = "NA";
     std::string Symbol = "NA";
     std::string Type = "NA"; // is it estimated or fixed
     double Value = 0.0;
     double Uncertainty = 0.0;
     double Boundaries [2] = { 0.0 }; // lower or upper boundaries used during the estimation process;
     std::string Unit = "NA";
};

// A class to hold path to projections' input file
class ProjectionInputFiles
{
   public:
     std::string LongName = "NA";
     std::string ShortName = "NA";
     std::string Path = "NA";
};

// A class to stock-recruitment function information
// Note: parameters and uncertainty were evaluate on the log-scale
class Stock_Recruitment_Relationship
{
 public:
  std::string Name = "NA";
  std::vector<double> par;
};

#endif
