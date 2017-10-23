// CREATED  29 April 2014
// MODIFIED 31 Aug 2017

// AUTHOR Marco.Kienzle@gmail.com

// PURPOSE some useful functions

#include "UsefulClasses.h"

// von Mises distribution function
std::vector<double> vonMisesRecDist(double a, double b, const int& NIPY);

// Function to access vector of parameters
double GetParameterValueAccordingToSymbol(const std::vector<Parameter> &ParVector, const std::string &Symbol);
double GetParameterLowerLimitAccordingToSymbol(const std::vector<Parameter> &ParVector, const std::string &Symbol);
double GetParameterUpperLimitAccordingToSymbol(const std::vector<Parameter> &ParVector, const std::string &Symbol);

double GetParameterValueAccordingToShortName(const std::vector<Parameter> &ParVector, const std::string &ShortName);
double GetParameterLowerLimitAccordingToShortName(const std::vector<Parameter> &ParVector, const std::string &ShortName);
double GetParameterUpperLimitAccordingToShortName(const std::vector<Parameter> &ParVector, const std::string &ShortName);

// Function to set value in vector of parameters
void SetParameterValueAccordingToShortName(std::vector<Parameter> &ParVector, const std::string &ShortName, double Value);
void SetParameterUncertaintyAccordingToShortName(std::vector<Parameter> &ParVector, const std::string &ShortName, double Uncertainty);

// Output info about the delay difference
void banner(std::string version_nb, std::string version_date);

// Read variable from files
void display_file(std::string& filename);
double fill_from_file(std::string& filename, const std::string VariableDescription);
std::vector<Parameter> ReadParameterDescription(const std::string& filename);
std::vector<Parameter> ReadOutputParameterDescription(const std::string& filename);
std::vector<ProjectionInputFiles> ReadProjectionInputFileDescription(const std::string& filename);

std::vector<double> ReadParameterFromSingleColumnFile(const std::string& filename);

// Function to access the path of file
std::string GetProjectionFilePathAccordingToShortName(const std::vector<ProjectionInputFiles> &InputFileVector, const std::string &ShortName);

// Write parameters to file
void WriteParameterToCSV(const std::string& filename, std::vector<Parameter> ParVec);

