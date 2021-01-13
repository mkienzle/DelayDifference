// CREATED  29 Apr 2014
// MODIFIED 18 Nov 2020

// AUTHOR Marco.Kienzle@gmail.com

// PURPOSE some useful functions

// Read variable from files
#include <iostream>
#include <fstream>
#include <string>
#include "lib_facilities.h"
#include "UsefulClasses.h"
#include <boost/algorithm/string.hpp>


// Output info about the delay difference
void banner(std::string version_nb, std::string version_date)
{

  // Print a starting message
  std::cout << "\"Down under\" delay difference model version " << version_nb << " (" << version_date << ")\n";
  std::cout << "Copyright (C) 2013-17 Queensland Government, Department of Agriculture and Fisheries\n";
  std::cout << "This code is distributed under the GNU GPL v.3 license (http://www.gnu.org/licenses/)\n\n";

}

// Display file content
void display_file(std::string& filename)
{
  string line;
  unsigned int line_number=0;

    ifstream ist(filename.c_str());    // open file for reading
    if (!ist) error("can't open input file ",filename);
    cout << "hello\n";
    // ... use ist ...
    while(!ist.eof()){
      line_number++;
      getline(ist,line);
      if(line.length() > 0){
	cout << "Line " << line_number << ": " << line << endl;
      
      if(boost::find_first(line, "Natural mortality")){
	cout << "The line below contains the parameter you want" << endl;
      }
    }

    }
    cout << "bye\n";
    // the file is implicitly closed when we leave the function
}

//
// Get input parameters to the delay difference from descriptions in a CSV file
//

std::vector<Parameter> ReadParameterDescription(const std::string& filename)
{

  // file variable 
  std::string line;  // store lines in file
  std::vector<int> N; // Number of elements per line
  unsigned int line_number=0; // counter

  // Container for the parameters information in the file
  std::vector<Parameter> ParameterVector;

  // Open file for reading
  ifstream file(filename.c_str());    // open file for reading
  if (!file) error("can't open input file ",filename);

  // Read lines from file
  while(getline(file,line)){

    if(line.length() > 0){
      line_number++; // count the number of lines

      // Read lines content into a vector
      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
    
      while(std::getline(iss, token, ','))
	tokens.push_back(token);

      // How many element were read from each line 
      N.push_back(tokens.size());

      if(line_number == 1) {} // skip the headerdo nothing with the file header
    
      else{
    
	Parameter SingleParameter;

	// Assign entries in each line to a variable of type Parameter
	SingleParameter.LongName = tokens[0];
	SingleParameter.ShortName = tokens[1];
	SingleParameter.Symbol = tokens[2];
	SingleParameter.Type = tokens[3];
	SingleParameter.Value = atof(tokens[4].c_str());
	SingleParameter.Boundaries[0] = atof(tokens[5].c_str()); // NA converted in 0.0
	SingleParameter.Boundaries[1] = atof(tokens[6].c_str()); // NA converted in 0.0
	SingleParameter.Unit = tokens[7];

	ParameterVector.push_back(SingleParameter); // Add parameter to the vector of parameters

      } // End of else
    } // End of if(line.length() > 0
  } // End of while

  file.close();

  // Check that all lines contains the same number of elements separated by commas
  if (std::adjacent_find(N.begin(), N.end(), std::not_equal_to<int>() ) != N.end() )
    error("Inconsistent number of elements separated by commas in ",filename);
  
  // End of function execution message
  //std::cout << "Finished reading " << ParameterVector.size() << " parameters descriptions from file " << filename << "\n";
 
  return ParameterVector;
}

//
// Get parameters output by the delay difference into a CSV file
//

std::vector<Parameter> ReadOutputParameterDescription(const std::string& filename)
{

  // file variable 
  std::string line;  // store lines in file
  std::vector<int> N; // Number of elements per line
  unsigned int line_number=0; // counter

  // Container for the parameters information in the file
  std::vector<Parameter> ParameterVector;

  // Open file for reading
  ifstream file(filename.c_str());    // open file for reading
  if (!file) error("can't open input file ",filename);

  // Read lines from file
  while(getline(file,line)){

    if(line.length() > 0){
      line_number++; // count the number of lines

      // Read lines content into a vector
      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
    
      while(std::getline(iss, token, ','))
	tokens.push_back(token);

      // How many element were read from each line 
      N.push_back(tokens.size());

      if(line_number == 1) {} // skip the headerdo nothing with the file header
    
      else{
    
	Parameter SingleParameter;

	// Assign entries in each line to a variable of type Parameter
	SingleParameter.LongName = tokens[0];
	SingleParameter.ShortName = tokens[1];
	SingleParameter.Symbol = tokens[2];
	SingleParameter.Type = tokens[3];
	SingleParameter.Value = atof(tokens[4].c_str());
	SingleParameter.Uncertainty = atof(tokens[5].c_str()); // NA converted in 0.0
	SingleParameter.Unit = tokens[6];

	ParameterVector.push_back(SingleParameter); // Add parameter to the vector of parameters

      } // End of else
    } // End of if(line.length() > 0
  } // End of while

  file.close();

  // Check that all lines contains the same number of elements separated by commas
  if (std::adjacent_find(N.begin(), N.end(), std::not_equal_to<int>() ) != N.end() )
    error("Inconsistent number of elements separated by commas in ",filename);
  
  // End of function execution message
  //std::cout << "Finished reading " << ParameterVector.size() << " parameters descriptions from file " << filename << "\n";
 
  return ParameterVector;
}

//
// Write estimated and fixed parameters to the delay difference to a CSV file
//

void WriteParameterToCSV(const std::string& filename, std::vector<Parameter> ParameterVector)
{

  // file variable 
  //std::string line;  // store lines in file
  //std::vector<int> N; // Number of elements per line
  //unsigned int line_number=0; // counter

  // Container for the parameters information in the file
  //std::vector<Parameter> ParameterVector;

  // Open file for reading
  ofstream file(filename.c_str());    // open file for reading
  if (!file) error("can't open input file ",filename);

  // Write a header
  file << "LongName,ShortName,Symbol,Type,Value,Uncertainty,Unit\n";

  // Write the data
  for(unsigned int i=0; i < ParameterVector.size(); i++)
    {
      file << ParameterVector[i].LongName + "," 
	+ ParameterVector[i].ShortName + ","
	+ ParameterVector[i].Symbol + ","
	+ ParameterVector[i].Type + ","
	   << ParameterVector[i].Value << ","
	   << ParameterVector[i].Uncertainty << ","
	+ ParameterVector[i].Unit + "\n";
    }

  file.close();
}

// Get parameter value given its symbol
double GetParameterValueAccordingToSymbol(const std::vector<Parameter> &ParVector, const std::string &Symbol){

    for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].Symbol.compare(Symbol) == 0) 
	{
	  return ParVector[i].Value;
	}
      }
    // 	Return an error if the parameter looked for was not found
    error("Error in GetParameterValueAccordingToSymbol(): no value found for parameter " + Symbol + " in input parameter file.");
    return 1;
}

// Get parameter lower limit given its symbol
double GetParameterLowerLimitAccordingToSymbol(const std::vector<Parameter> &ParVector, const std::string &Symbol){

    for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].Symbol.compare(Symbol) == 0) 
	{
	  return ParVector[i].Boundaries[0];
	}
      }
    // 	Return an error if the parameter looked for was not found
    error("Error in GetParameterLowerLimitAccordingToSymbol(): no value found for parameter " + Symbol + " in input parameter file.");
    return 1;
}

// Get parameter upper limit given its symbol
double GetParameterUpperLimitAccordingToSymbol(const std::vector<Parameter> &ParVector, const std::string &Symbol){

    for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].Symbol.compare(Symbol) == 0) 
	{
	  return ParVector[i].Boundaries[1];
	}
      }
    // 	Return an error if the parameter looked for was not found
    error("Error in GetParameterUpperLimitAccordingToSymbol(): no value found for parameter " + Symbol + " in input parameter file.");
    return 1;
}

// Get parameter value given its short name
double GetParameterValueAccordingToShortName(const std::vector<Parameter> &ParVector, const std::string &ShortName){

    for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].ShortName.compare(ShortName) == 0) 
	{
	  return ParVector[i].Value;
	}
      }
    // 	Return an error if the parameter looked for was not found
    error("Error in GetParameterValueAccordingToShortName(): no value found for parameter " + ShortName + " in input parameter file.");
    return 1;
}

// Get parameter lower limit given its short name
double GetParameterLowerLimitAccordingToShortName(const std::vector<Parameter> &ParVector, const std::string &ShortName){

    for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].ShortName.compare(ShortName) == 0) 
	{
	  return ParVector[i].Boundaries[0];
	}
      }
    // 	Return an error if the parameter looked for was not found
    error("Error in GetParameterLowerLimitAccordingToShortName(): no value found for parameter " + ShortName + " in input parameter file.");
    return 1;
}

// Get parameter upper limit given its short name
double GetParameterUpperLimitAccordingToShortName(const std::vector<Parameter> &ParVector, const std::string &ShortName){

    for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].ShortName.compare(ShortName) == 0) 
	{
	  return ParVector[i].Boundaries[1];
	}
      }
    // 	Return an error if the parameter looked for was not found
    error("Error in GetParameterUpperLimitAccordingToShortName(): no value found for parameter " + ShortName + " in input parameter file.");
    return 1;
}

//
// Set parameter value given its symbol
//

void SetParameterValueAccordingToShortName(std::vector<Parameter> &ParVector, const std::string &ShortName, double Value){

    for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].ShortName.compare(ShortName) == 0) 
	{
	  ParVector[i].Value = Value;
	}
      }
}

// Set parameter uncertainty given its symbol
void SetParameterUncertaintyAccordingToShortName(std::vector<Parameter> &ParVector, const std::string &ShortName, double Uncertainty){

  for(unsigned int i = 0; i < ParVector.size(); ++i)
      {
	if( ParVector[i].ShortName.compare(ShortName) == 0) 
	{
	  ParVector[i].Uncertainty = Uncertainty;
	}
      }
}

// 
// Get the value of a variable from filename given a VariableDescription
//
double fill_from_file(std::string& filename, const std::string VariableDescription)
{
  string line;
  unsigned int line_number=0;
  double ParameterFromFile=0;

  bool FoundVariable = 0;

    ifstream ist(filename.c_str());    // open file for reading
    if (!ist) error("can't open input file ",filename);

    // ... use ist ...
    while(!ist.eof()){
      line_number++;
      getline(ist,line);
      
      // If the line in filename match the description, get the value in the next line
      if(boost::find_first(line, VariableDescription)){
	//cout << "The line below contains the parameter you want" << endl; 
	ist >> ParameterFromFile;
	ist.close();
	FoundVariable = 1;
      }

    }

    if(FoundVariable == 0) cout << "*** Warning *** : could not find a value for " << VariableDescription << " in file " << filename << endl;
    // the file is implicitly closed when we leave the function
    return ParameterFromFile;
}

// PURPOSE compute the probability from a von mises distribution with parameter a and b for every week of the year

// DEPENDENCY on the C++ library PROB by John Burkardt
//            available at /home/mkienzle/mystuff/Software/Prob_library/prob.cpp


// CHANGES 2014-05-27 input a is modified into a modulo M_PI so a can vary between -Inf and +Inf while the mean of von Mises is [-M_PI; M_PI]

// INPUT a is allowed to vary between -Inf and Inf
//       b must be positive

std::vector<double> vonMisesRecDist(double a, double b, const int& NIPY){

// Global variables
//extern int NIPY;


  //  cout << "\n";
  //  cout << "  PDF parameter A =      " << a << "\n";
  //  cout << "  PDF parameter B =      " << b << "\n";

  //  cout << " a=" << a << endl;
  //  cout << " modulo=" << fmod(a, M_PI) << endl;

  // The input for the mean (a) is tranformed to make sure it vary within [-pi; pi]
  a = fmod(a, M_PI);

  if ( !von_mises_check ( a, b ) )
  {
    std::cout << "\n";
    std::cout << "vonMisesRecDist - Fatal error!\n";
    std::cout << "  The parameters are not legal.\n";
    std::cout << " a=" << a << endl;
    throw;
  }

  // Assume a year is made of NIPY interval per years
  
  // divide the range [-pi;pi] into NIPY  intervals
  vector<double> boundaries(NIPY+1, 0.0);
  for(unsigned int i = 0; i < NIPY + 1; ++i){
    boundaries[i] = -M_PI + i * 2 * M_PI / NIPY;
      }

      // von_mises_pdf(X,A,B) requires
      ////    Input, double X, the argument of the CDF.
      //    A - PI <= X <= A + PI.
      //
      //    Input, double A, B, the parameters of the PDF.
      //    -PI <= A <= PI,
      //    0.0 < B.

  ////// Calculate von Mises probability associated with each interval

  // create a vector that will contain the proportion of recruit in each weeks
  vector<double> intervals(NIPY, 0.0);

  double total = 0.0;

    for(unsigned int i = 0; i < NIPY; ++i){

      ///// Calculate probability in each interval using the cumulative distribution function

      if(boundaries[i+1] <= (a - M_PI) ){
	intervals[i] = von_mises_cdf(boundaries[i+1] + 2 * M_PI, a, b) - von_mises_cdf(boundaries[i] + 2 * M_PI, a, b);
      }

      if(boundaries[i] <= (a - M_PI) && (a - M_PI) < boundaries[i+1]){
	intervals[i] = von_mises_cdf(a + M_PI, a, b) - von_mises_cdf(boundaries[i] + 2 * M_PI, a, b)
	  + von_mises_cdf(boundaries[i+1], a, b) - von_mises_cdf(a - M_PI, a, b);
}

      if(boundaries[i] > (a - M_PI) && (a + M_PI) > boundaries[i+1]){
	intervals[i] = von_mises_cdf(boundaries[i+1], a, b) - von_mises_cdf(boundaries[i], a, b);
}

      if(boundaries[i] < (a + M_PI) && (a + M_PI) <= boundaries[i+1]){
	intervals[i] = von_mises_cdf(a + M_PI, a, b) - von_mises_cdf(boundaries[i], a, b)
	  + von_mises_cdf( boundaries[i+1] - 2 * M_PI, a,b ) - von_mises_cdf( a - M_PI, a, b);
}

      if(boundaries[i] >= (a + M_PI)){
	intervals[i] = von_mises_cdf(boundaries[i+1] - 2 * M_PI, a, b) - von_mises_cdf(boundaries[i] - 2 * M_PI, a, b);
}

      // sometimes small negative probability appears 
      if(intervals[i] < 0.0){
	// if they are small, replace them by zero
	if(abs(intervals[i]) < 1e-5){ intervals[i] = 0.0;}
        // otherwise throw an error
	else{throw;}
      }
      
     total += intervals[i];
      	}
 
    //std::cout << "von mises density equal " << total << "\n";
    // throw exception if cumulative density 
    if(abs(total - 1) > 1e-2){
      std::cout << "von mises density equal " << total << "\n";
      throw;
    }
    //// std::cout << "And the sum of prob is " << total << "\n";
    return intervals;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read data from a specific column in a file. Columns of data in this file are separated by space
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//std::vector<double> ReadData(std::string &filename, std::int ColNb){}

//
// Get information about projections's input files
//

std::vector<ProjectionInputFiles> ReadProjectionInputFileDescription(const std::string& filename)
{

  // file variable 
  std::string line;  // store lines in file
  std::vector<int> N; // Number of elements per line
  unsigned int line_number=0; // counter

  // Container for the parameters information in the file
  std::vector<ProjectionInputFiles> ProjectionInputFilesVector;

  // Open file for reading
  ifstream file(filename.c_str());    // open file for reading
  if (!file) error("can't open input file ",filename);

  // Read lines from file
  while(getline(file,line)){

    if(line.length() > 0){
      line_number++; // count the number of lines

      // Read lines content into a vector
      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
    
      while(std::getline(iss, token, ','))
	tokens.push_back(token);

      // How many element were read from each line 
      N.push_back(tokens.size());

      if(line_number == 1) {} // skip the headerdo nothing with the file header
    
      else{
    
	ProjectionInputFiles SingleLine;

	// Assign entries in each line to a variable of type Parameter
	SingleLine.LongName = tokens[0];
	SingleLine.ShortName = tokens[1];
	SingleLine.Path = tokens[2];

	//std::cout << "Token[0] is " << tokens[0] << std::endl;
	ProjectionInputFilesVector.push_back(SingleLine); // Add parameter to the vector of parameters

      } // End of else
    } // End of if(line.length() > 0
  } // End of while

  file.close();

  // Check that all lines contains the same number of elements separated by commas
  if (std::adjacent_find(N.begin(), N.end(), std::not_equal_to<int>() ) != N.end() )
    error("Inconsistent number of elements separated by commas in ",filename);
  
  // End of function execution message
  //std::cout << "Finished reading " << ParameterVector.size() << " parameters descriptions from file " << filename << "\n";
 
  return ProjectionInputFilesVector;
}

// Get projection file path according to short name
std::string GetProjectionFilePathAccordingToShortName(const std::vector<ProjectionInputFiles> &InputFileVector, const std::string &ShortName){

    for(unsigned int i = 0; i < InputFileVector.size(); ++i)
      {
	if( InputFileVector[i].ShortName.compare(ShortName) == 0) 
	{
	  return InputFileVector[i].Path;
	}
      }
    // 	Return an error if the parameter looked for was not found
    error("Error in GetProjectionFilePathAccordingToShortName(): no value found for short name " + ShortName + " in input parameter file.");
    return "failed";
}

// A function to read a single raw of data in a file
// applies to Stock Recruitment parameters, weekly effort pattern, weekly availability, weekly percentage of maturity
std::vector<double> ReadParameterFromSingleColumnFile(const std::string& filename)
{

  // file variable 
  std::string line;  // store lines in file
  std::vector<int> N; // Number of elements per line
  unsigned int line_number=0; // counter

  // Container for the parameters information in the file
  std::vector<double> ValueVector;

  // Open file for reading
  ifstream file(filename.c_str());    // open file for reading
  if (!file) error("can't open input file ", filename);

  // Read lines from file
  while(getline(file,line)){

    if(line.length() > 0){
      line_number++; // count the number of lines

      // Read lines content into a vector
      std::istringstream iss(line);
      std::string token;
      std::vector<std::string> tokens;
    
      while(std::getline(iss, token, ','))
	tokens.push_back(token);

      // How many element were read from each line 
      N.push_back(tokens.size());

	ValueVector.push_back(atof(token.c_str())); // Add parameter to the vector of parameters

    } // End of if(line.length() > 0

  } // End of while

  file.close();

  // Check that all lines contains the same number of elements separated by commas
  if(std::adjacent_find(N.begin(), N.end(), std::not_equal_to<int>() ) != N.end() )
    error("Inconsistent number of elements separated by commas in ",filename);

  // End of function execution message
  //std::cout << "Finished reading " << ValueVector.size() << " parameters descriptions from file " << filename << "\n";

  return ValueVector;
}

