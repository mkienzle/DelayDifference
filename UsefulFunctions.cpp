// CREATED  29 April 2014
// MODIFIED 21 July  2014

// AUTHOR Marco.Kienzle@gmail.com

// PURPOSE some useful functions

// Read variable from files
#include <iostream>
#include <fstream>
#include <string>
#include "lib_facilities.h"
#include <boost/algorithm/string.hpp>

// Output info about the delay difference
void banner(string version_nb, string version_date)
{

  // Print a starting message
  std::cout << "\"Down under\" delay difference model version " << version_nb << " (" << version_date << ")\n";
  std::cout << "Copyright (C) 2013-15 Queensland Government, Department of Agriculture, Forestry and Fisheries\n";
  std::cout << "This code is distributed under the GNU GPL v.3 license (http://www.gnu.org/licenses/)\n\n";

}

// Display file content
void display_file(string& filename)
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
      cout << "Line " << line_number << ": " << line << endl;
      
      if(boost::find_first(line, "Natural mortality")){
	cout << "The line below contains the parameter you want" << endl;
      }

    }
    cout << "bye\n";
    // the file is implicitly closed when we leave the function
}

// 
// Get the value of a variable from filename given a VariableDescription
//
double fill_from_file(string& filename, const string VariableDescription)
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

std::vector<double> vonMisesRecDist(double a, double b){

// Global variables
extern int NIPY;


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
      if(intervals[i] < 0.0)
	// if they are small, replace them by zero
	if(abs(intervals[i]) < 1e-5){ intervals[i] = 0.0;}
      // otherwise throw an error
	else{
	  throw;}

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
