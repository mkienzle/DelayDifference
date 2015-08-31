// CREATED  29 April 2014
// MODIFIED 29 April 2014

// AUTHOR Marco.Kienzle@gmail.com

// PURPOSE some useful functions

// Read variable from files
#include <iostream>
#include <fstream>
#include <string>
#include "../std_lib_facilities.h"
#include <boost/algorithm/string.hpp>

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

    if(FoundVariable == 0) cout << "*** Warning *** : could find a value for " << VariableDescription << " in file " << filename << endl;
    // the file is implicitly closed when we leave the function
    return ParameterFromFile;
}

