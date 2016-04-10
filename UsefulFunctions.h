// CREATED  29 April 2014
// MODIFIED  5 April 2016

// AUTHOR Marco.Kienzle@gmail.com

// PURPOSE some useful functions

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

// Write parameters to file
void WriteParameterToCSV(const std::string& filename, std::vector<Parameter> ParVec);

