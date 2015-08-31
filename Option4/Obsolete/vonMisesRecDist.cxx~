// CREATED  12 November 2013
// MODIFIED 20 May 2014

// PURPOSE compute the probability from a von mises distribution with parameter a and b for every week of the year

// DEPENDENCY on the C++ library PROB by John Burkardt
//            available at /home/mkienzle/mystuff/Software/Prob_library/prob.cpp
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>

#include <prob.hpp>

//#include "FixPar.h"

using namespace std;

std::vector<double> vonMisesRecDist(double a, double b){

  // Global variables
  extern int NIPY;

  //  cout << "\n";
  //  cout << "  PDF parameter A =      " << a << "\n";
  //  cout << "  PDF parameter B =      " << b << "\n";

  if ( !von_mises_check ( a, b ) )
  {
    cout << "\n";
    cout << "vonMisesRecDist - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    throw;
  }

  // Assume a year is made of NIPY intervals per year
  
  // divide the range [-pi;pi] into NIPY intervals
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
 
    //cout << "von mises density equal " << total << "\n";
    // throw exception if cumulative density 
    if(abs(total - 1) > 1e-2){
      cout << "von mises density equal " << total << "\n";
      throw;
    }
    //// cout << "And the sum of prob is " << total << "\n";
    return intervals;

}
