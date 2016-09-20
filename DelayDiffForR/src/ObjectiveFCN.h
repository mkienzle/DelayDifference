// @(#)root/minuit2:$Id: GaussFcn.h 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

// REFERENCE C++ MINUIT User's guide by F. James and M. Winkler, CERN, Geneva, 16 Jun 2004

#ifndef MN_ObjectiveFCN_H_
#define MN_ObjectiveFCN_H_

#include "Minuit2/FCNBase.h"
#include <Rcpp.h>

#include <vector>

namespace ROOT {

   namespace Minuit2 {

     // This objective function inherit from the abstract class FCNBase (more at root.cern.ch/root/html/ROOT_Minuit2_FCNBase.html)
class ObjectiveFCN : public FCNBase {

public:

 ObjectiveFCN(const Rcpp::NumericVector &length, const Rcpp::NumericVector &weight, const Rcpp::NumericVector &sd) : fLength(length), fWeight(weight), fsd(sd),
	 // Error definition of the function is 0.5 for negative log likelihood (see Minuit manual for more infos)
	 // error definition (=1. for getting 1 sigma error for chi2 fits)
					      fErrorDef(1.) {}

  ~ObjectiveFCN() {}

  virtual double Up() const {return fErrorDef;}
  virtual double operator()(const std::vector<double>&) const;


  //std::vector<double> Xvalues() const {return fLength;}
  //std::vector<double> Yvalues() const {return fWeight;}
  //std::vector<double> Zvalues() const {return fsd;}
		  
  void SetErrorDef(double def) {fErrorDef = def;}

private:

  //std::vector<double> fLength;
  //std::vector<double> fWeight;
  //std::vector<double> fsd;

  Rcpp::NumericVector fLength;
  Rcpp::NumericVector fWeight;
  Rcpp::NumericVector fsd;
  
  double fErrorDef;
};

  }  // namespace Minuit2

}  // namespace ROOT

#endif //MN_ObjectiveFCN_H_
