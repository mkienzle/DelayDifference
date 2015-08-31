// @(#)root/minuit2:$Id: GaussFcn.h 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

// REFERENCE C++ MINUIT User's guide by F. James and M. Winkler, CERN, Geneva, 16 Jun 2004

#ifndef MN_LogLikelihoodFunction_H_
#define MN_LogLikelihoodFunction_H_

#include "Minuit2/FCNBase.h"

#include <vector>

namespace ROOT {

   namespace Minuit2 {

     // This objective function inherit from the abstract class FCNBase (more at root.cern.ch/root/html/ROOT_Minuit2_FCNBase.html)
class LogLikelihoodFunction : public FCNBase {

public:

     LogLikelihoodFunction(const std::vector<double>& catch_values,
			   const std::vector<double>& targetedeffort_values, const std::vector<double>& nontargetedeffort_values) : fCatchValues(catch_values), fTargetedEffortValues(targetedeffort_values), fNontargetedEffortValues(nontargetedeffort_values),
	 // Error definition of the function is 0.5 for negative log likelihood (see Minuit manual for more infos)
					      fErrorDef(0.5) {}

  ~LogLikelihoodFunction() {}

  virtual double Up() const {return fErrorDef;}
  virtual double operator()(const std::vector<double>&) const;

  std::vector<double> Zvalues() const {return fNontargetedEffortValues;}  
  std::vector<double> Yvalues() const {return fTargetedEffortValues;}
  std::vector<double> Xvalues() const {return fCatchValues;}

  void SetErrorDef(double def) {fErrorDef = def;}

private:

  std::vector<double> fNontargetedEffortValues;
    std::vector<double> fTargetedEffortValues;
  std::vector<double> fCatchValues;
  double fErrorDef;
};

  }  // namespace Minuit2

}  // namespace ROOT

#endif //MN_LogLikelihoodFunction_H_
