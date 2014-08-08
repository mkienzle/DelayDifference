// @(#)root/minuit2:$Id: GaussFcn.h 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_LogLikLinearReg_H_
#define MN_LogLikLinearReg_H_

#include "Minuit2/FCNBase.h"

#include <vector>


// linear function
     double line(const double *x, const double *par);

namespace ROOT {

   namespace Minuit2 {



class LogLikLinearReg : public FCNBase {

public:

     LogLikLinearReg(const std::vector<double>& x_values,
	      const std::vector<double>& y_values) : fXvalues(x_values),
	                                      fYvalues(y_values),
					      fErrorDef(0.5) {}

  ~LogLikLinearReg() {}

  virtual double Up() const {return fErrorDef;}
  virtual double operator()(const std::vector<double>&) const;
  
  std::vector<double> Yvalues() const {return fYvalues;}
  std::vector<double> Xvalues() const {return fXvalues;}

  void SetErrorDef(double def) {fErrorDef = def;}

private:

  
  std::vector<double> fYvalues;
  std::vector<double> fXvalues;
  double fErrorDef;
};

  }  // namespace Minuit2

}  // namespace ROOT

#endif //MN_LogLikLinearReg_H_
