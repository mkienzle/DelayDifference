// CREATED   9 Sep. 2013
// MODIFIED 10 Sep. 2013

// PURPOSE demonstrate the use of MINUIT for fitting a linear regression by maximum likelihood

// STATUS working

// MODIFIED FROM @(#)root/minuit2:$Id: DemoGaussSim.cxx 20880 2007-11-19 11:23:41Z rdm $
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  
// Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT 


// COMPILE USING g++ -g -o FitLinearRegression LogLikLinearReg.cxx FitLinearRegression.cxx `root-config --glibs --cflags` -lMinuit2
// OR COMPILE USING g++ -g -o FitLinearRegression LogLikLinearReg.cxx FitLinearRegression.cxx -pthread -m64 -I/usr/local/include/root -L/usr/local/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lMinuit2

// NOTE simulated data were create om R using
// dta <- cbind(rep(seq(1,50), 2), rnorm( 100, mean = -10 + 2.35 * rep(seq(1,50), 2), sd = 1.5))
// write.table(file = "SimData.txt", dta, row.names = FALSE, col.names = FALSE)

#include "LogLikLinearReg.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
# include "Riostream.h"
# include "TFile.h"

# include "TH1F.h"
# include "TCanvas.h"
# include "TStyle.h"
# include "TGraph.h"
# include "TROOT.h"
# include "TH2F.h"
# include "TLine.h"
# include "TF1.h"

using namespace ROOT::Minuit2;



int main() {

  // Read in data from file
  std::ifstream ifs( "SimData.txt" );

  std::vector< double > col1, col2;
  double a, b;
  int counter=0;

  while( ifs >> a >> b ){
    col1.push_back( a );
    col2.push_back( b );
   //   printf(" Reading %.2f from file \n", col1);
    counter++;
 }

  std::vector<double> xvalues = col1;
  std::vector<double> yvalues = col2;

  int x1[5] = {1,10,100,1000,10000};

  TCanvas *cvs = new TCanvas();
  TGraph* graph = new TGraph(counter, &xvalues[0], &yvalues[0]);
  //TH1F * h = new TH1F("h", "my first TH1F", 20, 0, 120);
  //TGraphErrors *gr = new TGraphErrors(); //counter, xvalues, yvalues, 0, 0);
  cvs->Divide(2,1);
  cvs->cd(1);
  graph->SetTitle("Linear regression by maximum likelihood");
  graph->SetMarkerStyle(20);
  graph->Draw("AP");
  //h->Draw();

  // create FCN function  
  LogLikLinearReg fFCN(xvalues, yvalues);


  // {
  //   // demonstrate minimal required interface for minimization
  //   // create Minuit parameters without names

  //   // starting values for parameters
  //   std::vector<double> init_par; 
  //   init_par.push_back(1.); 
  //   init_par.push_back(2.); 
  //   init_par.push_back(10.);

  //   // starting values for initial uncertainties
  //   std::vector<double> init_err; 
  //   init_err.push_back(0.1); 
  //   init_err.push_back(0.1); 
  //   init_err.push_back(1.0);

  //   printf("The value of the line is %.2f\n", line( &xvalues[1], &init_par[0]));
    
  //   // create minimizer (default constructor)
  //   VariableMetricMinimizer fMinimizer;
    
  //   // Minimize
  //   FunctionMinimum min = fMinimizer.Minimize(fFCN, init_par, init_err);

  //   // output
  //   std::cout<<"minimum: "<<min<<std::endl;

  //  }

  {
    // demonstrate standard minimization using MIGRAD
    // create Minuit parameters with names
    MnUserParameters upar;
    upar.Add("intercept", 1, 0.1);
    upar.Add("slope", 2, 0.1);
    upar.Add("sigma", 1, 0.2);

    // create MIGRAD minimizer
    MnMigrad migrad(fFCN, upar);

    // Minimize
    FunctionMinimum min = migrad();

    // output
    std::cout<<"minimum: "<<min<<std::endl;

    // Overlay regression line with data
    //TLine *line = new TLine(0,0, 20, 100);
    //line->SetLineColor(kRed);
    //line->Draw();
    TF1 *f1 = new TF1("f1", line, -100, 100, 2);
    f1->SetParameters(min.UserState().Value("intercept"),min.UserState().Value("slope"));
    f1->Draw("same");
    //    printf("Par 1 is %.f, par2 is %.f\n", min.Value("intercept"), upar.Value("slope"));
    //    m.getParameters("intercept")
  }

  // {
  //   // demonstrate full interaction with parameters over subsequent 
  //   // minimizations

  //   // create Minuit parameters with names
  //   MnUserParameters upar;
  //   upar.Add("mean", mean, 0.1);
  //   upar.Add("sigma", rms, 0.1);
  //   upar.Add("area", area, 0.1);

  //   // access Parameter by Name to set limits...
  //   upar.SetLimits("mean", mean-0.01, mean+0.01);

  //   // ... or access Parameter by Index
  //   upar.SetLimits(1, rms-0.1, rms+0.1);
    
  //   // create Migrad minimizer
  //   MnMigrad migrad(fFCN, upar);

  //   // Fix a Parameter...
  //   migrad.Fix("mean");

  //   // ... and Minimize
  //   FunctionMinimum min = migrad();

  //   // output
  //   std::cout<<"minimum: "<<min<<std::endl;

  //   // Release a Parameter...
  //   migrad.Release("mean");

  //   // ... and Fix another one
  //   migrad.Fix(1);

  //   // and Minimize again
  //   FunctionMinimum min1 = migrad();
 
  //   // output
  //   std::cout<<"minimum1: "<<min1<<std::endl;

  //   // Release the Parameter...
  //   migrad.Release(1);

  //   // ... and Minimize with all three parameters (still with limits!)
  //   FunctionMinimum min2 = migrad();
    
  //   // output
  //   std::cout<<"minimum2: "<<min2<<std::endl;

  //   // remove all limits on parameters...
  //   migrad.RemoveLimits("mean");
  //   migrad.RemoveLimits("sigma");

  //   // ... and Minimize again with all three parameters (now without limits!)
  //   FunctionMinimum min3 = migrad();

  //   // output
  //   std::cout<<"minimum3: "<<min3<<std::endl;
  // }

  // {
  //   // test single sided limits
  //   MnUserParameters upar;
  //   upar.Add("mean", mean, 0.1);
  //   upar.Add("sigma", rms-1., 0.1);
  //   upar.Add("area", area, 0.1);

  //   // test Lower limits
  //   upar.SetLowerLimit("mean", mean-0.01);

  //   // test Upper limits
  //   upar.SetUpperLimit("sigma", rms-0.5);

  //   // create MIGRAD minimizer
  //   MnMigrad migrad(fFCN, upar);

  //   // ... and Minimize
  //   FunctionMinimum min = migrad();
  //   std::cout<<"test Lower limit minimim= "<<min<<std::endl;
  // }

  // {
  //   // demonstrate MINOS Error analysis

  //   // create Minuit parameters with names
  //   MnUserParameters upar;
  //   upar.Add("mean", mean, 0.1);
  //   upar.Add("sigma", rms, 0.1);
  //   upar.Add("area", area, 0.1);

  //   // create Migrad minimizer
  //   MnMigrad migrad(fFCN, upar);

  //   // Minimize
  //   FunctionMinimum min = migrad();

  //   // create MINOS Error factory
  //   MnMinos Minos(fFCN, min);

  //   {
  //     // 1-sigma MINOS errors (minimal interface)
  //     std::pair<double,double> e0 = Minos(0);
  //     std::pair<double,double> e1 = Minos(1);
  //     std::pair<double,double> e2 = Minos(2);
      
  //     // output
  //     std::cout<<"1-sigma Minos errors: "<<std::endl;
  //     std::cout<<"par0: "<<min.UserState().Value("mean")<<" "<<e0.first<<" "<<e0.second<<std::endl;
  //     std::cout<<"par1: "<<min.UserState().Value(1)<<" "<<e1.first<<" "<<e1.second<<std::endl;
  //     std::cout<<"par2: "<<min.UserState().Value("area")<<" "<<e2.first<<" "<<e2.second<<std::endl;
  //   }

  //   {
  //     // 2-sigma MINOS errors (rich interface)
  //     fFCN.SetErrorDef(4.);
  //     MinosError e0 = Minos.Minos(0);
  //     MinosError e1 = Minos.Minos(1);
  //     MinosError e2 = Minos.Minos(2);
      
  //     // output
  //     std::cout<<"2-sigma Minos errors: "<<std::endl;
  //     std::cout<<e0<<std::endl;
  //     std::cout<<e1<<std::endl;
  //     std::cout<<e2<<std::endl;
  //   }
  // }

  // {
  //   // demostrate MINOS Error analysis with limits

  //   // create Minuit parameters with names
  //   MnUserParameters upar;
  //   upar.Add("mean", mean, 0.1);
  //   upar.Add("sigma", rms, 0.1);
  //   upar.Add("area", area, 0.1);

  //   double meanLow = -50.03;
  //   double rmsUp = 1.55;
  //   std::cout << "sigma Limit: " << rmsUp << "\tmean limit: " << meanLow << std::endl;
  //   // test Lower limits
  //   upar.SetLowerLimit("mean", meanLow);
  //   // test Upper limits
  //   upar.SetUpperLimit("sigma", rmsUp);

  //   // create Migrad minimizer
  //   MnMigrad migrad(fFCN, upar);

  //   // Minimize
  //   FunctionMinimum min = migrad();

  //   // create MINOS Error factory
  //   MnMinos Minos(fFCN, min);

  //   {
  //     // 3-sigma MINOS errors (minimal interface)
  //     fFCN.SetErrorDef(9.);
  //     std::pair<double,double> e0 = Minos(0);
  //     std::pair<double,double> e1 = Minos(1);
  //     std::pair<double,double> e2 = Minos(2);

      
  //     // output
  //     std::cout<<"3-sigma Minos errors with limits: "<<std::endl;
  //     std::cout.precision(16);
  //     std::cout<<"par0: "<<min.UserState().Value("mean")<<" "<<e0.first<<" "<<e0.second<<std::endl;
  //     std::cout<<"par1: "<<min.UserState().Value(1)<<" "<<e1.first<<" "<<e1.second<<std::endl;
  //     std::cout<<"par2: "<<min.UserState().Value("area")<<" "<<e2.first<<" "<<e2.second<<std::endl;


  //   }

  // }

  // {
  //   // demonstrate how to use the CONTOURs

  //   // create Minuit parameters with names
  //   MnUserParameters upar;
  //   upar.Add("mean", mean, 0.1);
  //   upar.Add("sigma", rms, 0.1);
  //   upar.Add("area", area, 0.1);

  //   // create Migrad minimizer
  //   MnMigrad migrad(fFCN, upar);

  //   // Minimize
  //   FunctionMinimum min = migrad();

  //   // create contours factory with FCN and Minimum
  //   MnContours contours(fFCN, min);
  
  //   //70% confidence level for 2 parameters Contour around the Minimum
  //   // (minimal interface)
  //   fFCN.SetErrorDef(2.41);
  //   std::vector<std::pair<double,double> > cont = contours(0, 1, 20);

  //   //95% confidence level for 2 parameters Contour
  //   // (rich interface)
  //   fFCN.SetErrorDef(5.99);
  //   ContoursError cont4 = contours.Contour(0, 1, 20);
    
  //   // plot the contours
  //   MnPlot plot;
  //   cont.insert(cont.end(), cont4().begin(), cont4().end());
  //   plot(min.UserState().Value("mean"), min.UserState().Value("sigma"), cont);

  //   // print out one Contour
  //   std::cout<<cont4<<std::endl;
  // }

  // Output graphics
  cvs->Print("Results/Graphics/LinearRegression.pdf");
  gStyle->SetPaperSize(10.,10.);
  gPad->Print("Results/Graphics/LinearRegression.tex");

  return 0;
}
