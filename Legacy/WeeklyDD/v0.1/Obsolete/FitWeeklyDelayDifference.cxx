// CREATED  10 Sep 2013
// MODIFIED 11 Sep 2013

// PURPOSE fit a delay difference model (with weekly time-steps) to a series of catch and effort data
//         by maximum likelihood

// COMPILE g++ -g -o FitWeeklyDelayDifference FitWeeklyDelayDifference.cxx LogLikelihoodFunction.cxx WeeklyDelayDifference.cxx `root-config --glibs --cflags` -lMinuit2;

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"

#include <stdio.h>
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
# include "TTree.h"

# include "WeeklyDelayDifference.h"
# include "LogLikelihoodFunction.h"
# include "FixPar.h"

using namespace ROOT::Minuit2;

int main() {

  // Read catch and effort data from file
  TTree *TigerPrawn = new TTree("TigerPrawnData", "Moreton Bay tiger prawn");
  TigerPrawn->ReadFile("Data/SimData2.txt", "TimeStep/D:Catch/D:Effort/D");

  int counter = (int) TigerPrawn->GetEntries();

 //  std::ifstream ifs( "Data/SimData2.txt" );
  std::vector< double > timestep;
 std::vector< double > TimeStep, FishCatch, Effort;
 //  std::vector< double > timestep, FishCatch, Effort;
 //  double a, b, c;
 //  int counter=0;

 //  while( ifs >> a >> b >> c){
 //    timestep.push_back( a ); FishCatch.push_back( b ); Effort.push_back( c);
 //    //printf(" Reading %.2f %.2f %.2f from file \n", a, b, c);
 //    counter++;
 // }

  // Plot a timeseries of catch
  //  Int_t TimeStep;
  //  Double_t TimeStep, FishCatch;
  
  //  TigerPrawn->SetBranchAddress("Catch", &FishCatch);
  //  TigerPrawn->SetBranchAddress("TimeStep", &TimeStep);

  TCanvas *cvs = new TCanvas();
  TigerPrawn->Draw("Catch:TimeStep");

  //TGraph* graph = new TGraph(counter, &TimeStep[0], &FishCatch[0]);
  //  TGraph* graph = new TGraph(counter, &timestep[0], &FishCatch[0]);

  // TGraph* graph = new TGraph((int) TigerPrawn->GetEntries(), TigerPrawn->GetV1(), TigerPrawn->GetV2());

  //  graph->SetTitle("");
  //  graph->SetMarkerStyle(20);
  //  graph->SetMarkerSize(0.75);
  //  graph->Draw("ALP");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Standard minimization using MIGRAD
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Declare the objective function (also called FCN function in documentation)  
  TigerPrawn->Branch("Catch", "vector<double>", &FishCatch);
  TigerPrawn->Branch("Effort", "vector<double>", &Effort);

    LogLikelihoodFunction fFCN(FishCatch, Effort);
  //  LogLikelihoodFunction fFCN(TigerPrawn->GetV2(), TigerPrawn->GetV3());

  // create Minuit parameters with names
  MnUserParameters upar;
  upar.Add("catchability", 1, 2);
  upar.Add("sigma", 10, 0.1);
 
  // Assert parameters domain
  upar.SetLimits("catchability", 0, 1e3);
  upar.SetLimits("sigma", 0, 1e3);

  // create MIGRAD minimizer
    MnMigrad migrad(fFCN, upar);

    // Minimize
    FunctionMinimum min = migrad();

    // output
    std::cout<<"minimum: "<<min<<std::endl;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plot the results
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
    // Calculate estimated catch
    double EstimatedCatches[counter], EstimatedBiomass[counter], EstimatedFishingMort[counter];
    double estPar;
    estPar = min.UserState().Value("catchability");

    for(unsigned int iii=0; iii < counter; iii++)
      {
	EstimatedFishingMort[iii] = min.UserState().Value("catchability") * CatchabilityScalingFactor * Effort[iii];
	EstimatedBiomass[iii] = WeeklyDD(&iii, &Effort[0], &estPar);
	EstimatedCatches[iii] = EstimatedFishingMort[iii] / (M + EstimatedFishingMort[iii]) * EstimatedBiomass[iii] * \
	  (1 - exp(- (M + EstimatedFishingMort[iii])));
      }

    TGraph* gr2 = new TGraph(counter, &timestep[0], &EstimatedCatches[0]);
    gr2->SetLineWidth(2);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(0.01);
    gr2->SetLineColor(4); // blue
    gr2->Draw("CP");

    // Save graph into a PDF
    cvs->Print("Results/Graphics/WeeklyDDFittedToData.pdf");

    // Plot an histogram of residuals of the fit
    TCanvas *c2 = new TCanvas();
    c2->Divide(1,2);
    c2->cd(1);

    TH1F *h1 = new TH1F("Summary", "Residuals of weekly delay difference model fit",50,-5,5);

    // Fill histogram with residuals
    for(unsigned int iii=0; iii < counter; iii++)
      h1->Fill( (sqrt(EstimatedCatches[iii]) - sqrt(FishCatch[iii])) / min.UserState().Value("sigma"));

    h1->SetFillColor(4);
    h1->Draw();

    h1->Fit("gaus");

    ///////////////////////////////////////////////
    /// Below, I would like to plot the shape of the likelihood function
    // see rootUsers_Guide_5_08.pdf p.28
    // NOT WORKING FOR THE MOMENT
    /////////////////////////////////////////////////////


// //make it:
// TH2D* h2=new TH2D("Likelihood surface","summary",5,2,8,5,10,30);

// //fill it:
// // std::vector< double > ParVect;

// for (int i=2000;i<8000;i++){for (int j=10000;j<30000;j++){

 
//  std::vector< double > ParVect;
//  ParVect.push_back( (double) i/1e3 );
//  ParVect.push_back( (double) j/1e3 );

//  h2->Fill((double) i/1e3, (double) j/1e3, LogLikelihoodFunction::operator()(&ParVect));
//          }}

//      c2->cd(2);

//        h2->Draw("lego2");

//     c2->Print("Results/Graphics/test.pdf");

  //  // Wec generate a 2-D function
  //  TF2 *f2 = new TF2("f2","x**2 + y**2 - x**3 -8*x*y**4",-1,1.2,-1.5,1.5);
  // TF1 *f1 = new TF1("f1", line, -100, 100, 2);

  //  f2->SetContour(48);
  //  f2->SetFillColor(45);

  //  // Draw this function in pad1 with Gouraud shading option
  //  pad1->cd();
  //  pad1->SetPhi(-80);
  //  pad1->SetLogz();
  //  f2->Draw("surf4");

  //  // Draw this function in pad2 with color mesh option
  //  pad2->cd();
  //  pad2->SetTheta(25);
  //  pad2->SetPhi(-110);
  //  pad2->SetLogz();
  //  f2->Draw("surf1");

  return 0;
}
