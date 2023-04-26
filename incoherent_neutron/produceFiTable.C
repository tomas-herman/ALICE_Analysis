
//
// produces a table with the infoherent fractions values  
//

// -----------------------------------------------------------------
// all headers are defined here

// c++ headers
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

// root headers
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

// RooFit headers
#include <RooGlobalFunc.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooAbsReal.h"
#include "RooGenericPdf.h"
#include "RooCrystalBall.h"

#include "utilities.h"

using namespace RooFit;

void produceFiTable(int iData = 0, int iSel = 1)
{
  // set mass ranges for the fit
  float mMin = 2.85; // 2.85;
  float mMax = 3.35; // 3.55;
  if (iData == 12) { // psi2s
    mMin = 3.3;
    mMax = 4;
  }
  // get the input tree
  TTree *dataTree = getTree(iData, iSel);
  // define kinematic ranges for data to fit 
  float minRap[4] = {-4.0, -4.00, -3.50, -3.00};
  float maxRap[4] = {-2.5, -3.50, -3.00, -2.50};

  auto znSelections = {0,1,2,3};
  
  double fI[2] = {0,0};
  // ---------
  // make a fit
  for (int i = 1; i < 4; i++) {
    cout << abs(maxRap[i]) << " < |y| < " << abs(minRap[i]) << endl;
    for ( auto& znSelection : znSelections ) {
      fI[0] = getPtFitResults(iSel,znSelection,"val","fiParam",mMin,mMax,minRap[i],maxRap[i]);
      fI[1] = getPtFitResults(iSel,znSelection,"err","fiParam",mMin,mMax,minRap[i],maxRap[i]);
      TString znClass;
      if (znSelection == 0) znClass = "0n0n";
      if (znSelection == 1) znClass = "0nXn";
      if (znSelection == 2) znClass = "Xn0n";
      if (znSelection == 3) znClass = "XnXn";
      cout << znClass.Data() << "   " << fI[0] << " +/- " << fI[1] << endl;
    }
    cout << "------------------------------------------------------" << endl;
  }
}

