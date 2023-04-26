//
// This program makes pt histograms from MC data passing given selection
// criteria. Those are then used to construct PDF templates for fitting 
//

// -----------------------------------------------------------------
// all headers are defined here

// c++ headers
#include <iostream>
#include <fstream>
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

#include "utilities.h"

//-----------------------------------------------
// Make a histo with pt template
void makeOneMcPtTemplate(TTree *dataTree, int iData, int iSel, 
                                         double minMass, double maxMass,
                                         double minRap, double maxRap)
{
  // TCanvas *c1 = new TCanvas(Form("c%i_%f_%f",iData, minRap, maxRap),"pT template PDFs",1200,800);

  // c1->cd();
  // // Set general options
  // gStyle->SetOptTitle(0);
  // gStyle->SetOptStat(0);
  // gStyle->SetPalette(kBird);
  // gStyle->SetPaintTextFormat("4.3f");
  // gStyle->SetFrameLineWidth(1);
  // gStyle->SetLabelSize(0.045,"xyz");
  // gStyle->SetTitleSize(0.05,"xyz");
  // gStyle->SetTextSize(0.04); 
 
  // Set pt binning
  int nBinsPt= 120;
  double ptRangeMin = 0.; 
  double ptRangeMax = 3.; 

  // set the input variables
  double pt;
  double mass;
  double rap;
  dataTree->SetBranchAddress("pt",&pt);
  dataTree->SetBranchAddress("mass",&mass);
  dataTree->SetBranchAddress("rap",&rap);

  // define histogram
  TH1D *h = new TH1D("h","pt template", nBinsPt, ptRangeMin, ptRangeMax);; 

  int nEntries = dataTree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    dataTree->GetEntry(i);
    if (pt < ptRangeMin) continue;
    if (pt > ptRangeMax) continue;
    if (mass < minMass) continue;
    if (mass > maxMass) continue;
    if (rap < minRap) continue;
    if (rap > maxRap) continue;
    h->Fill(pt);
  }

  h->SetLineColor(kBlack);
  h->SetLineWidth(1);
  h->GetXaxis()->SetTitle("pt [GeV/c]");
  h->GetYaxis()->SetTitle("#Entries");
  // h->Draw("");

  // c1->Draw();

  // set up the output file to save the pt hist
  TFile *fOut = new TFile(Form("ptHistos/ptHisto_data_%i_%.2f_%.2f.root", iData, abs(minRap), abs(maxRap)),"recreate");

  // write the histo to the output file
  fOut->cd();
  h->Write();
  fOut->Close();

  h->Delete();
}

void makeMcPtTemplates(int iData = 10, int iSel = 1)
{
  // set mass ranges for the pt template
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
  // float minRap[7] = {-4.0, -4.00, -3.75, -3.50, -3.25, -3.00, -2.75};
  // float maxRap[7] = {-2.5, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50};

  for (int i = 0; i < 4; i++){
    makeOneMcPtTemplate(dataTree, iData, iSel, mMin, mMax, minRap[i], maxRap[i]);
  }
}