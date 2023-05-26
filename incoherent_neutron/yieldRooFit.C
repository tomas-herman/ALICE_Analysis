
//
// makes a fit to the invariant mass distribution to data from produceFitTree.C 
// and produces gamma gamma yield values for pt template PDF 
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

// -----------------------------------------------------------------
// do one mass fit for MC
void doOneSimpleDscbFit(TTree *dataTree, int iSel, int znSelection,
        double minPt, double maxPt,
        double minMass, double maxMass,
        double minRap, double maxRap)
{
  // define number of bins
  auto nBins = 80;
  // define variables for the tree
  RooRealVar pt("pt","pt",minPt,maxPt);
  RooRealVar mass("mass","mass",2,5);
  mass.setRange("m_fitRange",minMass, maxMass);
  mass.setRange("m_JPsiRange",2.85, 3.35);
  RooRealVar rap("rap","rap",minRap,maxRap);
  RooRealVar znclass("znclass","znclass",znSelection,znSelection);

  // import the tree
  RooDataSet inData("inData", "inData", RooArgSet(pt, mass, rap, znclass), Import(*dataTree));
  int nEvents = inData.numEntries();

  // build pdfs
  // ---double sided CB
  RooRealVar m0("m0", "m0", 3.1,minMass,maxMass);
  RooRealVar sigmaL("sigmaL", "sigmaL",0.08, 0.04, 0.15);
  RooRealVar sigmaR("sigmaR", "sigmaR",0.08, 0.04, 0.15);  
  RooRealVar nL("nL", "nL", 10,10,10);
  nL.setConstant();    
  RooRealVar nR("nR", "nR", 10,10,10);
  nR.setConstant();    
  RooRealVar alphaL("alphaL", "alphaL", 1.27, 0.1, 3);
  RooRealVar alphaR("alphaR", "alphaR", 1.50, 0.1, 3);  
  RooCrystalBall crystalBall("crystalBall", "crystalBall", mass, m0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);

  // ---combine the PDFs
  RooRealVar N_CB("N_{J/Psi}","Number of CB events",1.00*nEvents,0.8*nEvents,1.2*nEvents);

  RooAddPdf fitFunc("fitFunc","DS CB PDF", RooArgList(crystalBall), RooArgList(N_CB));

   // fit
  RooFitResult* r = fitFunc.fitTo(inData,Extended(kTRUE),Range("m_fitRange"),Save());

  // Compute number of events from the mass fit
  // ---Define variabels
  double N_JPsi[2];
  // ---Compute the number of events
  N_JPsi[0] = N_CB.getVal();
  N_JPsi[1] =  N_CB.getError();

  //-----------------------------------------------------------
  // Draw mass histogram
  // -------------------------------------------------------------------------------- 
  
  // Set general options
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTextSize(0.04);
  
    //-----------------------------------------------------------
  // plot mass distributions

  TCanvas *c = new TCanvas(Form("pt in (%.2f,%.2f), ZN class = %d",minPt,maxPt,znSelection),
         Form("pt in (%.2f,%.2f), ZN class = %d",minPt,maxPt,znSelection),
         1600, 800);
  c->Divide(2);      
  
  //----------------------------
  //Draw correlation matrix
  c->cd(2);
  TPad *c2 = new TPad("c2","c2",0.001,0.001,0.999,0.999);
  // --- pad option
  c2->SetRightMargin(0.15);
  c2->SetLeftMargin(0.15);
  c2->Draw();
  c2->cd();
  // --- TH2D correlation matrix
  TH2* h_CorrM = r->correlationHist();
  h_CorrM->SetMarkerSize(1.2);
  h_CorrM->Draw("zcol,text");

  //------------------------------
  // Draw mass fit
  c->cd(1);
  TPad *c1 = new TPad("c1","c1",0.001,0.001,0.999,0.999);
   // --- pad option
  c1->SetRightMargin(0.03);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.12);
  c1->Draw();
  c1->cd();

  // plot mass distributions
  RooPlot *frame = mass.frame(Title(Form("pt in (%.2f,%.2f)",minPt,maxPt)));
  inData.plotOn(frame,Binning(nBins),MarkerStyle(20),MarkerSize(0.5));
  fitFunc.plotOn(frame, LineColor(kBlack));
  crystalBall.paramOn(frame,Layout(0.65, 0.9, 0.65));
  frame->getAttText()->SetTextSize(0.02);
  frame->GetXaxis()->SetTitle("#it{m_{#mu^{+}#mu^{-}}} (GeV/#it{c}^{2})");
  frame->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c^{2}} ", (maxMass-minMass)/nBins*1000));
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->Draw();

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.038);
  leg->AddEntry((TObject*)0,Form("%3.2f < p_{T} < %3.2f",minPt,maxPt),"");
  leg->AddEntry((TObject*)0,Form("N_{J/#psi} = %3.0f #pm %3.0f",N_JPsi[0],N_JPsi[1]),"");
  leg->Draw();
}

// -----------------------------------------------------------------
// do one mass fit 
void doOneFit(TTree *dataTree, int iSel, int znSelection,
		    double minPt, double maxPt,
		    double minMass, double maxMass,
		    double minRap, double maxRap,
        TFile *fOut, int bin)
{
  // define number of bins
  auto nBins = 80;
  // define variables for the tree
  RooRealVar pt("pt","pt",minPt,maxPt);
  RooRealVar mass("mass","mass",2,5);
  mass.setRange("m_fitRange",minMass, maxMass);
  mass.setRange("m_JPsiRange",2.85, 3.35);
  RooRealVar rap("rap","rap",minRap,maxRap);
  
  double znSelectionLow, znSelectionHigh;
  if (znSelection > -1.5 && znSelection < 3.5)
  {
    znSelectionLow  = znSelection-0.001;
    znSelectionHigh = znSelection+0.001;
  } else if (znSelection > 4.5 && znSelection < 5.5)
  {
    znSelectionLow  = -znSelection-0.001;
    znSelectionHigh =  znSelection+0.001;
  } else
  {
    cout << "Wrong znSelection!" << endl;
    throw;
  }

  RooRealVar znclass("znclass","znclass",znSelectionLow,znSelectionHigh);


  // import the tree
  RooDataSet inData("inData", "inData", RooArgSet(pt, mass, rap, znclass), Import(*dataTree));
  int nEvents = inData.numEntries();

  // build pdfs
  // ---double sided CB
  RooRealVar m0("m0", "m0", 3.1,minMass,maxMass);

  RooRealVar sigmaL("sigmaL", "sigmaL",0.08, 0.04, 0.15);
  if (minPt < 0.25) sigmaL.setVal(getMcMassParam(10,"val","sigmaLParam",0.00,3.00,-4.00,-2.50));
  if (minPt >= 0.25) sigmaL.setVal(getMcMassParam(11,"val","sigmaLParam",0.00,3.00,-4.00,-2.50));
  sigmaL.setConstant();

  RooRealVar sigmaR("sigmaR", "sigmaR",0.08, 0.04, 0.15);
  if (minPt < 0.25) sigmaR.setVal(getMcMassParam(10,"val","sigmaRParam",0.00,3.00,-4.00,-2.50));
  if (minPt >= 0.25) sigmaR.setVal(getMcMassParam(11,"val","sigmaRParam",0.00,3.00,-4.00,-2.50));
  // sigmaR.setConstant();

  RooRealVar nL("nL", "nL", 10);
  if (minPt < 0.25) nL.setVal(getMcMassParam(10,"val","nLParam",0.00,3.00,-4.00,-2.50));
  if (minPt >= 0.25) nL.setVal(getMcMassParam(11,"val","nLParam",0.00,3.00,-4.00,-2.50));
  nL.setConstant();    

  RooRealVar nR("nR", "nR", 10);
  if (minPt < 0.25) nR.setVal(getMcMassParam(10,"val","nRParam",0.00,3.00,-4.00,-2.50));
  if (minPt >= 0.25) nR.setVal(getMcMassParam(11,"val","nRParam",0.00,3.00,-4.00,-2.50));
  nR.setConstant();

  RooRealVar alphaL("alphaL", "alphaL", 1.27, 0.1, 3);
  if (minPt < 0.25) alphaL.setVal(getMcMassParam(10,"val","alphaLParam",0.00,3.00,-4.00,-2.50));
  if (minPt >= 0.25) alphaL.setVal(getMcMassParam(11,"val","alphaLParam",0.00,3.00,-4.00,-2.50));
  alphaL.setConstant();

  RooRealVar alphaR("alphaR", "alphaR", 1.50, 0.1, 3); 
  if (minPt < 0.25) alphaR.setVal(getMcMassParam(10,"val","alphaRParam",0.00,3.00,-4.00,-2.50));
  if (minPt >= 0.25) alphaR.setVal(getMcMassParam(11,"val","alphaRParam",0.00,3.00,-4.00,-2.50));
  alphaR.setConstant();

  RooCrystalBall crystalBall("crystalBall", "crystalBall", mass, m0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);

  // ---background gamma gamma
  // RooRealVar lambda("lambda","exponent",-0.9,-5.,-0.5);
  // RooRealVar a2("a2","parameter a2",0.72,0.01,1);
  // RooRealVar a3("a3","parameter a3",0.99,0.01,1);
  // RooRealVar a4("a4","parameter a4",0.26,0.01,1);

  // Polynom
  // RooGenericPdf Bkgd("Bkgd","(mass>4)? exp(mass*lambda) : exp(mass*lambda)*(1+a2*pow((mass-4),2)+a3*pow((mass-4),3)+a4*pow((mass-4),4))",RooArgSet(mass,lambda,a2,a3,a4));
  // Exp
  // RooExponential Bkgd("Bkgd", "Bkgd", mass, lambda);

  // Linear
  RooRealVar bkg_c1("bkg_c1", "slope of background", -0.21, -2.0, 0.0);
  bkg_c1.setConstant();
  RooPolynomial Bkgd("Bkgd", "linear function for background", mass, RooArgList(bkg_c1));
  // Flat
  // RooPolynomial Bkgd("Bkgd", "flat function", mass, RooArgList());

  // ---combine the PDFs
  RooRealVar N_CB("N_{J/Psi}","Number of CB events",0.8*nEvents,0.05*nEvents,nEvents);
  RooRealVar N_BG("N_{bg}","Number of BG events",0.1*nEvents,0,nEvents);

  RooAddPdf fitFunc("fitFunc","DS CB and Background PDF", RooArgList(crystalBall,Bkgd), RooArgList(N_CB,N_BG));

   // fit
  RooFitResult* r = fitFunc.fitTo(inData,Extended(kTRUE),Range("m_fitRange"),Save());

  // Compute number of events from the mass fit
  // ---Define variabels
  double N_JPsi[2];
  double N_bkgd_JPsi_mass_range[2];
  // ---Create integrals
  RooAbsReal *I_bkgd_JPsi_mass_range = Bkgd.createIntegral(mass,NormSet(mass),Range("m_JPsiRange"));
  // ---Compute the number of events
  N_JPsi[0] = N_CB.getVal();
  N_JPsi[1] =  N_CB.getError();
  N_bkgd_JPsi_mass_range[0] =  I_bkgd_JPsi_mass_range->getVal()*N_BG.getVal();
  N_bkgd_JPsi_mass_range[1] =  I_bkgd_JPsi_mass_range->getVal()*N_BG.getError();  

  // ---Save the J/Psi yields to a file
  TVector nCbYield(2);
  nCbYield[0] = N_JPsi[0];
  nCbYield[1] = N_JPsi[1];

  TVector nBkgdYield(2);
  nBkgdYield[0] = N_bkgd_JPsi_mass_range[0];
  nBkgdYield[1] = N_bkgd_JPsi_mass_range[1];

  // write the histo to the output file
  fOut->cd();
  nCbYield.Write(Form("nCbYield_%i",bin));
  nBkgdYield.Write(Form("nBkgdYield_%i",bin));

  //-----------------------------------------------------------
  // Draw mass histogram
  // -------------------------------------------------------------------------------- 
  
  // Set general options
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTextSize(0.04);
  
  //-----------------------------------------------------------
  // plot mass distributions

  TCanvas *c = new TCanvas(Form("pt in (%.2f,%.2f), ZN class = %d",minPt,maxPt,znSelection),
			   Form("pt in (%.2f,%.2f), ZN class = %d",minPt,maxPt,znSelection),
			   1600, 800);
	c->Divide(2);		   
  
  //----------------------------
  //Draw correlation matrix
  c->cd(2);
  TPad *c2 = new TPad("c2","c2",0.001,0.001,0.999,0.999);
  // --- pad option
  c2->SetRightMargin(0.15);
  c2->SetLeftMargin(0.15);
  c2->Draw();
  c2->cd();
  // --- TH2D correlation matrix
  TH2* h_CorrM = r->correlationHist();
  h_CorrM->SetMarkerSize(1.2);
  h_CorrM->Draw("zcol,text");

  //------------------------------
  // Draw mass fit
  c->cd(1);
  TPad *c1 = new TPad("c1","c1",0.001,0.001,0.999,0.999);
   // --- pad option
  c1->SetRightMargin(0.03);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.12);
  c1->Draw();
  c1->cd();

  // plot mass distributions
  RooPlot *frame = mass.frame(Title(Form("pt in (%.2f,%.2f)",minPt,maxPt)));
  inData.plotOn(frame,Binning(nBins),MarkerStyle(20),MarkerSize(0.5));
  fitFunc.plotOn(frame, LineColor(kBlack));
  fitFunc.plotOn(frame,Name("crystalBall"), Components(crystalBall), LineColor(kRed+1));
  fitFunc.plotOn(frame,Name("Bkgd"), Components(Bkgd), LineColor(kBlue+1));
  Bkgd.paramOn(frame,Layout(0.65, 0.9, 0.65));
  frame->getAttText()->SetTextSize(0.02);
  frame->GetXaxis()->SetTitle("#it{m_{#mu^{+}#mu^{-}}} (GeV/#it{c}^{2})");
  frame->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c^{2}} ", (maxMass-minMass)/nBins*1000));
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->Draw();

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.038);
  if (znSelection == 0) leg->AddEntry((TObject*)0, "0N0N", "");
  if (znSelection == 1) leg->AddEntry((TObject*)0, "0NXN", "");
  if (znSelection == 2) leg->AddEntry((TObject*)0, "XN0N", "");
  if (znSelection == 3) leg->AddEntry((TObject*)0, "XNXN", "");
  if (znSelection == 5) leg->AddEntry((TObject*)0, "any neutron", "");
  leg->AddEntry((TObject*)0,Form("%3.2f < p_{T} < %3.2f",minPt,maxPt),"");
  leg->AddEntry((TObject*)0,Form("N_{J/#psi} = %3.0f #pm %3.0f",N_JPsi[0],N_JPsi[1]),"");
  leg->AddEntry((TObject*)0,Form("N_{bg(2.85,3.35)} = %3.0f #pm %3.0f", N_bkgd_JPsi_mass_range[0],N_bkgd_JPsi_mass_range[1]),"");
  leg->Draw();
  //---Save as pdf
  TString filepath = Form("massPlots/Selection_%i/GammaMass_fit_%.2f_y_%.2f_%.2f_pt_%.2f_ZNclass_%d.png", iSel, minRap,maxRap,minPt,maxPt,znSelection);
  gSystem->mkdir(filepath, kTRUE);
  c1->SaveAs(filepath.Data());
  cout << "chi^2 = " << frame->chiSquare() << endl;
  cout << " Entries " << inData.numEntries() << endl;
}

void yieldRooFit(int iData = 0, int iSel = 1, int znSelection = 5)
{
  // set mass ranges for the fit
  float mMin = 2.55; // 2.85;
  float mMax = 3.45; // 3.35;
  if (iData == 12) { // psi2s
    mMin = 3.3;
    mMax = 4;
  }
  // get the input tree
  TTree *dataTree = getTree(iData, iSel);
  // define kinematic ranges for data to fit 
  float minRap[3] = {-4.0, -4.00, -3.25 };
  float maxRap[3] = {-2.5, -3.25, -2.50 };  
  // float minRap[4] = {-4.0, -4.00, -3.50, -3.00};
  // float maxRap[4] = {-2.5, -3.50, -3.00, -2.50};
  
  // Set pt binning
  setPtBinning(znSelection);

  // ---------
  // make a fit
  // doOneSimpleDscbFit(dataTree, iSel, znSelection, 0, 3, 2,5, minRap[0], maxRap[0]);
  // ---------

  // Create the output file for the yields
  TString fullpath = Form("zdcClassYields/Selection_%i/Yields_ZNclass%d.root",iSel,znSelection);
  gSystem->mkdir(fullpath, kTRUE);
  TFile *fOut = new TFile(fullpath,"recreate");

  // Do the fits
  for (auto bin = 0; bin < nPtBins; bin++) {
    doOneFit(dataTree, iSel, znSelection, ptBinBoundariesVec[bin], ptBinBoundariesVec[bin+1], mMin, mMax, minRap[0], maxRap[0], fOut, bin);
  }
  fOut->Close();
}

