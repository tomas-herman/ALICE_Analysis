
//
// makes a fit to the pt distribution to data from produceFitTree.C  
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
#include "/home/tomas/cernbox/work/ALICE_Analysis/old_coherent/LoadEff.C"

using namespace RooFit;

//-----------------------------------------------
// Get the histo with pt template
TH1 *getPtTemplateHisto(int iData, int znSelection, double minRap, double maxRap)
{
  TFile *inFile = nullptr;
  TH1 *h = nullptr;

  if (znSelection == -1){
    inFile = TFile::Open(Form("ptHistos/ptHisto_data_%i_%.2f_%.2f.root", iData, abs(minRap), abs(maxRap)));
    h = inFile->Get<TH1>("h");
  } else{
    inFile = TFile::Open(Form("ptHistos/ptHisto_gamma_%i.root",znSelection));
    h = inFile->Get<TH1>(Form("h%i",znSelection));
  }
  // return the histo
  return ((TH1 *) h);
}

//-----------------------------------------------
// compute fd value on the fly or make a function to compute, save and load them 
double computeFd(int ptCut, double minRap, double maxRap) 
{
  // Load efficiencies
  double effCohJpsiToMuPtCut;
  double effCohJpsiToMuPtAll;
  double effCohPsi2sToMuPtCut;
  double effCohPsi2sToMuPtAll;
  double effCohPsi2sToMuPiPtCut;
  double effCohPsi2sToMuPiPtAll;

  LoadEff(minRap, maxRap, 
          effCohJpsiToMuPtCut, effCohPsi2sToMuPtCut, effCohPsi2sToMuPiPtCut,
          effCohJpsiToMuPtAll, effCohPsi2sToMuPtAll, effCohPsi2sToMuPiPtAll, 
          "ADvetoOff", "CMUP6");

  // Psi' to J/Psi cross section ratio values taken from https://alice-publications.web.cern.ch/system/files/draft/5085/2019-08-02-paper_v10.pdf
  double ratio = 0.150;
  double ratio_err = 0.0285;

  // --------------------------------------------------------------------------------
  // Compute the feed down contribution
  double fdPt;
  double fdPtErr;

  if (ptCut == 0) {
    fdPt = ratio*effCohPsi2sToMuPiPtAll*0.61400/effCohJpsiToMuPtAll;
    fdPtErr = ratio/ratio_err*fdPt;
  }
  if (ptCut == 1) {
    fdPt = ratio*effCohPsi2sToMuPiPtCut*0.61400/effCohJpsiToMuPtCut;
    fdPtErr = ratio/ratio_err*fdPt;
  }
    return fdPt;
}


// -----------------------------------------------------------------
// do one mass fit for MC
void doOnePtFit(TTree *dataTree, int iSel, int znSelection,
                        		     double minPt, double maxPt,
                        		     double minMass, double maxMass,
                        		     double minRap, double maxRap)
{
  // define number of bins
  auto nBins = 120;
  // define variables for the tree
  RooRealVar pt("pt","pt",minPt,maxPt);
  pt.setRange("pt_fitRange",minPt, maxPt);
  RooRealVar mass("mass","mass",minMass,maxMass);
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
  // ---Create histograms
  RooDataHist dataHistCohJpsiToMu("dataHistCohJpsiToMu","dataHistCohJpsiToMu",pt,getPtTemplateHisto(10, -1, minRap, maxRap),1);
  RooDataHist dataHistIncohJpsiToMu("dataHistIncohJpsiToMu","dataHistIncohJpsiToMu",pt,getPtTemplateHisto(11, -1, minRap, maxRap),1);
  RooDataHist dataHistCohPsi2sToMuPi("dataHistCohPsi2sToMuPi","dataHistCohPsi2sToMuPi",pt,getPtTemplateHisto(15, -1, minRap, maxRap),1);
  RooDataHist dataHistIncohPsi2sToMuPi("dataHistIncohPsi2sToMuPi","dataHistIncohPsi2sToMuPi",pt,getPtTemplateHisto(16, -1, minRap, maxRap),1);
  // RooDataHist dataHistTwoGammaToMuMedium("dataHistTwoGammaToMuMedium","dataHistTwoGammaToMuMedium",pt,getPtTemplateHisto(13, -1, minRap, maxRap),1); // STARlight PDF
  RooDataHist dataHistTwoGammaToMuMedium("dataHistTwoGammaToMuMedium","dataHistTwoGammaToMuMedium",pt,getPtTemplateHisto(-1, znSelection, -1, -1),1); // Data driven PDF

  // return the PDF
  RooHistPdf pdfCohJpsiToMu("pdfCohJpsiToMu", "pdfCohJpsiToMu",pt,dataHistCohJpsiToMu,0);
  RooHistPdf pdfIncohJpsiToMu("pdfIncohJpsiToMu", "pdfIncohJpsiToMu",pt,dataHistIncohJpsiToMu,0);
  RooHistPdf pdfCohPsi2sToMuPi("pdfCohPsi2sToMuPi", "pdfCohPsi2sToMuPi",pt,dataHistCohPsi2sToMuPi,0);
  RooHistPdf pdfIncohPsi2sToMuPi("pdfIncohPsi2sToMuPi", "pdfIncohPsi2sToMuPi",pt,dataHistIncohPsi2sToMuPi,0);
  RooHistPdf pdfTwoGammaToMuMedium("pdfTwoGammaToMuMedium", "pdfTwoGammaToMuMedium",pt,dataHistTwoGammaToMuMedium,0);

  // ---Create Incoherent Disocitation PDF
  // Values taken from https://arxiv.org/abs/1304.5162
  RooRealVar b("b","b",1.79, 0.1, 6); //1.67, 1.91); 
  RooRealVar n("n","n",3.58, 1.5,8); //3.43,3.73); 

  // b.setConstant(kTRUE);
  n.setConstant(kTRUE);

  RooGenericPdf pdfIncohJpsiToX("pdfIncohJpsiToX","pt*pow((1+pow(pt,2)*b/n),-n)",RooArgSet(pt, b, n)); 

  // Get fD value
  double fdPtAll = computeFd(0, minRap, maxRap);
  RooRealVar rooFdPtAll("fdPtAll","fd for all pt range",fdPtAll);
  rooFdPtAll.setConstant(kTRUE);


  //Set the normalisation for the pt fit templates
  // ---JPsi
  RooRealVar nCohJpsiToMu("J/#psi_{coh}","number of coherent jpsi events",0.75*nEvents,0,nEvents);   
  RooRealVar nIncohJpsiToMu("J/#psi_{incoh}","number of incoherent jpsi events",0.1*nEvents,0,nEvents); 
  // ---Psi' feed down fixed by JPsi values * feed down coeficient
  RooFormulaVar nCohPsi2sToMuPi("psi'#to#pi_{coh}","@0*@1",RooArgList(nCohJpsiToMu,rooFdPtAll));
  RooFormulaVar nIncohPsi2sToMuPi("psi'#to#pi_{incoh}","@0*@1",RooArgList(nIncohJpsiToMu,rooFdPtAll));
  // ---Gamma gamma is fixed to background in the mass fit
  double nBkgd = getMassFitResults(iSel,znSelection,"val","nBkgdParam",minPt,maxPt,minRap,maxRap);
  RooRealVar nTwoGammaToMuMedium("#gamma#gamma","number of gg", nBkgd,0,nEvents);
  nTwoGammaToMuMedium.setConstant(kTRUE);
  // ---Disociative
  RooRealVar nIncohJpsiToX("N_{disoc}","number of incoherent jpsi disociated signal",0.1*nEvents,0,nEvents);

  // Create the model as the sum of the templates
  // --- For J/Psi
  RooAddPdf fitFunc("fitFunc","fitFunc",RooArgList(pdfCohJpsiToMu,pdfIncohJpsiToMu,pdfCohPsi2sToMuPi,pdfIncohPsi2sToMuPi,pdfTwoGammaToMuMedium,pdfIncohJpsiToX),
                                        RooArgList(nCohJpsiToMu,nIncohJpsiToMu,nCohPsi2sToMuPi,nIncohPsi2sToMuPi,nTwoGammaToMuMedium,nIncohJpsiToX));
  
  // fit
  RooFitResult* r = fitFunc.fitTo(inData,Extended(kTRUE),Range("pt_fitRange"),Save());

  // --------------------------------------------------------------------------------
  // Compute cohernt fraction in pt bins
  pt.setRange("inCoherentPtRange0",0.0, 3.0);
  pt.setRange("inCoherentPtRange1",0.3, 1.5);
  pt.setRange("inCoherentPtRange2",0.3, 0.5);
  pt.setRange("inCoherentPtRange3",0.5, 0.7);
  pt.setRange("inCoherentPtRange4",0.7, 0.9);
  pt.setRange("inCoherentPtRange5",0.9, 1.2);
  pt.setRange("inCoherentPtRange6",1.2, 1.5);
  // ---Create integrals for pt bins
  RooAbsReal *iCohJpsiToMu0 = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange0"));
  RooAbsReal *iIncohJpsiToMu0 = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange0"));
  RooAbsReal *iCohPsi2sToMuPi0 = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange0"));
  RooAbsReal *iIncohPsi2sToMuPi0 = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange0"));
  RooAbsReal *iTwoGammaToMuMedium0 = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange0"));
  RooAbsReal *iIncohJpsiToX0 = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange0"));

  RooAbsReal *iCohJpsiToMu1 = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange1"));
  RooAbsReal *iIncohJpsiToMu1 = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange1"));
  RooAbsReal *iCohPsi2sToMuPi1 = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange1"));
  RooAbsReal *iIncohPsi2sToMuPi1 = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange1"));
  RooAbsReal *iTwoGammaToMuMedium1 = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange1"));
  RooAbsReal *iIncohJpsiToX1 = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange1"));

  RooAbsReal *iCohJpsiToMu2 = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange2"));
  RooAbsReal *iIncohJpsiToMu2 = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange2"));
  RooAbsReal *iCohPsi2sToMuPi2 = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange2"));
  RooAbsReal *iIncohPsi2sToMuPi2 = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange2"));
  RooAbsReal *iTwoGammaToMuMedium2 = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange2"));
  RooAbsReal *iIncohJpsiToX2 = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange2"));

  RooAbsReal *iCohJpsiToMu3 = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange3"));
  RooAbsReal *iIncohJpsiToMu3 = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange3"));
  RooAbsReal *iCohPsi2sToMuPi3 = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange3"));
  RooAbsReal *iIncohPsi2sToMuPi3 = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange3"));
  RooAbsReal *iTwoGammaToMuMedium3 = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange3"));
  RooAbsReal *iIncohJpsiToX3 = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange3"));

  RooAbsReal *iCohJpsiToMu4 = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange4"));
  RooAbsReal *iIncohJpsiToMu4 = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange4"));
  RooAbsReal *iCohPsi2sToMuPi4 = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange4"));
  RooAbsReal *iIncohPsi2sToMuPi4 = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange4"));
  RooAbsReal *iTwoGammaToMuMedium4 = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange4"));
  RooAbsReal *iIncohJpsiToX4 = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange4"));

  RooAbsReal *iCohJpsiToMu5 = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange5"));
  RooAbsReal *iIncohJpsiToMu5 = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange5"));
  RooAbsReal *iCohPsi2sToMuPi5 = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange5"));
  RooAbsReal *iIncohPsi2sToMuPi5 = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange5"));
  RooAbsReal *iTwoGammaToMuMedium5 = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange5"));
  RooAbsReal *iIncohJpsiToX5 = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange5"));

  RooAbsReal *iCohJpsiToMu6 = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange6"));
  RooAbsReal *iIncohJpsiToMu6 = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange6"));
  RooAbsReal *iCohPsi2sToMuPi6 = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange6"));
  RooAbsReal *iIncohPsi2sToMuPi6 = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange6"));
  RooAbsReal *iTwoGammaToMuMedium6 = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange6"));
  RooAbsReal *iIncohJpsiToX6 = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("inCoherentPtRange6"));

  RooFormulaVar fc0("fc0","(@2*@5)/(@0*@3+@1*@4)",RooArgList(nIncohJpsiToMu,nIncohJpsiToX,nCohJpsiToMu,*iIncohJpsiToMu0,*iIncohJpsiToX0,*iCohJpsiToMu0));
  RooFormulaVar fc1("fc1","(@2*@5)/(@0*@3+@1*@4)",RooArgList(nIncohJpsiToMu,nIncohJpsiToX,nCohJpsiToMu,*iIncohJpsiToMu1,*iIncohJpsiToX1,*iCohJpsiToMu1));
  RooFormulaVar fc2("fc2","(@2*@5)/(@0*@3+@1*@4)",RooArgList(nIncohJpsiToMu,nIncohJpsiToX,nCohJpsiToMu,*iIncohJpsiToMu2,*iIncohJpsiToX2,*iCohJpsiToMu2));
  RooFormulaVar fc3("fc3","(@2*@5)/(@0*@3+@1*@4)",RooArgList(nIncohJpsiToMu,nIncohJpsiToX,nCohJpsiToMu,*iIncohJpsiToMu3,*iIncohJpsiToX3,*iCohJpsiToMu3));
  RooFormulaVar fc4("fc4","(@2*@5)/(@0*@3+@1*@4)",RooArgList(nIncohJpsiToMu,nIncohJpsiToX,nCohJpsiToMu,*iIncohJpsiToMu4,*iIncohJpsiToX4,*iCohJpsiToMu4));
  RooFormulaVar fc5("fc5","(@2*@5)/(@0*@3+@1*@4)",RooArgList(nIncohJpsiToMu,nIncohJpsiToX,nCohJpsiToMu,*iIncohJpsiToMu5,*iIncohJpsiToX5,*iCohJpsiToMu5));
  RooFormulaVar fc6("fc6","(@2*@5)/(@0*@3+@1*@4)",RooArgList(nIncohJpsiToMu,nIncohJpsiToX,nCohJpsiToMu,*iIncohJpsiToMu6,*iIncohJpsiToX6,*iCohJpsiToMu6));  

  //-----------------------------------------------------------
  // Draw pt histogram
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
  // Draw pt fit
  c->cd(1);
  TPad *c1 = new TPad("c1","c1",0.001,0.001,0.999,0.999);
   // --- pad option
  c1->SetRightMargin(0.03);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.12);
  c1->Draw();
  c1->cd();
  c1->SetLogy();

  // plot mass distributions
  RooPlot *frame = pt.frame(Title(Form("pt in (%.2f,%.2f)",minPt,maxPt)));
  inData.plotOn(frame,Binning(nBins),MarkerStyle(20),MarkerSize(0.5));
  fitFunc.plotOn(frame, LineColor(kBlack));
  fitFunc.plotOn(frame,Name("pdfCohJpsiToMu"), Components(pdfCohJpsiToMu), LineColor(kBlue+1), LineWidth(2));
  fitFunc.plotOn(frame,Name("pdfIncohJpsiToMu"),Components(pdfIncohJpsiToMu), LineColor(kRed+1), LineWidth(2));
  fitFunc.plotOn(frame,Name("pdfCohPsi2sToMuPi"),Components(pdfCohPsi2sToMuPi), LineColor(kCyan+1), LineWidth(2));
  fitFunc.plotOn(frame,Name("pdfIncohPsi2sToMuPi"),Components(pdfIncohPsi2sToMuPi), LineColor(kOrange), LineWidth(2));
  fitFunc.plotOn(frame,Name("pdfTwoGammaToMuMedium"),Components(pdfTwoGammaToMuMedium), LineColor(kGreen+2), LineWidth(2));
  fitFunc.plotOn(frame,Name("pdfIncohJpsiToX"),Components(pdfIncohJpsiToX), LineColor(kMagenta+1), LineWidth(2));
  frame->GetXaxis()->SetTitle("#it{p_{T}_{#mu^{+}#mu^{-}}} (GeV/#it{c})");
  frame->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c} ", (maxPt-minPt)/nBins*1000));
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetAxisRange(0.025,10*nEvents,"Y");
  frame->Draw();

  TLegend *leg2 = new TLegend(0.10,0.70,0.55,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.036);
  if (znSelection == 0) leg2->AddEntry((TObject*)0, "0N0N", "");
  if (znSelection == 1) leg2->AddEntry((TObject*)0, "0NXN", "");
  if (znSelection == 2) leg2->AddEntry((TObject*)0, "XN0N", "");
  if (znSelection == 3) leg2->AddEntry((TObject*)0, "XNXN", "");
  if (znSelection == 5) leg2->AddEntry((TObject*)0, "any neutron", "");
  leg2->AddEntry((TObject*)0,Form("%3.2f < m_{#mu^{+}#mu^{-}} < %3.2f",minMass,maxMass),"");
  leg2->AddEntry((TObject*)0,Form("%3.2f < y < %3.2f",minRap,maxRap),"");
  leg2->Draw();

  TLegend *leg = new TLegend(0.50,0.52,0.93,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.036);
  leg->AddEntry((TObject*)0,Form("f^{0.0,3.0}_{C} = %3.2f #pm %3.2f %%",fc0.getVal()*100, fc0.getPropagatedError(*r)*100),"");
  leg->AddEntry((TObject*)0,Form("f^{0.3,1.5}_{C} = %3.2f #pm %3.2f %%",fc1.getVal()*100, fc1.getPropagatedError(*r)*100),"");
  leg->AddEntry((TObject*)0,Form("f^{0.3,0.5}_{C} = %3.2f #pm %3.2f %%",fc2.getVal()*100, fc2.getPropagatedError(*r)*100),"");
  leg->AddEntry((TObject*)0,Form("f^{0.5,0.7}_{C} = %3.2f #pm %3.2f %%",fc3.getVal()*100, fc3.getPropagatedError(*r)*100),"");
  leg->AddEntry((TObject*)0,Form("f^{0.7,0.9}_{C} = %3.2f #pm %3.2f %%",fc4.getVal()*100, fc4.getPropagatedError(*r)*100),"");
  leg->AddEntry((TObject*)0,Form("f^{0.9,1.2}_{C} = %3.2f #pm %3.2f %%",fc5.getVal()*100, fc5.getPropagatedError(*r)*100),"");
  leg->AddEntry((TObject*)0,Form("f^{1.2,1.5}_{C} = %3.2f #pm %3.2f %%",fc6.getVal()*100, fc6.getPropagatedError(*r)*100),"");
  // leg->AddEntry("pdfCohJpsiToMu","Coherent J/#psi", "L");
  // leg->AddEntry("pdfIncohJpsiToMu","Incoherent J/#psi", "L");
  // leg->AddEntry("pdfCohPsi2sToMuPi","Coh. #psi' #rightarrow  J/#psi", "L");
  // leg->AddEntry("pdfIncohPsi2sToMuPi","Incoh. #psi' #rightarrow  J/#psi", "L");
  // leg->AddEntry("pdfTwoGammaToMuMedium","#gamma#gamma #rightarrow #mu#mu", "L");
  // leg->AddEntry("pdfIncohJpsiToX","Nucleon disoc.", "L");
  leg->Draw();

  c1->Draw();
  //---Save as pdf
  char ptFit[120];
  sprintf(ptFit,"ptPlots/Selection_%i/pt_fit_%.2f_y_%.2f_%.2f_mass_%.2f_ZNclass_%d.png", iSel, minRap,maxRap,minMass,maxMass,znSelection);
  c1->SaveAs(ptFit);
  cout << "chi^2 = " << frame->chiSquare() << endl;
  cout << " Entries " << inData.numEntries() << endl;

  // set up the output file to save the fit results
  TFile *fOut = new TFile(Form("ptFitResults/Selection_%i/Data_ZNclass%d_%.2f_y_%.2f_%.2f_mass_%.2f.root",iSel,znSelection,abs(minRap),abs(maxRap),minMass,maxMass),"recreate");

  TVector fcParam0(2), fcParam1(2), fcParam2(2), fcParam3(2), fcParam4(2), fcParam5(2), fcParam6(2);
  fcParam0[0] = fc0.getVal();
  fcParam0[1] = fc0.getPropagatedError(*r);
  fcParam1[0] = fc1.getVal();
  fcParam1[1] = fc1.getPropagatedError(*r);
  fcParam2[0] = fc2.getVal();
  fcParam2[1] = fc2.getPropagatedError(*r);
  fcParam3[0] = fc3.getVal();
  fcParam3[1] = fc3.getPropagatedError(*r);
  fcParam4[0] = fc4.getVal();
  fcParam4[1] = fc4.getPropagatedError(*r);
  fcParam5[0] = fc5.getVal();
  fcParam5[1] = fc5.getPropagatedError(*r);
  fcParam6[0] = fc6.getVal();
  fcParam6[1] = fc6.getPropagatedError(*r);

  // write the histo to the output file
  fOut->cd();
  fcParam0.Write("fcParam0");
  fcParam1.Write("fcParam1");
  fcParam2.Write("fcParam2");
  fcParam3.Write("fcParam3");
  fcParam4.Write("fcParam4");
  fcParam5.Write("fcParam5");
  fcParam6.Write("fcParam6");
  fOut->Close();

}

void ptRooFit(int iData = 0, int iSel = 1)
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
  float minRap[3] = {-4.0, -4.00, -3.25 };
  float maxRap[3] = {-2.5, -3.25, -2.50 };
  // float minRap[4] = {-4.0, -4.00, -3.50, -3.00};
  // float maxRap[4] = {-2.5, -3.50, -3.00, -2.50};
  
  auto znSelections = {5};
  // ---------
  // make a fit
  for ( auto& znSelection : znSelections ) {
    for (int i = 0; i < 3; i++) {
      doOnePtFit(dataTree, iSel, znSelection, 0., 3., mMin, mMax, minRap[i], maxRap[i]); 
    }
  }
}

