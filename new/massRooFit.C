
//
// makes a fit to the invariant mass distribution to data from produceFitTree.C  
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


//-----------------------------------------------
// Get a MC parameter for mass fit
double getMcMassParam(int iData,TString valOrErr,TString paramName,double minPt,double maxPt,double minRap,double maxRap)
{
  TFile *inFile = TFile::Open(Form("massMcParameters/mcParam_data_%i_%.2f_y_%.2f_%.2f_pt_%.2f.root",iData,abs(minRap),abs(maxRap),minPt,maxPt));

  TVector *ptrParam;
  if (paramName.Contains("m0Param")) { ptrParam = (TVector*)inFile->Get("m0Param"); }
  else if (paramName.Contains("sigmaLParam")) { ptrParam = (TVector*)inFile->Get("sigmaLParam"); }
  else if (paramName.Contains("sigmaRParam")) { ptrParam = (TVector*)inFile->Get("sigmaRParam"); }
  else if (paramName.Contains("nLParam")) { ptrParam = (TVector*)inFile->Get("nLParam"); }
  else if (paramName.Contains("nRParam")) { ptrParam = (TVector*)inFile->Get("nRParam"); }
  else if (paramName.Contains("alphaLParam")) { ptrParam = (TVector*)inFile->Get("alphaLParam"); }
  else if (paramName.Contains("alphaRParam")) { ptrParam = (TVector*)inFile->Get("alphaRParam"); }

  else if (paramName.Contains("lambdaParam")) { ptrParam = (TVector*)inFile->Get("lambdaParam"); }
  else if (paramName.Contains("a2Param")) { ptrParam = (TVector*)inFile->Get("a2Param"); }
  else if (paramName.Contains("a3Param")) { ptrParam = (TVector*)inFile->Get("a3Param"); }
  else if (paramName.Contains("a4Param")) { ptrParam = (TVector*)inFile->Get("a4Param"); }

  else {
    cout << "Parameter " << paramName.Data() << " not found! Exiting." << endl;
    throw;
  }
  
  TVector Param = *ptrParam;
  double par;
  // return the histo
  if (valOrErr.Contains("val")){ par = Param[0]; }
  else if (valOrErr.Contains("err")){ par = Param[1]; }
  else {
    cout << "Select val or err. " << valOrErr.Data() << " is not known! Exiting." << endl;
    throw;
  }
  return (par);  
}


// -----------------------------------------------------------------
// do one dscb fit for MC
void doOneDscbFit(TTree *dataTree, int iData, int iSel, int znSelection,
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
  RooRealVar nL("nL", "nL", 30,30,30);
  nL.setConstant();    
  RooRealVar nR("nR", "nR", 30,30,30);
  nR.setConstant();    
  RooRealVar alphaL("alphaL", "alphaL", 1.27, 0.1, 3);
  RooRealVar alphaR("alphaR", "alphaR", 1.50, 0.1, 3);  
  RooCrystalBall crystalBall("crystalBall", "crystalBall", mass, m0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);

  // ---combine the PDFs
  RooRealVar N_CB("N_{CB}","Number of CB events",1.00*nEvents,0.8*nEvents,1.2*nEvents);

  RooAddPdf fitFunc("fitFunc","DS CB PDF", RooArgList(crystalBall), RooArgList(N_CB));

   // fit
  RooFitResult* r = fitFunc.fitTo(inData,Extended(kTRUE),Range("m_fitRange"),Save());

  // Compute number of events from the mass fit
  // ---Define variabels
  double nCB[2];
  // ---Compute the number of events
  nCB[0] = N_CB.getVal();
  nCB[1] =  N_CB.getError();

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
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetAxisRange(0.025,nEvents/4,"Y");
  frame->Draw();

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.038);
  leg->AddEntry((TObject*)0,Form("%3.2f < p_{T} < %3.2f",minPt,maxPt),"");
  leg->AddEntry((TObject*)0,Form("%3.2f < y < %3.2f",minRap,maxRap),"");
  leg->AddEntry((TObject*)0,Form("N_{CB} = %3.0f #pm %3.0f",nCB[0],nCB[1]),"");
  leg->Draw();

  // set up the output file to save the fit result
  TFile *fOut = new TFile(Form("massMcParameters/mcParam_data_%i_%.2f_y_%.2f_%.2f_pt_%.2f.root",iData,abs(minRap),abs(maxRap),minPt,maxPt),"recreate");

  TVector m0Param(2);
  m0Param[0] = m0.getVal();
  m0Param[1] = m0.getError();

  TVector sigmaLParam(2);
  sigmaLParam[0] = sigmaL.getVal();
  sigmaLParam[1] = sigmaL.getError();

  TVector sigmaRParam(2);
  sigmaRParam[0] = sigmaR.getVal();
  sigmaRParam[1] = sigmaR.getError();

  TVector nLParam(2);
  nLParam[0] = nL.getVal();
  nLParam[1] = nL.getError();

  TVector nRParam(2);
  nRParam[0] = nR.getVal();
  nRParam[1] = nR.getError();

  TVector alphaLParam(2);
  alphaLParam[0] = alphaL.getVal();
  alphaLParam[1] = alphaL.getError();

  TVector alphaRParam(2);
  alphaRParam[0] = alphaR.getVal();
  alphaRParam[1] = alphaR.getError();

  // write the histo to the output file
  fOut->cd();
  m0Param.Write("m0Param");
  sigmaLParam.Write("sigmaLParam");
  sigmaRParam.Write("sigmaRParam");
  nLParam.Write("nLParam");
  nRParam.Write("nRParam");
  alphaLParam.Write("alphaLParam");
  alphaRParam.Write("alphaRParam");
  fOut->Close();
}

// -----------------------------------------------------------------
// do one bkgd fit for MC
void doOneBkgdFit(TTree *dataTree, int iData, int iSel, int znSelection,
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
  // ---background gamma gamma
  RooRealVar lambda("#lambda","exponent",-0.9,-5.,-0.5);
  RooRealVar a2("a_{2}","parameter a2",0.72,0.01,1);
  RooRealVar a3("a_{3}","parameter a3",0.99,0.01,1);
  RooRealVar a4("a_{4}","parameter a4",0.26,0.01,1);

  // Polynom
  RooGenericPdf Bkgd("Bkgd","(@0>4)? exp(@0*@1) : exp(@0*@1)*(1+@2*pow((@0-4),2)+@3*pow((@0-4),3)+@4*pow((@0-4),4))",RooArgSet(mass,lambda,a2,a3,a4));
  
  // ---combine the PDFs
  RooRealVar N_BG("N_{bg}","Number of BG events",1.0*nEvents,0.8*nEvents,1.2*nEvents);

  RooAddPdf fitFunc("fitFunc","bkgd PDF", RooArgList(Bkgd), RooArgList(N_BG));

   // fit
  RooFitResult* r = fitFunc.fitTo(inData,Extended(kTRUE),Range("m_fitRange"),Save());

  // Compute number of events from the mass fit
  // ---Define variabels
  double nBkgd[2];
  // ---Compute the number of events
  nBkgd[0] = N_BG.getVal();
  nBkgd[1] = N_BG.getError();

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
  Bkgd.paramOn(frame,Layout(0.65, 0.9, 0.65));
  frame->getAttText()->SetTextSize(0.02);
  frame->GetXaxis()->SetTitle("#it{m_{#mu^{+}#mu^{-}}} (GeV/#it{c}^{2})");
  frame->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c^{2}} ", (maxMass-minMass)/nBins*1000));
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetAxisRange(0.025,nEvents/24,"Y");
  frame->Draw();

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.038);
  leg->AddEntry((TObject*)0,Form("%3.2f < p_{T} < %3.2f",minPt,maxPt),"");
  leg->AddEntry((TObject*)0,Form("%3.2f < y < %3.2f",minRap,maxRap),"");
  leg->AddEntry((TObject*)0,Form("N_{bkgd} = %3.0f #pm %3.0f",nBkgd[0],nBkgd[1]),"");
  leg->Draw();

  // set up the output file to save the fit result
  TFile *fOut = new TFile(Form("massMcParameters/mcParam_data_%i_%.2f_y_%.2f_%.2f_pt_%.2f.root",iData,abs(minRap),abs(maxRap),minPt,maxPt),"recreate");

  TVector lambdaParam(2);
  lambdaParam[0] = lambda.getVal();
  lambdaParam[1] = lambda.getError();

  TVector a2Param(2);
  a2Param[0] = a2.getVal();
  a2Param[1] = a2.getError();

  TVector a3Param(2);
  a3Param[0] = a3.getVal();
  a3Param[1] = a3.getError();

  TVector a4Param(2);
  a4Param[0] = a4.getVal();
  a4Param[1] = a4.getError();

  // write the histo to the output file
  fOut->cd();
  lambdaParam.Write("lambdaParam");
  a2Param.Write("a2Param");
  a3Param.Write("a3Param");
  a4Param.Write("a4Param");
  fOut->Close();
}


// -----------------------------------------------------------------
// do one mass fit
void doOneFit(TTree *dataTree, int iSel, int znSelection,
		    double minPt, double maxPt,
		    double minMass, double maxMass,
		    double minRap, double maxRap)
{
  // define number of bins
  auto nBins = 80;
  // define variables for the tree
  RooRealVar pt("pt","pt",minPt,maxPt);
  RooRealVar mass("mass","mass",2,5);
  mass.setRange("m_fitRange",2.0, 5);
  mass.setRange("m_JPsiRange",2.85, 3.35);
  RooRealVar rap("rap","rap",minRap,maxRap);
  RooRealVar znclass("znclass","znclass",znSelection,znSelection);

  // import the tree
  RooDataSet inData("inData", "inData", RooArgSet(pt, mass, rap, znclass), Import(*dataTree));
  int nEvents = inData.numEntries();

  // build pdfs
  // ---double sided CB
  RooRealVar m0("m0", "m0", 3.1,minMass,maxMass);
  RooFormulaVar m02("m02","@0+3.686097-3.096900",RooArgList(m0));

  RooRealVar sigmaL("sigmaL", "sigmaL",0.08, 0.04, 0.15);
  sigmaL.setVal(getMcMassParam(10,"val","sigmaLParam",minPt,maxPt,minRap,maxRap));
  sigmaL.setConstant();
  RooRealVar sigmaL2("sigmaL2", "sigmaL2",0.08, 0.04, 0.15);
  sigmaL2.setVal(getMcMassParam(12,"val","sigmaLParam",minPt,maxPt,minRap,maxRap));
  sigmaL2.setConstant();

  RooRealVar sigmaR("sigmaR", "sigmaR",0.08, 0.04, 0.15);
  sigmaR.setVal(getMcMassParam(10,"val","sigmaRParam",minPt,maxPt,minRap,maxRap));
  sigmaR.setConstant();
  RooRealVar sigmaR2("sigmaR2", "sigmaR2",0.08, 0.04, 0.15);
  sigmaR2.setVal(getMcMassParam(12,"val","sigmaRParam",minPt,maxPt,minRap,maxRap));
  sigmaR2.setConstant();

  RooRealVar nL("nL", "nL", 30,30,30);
  nL.setConstant();    
  RooRealVar nR("nR", "nR", 30,30,30);
  nR.setConstant();

  RooRealVar alphaL("alphaL", "alphaL", 1.27, 0.1, 3);
  alphaL.setVal(getMcMassParam(10,"val","alphaLParam",minPt,maxPt,minRap,maxRap));
  alphaL.setConstant();
  RooRealVar alphaL2("alphaL2", "alphaL2", 1.27, 0.1, 3);
  alphaL2.setVal(getMcMassParam(12,"val","alphaLParam",minPt,maxPt,minRap,maxRap));
  alphaL2.setConstant();

  RooRealVar alphaR("alphaR", "alphaR", 1.50, 0.1, 3); 
  alphaR.setVal(getMcMassParam(10,"val","alphaRParam",minPt,maxPt,minRap,maxRap));
  alphaR.setConstant();
  RooRealVar alphaR2("alphaR2", "alphaR2", 1.50, 0.1, 3); 
  alphaR2.setVal(getMcMassParam(12,"val","alphaRParam",minPt,maxPt,minRap,maxRap));
  alphaR2.setConstant();

  RooCrystalBall crystalBall("crystalBall", "crystalBall", mass, m0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
  RooCrystalBall crystalBall2("crystalBall2", "crystalBall2", mass, m02, sigmaL2, sigmaR2, alphaL2, nL, alphaR2, nR);

  // ---background gamma gamma
  RooRealVar lambda("#lambda","exponent",-0.9,-5.,-0.5);
  lambda.setVal(getMcMassParam(13,"val","lambdaParam",minPt,maxPt,minRap,maxRap));
  lambda.setConstant();

  RooRealVar a2("a_{2}","parameter a2",0.72,0.01,1);
  a2.setVal(getMcMassParam(13,"val","a2Param",minPt,maxPt,minRap,maxRap));
  a2.setConstant();
   
  RooRealVar a3("a_{3}","parameter a3",0.99,0.01,1);
  a3.setVal(getMcMassParam(13,"val","a3Param",minPt,maxPt,minRap,maxRap));
  a3.setConstant();

  RooRealVar a4("a_{4}","parameter a4",0.26,0.01,1);
  a4.setVal(getMcMassParam(13,"val","a4Param",minPt,maxPt,minRap,maxRap));
  a4.setConstant();

  // Polynom
  RooGenericPdf Bkgd("Bkgd","(@0>4)? exp(@0*@1) : exp(@0*@1)*(1+@2*pow((@0-4),2)+@3*pow((@0-4),3)+@4*pow((@0-4),4))",RooArgSet(mass,lambda,a2,a3,a4));
  // ---combine the PDFs
  RooRealVar N_CB("N_{J/Psi}","Number of CB events",0.8*nEvents,0.05*nEvents,nEvents);
  RooRealVar N_CB2("N_{Psi'}","Number of CB2 events",0.05*nEvents,0,nEvents);
  RooRealVar N_BG("N_{bg}","Number of BG events",0.1*nEvents,0,nEvents);

  RooAddPdf fitFunc("fitFunc","DS CB and Background PDF", RooArgList(crystalBall,crystalBall2,Bkgd), RooArgList(N_CB,N_CB2,N_BG));

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
  fitFunc.plotOn(frame,Name("crystalBall2"), Components(crystalBall2), LineColor(kGreen+2));
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
  leg->AddEntry((TObject*)0,Form("%3.2f < p_{T} < %3.2f",minPt,maxPt),"");
  leg->AddEntry((TObject*)0,Form("%3.2f < y < %3.2f",minRap,maxRap),"");
  leg->AddEntry((TObject*)0,Form("N_{J/#psi} = %3.0f #pm %3.0f",N_JPsi[0],N_JPsi[1]),"");
  leg->AddEntry((TObject*)0,Form("N_{bg(2.85,3.35)} = %3.0f #pm %3.0f", N_bkgd_JPsi_mass_range[0],N_bkgd_JPsi_mass_range[1]),"");
  leg->Draw();
  //---Save as pdf
  char massFit[120];
  sprintf(massFit,"MassPlots/Selection_%i/Mass_fit_%.2f_y_%.2f_%.2f_pt_%.2f_ZNclass_%d.png", iSel, minRap,maxRap,minPt,maxPt,znSelection);
  c1->SaveAs(massFit);
  cout << "chi^2 = " << frame->chiSquare() << endl;
  cout << " Entries " << inData.numEntries() << endl;

  // set up the output file to save the fit results
  TFile *fOut = new TFile(Form("massFitResults/Selection_%i/Data_ZNclass%d_%.2f_y_%.2f_%.2f_pt_%.2f.root",iSel,znSelection,abs(minRap),abs(maxRap),minPt,maxPt),"recreate");

  TVector nCbParam(2);
  nCbParam[0] = N_JPsi[0];
  nCbParam[1] = N_JPsi[1];

  TVector nBkgdParam(2);
  nBkgdParam[0] = N_bkgd_JPsi_mass_range[0];
  nBkgdParam[1] = N_bkgd_JPsi_mass_range[1];

  // write the histo to the output file
  fOut->cd();
  nCbParam.Write("nCbParam");
  nBkgdParam.Write("nBkgdParam");
  fOut->Close();
}

void massRooFit(int iData = 0, int iSel = 1, int znSelection = 0)
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
  
  // ---------
  // make a fit
  for (int i = 0; i < 4; i++) {
    // doOneDscbFit(dataTree, iData, iSel, znSelection, 0, 3, 2,5, minRap[i], maxRap[i]);
    // doOneBkgdFit(dataTree, iData, iSel, znSelection, 0, 3, 2,5, minRap[i], maxRap[i]);
    doOneFit(dataTree, iSel, znSelection, 0, 3, mMin, mMax, minRap[i], maxRap[i]);
  }
}

