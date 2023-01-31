
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
#include "/home/tomas/cernbox/work/ALICE_Analysis/LoadEff.C"

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
  RooRealVar znclass("znclass","znclass",znSelection,znSelection);

  // import the tree
  RooDataSet inData("inData", "inData", RooArgSet(pt, mass, rap, znclass), Import(*dataTree));
  int nEvents = inData.numEntries();

  // build pdfs
  // ---Create histograms
  RooDataHist dataHistCohJpsiToMu("dataHistCohJpsiToMu","dataHistCohJpsiToMu",pt,getPtTemplateHisto(10, -1, minRap, maxRap),1);
  RooDataHist dataHistIncohJpsiToMu("dataHistIncohJpsiToMu","dataHistIncohJpsiToMu",pt,getPtTemplateHisto(11, -1, minRap, maxRap),1);
  RooDataHist dataHistCohPsi2sToMuPi("dataHistCohPsi2sToMuPi","dataHistCohPsi2sToMuPi",pt,getPtTemplateHisto(15, -1, minRap, maxRap),1);
  RooDataHist dataHistIncohPsi2sToMuPi("dataHistIncohPsi2sToMuPi","dataHistIncohPsi2sToMuPi",pt,getPtTemplateHisto(16, -1, minRap, maxRap),1);
  // RooDataHist dataHistTwoGammaToMuMedium("dataHistTwoGammaToMuMedium","dataHistTwoGammaToMuMedium",pt,getPtTemplateHisto(13, -1, minRap, maxRap),1);
  RooDataHist dataHistTwoGammaToMuMedium("dataHistTwoGammaToMuMedium","dataHistTwoGammaToMuMedium",pt,getPtTemplateHisto(-1, znSelection, -1, -1),1);

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
  // Compute incoherent fraction
  pt.setRange("coherentPtRange",0.0,0.25);
  // ---Create integrals for pt < 0.25
  RooAbsReal *iCohJpsiToMu = pdfCohJpsiToMu.createIntegral(pt,NormSet(pt),Range("coherentPtRange"));
  RooAbsReal *iIncohJpsiToMu = pdfIncohJpsiToMu.createIntegral(pt,NormSet(pt),Range("coherentPtRange"));
  RooAbsReal *iCohPsi2sToMuPi = pdfCohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("coherentPtRange"));
  RooAbsReal *iIncohPsi2sToMuPi = pdfIncohPsi2sToMuPi.createIntegral(pt,NormSet(pt),Range("coherentPtRange"));
  RooAbsReal *iTwoGammaToMuMedium = pdfTwoGammaToMuMedium.createIntegral(pt,NormSet(pt),Range("coherentPtRange"));
  RooAbsReal *iIncohJpsiToX = pdfIncohJpsiToX.createIntegral(pt,NormSet(pt),Range("coherentPtRange"));

  RooFormulaVar fi("fi","(@0*@5+@1*@6+@2*@7)/(@3*@8+@4*@9)",RooArgList(nIncohJpsiToMu,nIncohPsi2sToMuPi,nIncohJpsiToX,nCohJpsiToMu,nCohPsi2sToMuPi,*iIncohJpsiToMu,*iIncohPsi2sToMuPi,*iIncohJpsiToX,*iCohJpsiToMu,*iCohPsi2sToMuPi));  

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
  leg2->AddEntry((TObject*)0,Form("%3.2f < m_{#mu^{+}#mu^{-}} < %3.2f",minMass,maxMass),"");
  leg2->AddEntry((TObject*)0,Form("%3.2f < y < %3.2f",minRap,maxRap),"");
  leg2->Draw();

  TLegend *leg = new TLegend(0.55,0.55,0.95,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.036);
  leg->AddEntry((TObject*)0,Form("f_{I} = %3.2f #pm %3.2f %%",fi.getVal()*100, fi.getPropagatedError(*r)*100),"");
  leg->AddEntry("pdfCohJpsiToMu","Coherent J/#psi", "L");
  leg->AddEntry("pdfIncohJpsiToMu","Incoherent J/#psi", "L");
  leg->AddEntry("pdfCohPsi2sToMuPi","Coh. #psi' #rightarrow  J/#psi", "L");
  leg->AddEntry("pdfIncohPsi2sToMuPi","Incoh. #psi' #rightarrow  J/#psi", "L");
  leg->AddEntry("pdfTwoGammaToMuMedium","#gamma#gamma #rightarrow #mu#mu", "L");
  leg->AddEntry("pdfIncohJpsiToX","Nucleon disoc.", "L");
  leg->Draw();

  c1->Draw();
  //---Save as pdf
  char ptFit[120];
  sprintf(ptFit,"ptPlots/Selection_%i/Fixpt_fit_%.2f_y_%.2f_%.2f_mass_%.2f_ZNclass_%d.png", iSel, minRap,maxRap,minMass,maxMass,znSelection);
  c1->SaveAs(ptFit);
  cout << "chi^2 = " << frame->chiSquare() << endl;
  cout << " Entries " << inData.numEntries() << endl;

  // set up the output file to save the fit results
  TFile *fOut = new TFile(Form("ptFitResults/Selection_%i/Data_ZNclass%d_%.2f_y_%.2f_%.2f_mass_%.2f.root",iSel,znSelection,abs(minRap),abs(maxRap),minMass,maxMass),"recreate");

  TVector fiParam(2);
  fiParam[0] = fi.getVal();
  fiParam[1] = fi.getPropagatedError(*r);

  // write the histo to the output file
  fOut->cd();
  fiParam.Write("fiParam");
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
  float minRap[4] = {-4.0, -4.00, -3.50, -3.00};
  float maxRap[4] = {-2.5, -3.50, -3.00, -2.50};
  
  auto znSelections = {0,1,2,3};
  // ---------
  // make a fit
  for ( auto& znSelection : znSelections ) {
    for (int i = 0; i < 4; i++) {
      doOnePtFit(dataTree, iSel, znSelection, 0., 3., mMin, mMax, minRap[i], maxRap[i]); 
    }
  }
}

