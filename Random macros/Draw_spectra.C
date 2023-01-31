// Offline analysis of the forward two muon spectrum
//-----------------------------------------------------------------------------------

////////////////////////////////////////
// Including headers and functions 
////////////////////////////////////////

// c++ headers
#include <iostream>
#include <fstream>
#include <string>

// ROOT headers
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TString.h>

// RooFit headers
#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooArgList.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooBinning.h>
#include <RooGenericPdf.h>
#include <RooHistPdf.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooWorkspace.h>

// Setting namespace
using namespace RooFit;

// My headers
#include "/home/tomas/cernbox/work/ALICE_Analysis/Variables.h"

// My functions
#include "/home/tomas/cernbox/work/ALICE_Analysis/LoadData.C"


///////////////////////////////////////////////////////////
// Draw mass and pt spectra
///////////////////////////////////////////////////////////
void DrawFunction(Double_t y_min, Double_t y_max, Double_t m_min, Double_t m_max, Double_t pt_min, Double_t pt_max, TString ZDC_class, TString ADveto, TString trigger)
{
  // -------------------------------------------------------------------------------- 
  // Set general plotting options
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTextSize(0.04);
  // gStyle->SetPadTickX(1);
  // gStyle->SetPadTickY(1);

  // --------------------------------------------------------------------------------
  // Define binning
  const Double_t n_bins_pt= 108; 
  const Double_t n_bins_m = 80; 

  RooBinning bin_pt(n_bins_pt, pt_range_min, pt_range_max);
  RooBinning bin_m(n_bins_m, m_range_min, m_range_max);

  fMuMuPt.setBinning(bin_pt);
  fMuMuM.setBinning(bin_m);

  // --------------------------------------------------------------------------------
  // Define fit ranges
  fMuMuM.setRange("CB_fitRange",2.2,3.5); // 2.7-3.5
  fMuMuM.setRange("CB2_fitRange",2.5,4.1); // 
  fMuMuM.setRange("m_fitRange",2.1,6);
  // ---Define integral ranges
  fMuMuM.setRange("JPsiMassRange",m_min_jpsi,m_max_jpsi);
  fMuMuM.setRange("Psi(2S)MassRange",m_min_psi2,m_max_psi2);
  fMuMuPt.setRange("JPsiptRange",0,pt_cut);

  // --------------------------------------------------------------------------------
  // Cuts
  char cuts_pt[500];
  char cuts_m[500];

  if(ZDC_class.Contains("all")){ 
    sprintf(cuts_pt,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM<%f && "
                         "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                         "fIsZNAFired > -1 && fIsZNCFired > -1 "
                         ,y_min, y_max, 2.7, m_min, m_max);

    sprintf(cuts_m,"fMuMuY>%f && fMuMuY<%f && fMuMuPt>%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                        "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                        "fIsZNAFired > -1 && fIsZNCFired > -1 "
                        ,y_min, y_max, pt_min, pt_max, 2, 6);


  } else if(ZDC_class.Contains("0N0N")){ 
    sprintf(cuts_pt,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                         "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                         "fIsZNAFired == 0 && fIsZNCFired == 0 "    
                         ,y_min, y_max, 2.7, m_min, m_max);

    sprintf(cuts_m,"fMuMuY>%f && fMuMuY<%f && fMuMuPt>%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                        "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                        "fIsZNAFired == 0 && fIsZNCFired == 0 "  
                        ,y_min, y_max, pt_min, pt_max, 2, 6);

    
  } else if(ZDC_class.Contains("0NXN")){ 
    sprintf(cuts_pt,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                         "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                         "fIsZNAFired == 0 && fIsZNCFired == 1 "    
                         ,y_min, y_max, 2.7, m_min, m_max);

    sprintf(cuts_m,"fMuMuY>%f && fMuMuY<%f && fMuMuPt>%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                        "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                        "fIsZNAFired == 0 && fIsZNCFired == 1 "  
                        ,y_min, y_max, pt_min, pt_max, 2, 6);


  } else if(ZDC_class.Contains("XN0N")){ 
    sprintf(cuts_pt,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                         "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                         "fIsZNAFired == 1 && fIsZNCFired == 0 "    
                         ,y_min, y_max, 2.7, m_min, m_max);

    sprintf(cuts_m,"fMuMuY>%f && fMuMuY<%f && fMuMuPt>%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                        "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                        "fIsZNAFired == 1 && fIsZNCFired == 0 "  
                        ,y_min, y_max, pt_min, pt_max, 2, 6);


  } else if(ZDC_class.Contains("XNXN")){ 
    sprintf(cuts_pt,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                         "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                         "fIsZNAFired == 1 && fIsZNCFired == 1 "    
                         ,y_min, y_max, 2.7, m_min, m_max);

    sprintf(cuts_m,"fMuMuY>%f && fMuMuY<%f && fMuMuPt>%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                        "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                        "fIsZNAFired == 1 && fIsZNCFired == 1 "  
                        ,y_min, y_max, pt_min, pt_max, 2, 6);
      
  // --------------------------------------------------------------------------------
  // Load data
  RooArgSet Arguments(fMuMuM,fMuMuY,fMuMuPt,fRunNum,fADADecision,fADCDecision,fV0ADecision,fV0CDecision,fV0CFiredCells);
  Arguments.add( RooArgSet(fIsZNAFired,fIsZNCFired) ); // You can only add 9 arguments at a time
  // ---Real data sets
  RooDataSet *dataIN18q = new RooDataSet ("dataIN18q", "dataIN18q", Arguments);
  RooDataSet *dataIN18r = new RooDataSet ("dataIN18r", "dataIN18r", Arguments);
  RooDataSet *dataIN15o = new RooDataSet ("dataIN15o", "dataIN15o", Arguments);

  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18q_AnalysisResults.root", "DataRec" , "NanoMUONCMUP6", dataIN18q);
  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18r_AnalysisResults.root", "DataRec" , "NanoMUONCMUP6", dataIN18r);

  // LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18q_AnalysisResults.root", "DataRec", "NanoMUONCMUP11", dataIN18q);
  // LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18r_AnalysisResults.root", "DataRec", "NanoMUONCMUP11", dataIN18r);
  // LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC15o_AnalysisResults.root", "DataRec", "NanoMUONCMUP11", dataIN15o);
 
  // --------------------------------------------------------------------------------
  // Get real data
  // ---Merge data sets
  RooDataSet *dataAllIN = (RooDataSet*) dataIN18q;
  dataAllIN->append(*dataIN18r);
  // dataAllIN->append(*dataIN15o); 
  // ---Do cuts
  RooAbsData* data_m = dataAllIN->reduce(Cut(cuts_m));
  RooAbsData* data_pt = dataAllIN->reduce(Cut(cuts_pt));

  // // --------------------------------------------------------------------------------
  // // Create PDFs
  // // ---Prepare PDFs
  // RooAbsPdf* PDF_CohJpsiToMu; 
  // RooAbsPdf* PDF_CohPsi2sToMu; 
  // RooAbsPdf* PDF_CohPsi2sToMuPi; 
  // RooAbsPdf* PDF_IncohJpsiToMu; 
  // RooAbsPdf* PDF_IncohPsi2sToMu; 
  // RooAbsPdf* PDF_IncohPsi2sToMuPi;   
  // RooAbsPdf* PDF_TwoGammaToMuMedium;

  // char FileName_PDFs[120];
  // sprintf(FileName_PDFs,"/home/tomas/cernbox/work/ALICE_Analysis/Workspaces/PDFs/%s/%s_%.2f_%.2f_PDFs.root",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
  
  // // ---Open file
  // TFile *f_PDFs = new TFile(FileName_PDFs);
  // RooWorkspace* w_PDFs = (RooWorkspace*) f_PDFs->Get("w_PDFs");

  // // ---Load PDFs
  // PDF_CohJpsiToMu = w_PDFs->pdf("PDF_CohJpsiToMu");
  // PDF_CohPsi2sToMu = w_PDFs->pdf("PDF_CohPsi2sToMu");
  // PDF_CohPsi2sToMuPi = w_PDFs->pdf("PDF_CohPsi2sToMuPi");
  // PDF_IncohJpsiToMu = w_PDFs->pdf("PDF_IncohJpsiToMu");
  // PDF_IncohPsi2sToMu = w_PDFs->pdf("PDF_IncohPsi2sToMu");
  // PDF_IncohPsi2sToMuPi = w_PDFs->pdf("PDF_IncohPsi2sToMuPi");
  // PDF_TwoGammaToMuMedium = w_PDFs->pdf("PDF_TwoGammaToMuMedium");

  // // Values taken from https://arxiv.org/abs/1304.5162
  // RooRealVar b("b","b",1.79, 0.1, 6); //1.67, 1.91); 
  // RooRealVar n("n","n",3.58, 1.5,8); //3.43,3.73); 

  // // if(ADveto.Contains("ADvetoOn")) b.setConstant(kTRUE);
  // // n.setConstant(kTRUE);

  // RooGenericPdf PDF_IncohJpsiToX("PDF_IncohJpsiToX","fMuMuPt*pow((1+pow(fMuMuPt,2)*b/n),-n)",RooArgSet(fMuMuPt, b, n)); 
  

  // //Set the normalisation for the pt fit templates
  // // ---JPsi
  // RooRealVar N_CohJpsiToMu("J#psi_{coh}","number of coherent jpsi events",0.1*data_pt->numEntries(),0,data_pt->numEntries());   
  // // RooRealVar N_IncohJpsiToMu("J#psi_{incoh}","number of incoherent jpsi events",0.1*data_pt->numEntries(),0*data_pt->numEntries(),data_pt->numEntries()); 
  // // // Psi'
  // // RooRealVar N_CohPsi2sToMu("#psi'to#mu_{coh}","number of coherent psi2s to mu events",0.75*data_pt->numEntries(),0,data_pt->numEntries());   
  // // RooRealVar N_IncohPsi2sToMu("#psi'to#mu_{incoh}","number of incoherent psi2s to mu events",0.1*data_pt->numEntries(),0,data_pt->numEntries()); 
  // // // ---Psi' feed down fixed by JPsi values * feed down coeficient
  // // RooFormulaVar N_CohPsi2sToMuPi("#psi'to#pi_{coh}","J#psi_{coh}*f_D_PtAll",RooArgList(N_CohJpsiToMu, f_D_PtAll));
  // // RooFormulaVar N_IncohPsi2sToMuPi("#psi'to#pi_{incoh}","J#psi_{incoh}*f_D_PtAll",RooArgList(N_IncohJpsiToMu, f_D_PtAll));
  // // // ---Gamma gamma is fixed to background in the mass fit
  // RooRealVar N_TwoGammaToMuMedium("#gamma#gamma","number of gg",0.1*data_pt->numEntries(),0,data_pt->numEntries());
  // // N_TwoGammaToMuMedium.setConstant(kTRUE);
  // // ---Disociative
  // RooRealVar N_IncohJpsiToX("N_{disoc}","number of incoherent jpsi disociated signal",0.1*data_pt->numEntries(),0,data_pt->numEntries());

  // // Create the model as the sum of the templates
  // // --- For J/Psi
  // RooAddPdf Pt_fit_func("Pt_fit_func","Sum of templates",RooArgList(PDF_IncohJpsiToX),
  //                                                        RooArgList(N_IncohJpsiToX));

  // RooFitResult* fit_pt;
  // fit_pt = Pt_fit_func.fitTo(*data_pt,Extended(kTRUE),Save(),Range(0.25,pt_range_max));


// --------------------------------------------------------------------------------
// Plot the mass fit
  TCanvas *cM = new TCanvas("cM","cM",800,800);
  cM->cd();

  cM->SetLeftMargin(0.15);
  cM->SetRightMargin(0.01);
  cM->SetBottomMargin(0.12);

  RooPlot* frame_m = fMuMuM.frame(Title("Mass fit")) ; 
  data_m->plotOn(frame_m,Name("data_m"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
  // Pt_fit_func.plotOn(frame_pt,Name("Pt_fit_func"),LineColor(kBlack), LineWidth(2)) ;
  // Pt_fit_func.plotOn(frame_pt,Name("PDF_CohJpsiToMu"),Components(*PDF_CohJpsiToMu), LineColor(kBlue+1), LineWidth(2));
  // Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohJpsiToX"),Components(PDF_IncohJpsiToX), LineColor(kMagenta+1), LineWidth(2), Range(0,pt_range_max));
  frame_m->GetXaxis()->SetTitle("#it{m_{#mu^{+}#mu^{-}}} (GeV/#it{c}^{2})");
  frame_m->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c^{2}} ", (m_range_max-m_range_min)/n_bins_m*1000));
  frame_m->GetYaxis()->SetTitleOffset(1.6);
  frame_m->SetAxisRange(0,( data_m->numEntries() )/9,"Y");
  frame_m->Draw();

  TLatex * text_m3 = new TLatex (4.5,0.89*frame_m->GetMaximum(),Form("#bf{%.2f < #it{p}_{T} < %.2f}",pt_min, pt_max));
  text_m3->Draw();
  TLatex * text_m4 = new TLatex (4.5,0.79*frame_m->GetMaximum(),Form("#bf{ZDC: %s}",ZDC_class.Data()));
  text_m4->Draw();

  char name_cM[120];
  sprintf(name_cM,"plots/%s_Mass_fit_%.2f_pt_%.2f.png",ZDC_class.Data(),pt_min,pt_max);
  cM->SaveAs(name_cM);   

// Plot the pt fit
  TCanvas *cPt = new TCanvas("cPt","cPt",800,800);
  cPt->cd();
  cPt->SetLogy();

  cPt->SetLeftMargin(0.14);
  cPt->SetRightMargin(0.01);
  cPt->SetBottomMargin(0.12);

  RooPlot* frame_pt = fMuMuPt.frame(Title("Mass plot")) ; 
  data_pt->plotOn(frame_pt,Name("data_pt"),Binning(bin_pt),MarkerStyle(20),MarkerSize(0.5));
  // Pt_fit_func.plotOn(frame_pt,Name("Pt_fit_func"),LineColor(kBlack), LineWidth(2)) ;
  // Pt_fit_func.plotOn(frame_pt,Name("PDF_CohJpsiToMu"),Components(*PDF_CohJpsiToMu), LineColor(kBlue+1), LineWidth(2));
  // Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohJpsiToX"),Components(PDF_IncohJpsiToX), LineColor(kMagenta+1), LineWidth(2), Range(0,pt_range_max));
  frame_pt->SetAxisRange(0,3,"X");
  frame_pt->SetAxisRange(1,( data_pt->numEntries() )/2.5,"Y");
  frame_pt->GetXaxis()->SetTitle("Dimuon #it{p}_{T} (GeV/#it{c})");
  frame_pt->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c} ", (pt_range_max-pt_range_min)/n_bins_pt*1000));
  frame_pt->GetYaxis()->SetTitleOffset(1.3);
  frame_pt->Draw();

  TLatex * text_pt3 = new TLatex (1.5,0.30*frame_pt->GetMaximum(),Form("#bf{%.2f < #it{m_{#mu^{+}#mu^{-}}} < %.2f}", m_min, m_max));
  text_pt3->Draw();
  TLatex * text_pt4 = new TLatex (1.5,0.11*frame_pt->GetMaximum(),Form("#bf{ZDC: %s}",ZDC_class.Data()));
  text_pt4->Draw();

  char name_cPt[120];
  sprintf(name_cPt,"plots/%s_Pt_fit_%.2f_m_%.2f.png",ZDC_class.Data(),m_min,m_max);
  cPt->SaveAs(name_cPt);	 
}


//////////////////////////////////////////////////////////////////////
// Run the macro
//////////////////////////////////////////////////////////////////////
void Draw_spectra()
{

Double_t pt_min = 0.0; 
Double_t pt_max = 2.7;
Double_t m_min = 0;
Double_t m_max = 6;

 DrawFunction(-4.0, -2.5, 0, 6, pt_min, pt_max, "all", "ADvetoOff", "CMUP6");

 DrawFunction(-4.0, -2.5, 2.4, 2.7, pt_min, pt_max, "all", "ADvetoOff", "CMUP6");
 DrawFunction(-4.0, -2.5, 4, 6, pt_min, pt_max, "all", "ADvetoOff", "CMUP6");


 DrawFunction(-4.0, -2.5, m_min, m_max, 0, 2.7, "all", "ADvetoOff", "CMUP6");

 DrawFunction(-4.0, -2.5, m_min, m_max, 0, 0.2, "all", "ADvetoOff", "CMUP6");
 DrawFunction(-4.0, -2.5, m_min, m_max, 0.2, 0.4, "all", "ADvetoOff", "CMUP6");
 DrawFunction(-4.0, -2.5, m_min, m_max, 0.4, 0.7, "all", "ADvetoOff", "CMUP6");
 DrawFunction(-4.0, -2.5, m_min, m_max, 0.7, 1, "all", "ADvetoOff", "CMUP6");
 DrawFunction(-4.0, -2.5, m_min, m_max, 1, 2.7, "all", "ADvetoOff", "CMUP6");

}
