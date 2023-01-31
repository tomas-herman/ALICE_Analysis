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
#include <RooFormulaVar.h>
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
#include "Variables.h"
// #include "GoodRuns.h"

// My functions
// #include "LoadGoodRuns.C"
#include "LoadData.C"
#include "LoadEff.C"
#include "LoadRatio.C"
#include "LoadMCmassParam.C"
#include "LoadPtParam.C"
#include "SetCuts.C"

///////////////////////////////////////////////////////////
// Fit the invariant mass and pt spectra
///////////////////////////////////////////////////////////

void DoFit(Double_t y_min, Double_t y_max, TString ZDC_class, TString ADveto, TString trigger, 
          TString ratio_str, TString MC_fit_str, TString Disoc_fit_str, TString PDF_str)
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
  // Define cuts
  SetCuts(y_min, y_max, ZDC_class, ADveto);

  // --------------------------------------------------------------------------------
  // Load efficiencies
  double Eff_CohJpsiToMu_PtCut;
  double Eff_CohJpsiToMu_PtAll;
  double Eff_CohPsi2sToMu_PtCut;
  double Eff_CohPsi2sToMu_PtAll;
  double Eff_CohPsi2sToMuPi_PtCut;
  double Eff_CohPsi2sToMuPi_PtAll;

  LoadEff(y_min, y_max, 
          Eff_CohJpsiToMu_PtCut, Eff_CohPsi2sToMu_PtCut, Eff_CohPsi2sToMuPi_PtCut,
          Eff_CohJpsiToMu_PtAll, Eff_CohPsi2sToMu_PtAll, Eff_CohPsi2sToMuPi_PtAll, 
          ADveto, trigger);

  RooRealVar Var_Eff_CohJpsiToMu_PtCut("Var_Eff_CohJpsiToMu_PtCut","Var_Eff_CohJpsiToMu_PtCut",Eff_CohJpsiToMu_PtCut);
  RooRealVar Var_Eff_CohJpsiToMu_PtAll("Var_Eff_CohJpsiToMu_PtAll","Var_Eff_CohJpsiToMu_PtAll",Eff_CohJpsiToMu_PtAll);
  RooRealVar Var_Eff_CohPsi2sToMu_PtCut("Var_Eff_CohPsi2sToMu_PtCut","Var_Eff_CohPsi2sToMu_PtCut",Eff_CohPsi2sToMu_PtCut);
  RooRealVar Var_Eff_CohPsi2sToMu_PtAll("Var_Eff_CohPsi2sToMu_PtAll","Var_Eff_CohPsi2sToMu_PtAll",Eff_CohPsi2sToMu_PtAll);
  RooRealVar Var_Eff_CohPsi2sToMuPi_PtCut("Var_Eff_CohPsi2sToMuPi_PtCut","Var_Eff_CohPsi2sToMuPi_PtCut",Eff_CohPsi2sToMuPi_PtCut);
  RooRealVar Var_Eff_CohPsi2sToMuPi_PtAll("Var_Eff_CohPsi2sToMuPi_PtAll","Var_Eff_CohPsi2sToMuPi_PtAll",Eff_CohPsi2sToMuPi_PtAll);

  Var_Eff_CohJpsiToMu_PtCut.setConstant();
  Var_Eff_CohJpsiToMu_PtAll.setConstant();
  Var_Eff_CohPsi2sToMu_PtCut.setConstant();
  Var_Eff_CohPsi2sToMu_PtAll.setConstant();
  Var_Eff_CohPsi2sToMuPi_PtCut.setConstant();
  Var_Eff_CohPsi2sToMuPi_PtAll.setConstant();

  // --------------------------------------------------------------------------------
  // Load data
  RooArgSet Arguments(fMuMuM,fMuMuY,fMuMuPt,fRunNum,fADADecision,fADCDecision,fV0ADecision,fV0CDecision,fV0CFiredCells);
  Arguments.add( RooArgSet(fIsZNAFired,fIsZNCFired) ); // You can only add 9 arguments at a time
  // ---MC data sets
  // ------LHC18l7
  RooDataSet *dataIN18l7_CohJpsiToMu = new RooDataSet ("dataIN18l7_CohJpsiToMu", "dataIN18l7_CohJpsiToMu", Arguments);
  RooDataSet *dataIN18l7_CohPsi2sToMu = new RooDataSet ("dataIN18l7_CohPsi2sToMu", "dataIN18l7_CohPsi2sToMu", Arguments);
  RooDataSet *dataIN18l7_CohPsi2sToMuPi = new RooDataSet ("dataIN18l7_CohPsi2sToMuPi", "dataIN18l7_CohPsi2sToMuPi", Arguments);
  RooDataSet *dataIN18l7_IncohJpsiToMu = new RooDataSet ("dataIN18l7_IncohJpsiToMu", "dataIN18l7_IncohJpsiToMu", Arguments);
  RooDataSet *dataIN18l7_IncohPsi2sToMu = new RooDataSet ("dataIN18l7_IncohPsi2sToMu", "dataIN18l7_IncohPsi2sToMu", Arguments);
  RooDataSet *dataIN18l7_IncohPsi2sToMuPi = new RooDataSet ("dataIN18l7_IncohPsi2sToMuPi", "dataIN18l7_IncohPsi2sToMuPi", Arguments);
  RooDataSet *dataIN18l7_TwoGammaToMuMedium = new RooDataSet ("dataIN18l7_TwoGammaToMuMedium", "dataIN18l7_TwoGammaToMuMedium", Arguments);

  if (PDF_str.Contains("Create_PDFs"))
  {
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/CohJpsiToMu_AnalysisResults.root", "MCRec", "NanoMUONCMUP6", dataIN18l7_CohJpsiToMu);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/CohPsi2sToMu_AnalysisResults.root", "MCRec", "NanoMUONCMUP6", dataIN18l7_CohPsi2sToMu);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/CohPsi2sToMuPi_AnalysisResults.root", "MCRec", "NanoMUONCMUP6", dataIN18l7_CohPsi2sToMuPi);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/IncohJpsiToMu_AnalysisResults.root", "MCRec", "NanoMUONCMUP6", dataIN18l7_IncohJpsiToMu);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/IncohPsi2sToMu_AnalysisResults.root", "MCRec", "NanoMUONCMUP6", dataIN18l7_IncohPsi2sToMu);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/IncohPsi2sToMuPi_AnalysisResults.root", "MCRec", "NanoMUONCMUP6", dataIN18l7_IncohPsi2sToMuPi);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/TwoGammaToMuMedium_AnalysisResults.root", "MCRec", "NanoMUONCMUP6", dataIN18l7_TwoGammaToMuMedium);
  }

  // ------LHC16b2
  RooDataSet *dataIN16b2_CohJpsiToMu = new RooDataSet ("dataIN16b2_CohJpsiToMu", "dataIN16b2_CohJpsiToMu", Arguments);
  RooDataSet *dataIN16b2_CohPsi2sToMu = new RooDataSet ("dataIN16b2_CohPsi2sToMu", "dataIN16b2_CohPsi2sToMu", Arguments);
  RooDataSet *dataIN16b2_CohPsi2sToMuPi = new RooDataSet ("dataIN16b2_CohPsi2sToMuPi", "dataIN16b2_CohPsi2sToMuPi", Arguments);
  RooDataSet *dataIN16b2_IncohJpsiToMu = new RooDataSet ("dataIN16b2_IncohJpsiToMu", "dataIN16b2_IncohJpsiToMu", Arguments);
  RooDataSet *dataIN16b2_IncohPsi2sToMu = new RooDataSet ("dataIN16b2_IncohPsi2sToMu", "dataIN16b2_IncohPsi2sToMu", Arguments);
  RooDataSet *dataIN16b2_IncohPsi2sToMuPi = new RooDataSet ("dataIN16b2_IncohPsi2sToMuPi", "dataIN16b2_IncohPsi2sToMuPi", Arguments);
  RooDataSet *dataIN16b2_TwoGammaToMuMedium = new RooDataSet ("dataIN16b2_TwoGammaToMuMedium", "dataIN16b2_TwoGammaToMuMedium", Arguments);

  if (PDF_str.Contains("Create_PDFs"))
  {
    if (trigger.Contains("CMUP11"))
    {
      LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/578_20200702-1049-LHC16b2/CohJpsiToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_CohJpsiToMu);
      LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/578_20200702-1049-LHC16b2/CohPsi2sToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_CohPsi2sToMu);
      LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/578_20200702-1049-LHC16b2/CohPsi2sToMuPi_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_CohPsi2sToMuPi);
      LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/578_20200702-1049-LHC16b2/IncohJpsiToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_IncohJpsiToMu);
      LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/578_20200702-1049-LHC16b2/IncohPsi2sToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_IncohPsi2sToMu);
      LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/578_20200702-1049-LHC16b2/IncohPsi2sToMuPi_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_IncohPsi2sToMuPi);
      LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/578_20200702-1049-LHC16b2/TwoGammaToMuMedium_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_TwoGammaToMuMedium);
    }
  }

  // ---Real data sets
  RooDataSet *dataIN18q = new RooDataSet ("dataIN18q", "dataIN18q", Arguments);
  RooDataSet *dataIN18r = new RooDataSet ("dataIN18r", "dataIN18r", Arguments);
  RooDataSet *dataIN15o = new RooDataSet ("dataIN15o", "dataIN15o", Arguments);

  if (trigger.Contains("CMUP6"))
  {
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18q_AnalysisResults.root", "DataRec" , "NanoMUONCMUP6", dataIN18q);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18r_AnalysisResults.root", "DataRec" , "NanoMUONCMUP6", dataIN18r);
  }

  if (trigger.Contains("CMUP11"))
  {
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18q_AnalysisResults.root", "DataRec", "NanoMUONCMUP11", dataIN18q);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18r_AnalysisResults.root", "DataRec", "NanoMUONCMUP11", dataIN18r);
    LoadData("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC15o_AnalysisResults.root", "DataRec", "NanoMUONCMUP11", dataIN15o);
  }

  // --------------------------------------------------------------------------------
  // Merge MC data sets
  RooDataSet *dataINMC_CohJpsiToMu = (RooDataSet*) dataIN18l7_CohJpsiToMu;
  if (trigger.Contains("CMUP11")){ dataINMC_CohJpsiToMu->append(*dataIN16b2_CohJpsiToMu); }

  RooDataSet *dataINMC_CohPsi2sToMu = (RooDataSet*) dataIN18l7_CohPsi2sToMu;
  if (trigger.Contains("CMUP11")){ dataINMC_CohPsi2sToMu->append(*dataIN16b2_CohPsi2sToMu); }

  RooDataSet *dataINMC_CohPsi2sToMuPi = (RooDataSet*) dataIN18l7_CohPsi2sToMuPi;
  if (trigger.Contains("CMUP11")){ dataINMC_CohPsi2sToMuPi->append(*dataIN16b2_CohPsi2sToMuPi); }

  RooDataSet *dataINMC_IncohJpsiToMu = (RooDataSet*) dataIN18l7_IncohJpsiToMu;
  if (trigger.Contains("CMUP11")){ dataINMC_IncohJpsiToMu->append(*dataIN16b2_IncohJpsiToMu); }

  RooDataSet *dataINMC_IncohPsi2sToMu = (RooDataSet*) dataIN18l7_IncohPsi2sToMu;
  if (trigger.Contains("CMUP11")){ dataINMC_IncohPsi2sToMu->append(*dataIN16b2_IncohPsi2sToMu); }

  RooDataSet *dataINMC_IncohPsi2sToMuPi = (RooDataSet*) dataIN18l7_IncohPsi2sToMuPi;
  if (trigger.Contains("CMUP11")){ dataINMC_IncohPsi2sToMuPi->append(*dataIN16b2_IncohPsi2sToMuPi); }

  RooDataSet *dataINMC_TwoGammaToMuMedium = (RooDataSet*) dataIN18l7_TwoGammaToMuMedium;
  if (trigger.Contains("CMUP11")){ dataINMC_TwoGammaToMuMedium->append(*dataIN16b2_TwoGammaToMuMedium); }

  // --------------------------------------------------------------------------------
  // Create MC pt templates
  // ---Prepare PDFs
  RooAbsPdf* PDF_CohJpsiToMu; 
  RooAbsPdf* PDF_CohPsi2sToMu; 
  RooAbsPdf* PDF_CohPsi2sToMuPi; 
  RooAbsPdf* PDF_IncohJpsiToMu; 
  RooAbsPdf* PDF_IncohPsi2sToMu; 
  RooAbsPdf* PDF_IncohPsi2sToMuPi;   
  RooAbsPdf* PDF_TwoGammaToMuMedium;

  // ---Prepare integrals
  RooRealVar *ncoh = new RooRealVar("ncoh", "ncoh", 1, "" );
  RooRealVar *ncoh2 = new RooRealVar("ncoh2", "ncoh2", 1, "" );
  RooRealVar *nincoh = new RooRealVar("nincoh", "nincoh", 1, "" );
  RooRealVar *nincoh2 = new RooRealVar("nincoh2", "nincoh2", 1, "" ); 

  if (PDF_str.Contains("Create_PDFs"))
  {
    // ---Do kinematic cuts
    RooDataSet *dataMC_CohJpsiToMu_Template = (RooDataSet*) dataINMC_CohJpsiToMu->reduce(Cut(cuts_pt_mc));
    RooDataSet *dataMC_CohPsi2sToMu_Template = (RooDataSet*) dataINMC_CohPsi2sToMu->reduce(Cut(cuts_pt_mc));
    RooDataSet *dataMC_CohPsi2sToMuPi_Template = (RooDataSet*) dataINMC_CohPsi2sToMuPi->reduce(Cut(cuts_pt_mc));
    RooDataSet *dataMC_IncohJpsiToMu_Template = (RooDataSet*) dataINMC_IncohJpsiToMu->reduce(Cut(cuts_pt_mc));
    RooDataSet *dataMC_IncohPsi2sToMu_Template = (RooDataSet*) dataINMC_IncohPsi2sToMu->reduce(Cut(cuts_pt_mc));
    RooDataSet *dataMC_IncohPsi2sToMuPi_Template = (RooDataSet*) dataINMC_IncohPsi2sToMuPi->reduce(Cut(cuts_pt_mc));
    RooDataSet *dataMC_TwoGammaToMuMedium_Template = (RooDataSet*) dataINMC_TwoGammaToMuMedium->reduce(Cut(cuts_pt_mc));

    // ---Create histograms
    RooDataHist Hist_dataMC_CohJpsiToMu("Hist_dataMC_CohJpsiToMu","Hist_dataMC_CohJpsiToMu",fMuMuPt,*dataMC_CohJpsiToMu_Template,1);
    RooDataHist Hist_dataMC_CohPsi2sToMu("Hist_dataMC_CohPsi2sToMu","Hist_dataMC_CohPsi2sToMu",fMuMuPt,*dataMC_CohPsi2sToMu_Template,1);
    RooDataHist Hist_dataMC_CohPsi2sToMuPi("Hist_dataMC_CohPsi2sToMuPi","Hist_dataMC_CohPsi2sToMuPi",fMuMuPt,*dataMC_CohPsi2sToMuPi_Template,1);
    RooDataHist Hist_dataMC_IncohJpsiToMu("Hist_dataMC_IncohJpsiToMu","Hist_dataMC_IncohJpsiToMu",fMuMuPt,*dataMC_IncohJpsiToMu_Template,1);
    RooDataHist Hist_dataMC_IncohPsi2sToMu("Hist_dataMC_IncohPsi2sToMu","Hist_dataMC_IncohPsi2sToMu",fMuMuPt,*dataMC_IncohPsi2sToMu_Template,1);
    RooDataHist Hist_dataMC_IncohPsi2sToMuPi("Hist_dataMC_IncohPsi2sToMuPi","Hist_dataMC_IncohPsi2sToMuPi",fMuMuPt,*dataMC_IncohPsi2sToMuPi_Template,1);
    RooDataHist Hist_dataMC_TwoGammaToMuMedium("Hist_dataMC_TwoGammaToMuMedium","Hist_dataMC_TwoGammaToMuMedium",fMuMuPt,*dataMC_TwoGammaToMuMedium_Template,1);

    // ---Create PDFs
    RooHistPdf Hist_PDF_CohJpsiToMu("PDF_CohJpsiToMu", "PDF_CohJpsiToMu", fMuMuPt,Hist_dataMC_CohJpsiToMu,0); 
    RooHistPdf Hist_PDF_CohPsi2sToMu("PDF_CohPsi2sToMu", "PDF_CohPsi2sToMu", fMuMuPt,Hist_dataMC_CohPsi2sToMu,0); 
    RooHistPdf Hist_PDF_CohPsi2sToMuPi("PDF_CohPsi2sToMuPi", "PDF_CohPsi2sToMuPi", fMuMuPt,Hist_dataMC_CohPsi2sToMuPi,0); 
    RooHistPdf Hist_PDF_IncohJpsiToMu("PDF_IncohJpsiToMu", "PDF_IncohJpsiToMu", fMuMuPt,Hist_dataMC_IncohJpsiToMu,0); 
    RooHistPdf Hist_PDF_IncohPsi2sToMu("PDF_IncohPsi2sToMu", "PDF_IncohPsi2sToMu", fMuMuPt,Hist_dataMC_IncohPsi2sToMu,0); 
    RooHistPdf Hist_PDF_IncohPsi2sToMuPi("PDF_IncohPsi2sToMuPi", "PDF_IncohPsi2sToMuPi", fMuMuPt,Hist_dataMC_IncohPsi2sToMuPi,0);   
    RooHistPdf Hist_PDF_TwoGammaToMuMedium("PDF_TwoGammaToMuMedium", "PDF_TwoGammaToMuMedium", fMuMuPt,Hist_dataMC_TwoGammaToMuMedium,0);

    PDF_CohJpsiToMu = (RooAbsPdf*) &Hist_PDF_CohJpsiToMu; 
    PDF_CohPsi2sToMu = (RooAbsPdf*) &Hist_PDF_CohPsi2sToMu; 
    PDF_CohPsi2sToMuPi = (RooAbsPdf*) &Hist_PDF_CohPsi2sToMuPi; 
    PDF_IncohJpsiToMu = (RooAbsPdf*) &Hist_PDF_IncohJpsiToMu; 
    PDF_IncohPsi2sToMu = (RooAbsPdf*) &Hist_PDF_IncohPsi2sToMu; 
    PDF_IncohPsi2sToMuPi = (RooAbsPdf*) &Hist_PDF_IncohPsi2sToMuPi;   
    PDF_TwoGammaToMuMedium = (RooAbsPdf*) &Hist_PDF_TwoGammaToMuMedium;

    // ---Create integrals
    RooAbsReal *icoh = PDF_CohJpsiToMu->createIntegral(fMuMuPt,NormSet(fMuMuPt),Range("JPsiptRange"));
    RooAbsReal *icoh2 = PDF_CohPsi2sToMuPi->createIntegral(fMuMuPt,NormSet(fMuMuPt),Range("JPsiptRange"));
    RooAbsReal *iincoh = PDF_IncohJpsiToMu->createIntegral(fMuMuPt,NormSet(fMuMuPt),Range("JPsiptRange"));
    RooAbsReal *iincoh2 = PDF_IncohPsi2sToMuPi->createIntegral(fMuMuPt,NormSet(fMuMuPt),Range("JPsiptRange"));

    // ---Create variabels from integral values in given range
    ncoh->setVal(icoh->getVal());
    ncoh2->setVal(icoh2->getVal());
    nincoh->setVal(iincoh->getVal());
    nincoh2->setVal(iincoh2->getVal());

    // ---Safe PDFs and integrals to a file
    RooWorkspace *w_PDFs = new RooWorkspace("w_PDFs","workspace with PDFs") ;

    // ---PDFs
    w_PDFs->import(*PDF_CohJpsiToMu);
    w_PDFs->import(*PDF_CohPsi2sToMu);
    w_PDFs->import(*PDF_CohPsi2sToMuPi);
    w_PDFs->import(*PDF_IncohJpsiToMu);
    w_PDFs->import(*PDF_IncohPsi2sToMu);
    w_PDFs->import(*PDF_IncohPsi2sToMuPi);
    w_PDFs->import(*PDF_TwoGammaToMuMedium);

    // ---Integrals
    w_PDFs->import(*ncoh);
    w_PDFs->import(*ncoh2);
    w_PDFs->import(*nincoh);
    w_PDFs->import(*nincoh2);

    char FileName_PDFs[120];
    sprintf(FileName_PDFs,"Workspaces/PDFs/%s/%s_%.2f_%.2f_PDFs.root",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));

    w_PDFs->writeToFile(FileName_PDFs, kTRUE);
  }
  
  // ---Load PDFs from a file
  if (PDF_str.Contains("Load_PDFs"))
  {
    char FileName_PDFs[120];
    sprintf(FileName_PDFs,"Workspaces/PDFs/%s/%s_%.2f_%.2f_PDFs.root",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
    
    // ---Open file
    TFile *f_PDFs = new TFile(FileName_PDFs);
    RooWorkspace* w_PDFs = (RooWorkspace*) f_PDFs->Get("w_PDFs");

    // ---Load PDFs
    PDF_CohJpsiToMu = w_PDFs->pdf("PDF_CohJpsiToMu");
    PDF_CohPsi2sToMu = w_PDFs->pdf("PDF_CohPsi2sToMu");
    PDF_CohPsi2sToMuPi = w_PDFs->pdf("PDF_CohPsi2sToMuPi");
    PDF_IncohJpsiToMu = w_PDFs->pdf("PDF_IncohJpsiToMu");
    PDF_IncohPsi2sToMu = w_PDFs->pdf("PDF_IncohPsi2sToMu");
    PDF_IncohPsi2sToMuPi = w_PDFs->pdf("PDF_IncohPsi2sToMuPi");
    PDF_TwoGammaToMuMedium = w_PDFs->pdf("PDF_TwoGammaToMuMedium");

    // ---Load integrals
    ncoh = w_PDFs->var("ncoh");
    ncoh2 = w_PDFs->var("ncoh2");
    nincoh = w_PDFs->var("nincoh");
    nincoh2 = w_PDFs->var("nincoh2");
  }

  // ---Create Incoherent Disocitation PDF
  // Values taken from https://arxiv.org/abs/1304.5162
  RooRealVar b("b","b",1.79, 0.1, 6); //1.67, 1.91); 
  RooRealVar n("n","n",3.58, 1.5,8); //3.43,3.73); 

  if(ADveto.Contains("ADvetoOn")) n.setConstant(kTRUE);
  b.setConstant(kTRUE);

  RooGenericPdf PDF_IncohJpsiToX("PDF_IncohJpsiToX","fMuMuPt*pow((1+pow(fMuMuPt,2)*b/n),-n)",RooArgSet(fMuMuPt, b, n)); 
  
  // --------------------------------------------------------------------------------
  // Get real data
  // ---Merge data sets
  RooDataSet *dataAllIN = (RooDataSet*) dataIN18q;
  dataAllIN->append(*dataIN18r);
  if (trigger.Contains("CMUP11")){ dataAllIN->append(*dataIN15o); }
  // ---Do cuts
  RooAbsData* data_pt = dataAllIN->reduce(Cut(cuts_pt_data));
  RooAbsData* data_m = dataAllIN->reduce(Cut(cuts_m_data));

  // This is used if I only want to count number of entries
  // cout << data_m->numEntries() << endl;
  // gROOT->ProcessLine(".q");

  RooAbsData* data_m_PtAll = dataAllIN->reduce(Cut(cuts_m_data_PtAll));
  // --------------------------------------------------------------------------------
  // Get Toy MC data
  // ---Do cuts

  // RooDataHist *dataMC_CohJpsiToMu_m_PtCut_binned = new RooDataHist("dataMC_CohJpsiToMu_m_PtCut_binned", "binned version of dataMC_CohJpsiToMu_m_PtCut", RooArgSet(fMuMuM), *dataMC_CohJpsiToMu_m_PtCut);
  RooAbsData* dataMC_CohJpsiToMu_m_PtCut = dataINMC_CohJpsiToMu->reduce(Cut(cuts_m_mc));
  RooAbsData* dataMC_CohPsi2sToMu_m_PtCut = dataINMC_CohPsi2sToMu->reduce(Cut(cuts_m_mc));
  RooAbsData* dataMC_TwoGammaToMuMedium_m_PtCut = dataINMC_TwoGammaToMuMedium->reduce(Cut(cuts_m_mc));

  RooAbsData* dataMC_CohJpsiToMu_m_PtAll = dataINMC_CohJpsiToMu->reduce(Cut(cuts_m_mc_PtAll));
  RooAbsData* dataMC_CohPsi2sToMu_m_PtAll = dataINMC_CohPsi2sToMu->reduce(Cut(cuts_m_mc_PtAll));
  RooAbsData* dataMC_TwoGammaToMuMedium_m_PtAll = dataINMC_TwoGammaToMuMedium->reduce(Cut(cuts_m_mc_PtAll));

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Invariant mass fit
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Int_t N_Evnts_m_PtAll = data_m_PtAll->numEntries();
  Int_t N_Evnts_m = data_m->numEntries();

  Int_t N_Evnts_CohJpsiToMu_m = dataMC_CohJpsiToMu_m_PtCut->numEntries();
  Int_t N_Evnts_CohPsi2sToMu_m = dataMC_CohPsi2sToMu_m_PtCut->numEntries();
  Int_t N_Evnts_TwoGammaToMuMedium_m = dataMC_TwoGammaToMuMedium_m_PtCut->numEntries();

  // ---J/Psi parameters
  RooRealVar m0("m","mass",3.096,3.05,3.15) ;
  RooRealVar n_cb("n_{cb}","n_cb",8, 0, 120) ;
  RooRealVar sigma("sigma","J/Psi width",0.092,0.01,0.2) ;
  RooRealVar sigmaMC("sigmaMC","J/Psi width in MC",0.092,0.01,0.2) ;
  RooRealVar alpha("#alpha","alpha",0.97,0,20) ;
  // ---Psi' parameters
  RooFormulaVar m02("m02","m+3.686097-3.096900",RooArgList(m0));
  RooRealVar sigma2MC("sigma2MC","Psi' width in MC ",0.10,0.01,0.20) ;
  
  // ---Background gamma gamma parameters
  RooRealVar lambda("lambda","exponent",-0.9,-9.,-0.05);
  RooRealVar a2("a2","parameter a2",0.72,0.01,1);
  RooRealVar a3("a3","parameter a3",0.99,0.01,1);
  RooRealVar a4("a4","parameter a4",0.26,0.01,1);


  // --------------------------------------------------------------------------------
  // Create fit functions
  // ---To use Gauss Exp
    // RooRealVar cb("cb","Gauss Exp PDF",2,5) ;
    // RooFormulaVar cb("cb","(((fMuMuM-m0)/sigma)>=-alpha)? exp(-0.5*pow((fMuMuM-m0)/sigma,2)) : exp(0.5*pow(alpha,2)+alpha*(fMuMuM-m0)/sigma)",RooArgList(fMuMuM,m0,sigma,alpha));     
    // ---Hand written CB
    // RooGenericPdf cbMC("cbMC","(((fMuMuM-m0)/sigmaMC)>-alpha)? exp(-0.5*pow((fMuMuM-m0)/sigmaMC,2)) : pow(n_cb/abs(alpha),n_cb)*exp(-0.5*pow(alpha,2))*pow(n_cb/abs(alpha)-abs(alpha)-(fMuMuM-m0)/sigmaMC,-n_cb)",RooArgSet(fMuMuM,m0,sigmaMC,alpha,n_cb));

  // ---Crystal ball functions for J/Psi and Psi' in data and MC
  RooCBShape cb("cb","Crystal Ball 1 PDF",fMuMuM,m0,sigma,alpha,n_cb) ;
  RooCBShape cbMC("cbMC","Crystal Ball 1 MC PDF",fMuMuM,m0,sigmaMC,alpha,n_cb) ;
  RooCBShape cb2MC("cb2MC","Crystal Ball 2 MC PDF",fMuMuM,m02,sigma2MC,alpha,n_cb) ;

  //---PDF to describe the gamma gamma background
  RooGenericPdf Bkgd("Bkgd","(fMuMuM>4)? exp(fMuMuM*lambda) : exp(fMuMuM*lambda)*(1+a2*pow((fMuMuM-4),2)+a3*pow((fMuMuM-4),3)+a4*pow((fMuMuM-4),4))",RooArgSet(fMuMuM,lambda,a2,a3,a4));
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Pt All
  // Create PDF to fit MC
  // ---PDF to fit MC gamma gamma background 
  RooRealVar N_MC_BG_PtAll("N_MC_{bg}","Number of MC BG events",1.0*N_Evnts_TwoGammaToMuMedium_m,0,3*N_Evnts_TwoGammaToMuMedium_m);
  RooAddPdf Bkgd_mc_PtAll("Bkgd_mc_PtAll","Background PDF for MC fit", RooArgList(Bkgd), RooArgList(N_MC_BG_PtAll));

  // ---PDF to fit MC Coherent Psi'
  RooRealVar N_MC_CB2_PtAll("N_MC_{cb2}","Number of MC Psi' events",1.0*N_Evnts_CohPsi2sToMu_m,0,3*N_Evnts_CohPsi2sToMu_m);
  RooAddPdf cb2_mc_PtAll("cb2_mc_PtAll","CB2 PDF for MC fit", RooArgList(cb2MC), RooArgList(N_MC_CB2_PtAll));
  
  // ---PDF to fit MC Coherent J/Psi
  RooRealVar N_MC_CB_PtAll("N_MC_{cb}","Number of MC J/Psi events",1.0*N_Evnts_CohJpsiToMu_m,0,3*N_Evnts_CohJpsiToMu_m);
  RooAddPdf cb_mc_PtAll("cb_mc_PtAll","CB PDF for MC fit", RooArgList(cbMC), RooArgList(N_MC_CB_PtAll));

  // ---Get the CB and bkgd parameters by fitting MC
  if (MC_fit_str.Contains("Compute_MC_fit_param"))
  {
    // ---Fit MC data to fix gamma gamma background parameters
    RooFitResult* fit_bkgd_PtAll = Bkgd_mc_PtAll.fitTo(*dataMC_TwoGammaToMuMedium_m_PtAll,Extended(kTRUE),Range("m_fitRange"),Save());
    // ---Fit MC data to get the sigma ratio of the Psi' and J/Psi
    RooFitResult* fit_cb2_PtAll = cb2_mc_PtAll.fitTo(*dataMC_CohPsi2sToMu_m_PtAll,Extended(kTRUE),Range("CB2_fitRange"),Save());
    // ---Fit MC data to fix J/Psi crystall ball parameters
    RooFitResult* fit_cb_PtAll = cb_mc_PtAll.fitTo(*dataMC_CohJpsiToMu_m_PtAll,Extended(kTRUE),Range("CB_fitRange"),Save());
    
    // --------------------------------------------------------------------------------
    // ---Plotting the MC fits

    // ---gamma gamma background 
    TCanvas *cMC_PtAll_bkgd = new TCanvas("cMC_PtAll_bkgd","cMC_PtAll_bkgd",1600,800);
    cMC_PtAll_bkgd->cd();

    RooPlot* frame_MC_PtAll_bkgd = fMuMuM.frame(Title("Mass fit")) ; 
    dataMC_TwoGammaToMuMedium_m_PtAll->plotOn(frame_MC_PtAll_bkgd,Name("MC_bkgd_PtAll"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
    Bkgd_mc_PtAll.plotOn(frame_MC_PtAll_bkgd,Name("Bkgd_mc_PtAll"),LineColor(kBlack), LineWidth(1)) ;
    frame_MC_PtAll_bkgd->Draw();

    char name_cMC_PtAll_bkgd[120];
    sprintf(name_cMC_PtAll_bkgd,"Fit_param/%s/%s_%.2f_%.2f_MC_PtAll_bkgd.png",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
    cMC_PtAll_bkgd->SaveAs(name_cMC_PtAll_bkgd);

    // ---J/Psi crystall ball 
    TCanvas *cMC_PtAll_cb = new TCanvas("cMC_PtAll_cb","cMC_PtAll_cb",1600,800);
    cMC_PtAll_cb->SetLogy();

    RooPlot* frame_MC_PtAll_cb = fMuMuM.frame(Title("Mass fit")) ; 
    dataMC_CohJpsiToMu_m_PtAll->plotOn(frame_MC_PtAll_cb,Name("MC_cb_PtAll"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
    cb_mc_PtAll.plotOn(frame_MC_PtAll_cb,Range("m_fitRange"),Name("cb_mc_PtAll"),LineColor(kBlack), LineWidth(1)) ;
    frame_MC_PtAll_cb->SetAxisRange(0.025,10000000,"Y");
    frame_MC_PtAll_cb->Draw();

    char name_cMC_PtAll_cb[120];
    sprintf(name_cMC_PtAll_cb,"Fit_param/%s/%s_%.2f_%.2f_MC_PtAll_cb.png",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
    cMC_PtAll_cb->SaveAs(name_cMC_PtAll_cb);

    // ---Psi' crystall ball 
    TCanvas *cMC_PtAll_cb2 = new TCanvas("cMC_PtAll_cb2","cMC_PtAll_cb2",1600,800);
    cMC_PtAll_cb2->cd();
    cMC_PtAll_cb2->SetLogy();

    RooPlot* frame_MC_PtAll_cb2 = fMuMuM.frame(Title("Mass fit")) ; 
    dataMC_CohPsi2sToMu_m_PtAll->plotOn(frame_MC_PtAll_cb2,Name("MC_cb2_PtAll"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
    cb2_mc_PtAll.plotOn(frame_MC_PtAll_cb2,Range("m_fitRange"),Name("cb2_mc_PtAll"),LineColor(kBlack), LineWidth(1)) ;
    frame_MC_PtAll_cb2->SetAxisRange(0.025,10000000,"Y");
    frame_MC_PtAll_cb2->Draw();

    char name_cMC_PtAll_cb2[120];
    sprintf(name_cMC_PtAll_cb2,"Fit_param/%s/%s_%.2f_%.2f_MC_PtAll_cb2.png",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
    cMC_PtAll_cb2->SaveAs(name_cMC_PtAll_cb2);


    // --------------------------------------------------------------------------------
    // Write the values from the MC fits
    char FileName_MC_m_param_PtAll[120];
    sprintf(FileName_MC_m_param_PtAll,"Fit_param/%s/%s_%.2f_%.2f_MC_m_param_PtAll.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));

    ofstream MC_param_File_PtAll;
    MC_param_File_PtAll.open (FileName_MC_m_param_PtAll, ios::trunc);
    MC_param_File_PtAll << "-----------------------------------------------------------" << endl;
    MC_param_File_PtAll << "lambda_val = " << lambda.getVal() << endl;
    MC_param_File_PtAll << "a2_val = " << a2.getVal() << endl;
    MC_param_File_PtAll << "a3_val = " << a3.getVal() << endl;
    MC_param_File_PtAll << "a4_val = " << a4.getVal() << endl;
    MC_param_File_PtAll << "n_cb_val = " << n_cb.getVal() << endl;
    MC_param_File_PtAll << "alpha_val = " << alpha.getVal() << endl;
    MC_param_File_PtAll << "sigmaMC_val = " << sigmaMC.getVal() << endl;
    MC_param_File_PtAll << "sigma2MC_val = " << sigma2MC.getVal() << endl;
    MC_param_File_PtAll << "lambda_err = " << lambda.getError() << endl;
    MC_param_File_PtAll << "a2_err = " << a2.getError() << endl;
    MC_param_File_PtAll << "a3_err = " << a3.getError() << endl;
    MC_param_File_PtAll << "a4_err = " << a4.getError() << endl;
    MC_param_File_PtAll << "n_cb_err = " << n_cb.getError() << endl;
    MC_param_File_PtAll << "alpha_err = " << alpha.getError() << endl;
    MC_param_File_PtAll << "sigmaMC_err = " << sigmaMC.getError() << endl;
    MC_param_File_PtAll << "sigma2MC_err = " << sigma2MC.getError() << endl;
    MC_param_File_PtAll << "-----------------------------------------------------------" << endl;
    MC_param_File_PtAll << endl;  
    MC_param_File_PtAll.close();
  }

  double lambda_val, lambda_err;
  double a2_val, a2_err;
  double a3_val, a3_err;
  double a4_val, a4_err;
  double n_cb_val, n_cb_err;
  double alpha_val, alpha_err;
  double sigmaMC_val, sigmaMC_err;
  double sigma2MC_val, sigma2MC_err;

  // --------------------------------------------------------------------------------
  // Load mass fit parameters from MC
  if (MC_fit_str.Contains("Load_MC_fit_param"))
  {
    LoadMCmassParam(y_min, y_max, ADveto, trigger, "PtAll",
                    lambda_val, lambda_err,
                    a2_val, a2_err,
                    a3_val, a3_err,
                    a4_val, a4_err,
                    n_cb_val, n_cb_err,
                    alpha_val, alpha_err,
                    sigmaMC_val, sigmaMC_err,
                    sigma2MC_val, sigma2MC_err);

    lambda.setVal(lambda_val);
    a2.setVal(a2_val);
    a3.setVal(a3_val);
    a4.setVal(a4_val);
    n_cb.setVal(n_cb_val);
    alpha.setVal(alpha_val);
    sigmaMC.setVal(sigmaMC_val);
    sigma2MC.setVal(sigma2MC_val);

    lambda.setError(lambda_err);
    a2.setError(a2_err);
    a3.setError(a3_err);
    a4.setError(a4_err);
    n_cb.setError(n_cb_err);
    alpha.setError(alpha_err);
    sigmaMC.setError(sigmaMC_err);
    sigma2MC.setError(sigma2MC_err);
  }

  // lambda.setConstant(kTRUE);
  a2.setConstant(kTRUE);
  a3.setConstant(kTRUE);
  a4.setConstant(kTRUE);
  n_cb.setConstant(kTRUE);
  alpha.setConstant(kTRUE);
  sigmaMC.setConstant(kTRUE);
  sigma2MC.setConstant(kTRUE);
  
  RooFormulaVar sigma_ratio("sigma_ratio","sigma2MC/sigmaMC",RooArgList(sigma2MC,sigmaMC));

  // ---Crystal ball functions for Psi'
  RooFormulaVar sigma2("sigma2","sigma*sigma_ratio",RooArgList(sigma, sigma_ratio));
  RooCBShape cb2("cb2","Crystal Ball 2 PDF",fMuMuM,m02,sigma2,alpha,n_cb) ;

  // --------------------------------------------------------------------------------
  // Fit the real data invariant mass spectra
  // ---------------------------
  // Compute number background events in an ivariant mass fit without a pt cut
  // ---Invarant mass fit function with no pt cut
  RooRealVar N_CB_PtAll("N_Jpsi_PtAll","Number of CB1 events PtAll",0.5*N_Evnts_m_PtAll,0,N_Evnts_m_PtAll);
  RooRealVar N_CB2_PtAll("Npsi2PtAll","Number of CB2 events PtAll",0.1*N_Evnts_m_PtAll,0,N_Evnts_m_PtAll);
  RooRealVar N_BG_PtAll("N_{bg}_PtAll","Number of BG events PtAll",0.5*N_Evnts_m_PtAll,0,N_Evnts_m_PtAll);

  RooAddPdf M_fit_func_PtAll("M_fit_func_PtAll","Crystal Ball 1 and 2 plus Background PDF for Pt All", RooArgList(Bkgd,cb,cb2), RooArgList(N_BG_PtAll,N_CB_PtAll,N_CB2_PtAll));
  // ---Fit the pt all data
  RooFitResult* fit_m_PtAll = M_fit_func_PtAll.fitTo(*data_m_PtAll,Extended(kTRUE),Range("m_fitRange"),Save());

  // -------------------------------------------------------------------------------- 
  // Plot the pt all fit
  TCanvas *cM_PtAll = new TCanvas("cM_PtAll","cM_PtAll",1600,800);
  cM_PtAll->cd();

  RooPlot* frame_m_PtAll = fMuMuM.frame(Title("Mass fit")) ; 
  data_m_PtAll->plotOn(frame_m_PtAll,Name("data_m"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
  M_fit_func_PtAll.plotOn(frame_m_PtAll,Name("M_fit_func_PtAll"),LineColor(kBlack), LineWidth(1)) ;
  M_fit_func_PtAll.plotOn(frame_m_PtAll,Name("cb"), Components(cb), LineColor(kRed+1),LineWidth(1));
  M_fit_func_PtAll.plotOn(frame_m_PtAll,Name("cb2"), Components(cb2), LineColor(kGreen+2), LineWidth(1));
  M_fit_func_PtAll.plotOn(frame_m_PtAll,Name("Bkgd"), Components(Bkgd), LineColor(kBlue+1), LineWidth(1));
  frame_m_PtAll->Draw();

  char name_cM_PtAll[120];
  sprintf(name_cM_PtAll,"Fit_param/%s/Mass_data_%s_%s_%.2f_%.2f_PtAll.png",trigger.Data(),ADveto.Data(),ZDC_class.Data(),abs(y_min),abs(y_max));
  cM_PtAll->SaveAs(name_cM_PtAll);
  
  // -------------------------------------------------------------------------------- 
  // Creare variables coming from pt all fit
  N_CB_PtAll.setConstant(kTRUE);
  N_CB2_PtAll.setConstant(kTRUE);
  N_BG_PtAll.setConstant(kTRUE);
  // ---Define variabels
  Double_t N_bkgd_PtAll_JPsi_mass_range[2];
  // ---Create integrals
  RooAbsReal *I_bkgd_PtAll_JPsi_mass_range = Bkgd.createIntegral(fMuMuM,NormSet(fMuMuM),Range("JPsiMassRange"));
  // ---Compute the number of events
  N_bkgd_PtAll_JPsi_mass_range[0] =  I_bkgd_PtAll_JPsi_mass_range->getVal()*N_BG_PtAll.getVal();
  N_bkgd_PtAll_JPsi_mass_range[1] =  I_bkgd_PtAll_JPsi_mass_range->getVal()*N_BG_PtAll.getError();  

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Pt Cut
    // // Create PDF to fit MC
    // // ---PDF to fit MC gamma gamma background 
    // RooRealVar N_MC_BG("N_MC_{bg}","Number of MC BG events",1.0*N_Evnts_TwoGammaToMuMedium_m,0,3*N_Evnts_TwoGammaToMuMedium_m);
    // RooAddPdf Bkgd_mc("Bkgd_mc","Background PDF for MC fit", RooArgList(Bkgd), RooArgList(N_MC_BG));

    // // ---PDF to fit MC Coherent Psi'
    // RooRealVar N_MC_CB2("N_MC_{cb2}","Number of MC Psi' events",1.0*N_Evnts_CohPsi2sToMu_m,0,3*N_Evnts_CohPsi2sToMu_m);
    // RooAddPdf cb2_mc("cb2_mc","CB2 PDF for MC fit", RooArgList(cb2MC), RooArgList(N_MC_CB2));
    
    // // ---PDF to fit MC Coherent J/Psi
    // RooRealVar N_MC_CB("N_MC_{cb}","Number of MC J/Psi events",1.0*N_Evnts_CohJpsiToMu_m,0,3*N_Evnts_CohJpsiToMu_m);
    // RooAddPdf cb_mc("cb_mc","CB PDF for MC fit", RooArgList(cbMC), RooArgList(N_MC_CB));

    // lambda.setConstant(kFALSE);
    // a2.setConstant(kFALSE);
    // a3.setConstant(kFALSE);
    // a4.setConstant(kFALSE);
    // n_cb.setConstant(kFALSE);
    // alpha.setConstant(kFALSE);
    // sigmaMC.setConstant(kFALSE);
    // sigma2MC.setConstant(kFALSE);

    // if (MC_fit_str.Contains("Compute_MC_fit_param"))
    // {
    //   // ---Fit MC data to fix gamma gamma background parameters
    //   RooFitResult* fit_bkgd_PtCut = Bkgd_mc.fitTo(*dataMC_TwoGammaToMuMedium_m_PtCut,Extended(kTRUE),Range("m_fitRange"),Save());
    //   // ---Fit MC data to get the sigma ratio of tje Psi' and J/Psi
    //   RooFitResult* fit_cb2_PtCut = cb2_mc.fitTo(*dataMC_CohPsi2sToMu_m_PtCut,Extended(kTRUE),Range("CB2_fitRange"),Save());
    //   // ---Fit MC data to fix J/Psi crystall ball parameters
    //   RooFitResult* fit_cb_PtCut = cb_mc.fitTo(*dataMC_CohJpsiToMu_m_PtCut,Extended(kTRUE),Range("CB_fitRange"),Save());

    //   // --------------------------------------------------------------------------------
    //   // ---Plotting the MC fits

    //   // ---gamma gamma background 
    //   TCanvas *cMC_PtCut_bkgd = new TCanvas("cMC_PtCut_bkgd","cMC_PtCut_bkgd",1600,800);
    //   cMC_PtCut_bkgd->cd();

    //   RooPlot* frame_MC_PtCut_bkgd = fMuMuM.frame(Title("Mass fit")) ; 
    //   dataMC_TwoGammaToMuMedium_m_PtCut->plotOn(frame_MC_PtCut_bkgd,Name("MC_bkgd_PtCut"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
    //   Bkgd_mc.plotOn(frame_MC_PtCut_bkgd,Name("Bkgd_mc_PtCut"),LineColor(kBlack), LineWidth(1)) ;
    //   frame_MC_PtCut_bkgd->Draw();

    //   char name_cMC_PtCut_bkgd[120];
    //   sprintf(name_cMC_PtCut_bkgd,"Fit_param/%s/%s_%.2f_%.2f_MC_PtCut_bkgd.png",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
    //   cMC_PtCut_bkgd->SaveAs(name_cMC_PtCut_bkgd);

    //   // ---J/Psi crystall ball 
    //   TCanvas *cMC_PtCut_cb = new TCanvas("cMC_PtCut_cb","cMC_PtCut_cb",1600,800);
    //   cMC_PtCut_cb->SetLogy();

    //   RooPlot* frame_MC_PtCut_cb = fMuMuM.frame(Title("Mass fit")) ; 
    //   dataMC_CohJpsiToMu_m_PtCut->plotOn(frame_MC_PtCut_cb,Name("MC_cb_PtCut"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
    //   cb_mc.plotOn(frame_MC_PtCut_cb,Range("m_fitRange"),Name("cb_mc_PtCut"),LineColor(kBlack), LineWidth(1)) ;
    //   frame_MC_PtCut_cb->SetAxisRange(0.025,10000000,"Y");
    //   frame_MC_PtCut_cb->Draw();

    //   char name_cMC_PtCut_cb[120];
    //   sprintf(name_cMC_PtCut_cb,"Fit_param/%s/%s_%.2f_%.2f_MC_PtCut_cb.png",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
    //   cMC_PtCut_cb->SaveAs(name_cMC_PtCut_cb);

    //   // ---Psi' crystall ball 
    //   TCanvas *cMC_PtCut_cb2 = new TCanvas("cMC_PtCut_cb2","cMC_PtCut_cb2",1600,800);
    //   cMC_PtCut_cb2->cd();
    //   cMC_PtCut_cb2->SetLogy();

    //   RooPlot* frame_MC_PtCut_cb2 = fMuMuM.frame(Title("Mass fit")) ; 
    //   dataMC_CohPsi2sToMu_m_PtCut->plotOn(frame_MC_PtCut_cb2,Name("MC_cb2_PtCut"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
    //   cb2_mc.plotOn(frame_MC_PtCut_cb2,Range("m_fitRange"),Name("cb2_mc_PtCut"),LineColor(kBlack), LineWidth(1)) ;
    //   frame_MC_PtCut_cb2->SetAxisRange(0.025,10000000,"Y");
    //   frame_MC_PtCut_cb2->Draw();

    //   char name_cMC_PtCut_cb2[120];
    //   sprintf(name_cMC_PtCut_cb2,"Fit_param/%s/%s_%.2f_%.2f_MC_PtCut_cb2.png",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
    //   cMC_PtCut_cb2->SaveAs(name_cMC_PtCut_cb2);


    //   // --------------------------------------------------------------------------------
    //   // Write the values from the MC fits
    //   char FileName_MC_m_param_PtCut[120];
    //   sprintf(FileName_MC_m_param_PtCut,"Fit_param/%s/%s_%.2f_%.2f_MC_m_param_PtCut.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));

    //   ofstream MC_param_File_PtCut;
    //   MC_param_File_PtCut.open (FileName_MC_m_param_PtCut, ios::trunc);
    //   MC_param_File_PtCut << "-----------------------------------------------------------" << endl;
    //   MC_param_File_PtCut << "lambda_val = " << lambda.getVal() << endl;
    //   MC_param_File_PtCut << "a2_val = " << a2.getVal() << endl;
    //   MC_param_File_PtCut << "a3_val = " << a3.getVal() << endl;
    //   MC_param_File_PtCut << "a4_val = " << a4.getVal() << endl;
    //   MC_param_File_PtCut << "n_cb_val = " << n_cb.getVal() << endl;
    //   MC_param_File_PtCut << "alpha_val = " << alpha.getVal() << endl;
    //   MC_param_File_PtCut << "sigmaMC_val = " << sigmaMC.getVal() << endl;
    //   MC_param_File_PtCut << "sigma2MC_val = " << sigma2MC.getVal() << endl;
    //   MC_param_File_PtCut << "lambda_err = " << lambda.getError() << endl;
    //   MC_param_File_PtCut << "a2_err = " << a2.getError() << endl;
    //   MC_param_File_PtCut << "a3_err = " << a3.getError() << endl;
    //   MC_param_File_PtCut << "a4_err = " << a4.getError() << endl;
    //   MC_param_File_PtCut << "n_cb_err = " << n_cb.getError() << endl;
    //   MC_param_File_PtCut << "alpha_err = " << alpha.getError() << endl;
    //   MC_param_File_PtCut << "sigmaMC_err = " << sigmaMC.getError() << endl;
    //   MC_param_File_PtCut << "sigma2MC_err = " << sigma2MC.getError() << endl;
    //   MC_param_File_PtCut << "-----------------------------------------------------------" << endl;
    //   MC_param_File_PtCut << endl;  
    //   MC_param_File_PtCut.close();
    // }

    // // --------------------------------------------------------------------------------
    // // Load mass fit parameters from MC
    // if (MC_fit_str.Contains("Load_MC_fit_param"))
    // {
    //   LoadMCmassParam(y_min, y_max, ADveto, trigger, "PtCut",
    //                   lambda_val, lambda_err,
    //                   a2_val, a2_err,
    //                   a3_val, a3_err,
    //                   a4_val, a4_err,
    //                   n_cb_val, n_cb_err,
    //                   alpha_val, alpha_err,
    //                   sigmaMC_val, sigmaMC_err,
    //                   sigma2MC_val, sigma2MC_err);

    //   lambda.setVal(lambda_val);
    //   a2.setVal(a2_val);
    //   a3.setVal(a3_val);
    //   a4.setVal(a4_val);
    //   n_cb.setVal(n_cb_val);
    //   alpha.setVal(alpha_val);
    //   sigmaMC.setVal(sigmaMC_val);
    //   sigma2MC.setVal(sigma2MC_val);

    //   lambda.setError(lambda_err);
    //   a2.setError(a2_err);
    //   a3.setError(a3_err);
    //   a4.setError(a4_err);
    //   n_cb.setError(n_cb_err);
    //   alpha.setError(alpha_err);
    //   sigmaMC.setError(sigmaMC_err);
    //   sigma2MC.setError(sigma2MC_err);
    // }

    // lambda.setConstant(kTRUE);
    // a2.setConstant(kTRUE);
    // a3.setConstant(kTRUE);
    // a4.setConstant(kTRUE);
    // n_cb.setConstant(kTRUE);
    // alpha.setConstant(kTRUE);
    // sigmaMC.setConstant(kTRUE);
    // sigma2MC.setConstant(kTRUE);

  // ---Invarant mass fit function
  RooRealVar N_CB("N_Jpsi","Number of CB1 events",0.5*N_Evnts_m,0,N_Evnts_m);
  RooRealVar N_CB2("N_psi2","Number of CB2 events",0.1*N_Evnts_m,0,N_Evnts_m);
  RooRealVar N_BG("N_{bg}","Number of BG events",0.5*N_Evnts_m,0,N_Evnts_m);

  RooAddPdf M_fit_func("M_fit_func","Crystal Ball 1 and 2 plus Background PDF", RooArgList(Bkgd,cb,cb2), RooArgList(N_BG,N_CB,N_CB2));

  // ---Fit data with the pt cut
  RooFitResult* fit_m = M_fit_func.fitTo(*data_m,Extended(kTRUE),Range("m_fitRange"),Save());

  // -------------------------------------------------------------------------------- 
  // Plot the pt cut fit
  TCanvas *cM_PtCut = new TCanvas("cM_PtCut","cM_PtCut",1600,800);
  cM_PtCut->cd();

  RooPlot* frame_m_PtCut = fMuMuM.frame(Title("Mass fit")) ; 
  data_m->plotOn(frame_m_PtCut,Name("data_m"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
  M_fit_func.plotOn(frame_m_PtCut,Name("M_fit_func"),LineColor(kBlack), LineWidth(1)) ;
  M_fit_func.plotOn(frame_m_PtCut,Name("cb"), Components(cb), LineColor(kRed+1),LineWidth(1));
  M_fit_func.plotOn(frame_m_PtCut,Name("cb2"), Components(cb2), LineColor(kGreen+2), LineWidth(1));
  M_fit_func.plotOn(frame_m_PtCut,Name("Bkgd"), Components(Bkgd), LineColor(kBlue+1), LineWidth(1));
  frame_m_PtCut->Draw();

  char name_cM_PtCut[120];
  sprintf(name_cM_PtCut,"Fit_param/%s/Mass_data_%s_%s_%.2f_%.2f_PtCut.png",trigger.Data(),ADveto.Data(),ZDC_class.Data(),abs(y_min),abs(y_max));
  cM_PtCut->SaveAs(name_cM_PtCut);
  
  if (MC_fit_str.Contains("Compute_MC_fit_param"))
  {
    cout<< "##########################################################################" << endl;
    cout<< "MC fit parameres computed and PDF templates created for " << ZDC_class << ", " << trigger << ", " << ADveto << ", " << y_min << " < y < " << y_max << endl;
    cout<< "##########################################################################" << endl;
    gROOT->ProcessLine(".q");
  }

  // -------------------------------------------------------------------------------- 
  // Creare variables coming from pt cut fit
  N_CB.setConstant(kTRUE);
  N_CB2.setConstant(kTRUE);
  N_BG.setConstant(kTRUE);

  // --------------------------------------------------------------------------------
  // Compute cross section ratio for Psi' and J/Psi
  if (ratio_str.Contains("RatioOnly")){ 
    RooFormulaVar ratio_PtAll("ratio_PtAll","Npsi2PtAll/N_Jpsi_PtAll*0.05961*Var_Eff_CohJpsiToMu_PtAll"
                                "/ ( 0.00800*Var_Eff_CohPsi2sToMu_PtAll-Npsi2PtAll/N_Jpsi_PtAll*0.61400*Var_Eff_CohPsi2sToMuPi_PtAll*0.05961 )",
                        RooArgList(N_CB2_PtAll, N_CB_PtAll, Var_Eff_CohJpsiToMu_PtAll, Var_Eff_CohPsi2sToMu_PtAll, Var_Eff_CohPsi2sToMuPi_PtAll));
    
    // RooFormulaVar ratio_PtCut("ratio_PtCut","N_{#psi'}/N_{J#psi}*0.05961*Var_Eff_CohJpsiToMu_PtCut"
    //                             "/ ( 0.00800*Var_Eff_CohPsi2sToMu_PtCut-N_{#psi'}/N_{J#psi}*0.61400*Var_Eff_CohPsi2sToMuPi_PtCut*0.05961 )",
    //                     RooArgList(N_CB2, N_CB, Var_Eff_CohJpsiToMu_PtCut, Var_Eff_CohPsi2sToMu_PtCut, Var_Eff_CohPsi2sToMuPi_PtCut));

    // ------------------------------------------------------------------------------------------------------------------------------------
    // Write the values needed for feed down correction computation to a file
    char FileName_ratio[120];
    sprintf(FileName_ratio,"Ratio/%s/%s_%.2f_%.2f_%s.txt",trigger.Data(), ADveto.Data(), abs(y_min),abs(y_max),ZDC_class.Data());

    cout << "##########################################################################" << endl;
    cout << N_CB2_PtAll.getVal() << endl;
    cout << N_CB_PtAll.getVal() << endl;
    cout << Var_Eff_CohJpsiToMu_PtAll.getVal() << endl;
    cout << Var_Eff_CohPsi2sToMu_PtAll.getVal() << endl;
    cout << Var_Eff_CohPsi2sToMuPi_PtAll.getVal() << endl;
    cout << ratio_PtAll.getVal() << endl;
    cout << ratio_PtAll.getPropagatedError(*fit_m_PtAll) << endl;
    cout << "##########################################################################" << endl;

    ofstream RatioFile;
    RatioFile.open (FileName_ratio, ios::trunc);
    RatioFile << "-----------------------------------------------------------" << endl;
    RatioFile << "ratio = " << ratio_PtAll.getVal() <<  endl;
    RatioFile << "ratio_err = " << ratio_PtAll.getPropagatedError(*fit_m_PtAll) <<  endl;
    RatioFile << "-----------------------------------------------------------" << endl;
    RatioFile << endl;  
    RatioFile.close();

    cout<< "##########################################################################" << endl;
    cout<< "Ratio Computed for " << ZDC_class << ", " << trigger << ", " << ADveto << ", " << y_min << " < y < " << y_max << endl;
    cout<< "##########################################################################" << endl;
    gROOT->ProcessLine(".q");
  }
  // --------------------------------------------------------------------------------
  // Load Psi' J/Psi cross section ratio
  double ratio;
  double ratio_err;

  // Values taken from https://alice-publications.web.cern.ch/system/files/draft/5085/2019-08-02-paper_v10.pdf
  ratio = 0.150;
  ratio_err = 0.0285;

  // LoadRatio(4.0, 2.5, "all", trigger, ADveto, ratio, ratio_err);

  RooRealVar var_ratio("var_ratio","var_ratio",ratio);
  var_ratio.setError(ratio_err);

  var_ratio.setConstant();
  // --------------------------------------------------------------------------------
  // Compute the feed down contribution
  RooFormulaVar f_D_PtAll("f_D_PtAll","var_ratio"
                                      "*Var_Eff_CohPsi2sToMuPi_PtAll*0.61400"
                                      "/Var_Eff_CohJpsiToMu_PtAll",
                          RooArgList(var_ratio, Var_Eff_CohJpsiToMu_PtAll, Var_Eff_CohPsi2sToMuPi_PtAll));

  RooFormulaVar f_D_PtCut("f_D_PtCut","var_ratio"
                                      "*Var_Eff_CohPsi2sToMuPi_PtCut*0.61400"
                                      "/Var_Eff_CohJpsiToMu_PtCut",
                          RooArgList(var_ratio, Var_Eff_CohJpsiToMu_PtCut, Var_Eff_CohPsi2sToMuPi_PtCut));

  // --------------------------------------------------------------------------------
  // Compute the Psi'/JPsi ratio
  RooFormulaVar R_Psi2_to_JPsi("R_Psi2_to_JPsi","N_psi2/N_Jpsi", RooArgList(N_CB2, N_CB));

  // --------------------------------------------------------------------------------
  // Compute number of events from the mass fit for different processes in given range
  // ---Define variabels
  Double_t N_bkgd_JPsi_mass_range[2];
  Double_t N_bkgd_Psi2S_mass_range[2];
  Double_t N_cb_JPsi[2];
  Double_t N_cb2_Psi2S[2];
  // ---Create integrals
  RooAbsReal *I_bkgd_JPsi_mass_range = Bkgd.createIntegral(fMuMuM,NormSet(fMuMuM),Range("JPsiMassRange"));
  RooAbsReal *I_bkgd_Psi2S_mass_range = Bkgd.createIntegral(fMuMuM,NormSet(fMuMuM),Range("Psi(2S)MassRange"));
  // RooAbsReal *I_cb_JPsi = cb.createIntegral(fMuMuM,NormSet(fMuMuM),Range("m_fitRange"));
  // RooAbsReal *I_cb2_Psi2S = cb2.createIntegral(fMuMuM,NormSet(fMuMuM),Range("m_fitRange"));
  // ---Compute the number of events
  N_cb_JPsi[0] = N_CB.getVal();
  N_cb_JPsi[1] =  N_CB.getError();
  N_cb2_Psi2S[0] = N_CB2.getVal();
  N_cb2_Psi2S[1] =  N_CB2.getError();
  N_bkgd_JPsi_mass_range[0] =  I_bkgd_JPsi_mass_range->getVal()*N_BG.getVal();
  N_bkgd_JPsi_mass_range[1] =  I_bkgd_JPsi_mass_range->getVal()*N_BG.getError();  
  N_bkgd_Psi2S_mass_range[0] =  I_bkgd_Psi2S_mass_range->getVal()*N_BG.getVal();
  N_bkgd_Psi2S_mass_range[1] =  I_bkgd_Psi2S_mass_range->getVal()*N_BG.getError();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Transverse momentum fit
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // --------------------------------------------------------------------------------
  //Set the normalisation for the pt fit templates
  // ---JPsi
  RooRealVar N_CohJpsiToMu("Jpsi_coh","number of coherent jpsi events",0.75*data_pt->numEntries(),0,data_pt->numEntries());   
  RooRealVar N_IncohJpsiToMu("Jpsi_incoh","number of incoherent jpsi events",0.1*data_pt->numEntries(),0*data_pt->numEntries(),data_pt->numEntries()); 
  // Psi'
  RooRealVar N_CohPsi2sToMu("#psi'to#mu_{coh}","number of coherent psi2s to mu events",0.75*data_pt->numEntries(),0,data_pt->numEntries());   
  RooRealVar N_IncohPsi2sToMu("#psi'to#mu_{incoh}","number of incoherent psi2s to mu events",0.1*data_pt->numEntries(),0,data_pt->numEntries()); 
  // ---Psi' feed down fixed by JPsi values * feed down coeficient
  RooFormulaVar N_CohPsi2sToMuPi("psi2topi_coh","Jpsi_coh*f_D_PtAll",RooArgList(N_CohJpsiToMu, f_D_PtAll));
  RooFormulaVar N_IncohPsi2sToMuPi("psi2topi_incoh","Jpsi_incoh*f_D_PtAll",RooArgList(N_IncohJpsiToMu, f_D_PtAll));
  // ---Gamma gamma is fixed to background in the mass fit
  RooRealVar N_TwoGammaToMuMedium("#gamma#gamma","number of gg", N_bkgd_JPsi_mass_range[0]);
  N_TwoGammaToMuMedium.setConstant(kTRUE);
  // ---Disociative
  RooRealVar N_IncohJpsiToX("N_disoc","number of incoherent jpsi disociated signal",0.1*data_pt->numEntries(),0,data_pt->numEntries());

  // Create the model as the sum of the templates
  // --- For J/Psi
  RooAddPdf Pt_fit_func("Pt_fit_func","Sum of templates",RooArgList(*PDF_CohJpsiToMu,*PDF_IncohJpsiToMu, *PDF_CohPsi2sToMuPi, *PDF_IncohPsi2sToMuPi, *PDF_TwoGammaToMuMedium, PDF_IncohJpsiToX),
                                                         RooArgList(N_CohJpsiToMu,N_IncohJpsiToMu, N_CohPsi2sToMuPi, N_IncohPsi2sToMuPi, N_TwoGammaToMuMedium, N_IncohJpsiToX));
  // --- For Psi(2S)
  // RooAddPdf Pt_fit_func("Pt_fit_func","Sum of templates",RooArgList(*PDF_CohPsi2sToMu, *PDF_IncohPsi2sToMu, *PDF_TwoGammaToMuMedium, PDF_IncohJpsiToX), RooArgList(N_CohPsi2sToMu, N_IncohPsi2sToMu, N_TwoGammaToMuMedium, N_IncohJpsiToX));

  // --------------------------------------------------------------------------------
  // Fit the real data pt spectra
  RooFitResult* fit_pt;

  if (Disoc_fit_str.Contains("Compute_disoc_fit_param"))
  {
    fit_pt = Pt_fit_func.fitTo(*data_pt,Extended(kTRUE),Save(),Range(0,pt_range_max));
    // --------------------------------------------------------------------------------
    // Plot the pt fit
    TCanvas *Compute_MC_fit_param = new TCanvas("Compute_MC_fit_param","Compute_MC_fit_param",1600,800);
    Compute_MC_fit_param->cd();

    Compute_MC_fit_param->SetLogy();
    Compute_MC_fit_param->SetLeftMargin(0.14);

    RooPlot* frame_pt_param = fMuMuPt.frame(Title("Pt fit")) ; 
    data_pt->plotOn(frame_pt_param,Name("data_pt"),Binning(bin_pt),MarkerStyle(20),MarkerSize(0.9));
    Pt_fit_func.plotOn(frame_pt_param,Name("Pt_fit_func"),LineColor(kBlack), LineWidth(2)) ;
    Pt_fit_func.plotOn(frame_pt_param,Name("PDF_CohJpsiToMu"), Components(*PDF_CohJpsiToMu), LineColor(kBlue+1), LineWidth(2));
    Pt_fit_func.plotOn(frame_pt_param,Name("PDF_IncohJpsiToMu"),Components(*PDF_IncohJpsiToMu), LineColor(kRed+1), LineWidth(2));
    Pt_fit_func.plotOn(frame_pt_param,Name("PDF_CohPsi2sToMuPi"),Components(*PDF_CohPsi2sToMuPi), LineColor(kCyan+1), LineWidth(2));
    Pt_fit_func.plotOn(frame_pt_param,Name("PDF_IncohPsi2sToMuPi"),Components(*PDF_IncohPsi2sToMuPi), LineColor(kOrange), LineWidth(2));
    Pt_fit_func.plotOn(frame_pt_param,Name("PDF_IncohJpsiToX"),Components(PDF_IncohJpsiToX), LineColor(kMagenta+1), LineWidth(2));
    Pt_fit_func.plotOn(frame_pt_param,Name("PDF_TwoGammaToMuMedium"),Components(*PDF_TwoGammaToMuMedium), LineColor(kGreen+2), LineWidth(2));
    frame_pt_param->SetAxisRange(0,3,"X");
    frame_pt_param->SetAxisRange(0.025,11000000,"Y");
    frame_pt_param->Draw();


    // ---Legend
    TLegend *leg_pt_param = new TLegend(0.7,0.7,0.9,0.9);
    leg_pt_param->SetFillStyle(0);
    leg_pt_param->SetBorderSize(0);
    leg_pt_param->AddEntry("PDF_CohJpsiToMu","Coherent J/#psi", "L");
    leg_pt_param->AddEntry("PDF_IncohJpsiToMu","Incoherent J/#psi", "L");
    leg_pt_param->AddEntry("PDF_CohPsi2sToMuPi","Coh. #psi' #rightarrow  J/#psi", "L");
    leg_pt_param->AddEntry("PDF_IncohPsi2sToMuPi","Incoh. #psi' #rightarrow  J/#psi", "L");
    leg_pt_param->AddEntry("PDF_IncohJpsiToX","Nucleon disoc.", "L");
    leg_pt_param->AddEntry("PDF_TwoGammaToMuMedium","#gamma#gamma #rightarrow #mu#mu", "L");
    leg_pt_param->Draw();

    char name_cPt_param[120];
    sprintf(name_cPt_param,"Fit_param/%s/Pt_data_param_%s_%s_%.2f_%.2f.png",trigger.Data(),ADveto.Data(),ZDC_class.Data(),abs(y_min),abs(y_max));
    Compute_MC_fit_param->SaveAs(name_cPt_param);

    // --------------------------------------------------------------------------------
    // Write the values of paameters from pt fit
    char FileName_pt_param[120];
    sprintf(FileName_pt_param,"Fit_param/%s/%s_%s_%.2f_%.2f_pt_param.txt",trigger.Data(),ADveto.Data(),ZDC_class.Data(),abs(y_min),abs(y_max));

    ofstream pt_param_File;
    pt_param_File.open (FileName_pt_param, ios::trunc);
    pt_param_File << "-----------------------------------------------------------" << endl;
    pt_param_File << "b_val = " << b.getVal() << endl;
    pt_param_File << "n_val = " << n.getVal() << endl;
    pt_param_File << "b_err = " << b.getError() << endl;
    pt_param_File << "n_err = " << n.getError() << endl;
    pt_param_File << "-----------------------------------------------------------" << endl;
    pt_param_File << endl;  
    pt_param_File.close();

    cout<< "##########################################################################" << endl;
    cout<< "pt fit parameres computed for " << ZDC_class << ", " << trigger << ", " << ADveto << ", " << y_min << " < y < " << y_max << endl;
    cout<< "##########################################################################" << endl;

    gROOT->ProcessLine(".q");
  }

  // --------------------------------------------------------------------------------
  // Load pt fit parameters for disociative background
  double b_val, b_err;
  double n_val, n_err;

  if (Disoc_fit_str.Contains("Load_disoc_fit_param"))
  {
    LoadPtParam(y_min, y_max, ADveto, "all", trigger,
                b_val, b_err,
                n_val, n_err);
    // 1.79 3.58
    if(ADveto.Contains("ADvetoOff")){
      // b.setVal(b_val);
      // b.setError(b_err);

      n.setVal(n_val);
      n.setError(n_err);
    }
  }

  if(ADveto.Contains("ADvetoOff")) n.setConstant(kTRUE);

  fit_pt = Pt_fit_func.fitTo(*data_pt,Extended(kTRUE),Save(),Range(0,pt_range_max));

  N_CohJpsiToMu.setConstant(kTRUE);
  N_IncohJpsiToMu.setConstant(kTRUE);
  N_CohPsi2sToMu.setConstant(kTRUE);
  N_IncohPsi2sToMu.setConstant(kTRUE);
  N_IncohJpsiToX.setConstant(kTRUE);

  // --------------------------------------------------------------------------------
  // Compute number of events from the pt fit for different processes in given range 
  // ---Create integrals
  RooAbsReal *idisoc = PDF_IncohJpsiToX.createIntegral(fMuMuPt,NormSet(fMuMuPt),Range("JPsiptRange"));
  // ---Create variabels from integral values in given range
  RooRealVar ndisoc("ndisoc", "ndisoc", idisoc->getVal(), "" );

  // --------------------------------------------------------------------------------
  // Compute the incoherent contribution
  RooFormulaVar var_ncoh("var_ncoh","Jpsi_coh*ncoh",RooArgList(N_CohJpsiToMu,*ncoh) );
  RooFormulaVar var_ncoh2("var_ncoh2","psi2topi_coh*ncoh2",RooArgList(N_CohPsi2sToMuPi,*ncoh2) );
  RooFormulaVar var_nincoh("var_nincoh","Jpsi_incoh*nincoh",RooArgList(N_IncohJpsiToMu,*nincoh) );
  RooFormulaVar var_nincoh2("var_nincoh2","psi2topi_incoh*nincoh2",RooArgList(N_IncohPsi2sToMuPi,*nincoh2) );
  RooFormulaVar var_ndisoc("var_ndisoc","N_disoc*ndisoc",RooArgList(N_IncohJpsiToX,ndisoc) );

  RooFormulaVar f_I("f_I","(var_nincoh+var_nincoh2+var_ndisoc)/(var_ncoh+var_ncoh2)", RooArgList(var_nincoh, var_nincoh2, var_ndisoc, var_ncoh, var_ncoh2));  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plotting
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // ------------------------------------------------------------------------------------------------------------------------------------
  // Mass plot
  //---Create mass canvas
  TCanvas *cM = new TCanvas("cM","cM",1600,800);
  cM->Divide(2);  

  //-----------------------------------------------------------------------------------
  // Draw Correlation Matrix
  cM->cd(2); 
  TPad *cM_2 = new TPad("cM_2", "cM_2",0.001,0.001,0.999,0.999);
  // ---Pad option
  cM_2->SetRightMargin(0.15);
  cM_2->SetLeftMargin(0.15);
  cM_2->Draw();
  cM_2->cd();
  // ---TH2D correlation matrix
  TH2* hCorr_m = fit_m->correlationHist();
  hCorr_m->SetMarkerSize(1.2);
  hCorr_m->GetXaxis()->SetLabelSize(0.045);
  hCorr_m->GetYaxis()->SetLabelSize(0.045);
  hCorr_m->Draw("zcol,text");

  //-----------------------------------------------------------------------------------
  // Draw mass Histogram  
  cM->cd(1);
  TPad *cM_1 = new TPad("cM_1", "cM_1",0.001,0.001,0.999,0.999);
  // ---Pad option
  cM_1->SetLeftMargin(0.15);
  cM_1->SetRightMargin(0.01);
  cM_1->SetBottomMargin(0.12);
  cM_1->Draw(); 
  cM_1->cd();
  // ---TH1 mass spectrum
  RooPlot* frame_m = fMuMuM.frame(Title("Mass fit")) ; 
  data_m->plotOn(frame_m,Name("data_m"),Binning(bin_m),MarkerStyle(20),MarkerSize(0.5));
  M_fit_func.plotOn(frame_m,Name("M_fit_func"),LineColor(kBlack), LineWidth(1)) ;
  M_fit_func.plotOn(frame_m,Name("cb"), Components(cb), LineColor(kRed+1),LineWidth(1));
  M_fit_func.plotOn(frame_m,Name("cb2"), Components(cb2), LineColor(kGreen+2), LineWidth(1));
  M_fit_func.plotOn(frame_m,Name("Bkgd"), Components(Bkgd), LineColor(kBlue+1), LineWidth(1));
  // ---Frame mass option
  frame_m->GetXaxis()->SetTitle("#it{m_{#mu^{+}#mu^{-}}} (GeV/#it{c}^{2})");
  frame_m->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c^{2}} ", (m_range_max-m_range_min)/n_bins_m*1000));
  frame_m->GetYaxis()->SetTitleOffset(1.6);
  frame_m->SetAxisRange(0,N_cb_JPsi[0]/3,"Y");
  frame_m->Draw();
  //---Calculate chi2 mass
  Double_t chi2_m = frame_m->chiSquare("M_fit_func","data_m",fit_m->floatParsFinal().getSize());
  // ---mass plot description
  TLatex * text_m1 = new TLatex (2.2,0.95*frame_m->GetMaximum(),"This thesis");
  text_m1->Draw();
  TLatex * text_m2 = new TLatex (3.5,0.95*frame_m->GetMaximum(),"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
  text_m2->Draw();
  TLatex * text_m3 = new TLatex (3.5,0.89*frame_m->GetMaximum(),Form("#bf{%.2f < y < %.2f, #it{p}_{T} < %.2f}",y_min,y_max, pt_cut));
  text_m3->Draw();
  TLatex * text_m4 = new TLatex (2.2,0.89*frame_m->GetMaximum(),Form("#bf{ZDC: %s}",ZDC_class.Data()));
  text_m4->Draw();
  // TLatex * text_m5 = new TLatex (2.2,0.81*frame_m->GetMaximum(),Form("#bf{%s}",ADveto.Data()));
  // text_m5->Draw();
  // TLatex * text_m6 = new TLatex (2.2,0.81*frame_m->GetMaximum(),Form("#bf{f_{D} = %.1f%% #pm %.1f%%}", 100*f_D_PtCut.getVal(), 100*f_D_PtCut.getVal()*ratio_err/ratio ));
  // text_m6->Draw();
  // ---Legend mass
  TLegend *leg_m = new TLegend(0.41,0.55,0.85,0.80);
  leg_m->SetFillStyle(0);
  leg_m->SetBorderSize(0);
  leg_m->SetTextSize(0.038);
  // leg_m->AddEntry("data_m","Data", "P");
  // leg_m->AddEntry("M_fit_func","Pt_fit_func","L");
  // leg_m->AddEntry("cb","Signal J/#psi","L");
  // leg_m->AddEntry("cb2","Signal #psi(2S)","L");
  // leg_m->AddEntry("Bkgd","Background","L");
  leg_m->AddEntry((TObject*)0,Form("#chi^{2}/dof: %.3f",chi2_m),"");
  leg_m->AddEntry((TObject*)0,Form("N_{J/#psi} = %3.0f #pm %3.0f",N_cb_JPsi[0],N_cb_JPsi[1]),"");
  leg_m->AddEntry((TObject*)0,Form("N_{#psi'} = %3.0f #pm %3.0f",N_cb2_Psi2S[0],N_cb2_Psi2S[1]),"");
  // leg_m->AddEntry((TObject*)0,Form("N_{#psi'}/N_{J/#psi} = %3.4f #pm %3.4f",R_Psi2_to_JPsi.getVal(),R_Psi2_to_JPsi.getPropagatedError(*fit_m)),"");
  // leg_m->AddEntry((TObject*)0,Form("N_{bg}(2.1,6.0) = %3.0f #pm %3.0f", N_BG.getVal(),N_BG.getError()),"");
  // leg_m->AddEntry((TObject*)0,Form("N_{bg}(%3.2f,%3.2f) = %3.0f #pm %3.0f", m_min_jpsi, m_max_jpsi, N_bkgd_JPsi_mass_range[0],N_bkgd_JPsi_mass_range[1]),"");
  // leg_m->AddEntry((TObject*)0,Form("N_{bg}(%3.2f,%3.2f) = %3.0f #pm %3.0f", m_min_psi2, m_max_psi2, N_bkgd_Psi2S_mass_range[0],N_bkgd_Psi2S_mass_range[1]),"");
  leg_m->AddEntry((TObject*)0,Form("# Entries = %.0i",data_m->numEntries()),"");
  leg_m->Draw();
  //---Save as pdf
  char massFit[120];
  sprintf(massFit,"Results/%s/%s_%s_Mass_fit_%.2f_y_%.2f.pdf",trigger.Data(),ADveto.Data(),ZDC_class.Data(),y_min,y_max);
  cM->SaveAs(massFit);

  // ------------------------------------------------------------------------------------------------------------------------------------
  // Pt plot
  //---Create pt canvas
  TCanvas *cPt = new TCanvas("cPt","cPt",1600,800);
  cPt->Divide(2);  

  //-----------------------------------------------------------------------------------
  // Draw Correlation Matrix
  cPt->cd(2); 
  TPad *cPt2 = new TPad("cPt2", "cPt2",0.001,0.001,0.999,0.999);
  // ---Pad option
  cPt2->SetRightMargin(0.15);
  cPt2->SetLeftMargin(0.15);
  cPt2->Draw();
  cPt2->cd();
  // ---TH2D correlation matrix
  TH2* hCorr_pt = fit_pt->correlationHist();
  hCorr_pt->SetMarkerSize(1.6);
  hCorr_pt->GetXaxis()->SetLabelSize(0.045);
  hCorr_pt->GetYaxis()->SetLabelSize(0.045);
  hCorr_pt->Draw("zcol,text");

  //-----------------------------------------------------------------------------------
  // Draw pt Histogram  
  cPt->cd(1);
  TPad *cPt1 = new TPad("cPt1", "cPt1",0.001,0.001,0.999,0.999);
  // ---Pad option
  cPt1->SetLogy();
  cPt1->SetLeftMargin(0.14);
  cPt1->SetRightMargin(0.01);
  cPt1->SetBottomMargin(0.12);
  cPt1->Draw(); 
  cPt1->cd();
  // ---TH1 pt spectrum
  RooPlot* frame_pt = fMuMuPt.frame(Title("Pt fit")) ;  
  data_pt->plotOn(frame_pt,Name("data_pt"),Binning(bin_pt),MarkerStyle(20),MarkerSize(0.9));
  Pt_fit_func.plotOn(frame_pt,Name("Pt_fit_func"),LineColor(kBlack), LineWidth(2)) ;
  // ---Psi(2S) 
  // Pt_fit_func.plotOn(frame_pt,Name("PDF_CohPsi2sToMu"), Components(PDF_CohPsi2sToMu), LineColor(kBlue+1), LineWidth(2));
  // Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohPsi2sToMu"),Components(PDF_IncohPsi2sToMu), LineColor(kRed+1), LineWidth(2));
  // ---J/Psi
  Pt_fit_func.plotOn(frame_pt,Name("PDF_CohJpsiToMu"), Components(*PDF_CohJpsiToMu), LineColor(kBlue+1), LineWidth(2));
  Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohJpsiToMu"),Components(*PDF_IncohJpsiToMu), LineColor(kRed+1), LineWidth(2));
  Pt_fit_func.plotOn(frame_pt,Name("PDF_CohPsi2sToMuPi"),Components(*PDF_CohPsi2sToMuPi), LineColor(kCyan+1), LineWidth(2));
  Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohPsi2sToMuPi"),Components(*PDF_IncohPsi2sToMuPi), LineColor(kOrange), LineWidth(2));
  Pt_fit_func.plotOn(frame_pt,Name("PDF_IncohJpsiToX"),Components(PDF_IncohJpsiToX), LineColor(kMagenta+1), LineWidth(2));
  Pt_fit_func.plotOn(frame_pt,Name("PDF_TwoGammaToMuMedium"),Components(*PDF_TwoGammaToMuMedium), LineColor(kGreen+2), LineWidth(2));
  // ---Frame pt option
  frame_pt->SetAxisRange(0,3,"X");
  frame_pt->SetAxisRange(0.025,11000000,"Y");
  frame_pt->GetXaxis()->SetTitle("Dimuon #it{p}_{T} (GeV/#it{c})");
  frame_pt->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c} ", (pt_range_max-pt_range_min)/n_bins_pt*1000));
  frame_pt->GetYaxis()->SetTitleOffset(1.3);
  frame_pt->Draw();
  //---Calculate chi2 pt
  Double_t chi2_pt = frame_pt->chiSquare("Pt_fit_func","data_pt",fit_pt->floatParsFinal().getSize()); 
  // ---pt plot description
  TLatex * text_pt1 = new TLatex (0.1,0.38*frame_pt->GetMaximum(),"This thesis");
  text_pt1->Draw();
  TLatex * text_pt2 = new TLatex (0.77,0.38*frame_pt->GetMaximum(),"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
  text_pt2->Draw();
  TLatex * text_pt3 = new TLatex (0.77,0.11*frame_pt->GetMaximum(),Form("#bf{%.2f < y < %.2f, %.2f < #it{m_{#mu^{+}#mu^{-}}} < %.2f}",y_min,y_max, m_min_jpsi, m_max_jpsi));
  text_pt3->Draw();
  TLatex * text_pt4 = new TLatex (0.1,0.11*frame_pt->GetMaximum(),Form("#bf{ZDC: %s}",ZDC_class.Data()));
  text_pt4->Draw();
  // TLatex * text_pt5 = new TLatex (0.1,0.03*frame_pt->GetMaximum(),Form("#bf{%s}",ADveto.Data()));
  // text_pt5->Draw();
  // TLatex * text_pt6 = new TLatex (0.1,0.03*frame_pt->GetMaximum(),Form("#bf{f_{I} = %.2f%% #pm %.2f%%}",100*f_I.getVal(), 100*f_I.getPropagatedError(*fit_pt) ));
  // text_pt6->Draw();
  // ---Legend pt
  TLegend *leg_pt = new TLegend(0.6,0.45,0.95,0.79);
  leg_pt->SetFillStyle(0);
  leg_pt->SetBorderSize(0);
  leg_pt->SetTextSize(0.04);
  leg_pt->AddEntry("data_pt","ALICE Data", "P");
  leg_pt->AddEntry("Pt_fit_func",Form("Fit: #chi^{2}/dof= %.3f",chi2_pt),"L");
  // ------Psi(2S)
  // leg_pt->AddEntry("PDF_CohPsi2sToMu","Coherent #psi(2S)", "L");
  // leg_pt->AddEntry("PDF_IncohPsi2sToMu","Incoherent #psi(2S)", "L");
  // ------J/Psi
  leg_pt->AddEntry("PDF_CohJpsiToMu","Coherent J/#psi", "L");
  leg_pt->AddEntry("PDF_IncohJpsiToMu","Incoherent J/#psi", "L");
  leg_pt->AddEntry("PDF_CohPsi2sToMuPi","Coh. #psi' #rightarrow  J/#psi", "L");
  leg_pt->AddEntry("PDF_IncohPsi2sToMuPi","Incoh. #psi' #rightarrow  J/#psi", "L");
  leg_pt->AddEntry("PDF_IncohJpsiToX","Nucleon disoc.", "L");
  leg_pt->AddEntry("PDF_TwoGammaToMuMedium","#gamma#gamma #rightarrow #mu#mu", "L");
  // leg_pt->AddEntry((TObject*)0,Form("#chi^{2}/NDF: %.3f",chi2_pt),"");
  leg_pt->Draw();
  //---Save as pdf
  char ptFit[120];
  sprintf(ptFit,"Results/%s/%s_%s_Pt_fit_%.2f_y_%.2f.pdf",trigger.Data(),ADveto.Data(),ZDC_class.Data(),y_min,y_max);
  cPt->SaveAs(ptFit);

  // ------------------------------------------------------------------------------------------------------------------------------------
  // Write the values needed for cross section computation to a file
  char FileName_CS[120];
  sprintf(FileName_CS,"Cross_section/%s/%s_%.2f_%.2f_%s.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max),ZDC_class.Data());

  ofstream CSFile;
  CSFile.open (FileName_CS, ios::trunc);
  CSFile << "-----------------------------------------------------------" << endl;
  CSFile << "N_JPsi = " << N_cb_JPsi[0] <<  endl;
  CSFile << "N_JPsi_err = " << N_cb_JPsi[1] <<  endl;
  CSFile << "f_D = " << f_D_PtCut.getVal() <<  endl;
  CSFile << "f_D_err = " << f_D_PtCut.getVal()*ratio_err/ratio <<  endl;
  CSFile << "f_I = " << f_I.getVal() <<  endl;
  CSFile << "f_I_err = " << f_I.getPropagatedError(*fit_pt) <<  endl;
  CSFile << "eff = " << Eff_CohJpsiToMu_PtCut <<  endl;
  CSFile << "-----------------------------------------------------------" << endl;
  CSFile << endl;  
  CSFile.close();

  cout<< "##########################################################################" << endl;
  cout << "Yield computed for " << ZDC_class << ", " << trigger << ", " << ADveto << ", " << y_min << " < y < " << y_max << endl;
  cout<< "##########################################################################" << endl;

  // cout << "I_bkgd_PtAll_JPsi_mass_range:" <<  I_bkgd_PtAll_JPsi_mass_range->getVal() << endl;
  // cout << "N_BG_PtAll:" <<  N_BG_PtAll.getVal() << endl;
  // cout << "N_bkgd_PtAll_JPsi_mass_range[0]:" <<  N_bkgd_PtAll_JPsi_mass_range[0] << endl;

  // cout << "I_bkgd_JPsi_mass_range:" <<  I_bkgd_JPsi_mass_range->getVal() << endl;
  // cout << "N_BG:" <<  N_BG.getVal() << endl;
  // cout << "N_bkgd_JPsi_mass_range[0]:" <<  N_bkgd_JPsi_mass_range[0] << endl;

  // cout << "N_CohJpsiToMu: " << N_CohJpsiToMu.getVal() << endl;
  // cout << "N_IncohJpsiToMu: " << N_IncohJpsiToMu.getVal() << endl;
  // cout << "N_CohPsi2sToMu: " << N_CohPsi2sToMu.getVal() << endl;
  // cout << "N_IncohPsi2sToMu: " << N_IncohPsi2sToMu.getVal() << endl;
  // cout << "N_CohPsi2sToMuPi: " << N_CohPsi2sToMuPi.getVal() << endl;
  // cout << "N_IncohPsi2sToMuPi: " << N_IncohPsi2sToMuPi.getVal() << endl;
  // cout << "N_TwoGammaToMuMedium: " << N_TwoGammaToMuMedium.getVal() << endl;
  // cout << "N_IncohJpsiToX: " << N_IncohJpsiToX.getVal() << endl;

  // ------------------------------------------------------------------------------------------------------------------------------------
  // Deleting
  if (dataIN18r) {delete dataIN18r;}
  if (dataIN18q) {delete dataIN18q;}
  if (dataIN15o) {delete dataIN15o;}
  if (dataIN18l7_CohJpsiToMu) {delete dataIN18l7_CohJpsiToMu;}
  if (dataIN18l7_CohPsi2sToMu) {delete dataIN18l7_CohPsi2sToMu;}
  if (dataIN18l7_CohPsi2sToMuPi) {delete dataIN18l7_CohPsi2sToMuPi;}
  if (dataIN18l7_IncohJpsiToMu) {delete dataIN18l7_IncohJpsiToMu;}
  if (dataIN18l7_IncohPsi2sToMu) {delete dataIN18l7_IncohPsi2sToMu;}
  if (dataIN18l7_IncohPsi2sToMuPi) {delete dataIN18l7_IncohPsi2sToMuPi;}
  if (dataIN18l7_TwoGammaToMuMedium) {delete dataIN18l7_TwoGammaToMuMedium;}
  if (dataIN16b2_CohJpsiToMu) {delete dataIN16b2_CohJpsiToMu;}
  if (dataIN16b2_CohPsi2sToMu) {delete dataIN16b2_CohPsi2sToMu;}
  if (dataIN16b2_CohPsi2sToMuPi) {delete dataIN16b2_CohPsi2sToMuPi;}
  if (dataIN16b2_IncohJpsiToMu) {delete dataIN16b2_IncohJpsiToMu;}
  if (dataIN16b2_IncohPsi2sToMu) {delete dataIN16b2_IncohPsi2sToMu;}
  if (dataIN16b2_IncohPsi2sToMuPi) {delete dataIN16b2_IncohPsi2sToMuPi;}
  if (dataIN16b2_TwoGammaToMuMedium) {delete dataIN16b2_TwoGammaToMuMedium;}
  if (cM) {delete cM;}
  if (text_m1) {delete text_m1;}
  if (text_m2) {delete text_m2;}
  if (text_m3) {delete text_m3;}
  // if (text_m4) {delete text_m4;}
  if (leg_m) {delete leg_m;}
  if (cPt) {delete cPt;}
  if (text_pt1) {delete text_pt1;}
  if (text_pt2) {delete text_pt2;}
  if (text_pt3) {delete text_pt3;}
  // if (text_pt4) {delete text_pt4;}
  if (leg_pt) {delete leg_pt;}
}

//////////////////////////////////////////////////////////////////////
// Do J/Psi analysis
//////////////////////////////////////////////////////////////////////

void AnalysisJPsi(Int_t y_bin, TString ZDC_class, TString ADveto, TString trigger, TString ratio_str, TString MC_fit_str, TString Disoc_fit_str, TString PDF_str)
{
  DoFit(y_bins_min[y_bin], y_bins_max[y_bin], ZDC_class, ADveto, trigger, ratio_str, MC_fit_str, Disoc_fit_str, PDF_str);
  // gROOT->ProcessLine(".q");
}
