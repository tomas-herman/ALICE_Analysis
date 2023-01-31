// Efficiency calculation
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
#include <TEfficiency.h>

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
#include "SetCuts.C"

///////////////////////////////////////////////////////////
// Compute the efficency
///////////////////////////////////////////////////////////

void ComputeEfficiency(Double_t y_min, Double_t y_max, TString ADveto, TString trigger)
{
  // --------------------------------------------------------------------------------
  // Define cuts
  SetCuts(y_min, y_max, "all", ADveto);

  // --------------------------------------------------------------------------------
  // Load MC data
  RooArgSet Arguments(fMuMuM,fMuMuY,fMuMuPt,fRunNum,fADADecision,fADCDecision,fV0ADecision,fV0CDecision,fV0CFiredCells);
  Arguments.add( RooArgSet(fIsZNAFired,fIsZNCFired) ); // You can only add 9 arguments at a time

  RooArgSet ArgumentsMC(fMCMuMuM,fMCMuMuY,fMCMuMuPt,fMCRunNum);

  // ---LHC18l7
  // ------Reconstructed
  RooDataSet *dataIN18l7_Rec_CohJpsiToMu = new RooDataSet ("dataIN18l7_Rec_CohJpsiToMu", "dataIN18l7_Rec_CohJpsiToMu", Arguments);
  RooDataSet *dataIN18l7_Rec_CohPsi2sToMu = new RooDataSet ("dataIN18l7_Rec_CohPsi2sToMu", "dataIN18l7_Rec_CohPsi2sToMu", Arguments);
  RooDataSet *dataIN18l7_Rec_CohPsi2sToMuPi = new RooDataSet ("dataIN18l7_Rec_CohPsi2sToMuPi", "dataIN18l7_Rec_CohPsi2sToMuPi", Arguments);

  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/569_20200509-1838-LHC18l7/CohJpsiToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN18l7_Rec_CohJpsiToMu);
  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/569_20200509-1838-LHC18l7/CohPsi2sToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN18l7_Rec_CohPsi2sToMu);
  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/569_20200509-1838-LHC18l7/CohPsi2sToMuPi_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN18l7_Rec_CohPsi2sToMuPi);
  
  // ------Generated
  RooDataSet *dataIN18l7_Gen_CohJpsiToMu = new RooDataSet ("dataIN18l7_Gen_CohJpsiToMu", "dataIN18l7_Gen_CohJpsiToMu", ArgumentsMC);
  RooDataSet *dataIN18l7_Gen_CohPsi2sToMu = new RooDataSet ("dataIN18l7_Gen_CohPsi2sToMu", "dataIN18l7_Gen_CohPsi2sToMu", ArgumentsMC);
  RooDataSet *dataIN18l7_Gen_CohPsi2sToMuPi = new RooDataSet ("dataIN18l7_Gen_CohPsi2sToMuPi", "dataIN18l7_Gen_CohPsi2sToMuPi", ArgumentsMC);

  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/569_20200509-1838-LHC18l7/CohJpsiToMu_AnalysisResults.root", "Gen", "NanoMUONScalingOn", dataIN18l7_Gen_CohJpsiToMu);
  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/569_20200509-1838-LHC18l7/CohPsi2sToMu_AnalysisResults.root", "Gen", "NanoMUONScalingOn", dataIN18l7_Gen_CohPsi2sToMu);
  LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/569_20200509-1838-LHC18l7/CohPsi2sToMuPi_AnalysisResults.root", "Gen", "NanoMUONScalingOn", dataIN18l7_Gen_CohPsi2sToMuPi);

  // ---LHC16b2
  // ------Reconstructed
  RooDataSet *dataIN16b2_Rec_CohJpsiToMu = new RooDataSet ("dataIN16b2_Rec_CohJpsiToMu", "dataIN16b2_Rec_CohJpsiToMu", Arguments);
  RooDataSet *dataIN16b2_Rec_CohPsi2sToMu = new RooDataSet ("dataIN16b2_Rec_CohPsi2sToMu", "dataIN16b2_Rec_CohPsi2sToMu", Arguments);
  RooDataSet *dataIN16b2_Rec_CohPsi2sToMuPi = new RooDataSet ("dataIN16b2_Rec_CohPsi2sToMuPi", "dataIN16b2_Rec_CohPsi2sToMuPi", Arguments);

  if (trigger.Contains("CMUP11")){
    LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/568_20200509-1837-LHC16b2/CohJpsiToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_Rec_CohJpsiToMu);
    LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/568_20200509-1837-LHC16b2/CohPsi2sToMu_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_Rec_CohPsi2sToMu);
    LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/568_20200509-1837-LHC16b2/CohPsi2sToMuPi_AnalysisResults.root", "MCRec", "NanoMUONScalingOn", dataIN16b2_Rec_CohPsi2sToMuPi);
  }
  
  // ------Generated
  RooDataSet *dataIN16b2_Gen_CohJpsiToMu = new RooDataSet ("dataIN16b2_Gen_CohJpsiToMu", "dataIN16b2_Gen_CohJpsiToMu", ArgumentsMC);
  RooDataSet *dataIN16b2_Gen_CohPsi2sToMu = new RooDataSet ("dataIN16b2_Gen_CohPsi2sToMu", "dataIN16b2_Gen_CohPsi2sToMu", ArgumentsMC);
  RooDataSet *dataIN16b2_Gen_CohPsi2sToMuPi = new RooDataSet ("dataIN16b2_Gen_CohPsi2sToMuPi", "dataIN16b2_Gen_CohPsi2sToMuPi", ArgumentsMC);

  if (trigger.Contains("CMUP11")){
    LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/568_20200509-1837-LHC16b2/CohJpsiToMu_AnalysisResults.root", "Gen", "NanoMUONScalingOn", dataIN16b2_Gen_CohJpsiToMu);
    LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/568_20200509-1837-LHC16b2/CohPsi2sToMu_AnalysisResults.root", "Gen", "NanoMUONScalingOn", dataIN16b2_Gen_CohPsi2sToMu);
    LoadData("/mnt/data/Data_processing/Processed_data/TrainResults/568_20200509-1837-LHC16b2/CohPsi2sToMuPi_AnalysisResults.root", "Gen", "NanoMUONScalingOn", dataIN16b2_Gen_CohPsi2sToMuPi);
  }

  // --------------------------------------------------------------------------------
  // Merge MC data sets
  // ------Reconstructed
  RooDataSet *dataINMC_Rec_CohJpsiToMu = (RooDataSet*) dataIN18l7_Rec_CohJpsiToMu;
    if (trigger.Contains("CMUP11")){ dataINMC_Rec_CohJpsiToMu->append(*dataIN16b2_Rec_CohJpsiToMu); }

  RooDataSet *dataINMC_Rec_CohPsi2sToMu = (RooDataSet*) dataIN18l7_Rec_CohPsi2sToMu;
    if (trigger.Contains("CMUP11")){ dataINMC_Rec_CohPsi2sToMu->append(*dataIN16b2_Rec_CohPsi2sToMu); }

  RooDataSet *dataINMC_Rec_CohPsi2sToMuPi = (RooDataSet*) dataIN18l7_Rec_CohPsi2sToMuPi;
    if (trigger.Contains("CMUP11")){ dataINMC_Rec_CohPsi2sToMuPi->append(*dataIN16b2_Rec_CohPsi2sToMuPi); }

  // ------Generated
  RooDataSet *dataINMC_Gen_CohJpsiToMu = (RooDataSet*) dataIN18l7_Gen_CohJpsiToMu;
    if (trigger.Contains("CMUP11")){ dataINMC_Gen_CohJpsiToMu->append(*dataIN16b2_Gen_CohJpsiToMu); }

  RooDataSet *dataINMC_Gen_CohPsi2sToMu = (RooDataSet*) dataIN18l7_Gen_CohPsi2sToMu;
    if (trigger.Contains("CMUP11")){ dataINMC_Gen_CohPsi2sToMu->append(*dataIN16b2_Gen_CohPsi2sToMu); }

  RooDataSet *dataINMC_Gen_CohPsi2sToMuPi = (RooDataSet*) dataIN18l7_Gen_CohPsi2sToMuPi;
    if (trigger.Contains("CMUP11")){ dataINMC_Gen_CohPsi2sToMuPi->append(*dataIN16b2_Gen_CohPsi2sToMuPi); }

  // --------------------------------------------------------------------------------
  // Do kinematic cuts
  // // ---Reconstructed
  RooAbsData *dataMC_Rec_CohJpsiToMu_PtCut = dataINMC_Rec_CohJpsiToMu->reduce(Cut(cuts_m_mc_Rec_PtCut));
  RooAbsData *dataMC_Rec_CohPsi2sToMu_PtCut = dataINMC_Rec_CohPsi2sToMu->reduce(Cut(cuts_m_mc_Rec_PtCut));
  RooAbsData *dataMC_Rec_CohPsi2sToMuPi_PtCut = dataINMC_Rec_CohPsi2sToMuPi->reduce(Cut(cuts_m_mc_Rec_PtCut));

  RooAbsData *dataMC_Rec_CohJpsiToMu_PtAll = dataINMC_Rec_CohJpsiToMu->reduce(Cut(cuts_m_mc_Rec_PtAll));
  RooAbsData *dataMC_Rec_CohPsi2sToMu_PtAll = dataINMC_Rec_CohPsi2sToMu->reduce(Cut(cuts_m_mc_Rec_PtAll));
  RooAbsData *dataMC_Rec_CohPsi2sToMuPi_PtAll = dataINMC_Rec_CohPsi2sToMuPi->reduce(Cut(cuts_m_mc_Rec_PtAll));

  // ---Generated
  RooAbsData *dataMC_Gen_CohJpsiToMu = dataINMC_Gen_CohJpsiToMu->reduce(Cut(cuts_m_mc_Gen));
  RooAbsData *dataMC_Gen_CohPsi2sToMu = dataINMC_Gen_CohPsi2sToMu->reduce(Cut(cuts_m_mc_Gen));
  RooAbsData *dataMC_Gen_CohPsi2sToMuPi = dataINMC_Gen_CohPsi2sToMuPi->reduce(Cut(cuts_m_mc_Gen));

  // --------------------------------------------------------------------------------
  // Count the number of events
  // ---Reconstructed
  Int_t N_Rec_CohJpsiToMu_PtCut = dataMC_Rec_CohJpsiToMu_PtCut->numEntries();
  Int_t N_Rec_CohPsi2sToMu_PtCut = dataMC_Rec_CohPsi2sToMu_PtCut->numEntries();
  Int_t N_Rec_CohPsi2sToMuPi_PtCut = dataMC_Rec_CohPsi2sToMuPi_PtCut->numEntries();

  Int_t N_Rec_CohJpsiToMu_PtAll = dataMC_Rec_CohJpsiToMu_PtAll->numEntries();
  Int_t N_Rec_CohPsi2sToMu_PtAll = dataMC_Rec_CohPsi2sToMu_PtAll->numEntries();  
  Int_t N_Rec_CohPsi2sToMuPi_PtAll = dataMC_Rec_CohPsi2sToMuPi_PtAll->numEntries();

  // ---Generated
  Int_t N_Gen_CohJpsiToMu = dataMC_Gen_CohJpsiToMu->numEntries();
  Int_t N_Gen_CohPsi2sToMu = dataMC_Gen_CohPsi2sToMu->numEntries();
  Int_t N_Gen_CohPsi2sToMuPi = dataMC_Gen_CohPsi2sToMuPi->numEntries();

  // // --------------------------------------------------------------------------------
  // // Create efficiency histogram
  // TEfficiency *Hist_Eff_CohJpsiToMu_PtCut  = new TEfficiency("CohJpsiToMu_PtCut","CohJpsiToMu_PtCut",1,0.,1.);
  // TEfficiency *Hist_Eff_CohPsi2sToMu_PtCut  = new TEfficiency("CohPsi2sToMu_PtCut","CohPsi2sToMu_PtCut",1,0.,1.);
  // TEfficiency *Hist_Eff_CohPsi2sToMuPi_PtCut  = new TEfficiency("CohPsi2sToMuPi_PtCut","CohPsi2sToMuPi_PtCut",1,0.,1.);

  // TEfficiency *Hist_Eff_CohJpsiToMu_PtAll  = new TEfficiency("CohJpsiToMu_PtAll","CohJpsiToMu_PtAll",1,0.,1.);
  // TEfficiency *Hist_Eff_CohPsi2sToMu_PtAll  = new TEfficiency("CohPsi2sToMu_PtAll","CohPsi2sToMu_PtAll",1,0.,1.);  
  // TEfficiency *Hist_Eff_CohPsi2sToMuPi_PtAll = new TEfficiency("CohPsi2sToMuPi_PtAll","CohPsi2sToMuPi_PtAll",1,0.,1.);

  // // ---Fill histograms
  // for (auto i=0; i<N_Rec_CohJpsiToMu_PtCut; i++){
  //     Hist_Eff_CohJpsiToMu_PtCut->Fill(kTRUE,0.5);
  // }
  // for (auto i=0; i<N_Rec_CohJpsiToMu_PtAll; i++){
  //     Hist_Eff_CohJpsiToMu_PtAll->Fill(kTRUE,0.5);
  // }  
  // for (auto i=0; i<(N_Gen_CohJpsiToMu-N_Rec_CohJpsiToMu_PtCut); i++){
  //     Hist_Eff_CohJpsiToMu_PtCut->Fill(kFALSE,0.5);
  // }
  // for (auto i=0; i<(N_Gen_CohJpsiToMu-N_Rec_CohJpsiToMu_PtAll); i++){
  //     Hist_Eff_CohJpsiToMu_PtAll->Fill(kFALSE,0.5);
  // }

  // for (auto i=0; i<N_Rec_CohPsi2sToMu_PtCut; i++){
  //     Hist_Eff_CohPsi2sToMu_PtCut->Fill(kTRUE,0.5);
  // }
  // for (auto i=0; i<N_Rec_CohPsi2sToMu_PtAll; i++){
  //     Hist_Eff_CohPsi2sToMu_PtAll->Fill(kTRUE,0.5);
  // }  
  // for (auto i=0; i<(N_Gen_CohPsi2sToMu-N_Rec_CohPsi2sToMu_PtCut); i++){
  //     Hist_Eff_CohPsi2sToMu_PtCut->Fill(kFALSE,0.5);
  // }
  // for (auto i=0; i<(N_Gen_CohPsi2sToMu-N_Rec_CohPsi2sToMu_PtAll); i++){
  //     Hist_Eff_CohPsi2sToMu_PtAll->Fill(kFALSE,0.5);
  // }
  
  // for (auto i=0; i<N_Rec_CohPsi2sToMuPi_PtCut; i++){
  //     Hist_Eff_CohPsi2sToMuPi_PtCut->Fill(kTRUE,0.5);
  // }
  // for (auto i=0; i<N_Rec_CohPsi2sToMuPi_PtAll; i++){
  //   Hist_Eff_CohPsi2sToMuPi_PtAll->Fill(kTRUE,0.5);
  // }  
  // for (auto i=0; i<(N_Gen_CohPsi2sToMuPi-N_Rec_CohPsi2sToMuPi_PtCut); i++){
  //     Hist_Eff_CohPsi2sToMuPi_PtCut->Fill(kFALSE,0.5);
  // }
  // for (auto i=0; i<(N_Gen_CohPsi2sToMuPi-N_Rec_CohPsi2sToMuPi_PtAll); i++){
  //     Hist_Eff_CohPsi2sToMuPi_PtAll->Fill(kFALSE,0.5);
  // }

  // // --------------------------------------------------------------------------------
  // // Compute efficiency
  // auto Eff_CohJpsiToMu_PtCut = Hist_Eff_CohJpsiToMu_PtCut->GetEfficiency(1);
  // auto Err_Up_Eff_CohJpsiToMu_PtCut = Hist_Eff_CohJpsiToMu_PtCut->GetEfficiencyErrorUp(1);
  // auto Err_Down_Eff_CohJpsiToMu_PtCut = Hist_Eff_CohJpsiToMu_PtCut->GetEfficiencyErrorLow(1);

  // auto Eff_CohJpsiToMu_PtAll = Hist_Eff_CohJpsiToMu_PtAll->GetEfficiency(1);
  // auto Err_Up_Eff_CohJpsiToMu_PtAll = Hist_Eff_CohJpsiToMu_PtAll->GetEfficiencyErrorUp(1);
  // auto Err_Down_Eff_CohJpsiToMu_PtAll = Hist_Eff_CohJpsiToMu_PtAll->GetEfficiencyErrorLow(1);
  // // -----
  // auto Eff_CohPsi2sToMu_PtCut = Hist_Eff_CohPsi2sToMu_PtCut->GetEfficiency(1);
  // auto Err_Up_Eff_CohPsi2sToMu_PtCut = Hist_Eff_CohPsi2sToMu_PtCut->GetEfficiencyErrorUp(1);
  // auto Err_Down_Eff_CohPsi2sToMu_PtCut = Hist_Eff_CohPsi2sToMu_PtCut->GetEfficiencyErrorLow(1);

  // auto Eff_CohPsi2sToMu_PtAll = Hist_Eff_CohPsi2sToMu_PtAll->GetEfficiency(1);
  // auto Err_Up_Eff_CohPsi2sToMu_PtAll = Hist_Eff_CohPsi2sToMu_PtAll->GetEfficiencyErrorUp(1);
  // auto Err_Down_Eff_CohPsi2sToMu_PtAll = Hist_Eff_CohPsi2sToMu_PtAll->GetEfficiencyErrorLow(1);
  // // -----
  // auto Eff_CohPsi2sToMuPi_PtCut = Hist_Eff_CohPsi2sToMuPi_PtCut->GetEfficiency(1);
  // auto Err_Up_Eff_CohPsi2sToMuPi_PtCut = Hist_Eff_CohPsi2sToMuPi_PtCut->GetEfficiencyErrorUp(1);
  // auto Err_Down_Eff_CohPsi2sToMuPi_PtCut = Hist_Eff_CohPsi2sToMuPi_PtCut->GetEfficiencyErrorLow(1);

  // auto Eff_CohPsi2sToMuPi_PtAll = Hist_Eff_CohPsi2sToMuPi_PtAll->GetEfficiency(1);
  // auto Err_Up_Eff_CohPsi2sToMuPi_PtAll = Hist_Eff_CohPsi2sToMuPi_PtAll->GetEfficiencyErrorUp(1);
  // auto Err_Down_Eff_CohPsi2sToMuPi_PtAll = Hist_Eff_CohPsi2sToMuPi_PtAll->GetEfficiencyErrorLow(1);

  // --------------------------------------------------------------------------------
  // Compute efficiency
  Double_t Eff_CohJpsiToMu_PtCut = (Double_t) N_Rec_CohJpsiToMu_PtCut/N_Gen_CohJpsiToMu;

  Double_t Eff_CohJpsiToMu_PtAll = (Double_t) N_Rec_CohJpsiToMu_PtAll/N_Gen_CohJpsiToMu;
  // -----
  Double_t Eff_CohPsi2sToMu_PtCut = (Double_t) N_Rec_CohPsi2sToMu_PtCut/N_Gen_CohPsi2sToMu;

  Double_t Eff_CohPsi2sToMu_PtAll = (Double_t) N_Rec_CohPsi2sToMu_PtAll/N_Gen_CohPsi2sToMu;
  // -----
  Double_t Eff_CohPsi2sToMuPi_PtCut = (Double_t) N_Rec_CohPsi2sToMuPi_PtCut/N_Gen_CohPsi2sToMuPi;

  Double_t Eff_CohPsi2sToMuPi_PtAll = (Double_t) N_Rec_CohPsi2sToMuPi_PtAll/N_Gen_CohPsi2sToMuPi;

  cout << N_Rec_CohJpsiToMu_PtCut << endl;
  cout << N_Gen_CohJpsiToMu << endl;
  cout << N_Rec_CohJpsiToMu_PtCut / N_Gen_CohJpsiToMu << endl;
  cout << Eff_CohJpsiToMu_PtCut << endl;
  // --------------------------------------------------------------------------------
  // Write the values to a file
  char FileName[120];
  sprintf(FileName,"Efficiency/%s/%s_%.2f_%.2f.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));


  // ofstream EffFile;
  // EffFile.open (FileName, ios::trunc);
  // EffFile << "-----------------------------------------------------------" << endl;
  // EffFile << "Eff_CohJpsiToMu_PtCut = " << Eff_CohJpsiToMu_PtCut << " + " << Err_Up_Eff_CohJpsiToMu_PtCut << " - " << Err_Down_Eff_CohJpsiToMu_PtCut << endl;
  // EffFile << "Eff_CohJpsiToMu_PtAll = " << Eff_CohJpsiToMu_PtAll << " + " << Err_Up_Eff_CohJpsiToMu_PtAll << " - " << Err_Down_Eff_CohJpsiToMu_PtAll << endl;
  // EffFile << "Eff_CohPsi2sToMu_PtCut = " << Eff_CohPsi2sToMu_PtCut << " + " << Err_Up_Eff_CohPsi2sToMu_PtCut << " - " << Err_Down_Eff_CohPsi2sToMu_PtCut << endl;
  // EffFile << "Eff_CohPsi2sToMu_PtAll = " << Eff_CohPsi2sToMu_PtAll << " + " << Err_Up_Eff_CohPsi2sToMu_PtAll << " - " << Err_Down_Eff_CohPsi2sToMu_PtAll << endl;
  // EffFile << "Eff_CohPsi2sToMuPi_PtCut = " << Eff_CohPsi2sToMuPi_PtCut << " + " << Err_Up_Eff_CohPsi2sToMuPi_PtCut << " - " << Err_Down_Eff_CohPsi2sToMuPi_PtCut << endl;
  // EffFile << "Eff_CohPsi2sToMuPi_PtAll = " << Eff_CohPsi2sToMuPi_PtAll << " + " << Err_Up_Eff_CohPsi2sToMuPi_PtAll << " - " << Err_Down_Eff_CohPsi2sToMuPi_PtAll << endl;
  // EffFile << "-----------------------------------------------------------" << endl;
  // EffFile << endl;  
  // EffFile.close();

  ofstream EffFile;
  EffFile.open (FileName, ios::trunc);
  EffFile << "-----------------------------------------------------------" << endl;
  EffFile << "Eff_CohJpsiToMu_PtCut = " << Eff_CohJpsiToMu_PtCut << endl;
  EffFile << "Eff_CohJpsiToMu_PtAll = " << Eff_CohJpsiToMu_PtAll << endl;
  EffFile << "Eff_CohPsi2sToMu_PtCut = " << Eff_CohPsi2sToMu_PtCut << endl;
  EffFile << "Eff_CohPsi2sToMu_PtAll = " << Eff_CohPsi2sToMu_PtAll << endl;
  EffFile << "Eff_CohPsi2sToMuPi_PtCut = " << Eff_CohPsi2sToMuPi_PtCut << endl;
  EffFile << "Eff_CohPsi2sToMuPi_PtAll = " << Eff_CohPsi2sToMuPi_PtAll << endl;
  EffFile << "-----------------------------------------------------------" << endl;
  EffFile << endl;  
  EffFile.close();

  // TCanvas *cM = new TCanvas("cM","cM",1600,800);
  // cM->cd(1);
  // RooPlot* frame_m = fRunNum.frame(Title("Mass fit")) ; 
  // data18l7_CohJpsiToMu_Template2->plotOn(frame_m,Name("data_m"),Binning(2200),MarkerStyle(20),MarkerSize(0.5));
  // frame_m->Draw();
  // cM->SaveAs("test.png");

  // ------------------------------------------------------------------------------------------------------------------------------------
  // Deleting
  if (dataIN18l7_Rec_CohJpsiToMu) {delete dataIN18l7_Rec_CohJpsiToMu;}
  if (dataIN18l7_Rec_CohPsi2sToMu) {delete dataIN18l7_Rec_CohPsi2sToMu;}
  if (dataIN18l7_Rec_CohPsi2sToMuPi) {delete dataIN18l7_Rec_CohPsi2sToMuPi;}
  if (dataIN18l7_Gen_CohJpsiToMu) {delete dataIN18l7_Gen_CohJpsiToMu;}
  if (dataIN18l7_Gen_CohPsi2sToMu) {delete dataIN18l7_Gen_CohPsi2sToMu;}
  if (dataIN18l7_Gen_CohPsi2sToMuPi) {delete dataIN18l7_Gen_CohPsi2sToMuPi;}
  if (dataIN16b2_Rec_CohJpsiToMu) {delete dataIN16b2_Rec_CohJpsiToMu;}
  if (dataIN16b2_Rec_CohPsi2sToMu) {delete dataIN16b2_Rec_CohPsi2sToMu;}
  if (dataIN16b2_Rec_CohPsi2sToMuPi) {delete dataIN16b2_Rec_CohPsi2sToMuPi;}
  if (dataIN16b2_Gen_CohJpsiToMu) {delete dataIN16b2_Gen_CohJpsiToMu;}
  if (dataIN16b2_Gen_CohPsi2sToMu) {delete dataIN16b2_Gen_CohPsi2sToMu;}
  if (dataIN16b2_Gen_CohPsi2sToMuPi) {delete dataIN16b2_Gen_CohPsi2sToMuPi;}
  // if (Hist_Eff_CohJpsiToMu_PtCut) {delete Hist_Eff_CohJpsiToMu_PtCut;}
  // if (Hist_Eff_CohPsi2sToMu_PtCut) {delete Hist_Eff_CohPsi2sToMu_PtCut;}
  // if (Hist_Eff_CohPsi2sToMuPi_PtCut) {delete Hist_Eff_CohPsi2sToMuPi_PtCut;}
  // if (Hist_Eff_CohJpsiToMu_PtAll) {delete Hist_Eff_CohJpsiToMu_PtAll;}
  // if (Hist_Eff_CohPsi2sToMu_PtAll) {delete Hist_Eff_CohPsi2sToMu_PtAll;}
  // if (Hist_Eff_CohPsi2sToMuPi_PtAll) {delete Hist_Eff_CohPsi2sToMuPi_PtAll;}
}

//////////////////////////////////////////////////////////////////////
// Do the efficiency  analysis
//////////////////////////////////////////////////////////////////////

void Efficiency(Int_t y_bin, TString ADveto, TString trigger)
{
  ComputeEfficiency(y_bins_min[y_bin], y_bins_max[y_bin], ADveto, trigger);
  // gROOT->ProcessLine(".q");
}

