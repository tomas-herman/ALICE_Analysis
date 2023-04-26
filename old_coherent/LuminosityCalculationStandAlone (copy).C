// Luminosity calculation
//-----------------------------------------------------------------------------------

// C/C++
#include "iostream"
#include "iomanip"
#include "fstream"
#include "sstream"
#include "stdio.h"
#include "string.h"

// ROOT
#include "TCanvas.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TEfficiency.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitter.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPad.h"
#include "TParticle.h"
#include "TPaveText.h"
#include "TPaveLabel.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector.h"

// AliROOT
#include "AliCDBManager.h"
#include "AliTriggerScalers.h"
#include "AliTriggerRunScalers.h"
#include "AliTimeStamp.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerConfiguration.h"
#include "AliLHCData.h"
#include "AliTriggerClass.h"
#include "AliTriggerBCMask.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

using namespace std;

// My headers
#include "GoodRuns.h"

///////////////////////////////////////////////////////////
// Set histogram look
///////////////////////////////////////////////////////////

void SetLumiHisto(TH1D* h, Color_t color = kBlue)
{
  gStyle->SetOptStat(0);
  h->SetTitleFont(43);
  h->SetTitleSize(25);
  h->GetYaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelSize(10);
  h->GetYaxis()->SetLabelSize(25);
  h->GetYaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetTitleOffset(0.6);
  h->GetYaxis()->SetDecimals(1);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetFillColor(color);
  h->SetMarkerSize(0.7);
  h->SetMarkerStyle(kFullCross);
  h->LabelsOption("v");
  h->SetMinimum(0);
  h->SetLineWidth(2);
  h->Sumw2(kFALSE);
}

///////////////////////////////////////////////////////////
// Plot Lumi per run
///////////////////////////////////////////////////////////
void SumLumi(TString period, TString trigger)
{

  // Open input file and read lumi per run
  TString LumiHistoName = "Lumi/LumiHisto";
  LumiHistoName += period;
  LumiHistoName += trigger;
  LumiHistoName += ".root";
  TFile *fInputLumiHisto = TFile::Open(LumiHistoName.Data(), "read");

  // Create plots
  TH1D *hLumi = dynamic_cast<TH1D*> (fInputLumiHisto->Get("hLumi"));
  TH1D *hLumiS = dynamic_cast<TH1D*> (fInputLumiHisto->Get("hLumiS"));
  Double_t integrated_lumi_ana(0.), integrated_lumi_rec(0.);

  char FileNamePDF[200];

  if (period.Contains("18q")){
    if (trigger.Contains("CMUP6")){
      hLumi->SetTitle("Luminosity 18q-CMUP6");
      sprintf(FileNamePDF,"Lumi/18q_lumi_CMUP6.pdf");
    } else if (trigger.Contains("CMUP11")){
      hLumi->SetTitle("Luminosity 18q-CMUP11");
      sprintf(FileNamePDF,"Lumi/18q_lumi_CMUP11.pdf");
    }
  } else if (period.Contains("18r")){
    if (trigger.Contains("CMUP6")){
      hLumi->SetTitle("Luminosity 18r_CMUP6");
      sprintf(FileNamePDF,"Lumi/18r_lumi_CMUP6.pdf");
    } else if (trigger.Contains("CMUP11")){
      hLumi->SetTitle("Luminosity 18r_CMUP11");
      sprintf(FileNamePDF,"Lumi/18r_lumi_CMUP11.pdf");
    }  
  } else if (period.Contains("15o")){
    if (trigger.Contains("CMUP10")){
      hLumi->SetTitle("Luminosity 15o_CMUP10");
      sprintf(FileNamePDF,"Lumi/15o_lumi_CMUP10.pdf");
    } else if (trigger.Contains("CMUP11")){
      hLumi->SetTitle("Luminosity 15o_CMUP11");
      sprintf(FileNamePDF,"Lumi/15o_lumi_CMUP11.pdf");
    } 
  }

  Double_t RangeMax = 0;
  // Get Lumi
  for(Int_t i(0);i<hLumiS->GetNbinsX();i++){
    integrated_lumi_ana += hLumiS->GetBinContent(i+1);
    integrated_lumi_rec += hLumi->GetBinContent(i+1);

    if(hLumi->GetBinContent(i+1) > RangeMax) RangeMax = hLumi->GetBinContent(i+1);
  }

  // Plotting
  SetLumiHisto(hLumi);
  SetLumiHisto(hLumiS,kRed);

  TCanvas *cLumi = new TCanvas("cLumi","cLumi",1500,500);

  hLumi->GetYaxis()->SetTitle("L_{int}[ub^{-1}]");
  hLumi->SetMaximum(RangeMax*1.2);
  hLumi->Draw();
  hLumiS->Draw("sameP0");
  TLegend* legLumi = new TLegend(0.7,0.8,0.9,0.9);
  legLumi->SetFillColor(kWhite);
  legLumi->AddEntry(hLumi,Form("Total lumi recorded: %.3f ub^{-1}",integrated_lumi_rec),"l");
  legLumi->AddEntry(hLumiS,Form("Total lumi analysed: %.3f ub^{-1}",integrated_lumi_ana),"l");
  legLumi->Draw();

  // Save plot   
  cLumi->SaveAs(FileNamePDF);

  cout << "Luminosity computed for " << trigger << ", " << period << endl;

}//end sumlumi

///////////////////////////////////////////////////////////
// Compute the luminosity
///////////////////////////////////////////////////////////
void LuminosityCalculationStandAlone(TString period, TString trigger)
{
  // Fill good runs
  vector<Int_t> runList;

  if (period.Contains("18q")){
    for(auto i : Vector_GoodRuns18q){
      runList.emplace_back(i);
    }
  } else if (period.Contains("18r")){
    for(auto i : Vector_GoodRuns18r){
      runList.emplace_back(i);
    }
  } else if (period.Contains("15o")){
    if (trigger.Contains("CMUP10")){
      for(auto i : Vector_GoodRuns15o_CMUP10){
        runList.emplace_back(i);
      }
    }  
    if (trigger.Contains("CMUP11")){
      for(auto i : Vector_GoodRuns15o_CMUP11){
        runList.emplace_back(i);
      }
    }  
  }

  Int_t nRunsInList = runList.size();

  // Here you put histogram with fired triggers per run in your dataset
  TH1D* hCMUPtrg = new TH1D("hCMUPtrg","",nRunsInList,0,nRunsInList);
  // Filling the histogram from my dataset and creating name for output file
  char FilePath[400];
  char FileNameTXT[200];
  if (period.Contains("18q")){
    if (trigger.Contains("CMUP6")){
      sprintf(FilePath,"/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18q_AnalysisResults.root");
      sprintf(FileNameTXT,"Lumi/18q_lumi_CMUP6.txt");
    } else if (trigger.Contains("CMUP11")){
      sprintf(FilePath,"/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18q_AnalysisResults.root");
      sprintf(FileNameTXT,"Lumi/18q_lumi_CMUP11.txt");
    }
  } else if (period.Contains("18r")){
    if (trigger.Contains("CMUP6")){
      sprintf(FilePath,"/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18r_AnalysisResults.root");
      sprintf(FileNameTXT,"Lumi/18r_lumi_CMUP6.txt");    
    } else if (trigger.Contains("CMUP11")){
      sprintf(FilePath,"/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC18r_AnalysisResults.root");
      sprintf(FileNameTXT,"Lumi/18r_lumi_CMUP11.txt");    
    }  
  } else if (period.Contains("15o")){
    if (trigger.Contains("CMUP10")){
      sprintf(FilePath,"/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC15o_AnalysisResults.root");
      sprintf(FileNameTXT,"Lumi/15o_lumi_CMUP10.txt");    
    } else if (trigger.Contains("CMUP11")){
      sprintf(FilePath,"/mnt/data/Data_processing/Processed_data/TrainResults/580_20200706-1520-Data/LHC15o_AnalysisResults.root");
      sprintf(FileNameTXT,"Lumi/15o_lumi_CMUP11.txt");    
    } 
  }

  TFile *TrgFile= new TFile(FilePath); 
  TTree *TrgTree;

  TString folder;
  if (period.Contains("15o")){
    if (trigger.Contains("CMUP")) {folder = "NanoMUONCMUP11";}
  }
  if (period.Contains("18")){
    if (trigger.Contains("CMUP6")) {folder = "NanoMUONCMUP6";}
    if (trigger.Contains("CMUP11")) {folder = "NanoMUONCMUP11";}
  }

  TrgFile->GetObject(folder+"/fTrgTree",TrgTree); 
  Int_t fTrgRunNum;
  Int_t fCMUP6;
  Int_t fCMUP10;
  Int_t fCMUP11;
  TrgTree->SetBranchAddress("fTrgRunNum",&fTrgRunNum);
  TrgTree->SetBranchAddress("fCMUP6",&fCMUP6);
  TrgTree->SetBranchAddress("fCMUP10",&fCMUP10);
  TrgTree->SetBranchAddress("fCMUP11",&fCMUP11);

  Int_t nEntries = TrgTree->GetEntries();
  for (Int_t irun=0;irun<nRunsInList;irun++){
    for (Int_t i=0;i<nEntries;i++){
      TrgTree->GetEntry(i);
      if (fTrgRunNum==runList[irun]){
        if (trigger.Contains("CMUP6")){
          if (fCMUP6==1) hCMUPtrg->Fill(irun);
        } else if (trigger.Contains("CMUP10")){  
          if (fCMUP10==1) hCMUPtrg->Fill(irun);      
        } else if (trigger.Contains("CMUP11")){  
          if (fCMUP11==1) hCMUPtrg->Fill(irun);  
        }       
      }
    }  
    cout << "\e[A\r\e[0K"<< "Computed " << irun+1 << " runs out of " << nRunsInList <<endl;
  }

  // Prepare objects for luminosity calculation
  const Int_t nrunsmax = 300;

  TString className;
  if (trigger.Contains("CMUP6")){
    className.Append("CMUP6-B-NOPF-MUFAST");
  } else if (trigger.Contains("CMUP10")){  
    className.Append("CMUP10-B-NOPF-MUFAST");
  } else if (trigger.Contains("CMUP11")){  
    className.Append("CMUP11-B-NOPF-MUFAST");
  }

  TChain* t = new TChain("trending");
  if (period.Contains("18q")){
    t->AddFile("Lumi/trending_merged_PbPb2018.root");
  } else if (period.Contains("18r")){
    t->AddFile("Lumi/trending_merged_PbPb2018.root");
  } else if (period.Contains("15o")){
    t->AddFile("Lumi/trending_merged_PbPb2015.root");
  }

  t->LoadTree(0);
  TObjArray* classes = new TObjArray();
  Double_t  lumi_seen[nrunsmax] = {0};
  Double_t  class_lumi[nrunsmax] = {0};
  Double_t  class_ds[nrunsmax] = {0};
  ULong64_t  class_l2a[nrunsmax] = {0};
  Int_t run;
  Double_t mu = 0;
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("classes",&classes);
  t->SetBranchAddress("lumi_seen",&lumi_seen);
  t->SetBranchAddress("class_lumi",&class_lumi);
  t->SetBranchAddress("class_ds",&class_ds);
  t->SetBranchAddress("class_l2a",&class_l2a);
  t->BuildIndex("run");

  TH1D* hLumi = new TH1D("hLumi","",nRunsInList,0,nRunsInList);
  TH1D* hLumiS = new TH1D("hLumiS","",nRunsInList,0,nRunsInList);
  TH1D* hScale = new TH1D("hScale","",nRunsInList,0,nRunsInList);
  TH1D* hDS = new TH1D("hDS","Trigger downscaling",nRunsInList,0,nRunsInList);

  ofstream LumiFile;
  LumiFile.open (FileNameTXT, ios::trunc);

  // Calculate seen luminosity for specific trigger class
  for (Int_t i=0;i<nRunsInList;i++){  // ordering of runs in histogram
  // for (Int_t i=nRunsInList-1;i>=0;i--){
    Int_t r = runList[i];
    char* srun = Form("%i",r);
    //Printf("%s %i %i",srun, i, r);
    t->GetEntryWithIndex(r);
    AliTriggerClass* cl = (AliTriggerClass*) classes->FindObject(className.Data());
    if (!cl) continue;
    Int_t iclass = classes->IndexOf(cl);
    // Printf("%i %f %i %s",run, mu, iclass, cl->GetName());
    Double_t l2a = (Double_t) class_l2a[iclass];
    //Printf("%.llu",class_l2a[iclass]);
    //Printf("%.10f",class_lumi[iclass]);
    Double_t scale = hCMUPtrg->GetBinContent(i+1)/l2a;
    hScale->Fill(srun,scale);
    hLumiS->Fill(srun,scale*class_lumi[iclass]);
    hLumi->Fill(srun,class_lumi[iclass]);
    hDS->Fill(srun,class_ds[iclass]);

    LumiFile << srun << ", " << scale*class_lumi[iclass] << "; " << mu <<  endl; 
  }
  LumiFile.close();

  TString LumiHistoName = "Lumi/LumiHisto";
  LumiHistoName += period;
  LumiHistoName += trigger;
  LumiHistoName += ".root";
  TFile* fOutputLumiHisto = new TFile(LumiHistoName.Data(),"recreate");
  hLumi->Write();
  hLumiS->Write();
  hScale->Write();
  hDS->Write();
  hCMUPtrg->Write();
  fOutputLumiHisto->Close();
  
  SumLumi(period, trigger);
}
