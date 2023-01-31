//
// this program computes the acceptance time efficiency of MC data
// for the forward j/psi analysis
//

// c++ headers
#include <iostream>

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMath.h"

// my headers
#include "GoodRuns.h"
#include "TreeVariables.h"
#include "Selection.h"

// ----------------------------------------------------

void GetRec(TTree *tr,  Double_t *Runs,
	    Double_t y_min, Double_t y_max,
	    Double_t m_min, Double_t m_max,
	    Double_t pt_min, Double_t pt_max)
{
  SetRecTree(tr);
  Int_t nevt = tr->GetEntries();
  for(Int_t i=0;i<nevt;i++) {
    tr->GetEntry(i);
    // apply selection
    if (!Sel_Pair_Charge()) continue;
    if (!Sel_Pair_Match()) continue;
    if (!Sel_Pair_pXdca()) continue;
    if (!Sel_PairMass(m_min,m_max)) continue;
    if (!Sel_PairPt(pt_min,pt_max)) continue;
    if (!Sel_PairRap(y_min,y_max)) continue;
    // event passes
    Runs[jRunN-FirstRun]+=1.0;
  }

}

// ----------------------------------------------------

void GetGen(TTree *tg, Double_t *Runs,
	    Double_t y_min, Double_t y_max,
	    Double_t m_min, Double_t m_max)
{
  // set variables
  tg->SetBranchAddress("jGenY",&jGenY);
  tg->SetBranchAddress("jGenM",&jGenM);  
  tg->SetBranchAddress("jRunN",&jRunN);
  Int_t nevt = tg->GetEntries();
  for(Int_t i=0;i<nevt;i++) {
    tg->GetEntry(i);
    // select
    if (jGenM>m_max) continue;
    if (jGenM<m_min) continue;
    if (jGenY<y_min) continue;
    if (jGenY>y_max) continue;
    // event passes
    Runs[jRunN-FirstRun]+=1.0;
  }
}

// ----------------------------------------------------

Double_t ErrDiv(Double_t x, Double_t y)
// error of x/y
{
  Double_t sigma = x*TMath::Sqrt((1./x)+(1./y))/y;
  return sigma;
}

// ----------------------------------------------------

void GetAxE(TH1D *AxE_h, Int_t nGoodRuns, Int_t *GoodRuns,
	    Double_t *GoodGen,Double_t *GoodRec,Double_t *tot16s)
{
  for(Int_t i=0;i<nGoodRuns;i++) {
    AxE_h->GetXaxis()->SetBinLabel(i+1,Form("%i",GoodRuns[i]));
    if (GoodGen[GoodRuns[i]-FirstRun] <0.1) continue;
    if (GoodRec[GoodRuns[i]-FirstRun] <0.1) continue;    
    Double_t axe_i = GoodRec[GoodRuns[i]-FirstRun]
      /GoodGen[GoodRuns[i]-FirstRun];
    Double_t axe_err_i = ErrDiv(GoodRec[GoodRuns[i]-FirstRun],
				GoodGen[GoodRuns[i]-FirstRun]);
    AxE_h->SetBinContent(i+1,axe_i);
    AxE_h->SetBinError(i+1,axe_err_i);
    tot16s[0] += GoodRec[GoodRuns[i]-FirstRun];
    tot16s[1] += GoodGen[GoodRuns[i]-FirstRun];    
  }
}

// ----------------------------------------------------

void MakeAxE(Int_t opt, Int_t trd)
// trd = 0 => all runs
// trd = 1 => runs with trd
// trd = 2 => runs without trd
{
  // define selection
  Double_t y_min = -4.0;
  Double_t y_max = -2.5;
  Double_t pt_min = 0.0;
  Double_t pt_max = 1.0;
  Double_t m_min = 2.8;
  Double_t m_max = 3.3;
  Double_t m_gen_min = 0;
  Double_t m_gen_max = 10;

  // define runs to study
  Double_t GoodRec[nRuns];
  Double_t GoodGen[nRuns];
  for (Int_t i=0;i<nRuns;i++) GoodRec[i]=GoodGen[i]=0.0;
  SetGoodRuns();

  // define histos to store the results
  TH1D *AxE16r_h = new TH1D("AxE16r_h","AxE16r_h",nGoodRunsLHC16r,0,nGoodRunsLHC16r);
  TH1D *AxE16s_h = new TH1D("AxE16s_h","AxE16s_h",nGoodRunsLHC16s,0,nGoodRunsLHC16r);
  
  // open the right file
  TFile *fIn = NULL;
  if (opt == 101) fIn = new TFile("ANA_MC_INCOH.root"); // 101 = incoherent
  if (opt == 102) fIn = new TFile("ANA_MC_gg.root"); // 102 = gg
  if (opt == 104) fIn = new TFile("ANA_MC_PSI2S.root"); // 104 = psi2s
  if (opt == 105) fIn = new TFile("ANA_MC_COH.root"); // 105 = gPb  
  
  if (opt == 102) { // gg => change mass range
    m_min = m_gen_min = 1.5;
    m_max = m_gen_max = 2.0;
  }
  
  // get trees
  TTree *tr = (TTree*) fIn->Get("jRecTree");
  TTree *tg = (TTree*) fIn->Get("jGenTree");  
  
  // get events that pass rec and gen selections
  GetRec(tr,GoodRec,y_min,y_max,m_min,m_max,pt_min,pt_max);
  GetGen(tg,GoodGen,y_min,y_max,m_gen_min,m_gen_max); 

  // get AxE for LHC16s
  Double_t AxE_tot_LHC16s = 0;
  Double_t AxE_tot_err_LHC16s = 0;
  Double_t AxE_tot_LHC16r = 0;
  Double_t AxE_tot_err_LHC16r = 0;
  if (trd == 0) {
    Double_t  tot16s[2];
    GetAxE(AxE16s_h,nGoodRunsLHC16s,GoodRunsLHC16s,GoodGen,GoodRec,tot16s);
    AxE_tot_LHC16s=tot16s[0]/tot16s[1];
    AxE_tot_err_LHC16s=ErrDiv(tot16s[0],tot16s[1]);  
    // get AxE for LHC16r
    Double_t  tot16r[2];
    GetAxE(AxE16r_h,nGoodRunsLHC16r,GoodRunsLHC16r,GoodGen,GoodRec,tot16r);
    AxE_tot_LHC16r=tot16r[0]/tot16r[1];
    AxE_tot_err_LHC16r=ErrDiv(tot16r[0],tot16r[1]);
    cout << " AxE 16r = " << AxE_tot_LHC16r << " +/- " << AxE_tot_err_LHC16r << endl;
    cout << " AxE 16s = " << AxE_tot_LHC16s << " +/- " << AxE_tot_err_LHC16s << endl;
  } else if (trd == 1) {
    Double_t  tot16rTRD[2];
    GetAxE(AxE16r_h,nTRDLHC16rRun,TRDLHC16rRun,GoodGen,GoodRec,tot16rTRD);
    Double_t AxE_tot_LHC16rTRD=tot16rTRD[0]/tot16rTRD[1];
    Double_t AxE_tot_err_LHC16rTRD=ErrDiv(tot16rTRD[0],tot16rTRD[1]);  
    Double_t  tot16sTRD[2];
    GetAxE(AxE16s_h,nTRDLHC16sRun,TRDLHC16sRun,GoodGen,GoodRec,tot16sTRD);
    Double_t AxE_tot_LHC16sTRD=tot16sTRD[0]/tot16sTRD[1];
    Double_t AxE_tot_err_LHC16sTRD=ErrDiv(tot16sTRD[0],tot16sTRD[1]);
    cout << " AxE 16r TRD = " << AxE_tot_LHC16rTRD << " +/- " << AxE_tot_err_LHC16rTRD << endl;
    cout << " AxE 16s TRD = " << AxE_tot_LHC16sTRD << " +/- " << AxE_tot_err_LHC16sTRD << endl;
  } else if (trd == 2) {
    // get AxE for LHC16r without TRD
    Double_t  tot16rNoTRD[2];
    GetAxE(AxE16r_h,nNoTRDLHC16rRun,NoTRDLHC16rRun,GoodGen,GoodRec,tot16rNoTRD);
    Double_t AxE_tot_LHC16rNoTRD=tot16rNoTRD[0]/tot16rNoTRD[1];
    Double_t AxE_tot_err_LHC16rNoTRD=ErrDiv(tot16rNoTRD[0],tot16rNoTRD[1]);  
    // get AxE for LHC16s without TRD
    Double_t  tot16sNoTRD[2];
    GetAxE(AxE16s_h,nNoTRDLHC16sRun,NoTRDLHC16sRun,GoodGen,GoodRec,tot16sNoTRD);
    Double_t AxE_tot_LHC16sNoTRD=tot16sNoTRD[0]/tot16sNoTRD[1];
    Double_t AxE_tot_err_LHC16sNoTRD=ErrDiv(tot16sNoTRD[0],tot16sNoTRD[1]);
    cout << " AxE 16r No TRD = " << AxE_tot_LHC16rNoTRD << " +/- " << AxE_tot_err_LHC16rNoTRD << endl;
    cout << " AxE 16s No TRD = " << AxE_tot_LHC16sNoTRD << " +/- " << AxE_tot_err_LHC16sNoTRD << endl;
  }

        
  if (trd == 0) {
    // plot histos
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas *cAxEr = new TCanvas("cAxEr","cAxEr",1200,500);
    AxE16r_h->Draw("e1");
    AxE16r_h->GetYaxis()->SetRangeUser(0,AxE_tot_LHC16r*1.5);
    TLatex *AxErLTX = new TLatex();
    AxErLTX->SetTextFont(42);
    AxErLTX->DrawLatex(2, AxE_tot_LHC16r*0.5,
		       Form("AxE in 16r = %7.6f #pm %7.6f",
			    AxE_tot_LHC16r,AxE_tot_err_LHC16r));
    AxErLTX->DrawLatex(2, AxE_tot_LHC16r*0.3,
		       Form("%3.2f < y < %3.2f",y_min,y_max));
    
    char cname[120];
    sprintf(cname,"AxE_LHC16r_%i_%3.0f_y_%3.0f.pdf",opt,TMath::Abs(y_max)*100,TMath::Abs(y_min)*100);
    cAxEr->Print(cname);
    
    TCanvas *cAxEs = new TCanvas("cAxEs","cAxEs",100,100,1200,500);
    AxE16s_h->Draw("e1");
    AxE16s_h->GetYaxis()->SetRangeUser(0,AxE_tot_LHC16s*1.5);
    TLatex *AxEsLTX = new TLatex();
    AxEsLTX->SetTextFont(42);
    AxEsLTX->DrawLatex(2, AxE_tot_LHC16s*0.5,
		       Form("AxE in 16s = %7.6f #pm %7.6f",
			    AxE_tot_LHC16s,AxE_tot_err_LHC16s));
    
    AxEsLTX->DrawLatex(2, AxE_tot_LHC16s*0.3,
		       Form("%3.2f < y < %3.2f",y_min,y_max));
    
    sprintf(cname,"AxE_LHC16s_%i_%3.0f_y_%3.0f.pdf",opt,TMath::Abs(y_max)*100,TMath::Abs(y_min)*100);
    cAxEs->Print(cname);
  }
}
