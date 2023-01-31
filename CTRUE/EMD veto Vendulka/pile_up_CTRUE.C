
// data from 2018 from
// http://aliqaevs.web.cern.ch/aliqaevs/data/2018/LHC18q/muon_calo_pass2/
// http://aliqaevs.web.cern.ch/aliqaevs/data/2018/LHC18r/muon_calo_pass2/

// analysis of CTRUE events to study veto inefficiency due to pile-up from different colisions
// file: /alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/543_20191009-1453/merge/AnalysisResults.root

// c++ headers
#include <iostream>
#include <fstream>

// root headers
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2D.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <TBits.h>
#include <TLatex.h>

// my headers
#include "GoodRuns_2018.h"
//#include "GoodRuns_2015.h"
#include "TreeVariables1.h"
#include "TriggerInputDefinitions.h"
// #include "LeadTailBC.h"

// global variable
Double_t Ctrue[nRuns];

// glabal bools
Bool_t f0VBA;
Bool_t f0VBC;
Bool_t f0UBA;
Bool_t f0UBC;
Bool_t f1ZED;
Bool_t f0MUL;
Bool_t f0STG;

// --------------------------------
// setup tree
// --------------------------------

void SetupTree(TFile *fIn)
{
  fAnaTree = (TTree*) fIn->Get("CtrueTask/fAnaTree");
  fAnaTree ->SetBranchAddress("fRunNum", &fRunNum);
  fAnaTree ->SetBranchAddress("fBCrossNum", &fBC);  
  fAnaTree ->SetBranchAddress("fTracklets", &fTracklets);	
  fAnaTree ->SetBranchAddress("fCtrue", &fCtrue);
  fAnaTree ->SetBranchAddress("fC1zed", &fC1zed);  
  fAnaTree ->SetBranchAddress("fL0inputs", &fL0inputs);
  fAnaTree ->SetBranchAddress("fL1inputs", &fL1inputs);
  fAnaTree ->SetBranchAddress("fZem1Energy", &fZem1Energy);
  fAnaTree ->SetBranchAddress("fZem2Energy", &fZem2Energy);  
  fAnaTree ->SetBranchAddress("fZNCEnergy", &fZNCEnergy);  
  fAnaTree ->SetBranchAddress("fZNAEnergy", &fZNAEnergy);
  fAnaTree ->SetBranchAddress("fZPCEnergy", &fZPCEnergy);
  fAnaTree ->SetBranchAddress("fZPAEnergy", &fZPAEnergy);  
  fAnaTree ->SetBranchAddress("fZNATDC", &fZNATDC[0]);
  fAnaTree ->SetBranchAddress("fZNCTDC", &fZNCTDC[0]);  
  fAnaTree ->SetBranchAddress("fZPATDC", &fZPATDC[0]);
  fAnaTree ->SetBranchAddress("fZPCTDC", &fZPCTDC[0]);  
  fAnaTree ->SetBranchAddress("fV0ADecision", &fV0ADecision);
  fAnaTree ->SetBranchAddress("fV0CDecision", &fV0CDecision);  
  fAnaTree ->SetBranchAddress("fADADecision", &fADADecision);
  fAnaTree ->SetBranchAddress("fADCDecision", &fADCDecision);  
  fAnaTree ->SetBranchAddress("fIR1Map", &fIR1Map);
  fAnaTree ->SetBranchAddress("fIR2Map", &fIR2Map);   

}

// ---------------------------------------------------------------------
Double_t berr(Double_t cut, Double_t ctrue)
// use Root formula for binomial error
// https://root-forum.cern.ch/t/how-to-calculate-binomial-efficiency-error-with-weights/3650/3
{
  Double_t eff = cut/ctrue;
  Double_t e1 = TMath::Sqrt(cut);
  Double_t e2 = TMath::Sqrt(ctrue);
  return TMath::Sqrt(TMath::Abs( ((1.0-2.0*eff)*cut + eff*eff*ctrue)/(ctrue*ctrue) ));
}

// ---------------------------------------------------------------------
Double_t fit_p1(Double_t *x, Double_t *p)
// fit model of a pol1 passing through 0
{
  return p[0]+p[1]*x[0];
}

// ---------------------------------------------------------------------
void DoEff(Int_t n, Double_t *weight, Double_t *mu,
	   Double_t p0, Double_t p0e, Double_t p1, Double_t p1e,
	   Double_t *eff)
{
  Double_t w_tot = 0.0;
  for(Int_t i=0;i<n;i++) {
    Double_t w = weight[i];
    w_tot += w;
    Double_t p = p0+p1*mu[i];
    eff[0] += (w*TMath::Exp(-p));
    Double_t pm = p0-p0e+(p1-p1e)*mu[i];	       
    eff[1] += (w*TMath::Exp(-pm));
    Double_t pp = p0+p0e+(p1+p1e)*mu[i];	           
    eff[2] += (w*TMath::Exp(-pp));            
    //  cout << i << " w " << w << " p " << p << " eff " <<  (w*TMath::Exp(-p)) << endl;
  }
  if (w_tot < 1e-12) return;
  for(Int_t i=0;i<3;i++) eff[i] = eff[i]/w_tot;
  
  cout << " efficiency  = " << eff[0]
       << " + " << (eff[1]-eff[0])
       << " - " << (eff[0]-eff[2])    
       << endl;
}

// ---------------------------------------------------------------------
void DoPlot(const char *title, Double_t *signal, Double_t *empty, Bool_t fix_p0 = kFALSE, Bool_t fix_p1 = kFALSE)
{
  // decide if constraint the origin of the line fit
  // to pass through the origin
  // Bool_t fix_p0 = kFALSE;

  // prepare for fractions to store weights to total efficiency
  Double_t *weight = new Double_t [nRuns];
  for (Int_t i=0;i<nRuns;i++) weight[i]=0.0;
  
  // prepare graphs
  Double_t ymin =10;
  Double_t ymax =-10;  
  Double_t all_signal = 0;
  Double_t all_empty = 0;  
  TGraphErrors *gr = new TGraphErrors();
  for(Int_t i=0;i<nRuns;i++) {
    // skip events without entries
    if (empty[i]<0.5) continue;
    if (signal[i]<0.5) continue;
    // increase counters
    all_signal+=signal[i];
    all_empty+=empty[i];  
    // get the point
    Double_t x = mu_all[i];
    Double_t y = signal[i]/empty[i];
    Double_t ye = berr(signal[i],empty[i]);
    if (ymin>y-ye) ymin = y-ye;
    if (ymax<y+ye) ymax = y+ye;
    // if (y < 0.01) continue;
    // fill graph
    Int_t n = gr->GetN();
    gr->SetPoint(n,x,y);
    gr->SetPointError(n,0,ye);

    weight[n]=w_all[i];
  }
  
  //plot graph
  // No stats and no titles
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1);  

  TCanvas *c = new TCanvas(title,title,1200,600);
  c->Divide(1,1);
  c->cd(1);
  gr->SetMarkerStyle(20);  gr->SetMarkerColor(kBlue+1); gr->SetLineColor(kBlue+1);
  // make some white space in the canvas
  ymin *= 0.5;
  if (ymin<0) ymin = 0.001;
  ymax *= 1.5;
  // define the histo
  Double_t mu_max = 0.0025;
  TH1F* frame1 = gPad->DrawFrame(0.,ymin,mu_max,ymax);
  frame1->GetYaxis()->SetTitleOffset(1.3);
  frame1->GetYaxis()->SetLabelSize(0.025);
  frame1->SetTitle(Form("%s;#mu;Probability",title));
  TF1 *pol = new TF1("pol",fit_p1, 0, 1,2);
  if (gr->GetN()>0) {
    gr->Draw("p,same");
    pol->SetParNames("p0","p1");
    pol->SetParameter(0,0.0);
    pol->SetParameter(1,1.0);
    pol->SetLineColor(kRed+1);
   // pol->SetParLimits(1,-10000,0);    
    if (fix_p0) pol->FixParameter(0,0.0);
    if (fix_p1) pol->FixParameter(1,0.0);
    gr->Fit("pol");
  } else { cout << " gr is empty " << endl;}
  gPad->Update();

  // compute efficiency
  Double_t eff[3]; // mean, minus, plus
  eff[0]=eff[1]=eff[2]=0.0;
  DoEff(gr->GetN(),weight,gr->GetX(),pol->GetParameter(0), pol->GetParError(0),
	pol->GetParameter(1), pol->GetParError(1), eff);

  // print eff in the canvas
  TLatex* l = new TLatex();
  l->SetTextFont(42);
  l->SetTextSize(0.045);
  l->DrawLatex(mu_max*0.05,ymax*0.9,Form("#varepsilon = %.5f + %.5f - %.5f",eff[0],(eff[1]-eff[0]),(eff[0]-eff[2])));
 
  // cout<<"eff = "<<eff[0]<<endl;
  // cout<<"err+_eff = "<<eff[1]-eff[0]<<endl;
  // cout<<"err-_eff = "<<eff[0]-eff[2]<<endl;

  // save
  // c->Print(Form("Canvas_%s.pdf",title));
  c->SaveAs(Form("%s.pdf",title));
  // clean up
  delete [] weight;
/*
 cout<<"eff = "<<eff[0]<<endl;
 cout<<"err+_eff = "<<eff[1]-eff[0]<<endl;
 cout<<"err-_eff = "<<eff[0]-eff[2]<<endl;
*/
}

// ---------------------------------------------------------------------
void pile_up_CTRUE(Int_t apply_period = 0, Int_t apply_pf=0, Int_t apply_bc=0)
{  
  // initialising
  SetGoodRuns(apply_period);
  SetGoodRunsIdx();
  Set2018PbPb();
  //Set2015PbPb();  
  Fill_mu();
  Fill_w();
  //  FillBCmaskInfo();

  // open files,  get trees
  TFile *InFile = new TFile("/mnt/data/Data_processing/Processed_data/CTRUE_AnalysisResults_PbPb_AOD_543.root"); //CTRUE 2018
  SetupTree(InFile);
 
  // ctrue and selection counters

  //ZNA pielup
  Double_t dZNA1a[nRuns] = {0};
  Double_t dZNA1b[nRuns] = {0};
  Double_t dZNA2a[nRuns] = {0};
  Double_t dZNA2b[nRuns] = {0};
  Double_t dZNA3a[nRuns] = {0};
  Double_t dZNA3b[nRuns] = {0};

  //ZNC pielup
  Double_t dZNC1a[nRuns] = {0};
  Double_t dZNC1b[nRuns] = {0};
  Double_t dZNC2a[nRuns] = {0};
  Double_t dZNC2b[nRuns] = {0};
  Double_t dZNC3a[nRuns] = {0};
  Double_t dZNC3b[nRuns] = {0};

  // All V0A
  Double_t E1[nRuns] = {0};

  Double_t VBA1[nRuns] = {0};
  Double_t VDA1a[nRuns] = {0};
  Double_t VDA1b[nRuns] = {0};

  //0nXn V0A
  Double_t d0nXn1a[nRuns] = {0};
  Double_t d0nXn1b[nRuns] = {0};
  Double_t d0nXn2a[nRuns] = {0};
  Double_t d0nXn2b[nRuns] = {0};
  Double_t d0nXn3a[nRuns] = {0};
  Double_t d0nXn3b[nRuns] = {0};

  //Xn0n V0A
  Double_t VBA_or_VDA_no_Tlets[nRuns] = {0};
  Double_t no_UDC_VDC_Tlets[nRuns] = {0};
  Double_t VBA_or_VDA_no_UDA[nRuns] = {0};
  Double_t no_UDA_UDC_VDC_Tlets[nRuns] = {0};
  Double_t VBA_or_VDA[nRuns] = {0};
  Double_t no_UDC_VDC[nRuns] = {0};

  //XnXn V0A
  Double_t VDA_or_VBA_no_Tlets[nRuns] = {0};
  Double_t no_UBA_UBC_Tlets[nRuns] = {0};
  Double_t VDA_or_VBA_no_AD_Tlets[nRuns] = {0};
  Double_t no_AD_Tlets[nRuns] = {0};
  Double_t VDA_or_VBA_ZNA_ZNC[nRuns] = {0};
  Double_t ZNA_ZNC[nRuns] = {0};


  // loop over events
  cout << " Entries " << fAnaTree->GetEntries() << endl;
  for(Int_t i=0;i<fAnaTree->GetEntries();i++) {
    // for(Int_t i=0;i<100;i++) {    
    // get event
    fAnaTree->GetEntry(i);

    // check that is an appropriate run
    Int_t run_idx = fRunNum-FirstRun;
 
    if (run_idx < 0 ) continue;
    if ((run_idx-1)>nRuns) continue;    
    if (GoodRuns[run_idx] == 0) continue;

    //only ctrue events
    if (fCtrue!=1) continue;
    Ctrue[run_idx] += 1.0;

   if (apply_period == 1){//LHC18 q
     if (fRunNum > 296625) continue;}
   if (apply_period == 2){//LHC18 r
     if (fRunNum < 296625) continue;}
    
   //    // check bc type
   //    if (apply_bc == 1) {
   //      if (fRunNum < 296625) { // 18q
  	// if(!Tail18q[GoodRunsIdx[run_idx]][fBC]) continue;
   //      } else {
  	// if(!Tail18r[GoodRunsIdx[run_idx]][fBC]) continue;
   //      }
   //    } else if (apply_bc == 2) {
   //      if (fRunNum < 296625) { // 18q
  	// if(!Lead18q[GoodRunsIdx[run_idx]][fBC]) continue;
   //      } else {
  	// if(!Lead18r[GoodRunsIdx[run_idx]][fBC]) continue;
   //      }
   //    }
 
    // evaluate past-future protection
    Bool_t PF_IR1 = (fIR1Map->CountBits()>0);
    Bool_t PF_IR2 = (fIR2Map->CountBits()>0);
    if (apply_pf == 1 && PF_IR1 ) continue;
    if (apply_pf == 2 && PF_IR2 ) continue;
    if (apply_pf == 3 && (PF_IR1 || PF_IR2) ) continue;        

    // get L0 trigger flags
    f0VBA = fL0inputs & 1 << (inputId_0VBA-1);
    f0VBC = fL0inputs & 1 << (inputId_0VBC-1);
    f0UBA = fL0inputs & 1 << (inputId_0UBA-1);
    f0UBC = fL0inputs & 1 << (inputId_0UBC-1);

    // set ZDC variables
    Bool_t ZNAhit = kFALSE;
    Bool_t ZNChit = kFALSE;
    if ((TMath::Abs(fZNATDC[0])<2) || (TMath::Abs(fZNATDC[1])<2) ||	(TMath::Abs(fZNATDC[2])<2) ||  (TMath::Abs(fZNATDC[3])<2))
      ZNAhit = kTRUE;
    
    if ((TMath::Abs(fZNCTDC[0])<2) || (TMath::Abs(fZNCTDC[1])<2) ||	(TMath::Abs(fZNCTDC[2])<2) ||  (TMath::Abs(fZNCTDC[3])<2))
      ZNChit = kTRUE;

    // ###################### ZNA ######################
    Bool_t EmptyZNA1 = (!ZNChit && fV0ADecision == 0 && fV0CDecision == 0 && fADADecision == 0 && fADCDecision == 0 && fTracklets == 0);
    Bool_t EmptyZNA2 = (!ZNChit && fV0ADecision == 0 && fV0CDecision == 0 && fADADecision == 0 && fADCDecision == 0);
    Bool_t EmptyZNA3 = (!ZNChit && fV0ADecision == 0 && fV0CDecision == 0 && fTracklets == 0);

    if (ZNAhit == 1 && EmptyZNA1) 
      dZNA1a[run_idx]+=1.0;
    if (EmptyZNA1)
      dZNA1b[run_idx]+=1.0;

    if (ZNAhit == 1 && EmptyZNA2) 
      dZNA2a[run_idx]+=1.0;
    if (EmptyZNA2)
      dZNA2b[run_idx]+=1.0;

    if (ZNAhit == 1 && EmptyZNA3) 
      dZNA3a[run_idx]+=1.0;
    if (EmptyZNA3)
      dZNA3b[run_idx]+=1.0;

    // ###################### ZNC ######################    
    Bool_t EmptyZNC1 = (!ZNAhit && fV0ADecision == 0 && fV0CDecision == 0 && fADADecision == 0 && fADCDecision == 0 && fTracklets == 0);
    Bool_t EmptyZNC2 = (!ZNAhit && fV0ADecision == 0 && fV0CDecision == 0 && fADADecision == 0 && fADCDecision == 0);
    Bool_t EmptyZNC3 = (!ZNAhit && fV0ADecision == 0 && fV0CDecision == 0 && fTracklets == 0);

    if (ZNChit == 1 && EmptyZNC1) 
      dZNC1a[run_idx]+=1.0;
    if (EmptyZNC1)
      dZNC1b[run_idx]+=1.0;

    if (ZNChit == 1 && EmptyZNC2) 
      dZNC2a[run_idx]+=1.0;
    if (EmptyZNC2)
      dZNC2b[run_idx]+=1.0;

    if (ZNChit == 1 && EmptyZNC3) 
      dZNC3a[run_idx]+=1.0;
    if (EmptyZNC3)
      dZNC3b[run_idx]+=1.0;

    // ###################### V0A ######################
    //--------All------------------------
    Bool_t Empty1 = (!ZNAhit && !ZNChit && fV0CDecision == 0 && fADADecision == 0 && fADCDecision == 0 && fTracklets == 0);

    if (Empty1)
      E1[run_idx]+=1.0;

    if (f0VBA && Empty1)
      VBA1[run_idx]+=1.0;
    if (fV0ADecision == 1 && Empty1)
      VDA1a[run_idx]+=1.0;
    if (fV0ADecision != 0 && Empty1)
      VDA1b[run_idx]+=1.0;

    //--------(0nXn)------------------------
    //V0A  Pile-up
     Bool_t EmptyA = (!ZNAhit && fADADecision == 0 && fTracklets == 0);

    //V0A Pile-up !! ------CAREFUL! each sample needs correction->for V0 I do not care about AD or ZN, for AD I want empty V0
       if (fV0ADecision == 1 && EmptyA) //no condition on AD, condition on 0 tracklets
        d0nXn1a[run_idx]+=1.0;
      if (EmptyA)
        d0nXn1b[run_idx]+=1.0;

       if (fV0ADecision == 1 && EmptyA) //empty AD, condition on 0 tracklets
        d0nXn2a[run_idx]+=1.0;
      if (fADADecision == 0 && EmptyA)
        d0nXn2b[run_idx]+=1.0;

       if (fV0ADecision == 1 && !ZNAhit && fADADecision == 0) //no condition on AD, no condition on tracklets
        d0nXn3a[run_idx]+=1.0;
      if (!ZNAhit && fADADecision == 0)
        d0nXn3b[run_idx]+=1.0;

    //--------(Xn0n)------------------------
    //V0A  Pile-up
     Bool_t EmptyC = (!ZNChit && fV0CDecision == 0 && fADCDecision == 0 && fTracklets == 0);

    //V0A Pile-up !! ------CAREFUL! each sample needs correction->for V0 I do not care about AD or ZN, for AD I want empty V0
       if (fV0ADecision == 1 && EmptyC) //no condition on AD, condition on 0 tracklets
        VBA_or_VDA_no_Tlets[run_idx]+=1.0;
      if (EmptyC)
        no_UDC_VDC_Tlets[run_idx]+=1.0;

       if (fV0ADecision == 1 && fADADecision == 0 && EmptyC) //empty AD, condition on 0 tracklets
        VBA_or_VDA_no_UDA[run_idx]+=1.0;
      if (fADADecision == 0 && EmptyC)
        no_UDA_UDC_VDC_Tlets[run_idx]+=1.0;

       if (fV0ADecision == 1 && !ZNChit && fV0CDecision == 0 && fADCDecision == 0) //no condition on AD, no condition on tracklets
        VBA_or_VDA[run_idx]+=1.0;
      if (!ZNChit && fV0CDecision == 0 && fADCDecision == 0)
        no_UDC_VDC[run_idx]+=1.0;

    //-------------(XnXn)-------------------------------
    Bool_t EmptyAC = (fTracklets == 0);
    
    //V0 Pile-up
    if(fV0ADecision == 1 && EmptyAC)
     VDA_or_VBA_no_Tlets[run_idx]+=1.0;
    if(EmptyAC)
     no_UBA_UBC_Tlets[run_idx]+=1.0;

    if(fV0ADecision == 1 && fADADecision == 0 && fADCDecision == 0 && EmptyAC)
     VDA_or_VBA_no_AD_Tlets[run_idx]+=1.0;
    if(EmptyAC  && fADADecision == 0 && fADCDecision == 0)
     no_AD_Tlets[run_idx]+=1.0;

    if(fV0ADecision == 1 &&  fADADecision == 0 && fADCDecision == 0)
     VDA_or_VBA_ZNA_ZNC[run_idx]+=1.0;
    if(fADADecision == 0 && fADCDecision == 0)
     ZNA_ZNC[run_idx]+=1.0;

  }

  // Do Plots
  // ###################### ZNA ###################### 
  DoPlot("CTRUE_ZNA1", dZNA1a, dZNA1b);
  DoPlot("CTRUE_ZNA2", dZNA2a, dZNA2b);
  DoPlot("CTRUE_ZNA3", dZNA3a, dZNA3b);

  // ###################### ZNC ###################### 
  DoPlot("CTRUE_ZNC1", dZNC1a, dZNC1b);
  DoPlot("CTRUE_ZNC2", dZNC2a, dZNC2b);
  DoPlot("CTRUE_ZNC3", dZNC3a, dZNC3b);

  // ###################### V0A ######################
  DoPlot("CTRUE_VBA1", VBA1, E1);
  DoPlot("CTRUE_VDA1a", VDA1a, E1);
  DoPlot("CTRUE_VDA1b", VDA1b, E1);

  DoPlot("CTRUE_VDA_0nXn1", d0nXn1a, d0nXn1b);
  DoPlot("CTRUE_VDA_0nXn2", d0nXn2a, d0nXn2b); //full requirements on empty events
  DoPlot("CTRUE_VDA_0nXn3", d0nXn3a, d0nXn3b);

  DoPlot("CTRUE_VDA_Xn0n1", VBA_or_VDA_no_Tlets, no_UDC_VDC_Tlets);
  DoPlot("CTRUE_VDA_Xn0n2", VBA_or_VDA_no_UDA, no_UDA_UDC_VDC_Tlets); //full requirements on empty events
  DoPlot("CTRUE_VDA_Xn0n3", VBA_or_VDA, no_UDC_VDC);

  DoPlot("CTRUE_VDA_XnXn1", VDA_or_VBA_no_Tlets, no_UBA_UBC_Tlets);//no tracklets
  DoPlot("CTRUE_VDA_XnXn2", VDA_or_VBA_no_AD_Tlets, no_AD_Tlets); //full requirements on empty events
  DoPlot("CTRUE_VDA_XnXn3", VDA_or_VBA_ZNA_ZNC, ZNA_ZNC);//no UDAC
}

