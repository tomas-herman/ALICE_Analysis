
// data from 2018 from
// http://aliqaevs.web.cern.ch/aliqaevs/data/2018/LHC18q/muon_calo_pass2/
// http://aliqaevs.web.cern.ch/aliqaevs/data/2018/LHC18r/muon_calo_pass2/

// analysis of C1ZED events to study veto inefficiency due to pileup cuased by EMD in the same collision
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
    pol->SetParameter(0,0.995);
    pol->SetParameter(1,-0.1);
    pol->SetLineColor(kRed+1);
   // pol->SetParLimits(1,-10000,0);    
    if (fix_p0) pol->FixParameter(0,0.0);
    if (fix_p1) pol->FixParameter(1,0.0);
    gr->Fit("pol");
  } else { cout << " gr is empty " << endl;}
  gPad->Update();

  if (!fix_p1) {
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
  } else {
    Double_t ratio = all_signal/all_empty;
    Double_t err = berr(all_signal,all_empty);
    TLatex* l = new TLatex();
    l->SetTextFont(42);
    l->SetTextSize(0.045);
    l->DrawLatex(mu_max*0.05,ymax*0.9,Form("ratio = %.5f #pm %.5f",ratio,err));

  }

  // save
  // c->Print(Form("Canvas_%s.pdf",title));
  c->SaveAs(Form("%s.pdf",title));
  // clean up
  delete [] weight;
}

// ---------------------------------------------------------------------
void Veto_CorrectionFactors(Int_t apply_period = 0)
{
  // apply_period controls the runs to be studied, either
  // = 0 all LHC18, = 1 LHC18q and = 2 LHC18r
  
  // initialising
  SetGoodRuns(apply_period);
  SetGoodRunsIdx();
  Set2018PbPb();
  //Set2015PbPb();  
  Fill_mu();
  Fill_w();

  // open file,  get tree
  TFile *InFile = new TFile("/mnt/data/Data_processing/Processed_data/CTRUE_AnalysisResults_PbPb_AOD_543.root"); //CTRUE 2018
  SetupTree(InFile);

  // ctrue and selection counters
  //(0n0n)
  Double_t d0n0n1a[nRuns] = {0};
  Double_t d0n0n1b[nRuns] = {0};
  Double_t d0n0n2a[nRuns] = {0};
  Double_t d0n0n2b[nRuns] = {0};
  Double_t d0n0n3a[nRuns] = {0};
  Double_t d0n0n3b[nRuns] = {0};
  //(0nXn)
  Double_t d0nXn1a[nRuns] = {0};
  Double_t d0nXn1b[nRuns] = {0};
  Double_t d0nXn2a[nRuns] = {0};
  Double_t d0nXn2b[nRuns] = {0};
  Double_t d0nXn3a[nRuns] = {0};
  Double_t d0nXn3b[nRuns] = {0};
  //(Xn0n)
  Double_t histo_EMDA_VDA_noVBA_emptyC1[nRuns] = {0};
  Double_t histo_EMDA_VDA_emptyC1[nRuns] = {0};
  Double_t histo_EMDA_VDA_noVBA_emptyC2[nRuns] = {0};
  Double_t histo_EMDA_VDA_emptyC2[nRuns] = {0};
  Double_t histo_EMDA_VDA_noVBA_emptyC3[nRuns] = {0};
  Double_t histo_EMDA_VDA_emptyC3[nRuns] = {0};
  Double_t histo_EMDA_VDA_noVBA_emptyC4[nRuns] = {0};
  Double_t histo_EMDA_VDA_emptyC4[nRuns] = {0};
  //(XnXn)
  Double_t histo_EMD_VDA_no_VBA_emptyAC1[nRuns] = {0};
  Double_t histo_EMD_VDA_emptyAC1[nRuns] = {0};
  Double_t histo_EMD_VDA_no_VBA_emptyAC2[nRuns] = {0};
  Double_t histo_EMD_VDA_emptyAC2[nRuns] = {0};
  Double_t histo_EMD_VDA_no_VBA_emptyAC3[nRuns] = {0};
  Double_t histo_EMD_VDA_emptyAC3[nRuns] = {0};

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

    //only 1ZED events
    if (fC1zed!=1) continue;

    if (apply_period == 1){//LHC18 q
     if (fRunNum > 296625) continue;}
    if (apply_period == 2){//LHC18 r
     if (fRunNum < 296625) continue;}

    // get L0 trigger flags
    f0VBA = fL0inputs & 1 << (inputId_0VBA-1);
    f0VBC = fL0inputs & 1 << (inputId_0VBC-1);
    f0UBA = fL0inputs & 1 << (inputId_0UBA-1);
    f0UBC = fL0inputs & 1 << (inputId_0UBC-1);

    // set ZDC variables
    Bool_t ZNAhit = kFALSE;
    Bool_t ZNChit = kFALSE;
    if ((TMath::Abs(fZNATDC[0])<2) || (TMath::Abs(fZNATDC[1])<2) || (TMath::Abs(fZNATDC[2])<2) ||  (TMath::Abs(fZNATDC[3])<2))
      ZNAhit = kTRUE;
    
    if ((TMath::Abs(fZNCTDC[0])<2) || (TMath::Abs(fZNCTDC[1])<2) || (TMath::Abs(fZNCTDC[2])<2) ||  (TMath::Abs(fZNCTDC[3])<2))
      ZNChit = kTRUE;

    //Define emptiness 
    Bool_t empty1 = (!ZNAhit && !ZNChit && fADADecision == 0 && fADCDecision == 0 && fV0CDecision == 0 && fTracklets == 0);
    Bool_t empty2 = (!ZNAhit && !ZNChit && fADADecision == 0 && fADCDecision == 0 && fV0CDecision == 0);
    Bool_t empty3 = (!ZNAhit && !ZNChit && fADADecision == 0 && fADCDecision == 0 && fTracklets == 0);

    Bool_t emptyA1 = (!ZNAhit && ZNChit && fADADecision == 0 && fTracklets == 0);
    Bool_t emptyA2 = (!ZNAhit && ZNChit && fADADecision == 0 && fADCDecision == 0 && fTracklets == 0);
    Bool_t emptyA3 = (!ZNAhit && ZNChit && fADADecision == 0 && fADCDecision == 0);

    Bool_t emptyC1 = (ZNAhit && !ZNChit && fADCDecision == 0 && fV0CDecision == 0 && fTracklets == 0);
    Bool_t emptyC2 = (ZNAhit && !ZNChit && fADCDecision == 0 && fV0CDecision == 0 && fADADecision == 0 && fTracklets == 0);
    Bool_t emptyC3 = (ZNAhit && !ZNChit && fADCDecision == 0 && fV0CDecision == 0 && fADADecision == 0);
    Bool_t emptyC4 = (ZNAhit && !ZNChit && fADCDecision == 0 && fADADecision == 0 && fTracklets == 0);

    Bool_t emptyAC1 = (ZNAhit && ZNChit && fADADecision == 0 && fADCDecision == 0 && fTracklets == 0); 
    Bool_t emptyAC2 = (ZNAhit && ZNChit && fTracklets == 0); 
    Bool_t emptyAC3 = (ZNAhit && ZNChit && fADADecision == 0 && fADCDecision == 0 ); 

    //(0n0n)
    Bool_t b0n0n1a = (fV0ADecision == 1 && !f0VBA && empty1);
    Bool_t b0n0n1b = (fV0ADecision == 1 && empty1);

    Bool_t b0n0n2a = (fV0ADecision == 1 && !f0VBA && empty2);
    Bool_t b0n0n2b = (fV0ADecision == 1 && empty2);

    Bool_t b0n0n3a = (fV0ADecision == 1 && !f0VBA && empty3);
    Bool_t b0n0n3b = (fV0ADecision == 1 && empty3);

    //(0nXn)
    Bool_t b0nXn1a = (fV0ADecision == 1 && !f0VBA && emptyA1);
    Bool_t b0nXn1b = (fV0ADecision == 1 && emptyA1);

    Bool_t b0nXn2a = (fV0ADecision == 1 && !f0VBA && emptyA2);
    Bool_t b0nXn2b = (fV0ADecision == 1 && emptyA2);

    Bool_t b0nXn3a = (fV0ADecision == 1 && !f0VBA && emptyA3);
    Bool_t b0nXn3b = (fV0ADecision == 1 && emptyA3);

    //(Xn0n)
    Bool_t EMDA_VDA_noVBA_emptyC1 = (fV0ADecision == 1 && !f0VBA && emptyC1);
    Bool_t EMDA_VDA_emptyC1 = (fV0ADecision == 1 && emptyC1);

    Bool_t EMDA_VDA_noVBA_emptyC2 = (fV0ADecision == 1 && !f0VBA && emptyC2);
    Bool_t EMDA_VDA_emptyC2 = (fV0ADecision == 1 && emptyC2);

    Bool_t EMDA_VDA_noVBA_emptyC3 = (fV0ADecision == 1 && !f0VBA && emptyC3);
    Bool_t EMDA_VDA_emptyC3 = (fV0ADecision == 1 && emptyC3);

    Bool_t EMDA_VDA_noVBA_emptyC4 = (fV0ADecision == 1 && !f0VBA && emptyC4);
    Bool_t EMDA_VDA_emptyC4 = (fV0ADecision == 1 && emptyC4);

    //(XnXn)
    Bool_t EMD_VDA_no_VBA_emptyAC1 = (emptyAC1 && fV0ADecision == 1 && !f0VBA);
    Bool_t EMD_VDA_emptyAC1 = (emptyAC1 && fV0ADecision == 1);

    Bool_t EMD_VDA_no_VBA_emptyAC2 = (emptyAC2 && fV0ADecision == 1 && !f0VBA);
    Bool_t EMDC_VDA_emptyAC2 = (emptyAC2 && fV0ADecision == 1);

    Bool_t EMD_VDA_no_VBA_emptyAC3 = (emptyAC3 && fV0ADecision == 1 && !f0VBA);
    Bool_t EMD_VDA_emptyAC3 = (emptyAC3 && fV0ADecision == 1);

    //----Filling the plots
    //--(0n0n)-----------------------------------
    if (b0n0n1a)
      d0n0n1a[run_idx]+=1.0;
    if (b0n0n1b)
      d0n0n1b[run_idx]+=1.0;

    if (b0n0n2a)
      d0n0n2a[run_idx]+=1.0;
    if (b0n0n2b)
      d0n0n2b[run_idx]+=1.0;

    if (b0n0n3a)
      d0n0n3a[run_idx]+=1.0;
    if (b0n0n3b)
      d0n0n3b[run_idx]+=1.0;

    //--(Xn0n)-----------------------------------
    if (b0nXn1a)
      d0nXn1a[run_idx]+=1.0;
    if (b0nXn1b)
      d0nXn1b[run_idx]+=1.0;

    if (b0nXn2a)
      d0nXn2a[run_idx]+=1.0;
    if (b0nXn2b)
      d0nXn2b[run_idx]+=1.0;

    if (b0nXn3a)
      d0nXn3a[run_idx]+=1.0;
    if (b0nXn3b)
      d0nXn3b[run_idx]+=1.0;

    //--(Xn0n)-----------------------------------
    if (EMDA_VDA_noVBA_emptyC1)
      histo_EMDA_VDA_noVBA_emptyC1[run_idx]+=1.0;
    if (EMDA_VDA_emptyC1)
      histo_EMDA_VDA_emptyC1[run_idx]+=1.0;

    if (EMDA_VDA_noVBA_emptyC2)
      histo_EMDA_VDA_noVBA_emptyC2[run_idx]+=1.0;
    if (EMDA_VDA_emptyC2)
      histo_EMDA_VDA_emptyC2[run_idx]+=1.0;

    if (EMDA_VDA_noVBA_emptyC3)
      histo_EMDA_VDA_noVBA_emptyC3[run_idx]+=1.0;
    if (EMDA_VDA_emptyC3)
      histo_EMDA_VDA_emptyC3[run_idx]+=1.0;

    if (EMDA_VDA_noVBA_emptyC4)
      histo_EMDA_VDA_noVBA_emptyC4[run_idx]+=1.0;
    if (EMDA_VDA_emptyC4)
      histo_EMDA_VDA_emptyC4[run_idx]+=1.0;

    //------(XnXn)-----------------------------------
    if (EMD_VDA_no_VBA_emptyAC1)
      histo_EMD_VDA_no_VBA_emptyAC1[run_idx]+=1.0;
    if (EMD_VDA_emptyAC1)
      histo_EMD_VDA_emptyAC1[run_idx]+=1.0;

    if (EMD_VDA_no_VBA_emptyAC2)
      histo_EMD_VDA_no_VBA_emptyAC2[run_idx]+=1.0;
    if (EMDC_VDA_emptyAC2)
      histo_EMD_VDA_emptyAC2[run_idx]+=1.0;

    if (EMD_VDA_no_VBA_emptyAC3)
      histo_EMD_VDA_no_VBA_emptyAC3[run_idx]+=1.0;
    if (EMD_VDA_emptyAC3)
      histo_EMD_VDA_emptyAC3[run_idx]+=1.0;
  }

  //----DoPlots-------
  //------(0n0n)--------------------------
  DoPlot("C1ZED_0n0n1",d0n0n1a,d0n0n1b);
  DoPlot("C1ZED_0n0n2",d0n0n2a,d0n0n2b);
  DoPlot("C1ZED_0n0n3",d0n0n3a,d0n0n3b);

  //------(0nXn)--------------------------
  DoPlot("C1ZED_0nXn1",d0nXn1a,d0nXn1b);
  DoPlot("C1ZED_0nXn2",d0nXn2a,d0nXn2b);
  DoPlot("C1ZED_0nXn3",d0nXn3a,d0nXn3b);

  //------(Xn0n)--------------------------
  DoPlot("C1ZED_Xn0n1",histo_EMDA_VDA_noVBA_emptyC1,histo_EMDA_VDA_emptyC1);
  DoPlot("C1ZED_Xn0n2",histo_EMDA_VDA_noVBA_emptyC2,histo_EMDA_VDA_emptyC2);
  DoPlot("C1ZED_Xn0n3",histo_EMDA_VDA_noVBA_emptyC3,histo_EMDA_VDA_emptyC3);
  DoPlot("C1ZED_Xn0n4",histo_EMDA_VDA_noVBA_emptyC4,histo_EMDA_VDA_emptyC4);

  //---(XnXn)-----------------------------------
  DoPlot("C1ZED_XnXn1",histo_EMD_VDA_no_VBA_emptyAC1,histo_EMD_VDA_emptyAC1);
  DoPlot("C1ZED_XnXn2",histo_EMD_VDA_no_VBA_emptyAC2,histo_EMD_VDA_emptyAC2);
  DoPlot("C1ZED_XnXn3",histo_EMD_VDA_no_VBA_emptyAC3,histo_EMD_VDA_emptyAC3);
}

