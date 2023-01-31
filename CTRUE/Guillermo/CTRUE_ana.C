
//
// code to analyse data from CTRUE events
//

//-----------------------------------------
// code to analyse data from CTRUE events
//-----------------------------------------

// c++ headers
#include <iostream>
using namespace std;

// root headers
#include "CandRoot.h"

// my headers
#include "TreeVariables.h"
#include "TriggerInputs.h"
#include "GoodRuns.h"
#include "mu.h"
#include "lumi_per_run_16rs.h"

//-----------------------------------------
// functions
//-----------------------------------------

//-----------------------------------------
// use Root formula for binomial error
Double_t berr(Double_t cut, Double_t ctrue)
{
  Double_t eff = cut/ctrue;
  Double_t e1 = TMath::Sqrt(cut);
  Double_t e2 = TMath::Sqrt(ctrue);
  return TMath::Sqrt(TMath::Abs( ((1.0-2.0*eff)*cut + eff*eff*ctrue)/(ctrue*ctrue) ));
}

//-----------------------------------------
// function to fit a pol1
Double_t fit_p1(Double_t *x, Double_t *p)
{
  return p[0]+p[1]*x[0];
}

//-----------------------------------------
// prepare a plot for the mu-dependence
void DoPlot(const char *title, Double_t *signal, Double_t *empty,
	    Bool_t fix_p0 = kFALSE, Bool_t fix_p1 = kFALSE)
{

  // prepare graphs
  Double_t ymin =10;
  Double_t ymax =-10;
  Double_t mumax = -1;
  
  TGraphErrors *gr_fit = new TGraphErrors();
  gr_fit->SetMarkerStyle(20);
  gr_fit->SetMarkerColor(kBlue);
  gr_fit->SetLineColor(kBlue);

  // fill graph
  for(Int_t i=0;i<nRuns;i++) {
    // skip events without entries
    if (empty[i]<0.5) continue;
    if (signal[i]<0.5) continue;
    // get the point
    Double_t x = mu_all[i];
    Double_t y = signal[i]/empty[i];
    Double_t ye = berr(signal[i],empty[i]);
    // reset boundries for plot
    if (ymin>y-ye) ymin = y-ye;
    if (ymax<y+ye) ymax = y+ye;
    if (mumax<x) mumax = x;
    // fill graph
    Int_t n_fit = gr_fit->GetN();
    gr_fit->SetPoint(n_fit,x,y);
    gr_fit->SetPointError(n_fit,0,ye);
  }

  // define and perform the fit
  TF1 *pol = new TF1("pol",fit_p1, 0, 1,2);
  pol->SetParNames("p0","p1");
  pol->SetParameter(0,1.1);
  pol->SetParameter(1,0.0);
  if (fix_p0) pol->FixParameter(0,0.0);
  if (fix_p1) pol->FixParameter(1,0.0);    
  gr_fit->Fit("pol");


  // compute efficiency
  Double_t eff = 0;
  Double_t eff_err = 0;
  Double_t tot_lumi = 0;  
  Double_t p1 = pol->GetParameter(1);
  Double_t p1e = pol->GetParError(1);  
   for(Int_t i=0;i<nRuns;i++) {
     if (GoodRunsIdx[i] < 0) continue;
     eff += (lumi_all[i]*mu_all[i]*p1);
     eff_err += (lumi_all[i]*mu_all[i]*(p1+p1e));
     tot_lumi += lumi_all[i];
   }

  //plot graph
  // No stats and no titles
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);  
  // canvas
  TCanvas *c = new TCanvas(title,title,1200,600);
  c->Divide(1,1);
  c->cd(1);
  // define the histo
  mumax *= 1.5;// make some white space in the canvas
  ymin *= 0.5; 
  if (ymin<0) ymin = 0.001;
  ymax *= 1.5;  
  TH1F* frame1 = gPad->DrawFrame(0.,ymin,mumax,ymax);
  frame1->GetYaxis()->SetTitleOffset(1.3);
  frame1->GetYaxis()->SetLabelSize(0.025);
  frame1->SetTitle(Form("%s;#mu;Prob.",title));
  gr_fit->Draw("p,same");
  gPad->Update();

 
  // print efficiency
  cout << "Eff for " << title << ": " << (eff/tot_lumi) << " +/- " << ((eff_err-eff)/tot_lumi) << endl;
  TLatex* l = new TLatex();
  l->SetTextFont(42);
  l->SetTextSize(0.045);
  l->DrawLatex(mumax*0.05,ymax*0.9,Form("pile-up = %.5f #pm %.5f",(eff/tot_lumi),((eff_err-eff)/tot_lumi)));

 // save
  c->Print(Form("Canvas_%s.pdf",title));

 
}


//-----------------------------------------
// check ZDC time
Bool_t isZDCtime(Double_t *zdctime)
{
  const Double_t max_time_zdc = 2;
  if ((TMath::Abs(zdctime[0])<max_time_zdc) ||
      (TMath::Abs(zdctime[1])<max_time_zdc) ||
      (TMath::Abs(zdctime[2])<max_time_zdc) ||
      (TMath::Abs(zdctime[3])<max_time_zdc))
    return kTRUE;
  return kFALSE;
}
   
//-----------------------------------------
// Printout message
void Printout(const char *name, Double_t den, Double_t num)
{
  cout << name << ": den = " << den
       << " num = " << num
       << " ratio = " << (num/den)
       << " +/- " << berr(num,den)
       << endl;
}

//-----------------------------------------
// main program
//-----------------------------------------

void CTRUE_ana(Int_t per)
// per = 0 => 16r; = 1 => 16s
{
  // open file,  get tree
  TFile *InFile = NULL;
  if (per == 0) InFile = new TFile("AnalysisResults_LHC16r_CTRUE.root");
  else if (per == 1) InFile = new TFile("AnalysisResults_LHC16s_CTRUE.root");
  else return;
  SetupTree(InFile);

  // correct input element 
  if (per == 0) inputId_0SH2 = 9;

  // set list of good runs
  SetGoodRuns_fwd(per);
  SetTRDRunsIdx(per);
  SetGoodRunsIdx_fwd(per);
  Set_mu_all();
  Set_lumi_all();
  
  // counters
  // --> zna
  Double_t count_zna_num = 0;
  Double_t count_zna_den = 0;
  Double_t histo_zna_num[nRuns];
  Double_t histo_zna_den[nRuns];
  // --> ada
  Double_t count_ada_num = 0;
  Double_t count_ada_den = 0;
  Double_t histo_ada_num[nRuns];
  Double_t histo_ada_den[nRuns];
  // --> v0a
  Double_t count_v0a_num = 0;
  Double_t count_v0a_den = 0;
  Double_t histo_v0a_num[nRuns];
  Double_t histo_v0a_den[nRuns];
  // --> zda||ada||v0a
  Double_t count_aside_num = 0;
  Double_t count_aside_den = 0;
  Double_t histo_aside_num[nRuns];
  Double_t histo_aside_den[nRuns];

  // --> znc
  Double_t count_znc_num = 0;
  Double_t count_znc_den = 0;
  Double_t histo_znc_num[nRuns];
  Double_t histo_znc_den[nRuns];
  // --> adc
  Double_t count_adc_num = 0;
  Double_t count_adc_den = 0;
  Double_t histo_adc_num[nRuns];
  Double_t histo_adc_den[nRuns];
  // --> v0c
  Double_t count_v0c_num = 0;
  Double_t count_v0c_den = 0;
  Double_t histo_v0c_num[nRuns];
  Double_t histo_v0c_den[nRuns];
  // --> zdc||adc||v0c
  Double_t count_cside_num = 0;
  Double_t count_cside_den = 0;
  Double_t histo_cside_num[nRuns];
  Double_t histo_cside_den[nRuns];


  // reset
  for(Int_t i=0;i<nRuns;i++) {
    histo_zna_num[i]=histo_zna_den[i]=0.0;
    histo_ada_num[i]=histo_ada_den[i]=0.0;    
    histo_v0a_num[i]=histo_v0a_den[i]=0.0;    
    histo_aside_num[i]=histo_aside_den[i]=0.0;    
    histo_znc_num[i]=histo_znc_den[i]=0.0;
    histo_adc_num[i]=histo_adc_den[i]=0.0;    
    histo_v0c_num[i]=histo_v0c_den[i]=0.0;    
    histo_cside_num[i]=histo_cside_den[i]=0.0;    
  }

  // loop over events
  cout << " Entries " << fAnaTree->GetEntries() << endl;
  for(Int_t i=0;i<fAnaTree->GetEntries();i++) {
    // get event
    fAnaTree->GetEntry(i);

    // only good runs
    Int_t run_idx = fRunNum-FirstRun;
    if (run_idx < 0 ) continue;
    if ((run_idx-1)>nRuns) continue;    
    if (GoodRunsIdx[run_idx] < 0) continue;
    if (TRDRunsIdx[run_idx] > -0.5) continue;
  
    //only ctrue-b events
    if (fCtrue!=1) continue;

    // get trigger flags
    Bool_t f0VBA = fL0inputs & 1 << (inputId_0VBA-1);
    Bool_t f0VBC = fL0inputs & 1 << (inputId_0VBC-1);
    Bool_t f0UBA = fL0inputs & 1 << (inputId_0UBA-1);
    Bool_t f0UBC = fL0inputs & 1 << (inputId_0UBC-1);
    Bool_t f0UGC = fL0inputs & 1 << (inputId_0UGC-1);
    Bool_t f0VGA = fL0inputs & 1 << (inputId_0VGA-1);
    Bool_t f0SH2 = fL0inputs & 1 << (inputId_0SH2-1);
    Bool_t f0VC5 = fL0inputs & 1 << (inputId_0VC5-1);

    // set ZDC variables
    Bool_t ZNAhit = isZDCtime(fZNATDC);
    Bool_t ZNChit = isZDCtime(fZNCTDC);

    // no activity at mid rap
    Bool_t noTrklts = (fTracklets==0);
    //noTrklts = kTRUE;
    
    // compute flags
    Bool_t c_side_empty = (!ZNChit && !f0UBC && !f0VBC && noTrklts);
    Bool_t a_side_signal = (ZNAhit || f0UBA || f0VBA || (fADADecision == 1)
			    || (fV0ADecision == 1) /*|| (fTracklets>0) */);

    Bool_t a_side_empty = (!ZNAhit && !f0UBA && !f0VBA && noTrklts);
    Bool_t c_side_signal = (ZNChit || f0UBC || f0VBC || (fADCDecision == 1)
			    || (fV0CDecision == 1) || f0UGC
			    || f0VGA || f0SH2 || f0VC5 /*|| (fTracklets>0) */);
    // update counters
    if (c_side_empty) {
      // zna
      count_zna_den+=1.0;
      histo_zna_den[fRunNum-FirstRun]+=1.0;
      if (ZNAhit) {
	count_zna_num+=1.0;
	histo_zna_num[fRunNum-FirstRun]+=1.0;
      }
      // ada  BB
      count_ada_den+=1.0;
      histo_ada_den[fRunNum-FirstRun]+=1.0;
     if (f0UBA || (fADADecision == 1)) {
    	count_ada_num+=1.0;
	histo_ada_num[fRunNum-FirstRun]+=1.0;
      }
      // v0a  BB
      count_v0a_den+=1.0;
      histo_v0a_den[fRunNum-FirstRun]+=1.0;
      if (f0VBA || (fV0ADecision == 1)) {
 	count_v0a_num+=1.0;
	histo_v0a_num[fRunNum-FirstRun]+=1.0;
      }
      // a_side signal
      count_aside_den+=1.0;
      histo_aside_den[fRunNum-FirstRun]+=1.0;
      if (a_side_signal) {
	count_aside_num+=1.0;
	histo_aside_num[fRunNum-FirstRun]+=1.0;
      }
    } // end of c_side_empty

    if (a_side_empty) {
      // znc
      count_znc_den+=1.0;
      histo_znc_den[fRunNum-FirstRun]+=1.0;
      if (ZNChit) {
	count_znc_num+=1.0;
	histo_znc_num[fRunNum-FirstRun]+=1.0;
      }
      // ada  BB
      count_adc_den+=1.0;
      histo_adc_den[fRunNum-FirstRun]+=1.0;
     if (f0UBC || (fADCDecision == 1)) {
    	count_adc_num+=1.0;
	histo_adc_num[fRunNum-FirstRun]+=1.0;
      }
      // v0c  BB
      count_v0c_den+=1.0;
      histo_v0c_den[fRunNum-FirstRun]+=1.0;
      if (f0VBC || (fV0CDecision == 1)) {
	count_v0c_num+=1.0;
	histo_v0c_num[fRunNum-FirstRun]+=1.0;
      }
      // c_side signal
      count_cside_den+=1.0;
      histo_cside_den[fRunNum-FirstRun]+=1.0;
      if (c_side_signal) {
	count_cside_num+=1.0;
	histo_cside_num[fRunNum-FirstRun]+=1.0;
      }
    } // end of a_side_empty
  } // end of loop over events

  // make plot
    DoPlot("zna",histo_zna_num,histo_zna_den);
    DoPlot("ada",histo_ada_num,histo_ada_den);
    DoPlot("v0a",histo_v0a_num,histo_v0a_den);
    DoPlot("aside",histo_aside_num,histo_aside_den);
    DoPlot("znc",histo_znc_num,histo_znc_den);
    DoPlot("adc",histo_adc_num,histo_adc_den);
    DoPlot("v0c",histo_v0c_num,histo_v0c_den);
    DoPlot("cside",histo_cside_num,histo_cside_den);

    /*
  // print out (not really needed, just for cross checks)
  if (per == 1) {
    Printout("zna", count_zna_den,count_zna_num);
    Printout("ada", count_ada_den,count_ada_num);
    Printout("v0a", count_v0a_den,count_v0a_num);
    Printout("aside", count_aside_den,count_aside_num);
  } else {
    Printout("znc", count_znc_den,count_znc_num);
    Printout("adc", count_adc_den,count_adc_num);
    Printout("v0c", count_v0c_den,count_v0c_num);
    Printout("cside", count_cside_den,count_cside_num);
  }
    */
}
