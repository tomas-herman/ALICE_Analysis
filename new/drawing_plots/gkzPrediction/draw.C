#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"

TGraph* GetShade(TGraph* gmin, TGraph* gmax, Color_t color){
  Int_t n = gmin->GetN();
  Double_t xmin = 0.0001;
  Double_t xmax = 0.1;
  Double_t dx = (xmax-xmin)/n;
  TGraph *gshade = new TGraph(2*n);
  for (int i=0;i<n;i++){
    double x    = gmin->GetX()[i];
    double ymin = gmin->GetY()[i];
    double ymax = gmax->GetY()[i];
    gshade->SetPoint(i      ,x,ymin);
    gshade->SetPoint(2*n-i-1,x,ymax);
  }
  gshade->SetFillColor(color);
  gshade->SetLineColor(color);
  return gshade;
}

void draw(){
  TFile* f = new TFile("gkz.root");

  TGraph* gltaw = (TGraph*) f->Get("cs_LTA_weak");
  TGraph* gltas = (TGraph*) f->Get("cs_LTA_strong");
  TGraph* gepsw = (TGraph*) f->Get("cs_EPS09_weak");
  TGraph* gepsc = (TGraph*) f->Get("cs_EPS09_central");
  TGraph* gepss = (TGraph*) f->Get("cs_EPS09_strong");

  TGraph* grltaw = (TGraph*) f->Get("r_LTA_weak");
  TGraph* grltas = (TGraph*) f->Get("r_LTA_strong");
  TGraph* grepsw = (TGraph*) f->Get("r_EPS09_weak");
  TGraph* grepsc = (TGraph*) f->Get("r_EPS09_central");
  TGraph* grepss = (TGraph*) f->Get("r_EPS09_strong");

  new TCanvas;
  TH1F* h1 = gPad->DrawFrame(10,1e-3,1000,0.2);
  h1->SetTitle(";W_{#gammaN};#sigma_{#gammaN}(W_{#gammaN}),mb");
  gPad->SetLogx();
  gPad->SetLogy();
  TGraph* g_lta_unc = GetShade(gltaw,gltas,kGreen-10);
  g_lta_unc->SetFillStyle(1001);
  g_lta_unc->Draw("f");
  gltaw->Draw("l");
  gltas->Draw("l");

  new TCanvas;
  TH1F* h2 = gPad->DrawFrame(10,1e-3,1000,0.2);
  h2->SetTitle(";W_{#gammaN};#sigma_{#gammaN}(W_{#gammaN}),mb");
  gPad->SetLogx();
  gPad->SetLogy();
  TGraph* g_eps_unc = GetShade(gepsw,gepss,kGreen-10);
  g_eps_unc->SetFillStyle(1001);
  g_eps_unc->Draw("f");
  gepsw->Draw("l");
  gepsc->Draw("l");
  gepss->Draw("l");

  new TCanvas;
  TH1F* h3 = gPad->DrawFrame(1e-5,0,1e-1,1.2);
  TGraph* g_rlta_unc = GetShade(grltaw,grltas,kGreen-10);
  g_rlta_unc->SetFillStyle(1001);
  g_rlta_unc->Draw("f");
  gPad->SetLogx();
  h3->SetTitle(";x;R_{Pb}(x,Q^{2})");
  grltaw->Draw("l");
  grltas->Draw("l");
 
  new TCanvas;
  TH1F* h4 = gPad->DrawFrame(1e-5,0,1e-1,1.2);
  gPad->SetLogx();
  h4->SetTitle(";x;R_{Pb}(x,Q^{2})");
  TGraph* g_reps_unc = GetShade(grepsw,grepss,kGreen-10);
  g_reps_unc->SetFillStyle(1001);
  g_reps_unc->Draw("f");
  grepsw->Draw("l");
  grepsc->Draw("l");
  grepss->Draw("l");

}