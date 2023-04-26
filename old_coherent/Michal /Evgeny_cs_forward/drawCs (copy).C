#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

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


void SetStyle(TH1D* h, Color_t color, Style_t style, Width_t width=3){
  h->SetLineColor(color);
  h->SetLineStyle(style);
  h->SetLineWidth(width);
}

void SetStyle(TGraph* g, Color_t color, Style_t style, Width_t width=3){
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
}

void drawCs(){
  gStyle->SetErrorX(0.01);
  gStyle->SetLineScalePS(2.3);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetTitleOffset(1.25,"XYZ");
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);

  TFile* fexp = new TFile("cs.root");
  TGraphErrors* gStat_data = (TGraphErrors*) fexp->Get("gStat_data_average");
  TGraphErrors* gSyst_data = (TGraphErrors*) fexp->Get("gSyst_data_average");
  
  gStat_data->SetMarkerColor(kBlack);
  gStat_data->SetMarkerStyle(kFullCircle);
  gStat_data->SetMarkerSize(0.7);
  gStat_data->SetLineColor(kBlack);
  gStat_data->SetLineWidth(1);

  gSyst_data->SetMarkerColor(kBlack);
  gSyst_data->SetMarkerStyle(kFullCircle);
  gSyst_data->SetMarkerSize(0.6);
  gSyst_data->SetLineColor(kBlack);
  gSyst_data->SetLineWidth(1);
  gSyst_data->SetFillStyle(1);

  TFile* fHepdata = new TFile("hepdata.root");
  fHepdata->cd("Table 1");
  fHepdata->ls();
  TGraphAsymmErrors* gHepData = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y1");

  TFile* fTheory = new TFile("predictions.root");
//  TGraph* g_lta_unc        = (TGraph*) fTheory->Get("lta_unc_cs_total");
//  TGraph* g_eps_unc        = (TGraph*) fTheory->Get("eps_unc_cs_total");
  TGraph* g_lta_l          = (TGraph*) fTheory->Get("glta_l");
  TGraph* g_lta_h          = (TGraph*) fTheory->Get("glta_h");
  TGraph* g_eps_l          = (TGraph*) fTheory->Get("geps_l");
  TGraph* g_eps_h          = (TGraph*) fTheory->Get("geps_h");
  TGraph* g_eps_val        = (TGraph*) fTheory->Get("geps_0");
  TGraph* g_lta_unc = GetShade(g_lta_l,g_lta_h,kWhite);
  TGraph* g_eps_unc = GetShade(g_eps_l,g_eps_h,kGreen-11);
  
  TH1D* hsljpsicoh         =  (TH1D*) fTheory->Get("hsljpsicoh");
  TH1D* hsljpsicoh_nonucl  =  (TH1D*) fTheory->Get("hsljpsicoh_nonucl");
  TH1D* hgmjpsicoh_iim_bg  =  (TH1D*) fTheory->Get("hgmjpsicoh_iim_bg");
  TH1D* hgmjpsicoh_bcgc_lc =  (TH1D*) fTheory->Get("hgmjpsicoh_bcgc_lc"); 
  TH1D* hipsat             =  (TH1D*) fTheory->Get("hhjpsicoh");
  TGraph* g_cck_gghs       = (TGraph*) fTheory->Get("g_cck_gghs");
  TGraph* g_cck_gshs       = (TGraph*) fTheory->Get("g_cck_gshs");
  TGraph* g_ls             = (TGraph*) fTheory->Get("g_ls");
  
  SetStyle(g_ls              ,kBlue-2  ,10);
  SetStyle(g_cck_gghs        ,kRed-4   ,7);
  SetStyle(g_cck_gshs        ,kRed+2   ,9);
  SetStyle(hsljpsicoh        ,kBlue    ,4);
  SetStyle(hsljpsicoh_nonucl ,kAzure+1 ,7);
  SetStyle(hgmjpsicoh_iim_bg ,kMagenta ,8);
  SetStyle(hgmjpsicoh_bcgc_lc,kViolet+1,6);
  SetStyle(hipsat            ,kOrange-6,9);
  SetStyle(g_lta_h           ,kGray+3  ,5);
  SetStyle(g_eps_val         ,kGreen+2 ,1);
  SetStyle(g_eps_unc         ,kGreen+2 ,1);
  SetStyle(gStat_data        ,kBlack   ,1);
  SetStyle(gSyst_data        ,kBlack   ,1);
  SetStyle(g_lta_unc         ,kGray+3  ,7);
  g_eps_unc->SetFillStyle(1001);
  g_eps_unc->SetFillColor(kGreen-10);


  TCanvas *c = new TCanvas("c","dsigma/dy",900,800);
  c->SetLeftMargin(0.10);
  c->SetBottomMargin(0.14);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.03);

  TH1F* frame = gPad->DrawFrame(-4.201,0,-2.199,7);
//  TH1F* frame = gPad->DrawFrame(-4.9,0,-1.9,8);
  frame->SetTitle(";y;d#sigma/dy (mb)");
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(0.7);
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetXaxis()->SetTitleFont(42);
  frame->GetYaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetLabelFont(42);
  frame->GetYaxis()->SetLabelFont(42);
  frame->GetXaxis()->SetNdivisions(210);
  frame->GetYaxis()->SetDecimals();
  
  gStat_data->SetLineWidth(2);
  gSyst_data->SetFillStyle(0);

  TGraphErrors* gSyst_legend = (TGraphErrors*) gSyst_data->Clone();
  gSyst_legend->Set(1);
  gSyst_legend->SetPoint(0,-4.035,6.17);
  gSyst_legend->SetPointError(0,0.1,0.15);

  gStyle->SetEndErrorSize(0);
  g_eps_unc->Draw("fsame");
  g_lta_h->Draw("same");
//  g_lta_unc->Draw("same");
  g_eps_val->Draw("same");
  hsljpsicoh->Draw("same c");
  hsljpsicoh_nonucl->Draw("same c");
  hgmjpsicoh_iim_bg->Draw("same c");
  printf("%f\n",hsljpsicoh_nonucl->GetBinContent(hsljpsicoh_nonucl->FindFixBin(2.5-0.125))); //3.875
//  hgmjpsicoh_bcgc_lc->Draw("same c");
  hipsat->Draw("same c"); 
  g_cck_gghs->Draw("same");
//  g_cck_gshs->Draw("same");
  g_ls->Draw("same");
  gSyst_legend->Draw("2 same");
//  gStat_legend->Draw("p");
  gSyst_data->Draw("2 same");
  gStat_data->Draw("p");
//  gHepData->Draw("p same");

  gPad->RedrawAxis();
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.55,0.94,"ALICE Pb+Pb #rightarrow Pb+Pb+J/#psi   #sqrt{#it{s}_{NN}} = 5.02 TeV");
  
  TLegend* l = new TLegend(0.13,0.55,0.60,0.90);
//  TLegend* l = new TLegend(0.14,0.52,0.55,0.92);
  l->SetMargin(0.18);
  l->SetFillStyle(0);
  l->AddEntry(gStat_data,"ALICE coherent J/#psi","p");
  l->AddEntry(hsljpsicoh_nonucl,"Impulse approximation","l");
  l->AddEntry(hsljpsicoh,"STARLIGHT","l");
  l->AddEntry(g_eps_val,"EPS09 LO (GKZ)","l");
  l->AddEntry(g_lta_h,"LTA (GKZ)","l");
//  l->AddEntry(g_cck_gshs,"GS-HS (CCK)","pl");
  l->AddEntry(hgmjpsicoh_iim_bg,"IIM BG (GM)","l");
  l->AddEntry(hipsat,"IPsat (LM)","l");
  l->AddEntry(g_ls,"BGK-I (LS)","l");
  l->AddEntry(g_cck_gghs,"GG-HS (CCK)","l");
//  l->AddEntry(hgmjpsicoh_bcgc_lc,"CGC (GM BCGC LC)","pl");
  l->Draw();
  gPad->Print("cs.png");
  gPad->Print("cs.pdf");
}
