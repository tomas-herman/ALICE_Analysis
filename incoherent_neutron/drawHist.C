//
// this program draws histograms of gamma gamma pt templates
//

// -----------------------------------------------------------------
// all headers are defined here

// c++ headers
#include <iostream>
#include <fstream>
#include <string> 

// root headers
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>

#include "utilities.h"

#include "gammaPDFextractionParam/valGammaGammaPDF_iSel_1_ZN_0.h"
#include "gammaPDFextractionParam/valGammaGammaPDF_iSel_1_ZN_1.h"
#include "gammaPDFextractionParam/valGammaGammaPDF_iSel_1_ZN_2.h"
#include "gammaPDFextractionParam/valGammaGammaPDF_iSel_1_ZN_3.h"
#include "gammaPDFextractionParam/valGammaGammaPDF_iSel_1_ZN_5.h"
#include "gammaPDFextractionParam/errGammaGammaPDF_iSel_1_ZN_0.h"
#include "gammaPDFextractionParam/errGammaGammaPDF_iSel_1_ZN_1.h"
#include "gammaPDFextractionParam/errGammaGammaPDF_iSel_1_ZN_2.h"
#include "gammaPDFextractionParam/errGammaGammaPDF_iSel_1_ZN_3.h"
#include "gammaPDFextractionParam/errGammaGammaPDF_iSel_1_ZN_5.h"
const int iSel = 1;

// #include "gammaPDFextractionParam/valGammaGammaPDF_iSel_0_ZN_0.h"
// #include "gammaPDFextractionParam/valGammaGammaPDF_iSel_0_ZN_1.h"
// #include "gammaPDFextractionParam/valGammaGammaPDF_iSel_0_ZN_2.h"
// #include "gammaPDFextractionParam/valGammaGammaPDF_iSel_0_ZN_3.h"
// #include "gammaPDFextractionParam/valGammaGammaPDF_iSel_0_ZN_5.h"
// #include "gammaPDFextractionParam/errGammaGammaPDF_iSel_0_ZN_0.h"
// #include "gammaPDFextractionParam/errGammaGammaPDF_iSel_0_ZN_1.h"
// #include "gammaPDFextractionParam/errGammaGammaPDF_iSel_0_ZN_2.h"
// #include "gammaPDFextractionParam/errGammaGammaPDF_iSel_0_ZN_3.h"
// #include "gammaPDFextractionParam/errGammaGammaPDF_iSel_0_ZN_5.h"
// const int iSel = 0;


//-----------------------------------------------
// Draw one plot
void drawOneHistNorm(int znSelection, TH1D *h, TLegend *legend)
{
  
  for (auto bin = 1; bin < nPtBins+1; bin++) {
    if (znSelection==0) {
      h->SetBinContent(bin,znGammaGammaVal0[bin-1]);
      h->SetBinError(bin,znGammaGammaErr0[bin-1]);
    }
    if (znSelection==1) {
      h->SetBinContent(bin,znGammaGammaVal1[bin-1]);
      h->SetBinError(bin,znGammaGammaErr1[bin-1]);
    }
    if (znSelection==2) {
      h->SetBinContent(bin,znGammaGammaVal2[bin-1]);
      h->SetBinError(bin,znGammaGammaErr2[bin-1]);
    }
    if (znSelection==3) {
      h->SetBinContent(bin,znGammaGammaVal3[bin-1]);
      h->SetBinError(bin,znGammaGammaErr3[bin-1]);
    }
    if (znSelection==5) {
      h->SetBinContent(bin,znGammaGammaVal5[bin-1]);
      h->SetBinError(bin,znGammaGammaErr5[bin-1]);
    }
  }

  // set up the output file to save the pt hist
  TFile *fOut = new TFile(Form("ptHistos/ptHisto_gamma_%i.root",znSelection),"recreate");

  double factor = 1.;
  if (h->GetSumw2N() == 0) h->Sumw2(kTRUE);
  h->Scale(factor/h->Integral(), "width");
  h->Draw("HIST, E, SAME");
  h->SetLineWidth(1);
  
  // write the histo to the output file
  fOut->cd();
  h->Write();
  fOut->Close();
  
  if (znSelection==0) {
    h->SetLineColor(kBlue+1);
    legend->AddEntry(h,"0N0N","f");
  }
  if (znSelection==1) {
    h->SetLineColor(kRed+1);
    legend->AddEntry(h,"0NXN","f");
  }
  if (znSelection==2) {
    h->SetLineColor(kCyan+1);
    legend->AddEntry(h,"XN0N","f");
  }
  if (znSelection==3) {
    h->SetLineColor(kOrange);
    legend->AddEntry(h,"XNXN","f");
  }
  if (znSelection==5) {
    h->SetLineColor(kBlack);
    legend->AddEntry(h,"any neutron","f");
  }
}

//-----------------------------------------------
// Draw the final comparison
void drawHist(int znSelection = 3)
{
  TCanvas *c1 = new TCanvas("c1","Gamm Gamma PDFs",1200,800);

  c1->cd();
  // Set general options
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTextSize(0.04); 
 
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  // pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pa

  auto legend = new TLegend(0.8,0.7,0.9,0.9);

  // Set pt binning
  setPtBinning(znSelection);
  double * ptBinBoundariesArr = &ptBinBoundariesVec[0];

  // get the input tree
  TTree *SLtree = getTree(13, iSel);
  double pt;
  SLtree->SetBranchAddress("pt",&pt);
  TH1D *SLh = new TH1D("h","gamma pt template",nPtBins,ptBinBoundariesArr);; 

  int nEntries = SLtree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    SLtree->GetEntry(i);
    SLh->Fill(pt);
  }

  SLh->SetLineColor(kMagenta+1);
  SLh->SetLineWidth(1);
  SLh->GetXaxis()->SetTitle("pt [GeV/c]");
  SLh->GetYaxis()->SetTitle("PDF norm to integral*bin width");
  double factor = 1.;
  if (SLh->GetSumw2N() == 0) SLh->Sumw2(kTRUE);
  SLh->Scale(factor/SLh->Integral(), "width");
  SLh->Draw("HIST, E, SAME");
  SLh->SetMinimum(0.0);  // Define Y minimum
  legend->AddEntry(SLh,"STARlight","f");

  TH1D *h0 = new TH1D("h0","gamma0 pt template",nPtBins,ptBinBoundariesArr);
  TH1D *h1 = new TH1D("h1","gamma1 pt template",nPtBins,ptBinBoundariesArr);
  TH1D *h2 = new TH1D("h2","gamma2 pt template",nPtBins,ptBinBoundariesArr);
  TH1D *h3 = new TH1D("h3","gamma3 pt template",nPtBins,ptBinBoundariesArr);
  TH1D *h5 = new TH1D("h5","gamma5 pt template",nPtBins,ptBinBoundariesArr);

  if (znSelection==0) drawOneHistNorm(0, h0, legend);
  if (znSelection==1) drawOneHistNorm(1, h1, legend);
  if (znSelection==2) drawOneHistNorm(2, h2, legend);
  if (znSelection==3) drawOneHistNorm(3, h3, legend);
  if (znSelection==5) drawOneHistNorm(5, h5, legend);

  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw("SAME");

  // lower Ratio plot will be in pad2
  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  // pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  gPad->SetLogy();

  TH1F *hr;
  // Define the ratio plot
  hr = (TH1F*)SLh->Clone("hr");

  hr->SetLineColor(kBlack);
  hr->SetMinimum(0.02);  // Define Y ..
  hr->SetMaximum(20.); // .. range
  hr->SetStats(0);      // No statistics on lower plot
  if (znSelection==0) hr->Divide(h0);
  if (znSelection==1) hr->Divide(h1);
  if (znSelection==2) hr->Divide(h2);
  if (znSelection==3) hr->Divide(h3);
  if (znSelection==5) hr->Divide(h5);
  // hr->SetMarkerStyle(21);
  hr->GetYaxis()->SetTitleOffset(0.45);
  hr->GetYaxis()->SetLabelSize(0.10);
  hr->GetXaxis()->SetLabelSize(0.10);
  hr->GetYaxis()->SetTitleSize(0.10);
  hr->GetXaxis()->SetTitleSize(0.10);
  hr->GetXaxis()->SetTitle("");
  hr->GetYaxis()->SetTitle("ratio SL/data");
  hr->Draw("");       // Draw the ratio plot

  //---Save as pdf
  char FileName_pdf[120];
  sprintf(FileName_pdf,"gammaTemplates_iSel_%i_ZN_%i.png", iSel, znSelection);

  c1->SaveAs(FileName_pdf);
}