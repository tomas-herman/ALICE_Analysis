//
// program to fit data on coherent j/psi photoproduction in PbPb
// in forward neutron classes to extract the gPb cross section
//
// It uses Eqs 14-16 of
// https://inspirehep.net/record/1491190?ln=en
// to obtain gPb cross sections and then a fit to a constant
//


//-------------------------------------------------------
// Header files from C++ and ROOT
//-------------------------------------------------------

// c++ headers
#include <iostream>
#include <fstream>
using std::cout;

// root headers
#include <TMath.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>

// -------------------------------------------------------------------------------
//  data (global variables)
// -------------------------------------------------------------------------------
#include "inputData.h"

// -------------------------------------------------------------------------------
//  other global variables
// -------------------------------------------------------------------------------
double totalUncorErr[nData];

// -------------------------------------------------------------------------------
// set the fluxes
// -------------------------------------------------------------------------------

void setFluxes(int n, double* fL, double *fL_src, double* fH, double* fH_src)
{
  for(int i=0;i<n;i++) {
    fL[i] = fL_src[i];
    fH[i] = fH_src[i];
  }
}

// -------------------------------------------------------------------------------
//  get low and high energy cross section from two UPC cross sections
// -------------------------------------------------------------------------------
void getXS(int i, int j, double *XS, double *XSerr)
{
  // get cross section
  double F = fluxL[i]*fluxH[j]-fluxL[j]*fluxH[i];
  XS[0] = (fluxH[j]*crossSection[i]-fluxH[i]*crossSection[j])/F;  
  XS[1] = (fluxL[i]*crossSection[j]-fluxL[j]*crossSection[i])/F;
  // compute errors (TBD)
  double e1 = fluxH[j]*totalUncorErr[i]*fluxH[j]*totalUncorErr[i];
  double e2 = fluxH[i]*totalUncorErr[j]*fluxH[i]*totalUncorErr[j];  
  XSerr[0] = TMath::Sqrt(e1+e2)/TMath::Abs(F);
  double e3 = fluxL[i]*totalUncorErr[j]*fluxL[i]*totalUncorErr[j];
  double e4 = fluxL[j]*totalUncorErr[i]*fluxL[j]*totalUncorErr[i];  
  XSerr[1] = TMath::Sqrt(e3+e4)/TMath::Abs(F);
  // print info
  cout << " ===>   fluxL[i] = " << fluxL[i] << " totalUncorErr[j] = " << totalUncorErr[j] << " fluxL[j] = " << fluxL[j] << " totalUncorErr[i] = " << totalUncorErr[i] << endl;
  cout << " ===>   F = " << F << " e1 = " << e1 << " e2 = " << e2 << " erL = " << XSerr[0] << endl;
  cout << " ===>   F = " << F << " e3 = " << e3 << " e4 = " << e4 << " erH = " << XSerr[1] << endl;  
  // change units
  XS[0]*=1000;
  XS[1]*=1000;  
  XSerr[0]*=1000;
  XSerr[1]*=1000;
}

// -------------------------------------------------------------------------------
//  fill histos
// -------------------------------------------------------------------------------
void fillHisto(int i, int j, int k, TH1F *hL, TH1F *hH)
{
  // to store current cross sections
  double gPbXS[2];
  double gPbXSerr[2];
  // get cross sections
  getXS(i,j,gPbXS,gPbXSerr);
  // update graphs
  hL->SetBinContent(k,gPbXS[0]);
  hL->SetBinError(k,gPbXSerr[0]);
  hH->SetBinContent(k,gPbXS[1]);
  hH->SetBinError(k,gPbXSerr[1]);
}

// -------------------------------------------------------------------------------
// get kinmatics for run2 PbPb collisions!
// -------------------------------------------------------------------------------
void getW(double y1, double y2, double *W)
{
  double mass = 3.096916; // GeV/c^2
  double Eb = 0.5*5020; // GeV
  double y = 0.5*(y1+y2);
  W[0] = TMath::Sqrt(mass*TMath::Exp(-y)*2.0*Eb);
  W[1] = TMath::Sqrt(mass*TMath::Exp(+y)*2.0*Eb);  
}

// -------------------------------------------------------------------------------
//  plot the result
// -------------------------------------------------------------------------------
void drawFit(const char *name, TH1F *histo)
{
  // draw
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas *c = new TCanvas(name,name,0,0,600,600);
  c->cd();
  histo->SetMarkerStyle(20);
  histo->SetYTitle("#sigma (mb)");
  histo->Draw("");
  // save
  c->SaveAs(Form("Fits/%s.pdf",name));
}

// -------------------------------------------------------------------------------
//  extract the cross sections at mid rapidity
// -------------------------------------------------------------------------------
void getXSmid()
{
  // fill the cross sections from the neutron classes
  TH1F *hXS = new TH1F("hXS","hXS",nData,-0.5,nData-0.5);
  const char *names[nData] = {"0n0n","0nXn+Xn0n","XnXn"};
  for (int i=0;i<nData;i++) {
    double xs = 1000*crossSection[i]/(2.0*fluxL[i]);
    double xse = 1000*totalUncorErr[i]/(2.0*fluxL[i]);    
    hXS->GetXaxis()->SetBinLabel(i+1,names[i]);
    hXS->SetBinContent(i+1,xs);
    hXS->SetBinError(i+1,xse);
    cout << names[i] << " xs = " << xs << " +/- " << xse << endl;
  }

  // fit
  hXS->Fit("pol0");

  // get energy
  double W[2];
  getW(yMin,yMax,W);
  hXS->SetTitle(Form("W=%4.0f GeV",W[0]));

  // draw
  drawFit(Form("W%.0fGeV",W[0]),hXS);
  
}

// -------------------------------------------------------------------------------
//  extract the cross sections pair wise
// -------------------------------------------------------------------------------
void getXS()
{
  // to store results for p0 fit 
  int nXS = 0;
  TH1F *hLowXS = NULL;
  TH1F *hHighXS = NULL;    
  
  // fill the graphs
  if (nData == 3) {
    const int nPairs = 3;
    int idx[3][2] = {{0,1},{0,2},{1,2}};
    const char *names[nPairs] = {"0n0n,0nXn","0n0n,XnXn","0nXn,XnXn"};
    hLowXS = new TH1F("hLowXS","hLowXS",nPairs,-0.5,nPairs-0.5);
    hHighXS = new TH1F("hHighXS","hHighXS",nPairs,-0.5,nPairs-0.5);      
    for (int i=0;i<nPairs;i++) {
      fillHisto(idx[i][0],idx[i][1],i+1,hLowXS,hHighXS);
      hLowXS->GetXaxis()->SetBinLabel(i+1,names[i]);
      hHighXS->GetXaxis()->SetBinLabel(i+1,names[i]);
      cout << names[i]
	   << " xsL = " << hLowXS->GetBinContent(i+1)
	   << " +/- " << hLowXS->GetBinError(i+1) << endl;
      cout << names[i]
	   << " xsH = " << hHighXS->GetBinContent(i+1)
	   << " +/- " << hHighXS->GetBinError(i+1) << endl;
    }
  } else if (nData == 4) {
    const int nPairs = 5;
    int idx[5][2] = {{0,1},{0,2},{0,3},{1,3},{2,3}};
    const char *names[nPairs] = {"0n0n,0nXn","0n0n,Xn0n","0n0n,XnXn","0nXn,XnXn","Xn0n,XnXn"};
    hLowXS = new TH1F("hLowXS","hLowXS",nPairs,-0.5,nPairs-0.5);
    hHighXS = new TH1F("hHighXS","hHighXS",nPairs,-0.5,nPairs-0.5);      
    for (int i=0;i<nPairs;i++) {
      fillHisto(idx[i][0],idx[i][1],i+1,hLowXS,hHighXS);
      hLowXS->GetXaxis()->SetBinLabel(i+1,names[i]);
      hHighXS->GetXaxis()->SetBinLabel(i+1,names[i]);
      cout << names[i]
	   << " xsL = " << hLowXS->GetBinContent(i+1)
	   << " +/- " << hLowXS->GetBinError(i+1) << endl;
      cout << names[i]
	   << " xsH = " << hHighXS->GetBinContent(i+1)
	   << " +/- " << hHighXS->GetBinError(i+1) << endl;
    }
  } else {
    cout << " number of data points not recognized " << endl;
    return;
  }
  
  // fit
  hLowXS->Fit("pol0");
  hHighXS->Fit("pol0");
  
  // get energy
  double W[2];
  getW(yMin,yMax,W);
  hLowXS->SetTitle(Form("W=%4.2f GeV",W[0]));
  hHighXS->SetTitle(Form("W=%4.2f GeV",W[1]));
  
  // draw
  drawFit(Form("W%.0fGeV",W[0]),hLowXS);
  drawFit(Form("W%.0fGeV",W[1]),hHighXS);
}

// -------------------------------------------------------------------------------
// Entry point
// -------------------------------------------------------------------------------

void fitgPb2()
{
  // set the total uncorrelated error to be used
  bool onlySTAT = kFALSE;
  for(int j=0;j<nData;j++) {
    if (onlySTAT) totalUncorErr[j]=statErr[j];
    else {
      double unc = 0;
      for(int i=0;i<nUncorErr;i++) unc+=TMath::Power(uncorErr[i][j],2);
      // convert to cross section
      unc = TMath::Sqrt(unc)*crossSection[j];
      totalUncorErr[j]=TMath::Sqrt((unc*unc)+(statErr[j]*statErr[j]));
    }
  }

  // set up the fluxes
  setFluxes(nData,fluxL,fluxL_noon_mn,fluxH,fluxH_noon_mn);
  
  // get cross section
  if (yMin>0) getXS();
  else getXSmid();

     
}
