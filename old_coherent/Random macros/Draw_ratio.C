// c++ headers
#include <iostream>
#include <fstream>

// ROOT headers
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TLatex.h"

// RooFit headers
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
using namespace RooFit;


// --------------------------------------------------------------------------------
//TTree selection parameters
const Double_t mMin = 2.85; //2.85 3.5
const Double_t mMax = 3.35; //3.35 3.9
const Double_t mrangeMin = 2.; //2.85
const Double_t mrangeMax = 6.; //3.35

const Double_t n_binsM = 80; //150

Double_t scalesMC[7] = {1, 0.0385, 0.1415, 0.2663, 0.2975, 0.2001, 0.0561};
Double_t scaleMC;

const Int_t nYBins = 7;
Float_t gYMin[nYBins] = {-4.0, -4.00, -3.75, -3.50, -3.25, -3.00, -2.75};
Float_t gYMax[nYBins] = {-2.5, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50};
const Int_t YBin = 0;

// void Draw_ratio()
// {

// RooRealVar jRecM("jRecM","jRecM",2,15) ;
// // --------------------------------------------------------------------------------
// // Create CB1 PDF

// RooRealVar mMC("mMC","massMC",3.09867,3.05,3.15) ;
// RooRealVar sigMC("#sigma_{MC}","resolution",0.0854761,0.01,0.2) ;
// RooRealVar n_cbMC("n_{cb,MC}","n_cbMC",93.1194,0,120) ;
// RooRealVar alphaMC("#alpha_{MC}","alphaMC",1.12159,0,10) ;

// mMC.setConstant(kTRUE);
// sigMC.setConstant(kTRUE);
// n_cbMC.setConstant(kTRUE);
// alphaMC.setConstant(kTRUE);

// RooCBShape cb1MC("cb1MC","Crystal Ball PDF",jRecM,mMC,sigMC,alphaMC,n_cbMC) ;

// // --------------------------------------------------------------------------------
// // Create CB2 PDF

// RooFormulaVar m2MC("m2MC","mMC+3.68609-3.096916",RooArgList(mMC));
// RooFormulaVar sig2MC("sig2MC","#sigma_{MC}*1.09",RooArgList(sigMC));

// RooCBShape cb2MC("cb2","Crystal Ball PDF",jRecM,m2MC,sig2MC,alphaMC,n_cbMC) ;

// // --------------------------------------------------------------------------------
// // Create BKGD PDF

// RooRealVar lambdaMC("#lambda_{MC}","exponent",-0.934756,-5.,0.);
// RooRealVar a1MC("a_{1,MC}","parameter a1",0.648354,0,2);
// RooRealVar a2MC("a_{2,MC}","parameter a2",0.900351,0,2);
// RooRealVar a3MC("a_{3,MC}","parameter a3",0.239378,0,2);

// lambdaMC.setConstant(kTRUE);
// a1MC.setConstant(kTRUE);
// a2MC.setConstant(kTRUE);
// a3MC.setConstant(kTRUE);

// RooGenericPdf BkgdMC("BkgdMC","jRecM>4? exp(jRecM*#lambda_{MC}) : exp(jRecM*#lambda_{MC})*(1+a_{1,MC}*pow((jRecM-4),2)+a_{2,MC}*pow((jRecM-4),3)+a_{3,MC}*pow((jRecM-4),4))",RooArgSet(jRecM,lambdaMC,a1MC,a2MC,a3MC));

// // plot results
// gStyle->SetOptTitle(0);
// gStyle->SetOptStat(0);
// // gStyle->SetPadTickX(1);
// // gStyle->SetPadTickY(1);
// gStyle->SetPalette(kBird);
// gStyle->SetPaintTextFormat("4.3f");
// gStyle->SetFrameLineWidth(1);

// gStyle->SetLabelSize(0.045,"xyz");
// gStyle->SetTitleSize(0.05,"xyz");
// gStyle->SetTextSize(0.04);

// TCanvas *c1 = new TCanvas("c1","c1",800,600);
// c1->cd();
// c1->SetLeftMargin(0.14);
// c1->SetBottomMargin(0.11);

// // Extract parameters from Evgeny's file
// TFile* fin = new TFile("EfitM.root");
// TH1D* hM_gl = (TH1D*) fin->Get(Form("gl_%i_%i_1",YBin,0));

// TF1 *fBG;
// fBG  = (TF1*) hM_gl->GetFunction("f2")->Clone("fBG");

// Double_t bE = fBG->GetParameter(1);
// Double_t a2E = fBG->GetParameter(3);
// Double_t a3E = fBG->GetParameter(4);
// Double_t a4E = fBG->GetParameter(5);


// RooRealVar m0("m","mass",3.096,3.05,3.15) ;
// RooRealVar n_cb("n_{cb}","n_cb",4,0,40) ;
// RooRealVar sig("#sigma","resolution",0.092,0.01,0.2) ;
// RooRealVar alpha("#alpha","alpha",0.97,0,20) ;
// RooFormulaVar m02("m02","m+3.68609-3.096916",RooArgList(m0));
// RooFormulaVar sig2("sig2","#sigma*1.09",RooArgList(sig));

// RooCBShape cb("cb","Crystal Ball PDF",jRecM,m0,sig,alpha,n_cb) ;
// RooCBShape cb2("cb2","Crystal Ball PDF",jRecM,m02,sig2,alpha,n_cb) ;

// // create a PDF to describe the background
// RooRealVar lambda("#lambda","exponent",-bE,-5.,0.);
// // lambda.setConstant(kTRUE);

// RooRealVar a1("a_{1}","parameter a1",a2E,0.0001,2);
// RooRealVar a2("a_{2}","parameter a2",a3E,0.0001,2);
// RooRealVar a3("a_{3}","parameter a3",a4E,0.0001,2);
// a1.setConstant(kTRUE);
// a2.setConstant(kTRUE);
// a3.setConstant(kTRUE);

// RooGenericPdf Bkgd("Bkgd","jRecM>4? exp(jRecM*#lambda) : exp(jRecM*#lambda)*(1+a_{1}*pow((jRecM-4),2)+a_{2}*pow((jRecM-4),3)+a_{3}*pow((jRecM-4),4))",RooArgSet(jRecM,lambda,a1,a2,a3));


// for(int i=1; i<7; i++){
//    scaleMC = scalesMC[i];

//    for(int j=0; j<50; j++){

//       RooDataSet *CB1MC = cb1MC.generate(jRecM,scaleMC*20998);
//       RooDataSet *CB2MC = cb2MC.generate(jRecM,scaleMC*512);
//       RooDataSet *bgMC = BkgdMC.generate(jRecM,scaleMC*29008);

//       CB1MC->append(*CB2MC);
//       CB1MC->append(*bgMC);
//       RooAbsData* dataM = CB1MC;

//       // #####################################################################################
//       // Invariant mass fit
//       // #####################################################################################

//       // ------------------------------------
//       Int_t nEvts = dataM->numEntries();

//       RooRealVar nCB("N_{J/#psi}","Number of CB events",0.5*nEvts,0,nEvts);
//       RooRealVar nCB2("N_{#psi'}","Number of CB events",0.1*nEvts,0,nEvts);
//       RooRealVar nBG("N_{bg}","Number of BG events",0.6*nEvts,0,nEvts);

//       RooAddPdf CbBG("CbBG","Crystal Ball plus Background PDF", RooArgList(Bkgd,cb,cb2), RooArgList(nBG,nCB,nCB2));

//       RooFitResult* rM = CbBG.fitTo(*dataM,Extended(kTRUE),Range(2,6),Save());

//       // compute background in range
//       Double_t nBkgd[2];
//       Double_t nCBM[2];
//       Double_t nCBM2[2];
//       jRecM.setRange("JPsiMassRange",mMin,mMax);
//       RooAbsReal *iBG = Bkgd.createIntegral(jRecM,NormSet(jRecM),Range("JPsiMassRange"));
//       RooAbsReal *iCB = cb.createIntegral(jRecM,NormSet(jRecM),Range(2.,6.));
//       RooAbsReal *iCB2 = cb2.createIntegral(jRecM,NormSet(jRecM),Range(2.,6.));

//       nCBM[0] = iCB->getVal()*nCB.getVal();
//       nCBM[1] =  iCB->getVal()*nCB.getError();
//       nCBM2[0] = iCB2->getVal()*nCB2.getVal();
//       nCBM2[1] =  iCB2->getVal()*nCB2.getError();
//       nBkgd[0] =  iBG->getVal()*nBG.getVal();
//       nBkgd[1] =  iBG->getVal()*nBG.getError();


//       TH1F *h = new TH1F("h1", "#psi'/J/#psi", 6,-4.0,-2.5);

//       h->SetBinContent(i,nCBM2[0]/nCBM[0]);
//       h->SetBinError(i,sqrt( ( (nCBM2[1]/nCBM[0])*(nCBM2[1]/nCBM[0]) ) +
//          ( (nCBM2[0]/nCBM[0]*nCBM[1]/nCBM[0])*(nCBM2[0]/nCBM[0]*nCBM[1]/nCBM[0]) ) -
//          ( (nCBM2[1]/nCBM[0])*(nCBM2[0]/nCBM[0]*nCBM[1]/nCBM[0]) ) ));

//       if(i==1 && j==0){
//          h->GetXaxis()->SetTitle("Number of events in given rapidity range");
//          h->GetYaxis()->SetTitle("N_{#psi'}/N_{J/#psi}");
//          h->SetAxisRange(0.,0.062, "Y");
//          // h->GetXaxis()->SetTitleOffset(0.9);
//          // h->GetYaxis()->SetTitleOffset(0.8);
//       }

//       h->SetMarkerColor(kBlue+1);
//       h->Draw("Same");

//       Float_t ymax = h->GetMaximum();

//       // CB1MC->Delete();
//       // CB2MC->Delete();
//       // bgMC->Delete();
//       // dataM->Delete();
//       // rM->Delete();
//       // h->Delete();
//    }
// }
   
//    TLine *line = new TLine(-4,0.0244,-2.5,0.0244);
//    line->SetLineColor(kRed+1);
//    line->SetLineWidth(2);
//    line->Draw();

//    c1->Print("ratioMC.pdf");
// }



void Draw_ratio()
{


// plot results
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
// gStyle->SetPadTickX(1);
// gStyle->SetPadTickY(1);
gStyle->SetPalette(kBird);
gStyle->SetPaintTextFormat("4.3f");
gStyle->SetFrameLineWidth(1);

gStyle->SetLabelSize(0.045,"xyz");
gStyle->SetTitleSize(0.05,"xyz");
gStyle->SetTextSize(0.04);

TCanvas *c1 = new TCanvas("c1","c1",800,600);
c1->cd();
c1->SetLeftMargin(0.14);
c1->SetBottomMargin(0.11);

TH1F *h = new TH1F("h1", "#psi'/J/#psi", 6,-4.0,-2.5);

h->GetXaxis()->SetTitle("Number of events in given rapidity range");
h->GetYaxis()->SetTitle("N_{#psi'}/N_{J/#psi}");
h->SetAxisRange(0.,0.062, "Y");

// ####### Data ####### 
h->SetBinContent(1,0.0527);
h->SetBinError(1,0.0114);

h->SetBinContent(2,0.0279);
h->SetBinError(2,0.0066);

h->SetBinContent(3,0.0181);
h->SetBinError(3,0.0056);

h->SetBinContent(4,0.0111);
h->SetBinError(4,0.0057);

h->SetBinContent(5,0.0357);
h->SetBinError(5,0.0072);

h->SetBinContent(6,0.0417);
h->SetBinError(6,0.0147);

h->SetMarkerColor(kBlue+1);
h->Draw("Same");

Float_t ymax = h->GetMaximum();

    

TLine *line = new TLine(-4,0.0242,-2.5,0.0242);
line->SetLineColor(kRed+1);
line->SetLineWidth(2);
line->Draw();

TBox *b = new TBox(-4,0.0242-0.0039,-2.5,0.0242+0.0039); 
b->SetFillColor(kRed+1); b->SetFillStyle(3002);
b->Draw();

c1-> Print("ratio_data_comment.pdf");
}


   // ####### Data ####### 
   // h->SetBinContent(1,0.0491);
   // h->SetBinError(1,0.0117);

   // h->SetBinContent(2,0.0281);
   // h->SetBinError(2,0.0067);

   // h->SetBinContent(3,0.0191);
   // h->SetBinError(3,0.0056);

   // h->SetBinContent(4,0.0111);
   // h->SetBinError(4,0.0057);

   // h->SetBinContent(5,0.0352);
   // h->SetBinError(5,0.0073);

   // h->SetBinContent(6,0.045);
   // h->SetBinError(6,0.015);

   // ####### MC #######
   // h->SetBinContent(1,0.0419);
   // h->SetBinError(1,0.0153);

   // h->SetBinContent(2,0.0345);
   // h->SetBinError(2,0.0080);

   // h->SetBinContent(3,0.0233);
   // h->SetBinError(3,0.0059);

   // h->SetBinContent(4,0.0268);
   // h->SetBinError(4,0.0055);

   // h->SetBinContent(5,0.0348);
   // h->SetBinError(5,0.0068);

   // h->SetBinContent(6,0.0136);
   // h->SetBinError(6,0.0128);