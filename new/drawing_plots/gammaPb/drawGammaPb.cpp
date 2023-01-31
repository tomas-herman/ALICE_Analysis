// c++ headers
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

// root headers
#include <Rtypes.h>
#include <TMath.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TPDF.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TVirtualFitter.h>
#include <TDatabasePDG.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TString.h>

using namespace std;

//_____________________________________________________________________________
// Global constants
const double mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
//_____________________________________________________________________________
// read model from PHYSICAL REVIEW C 97, 024901 (2018), GG-hs (red solid line)
void readHS(TGraph *gr)
{
  // square of jpsi mass
  
  double jpsiMass2 = mjpsi*mjpsi;
  
  // tmp storage of input data
  double x = 0;
  double sig = 0;

  // reading input file
  ifstream ifs("predictions/fig2_sigma_GGhs_xBj.txt");
  int nlines = 41;
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> x;
    ifs >> sig;
    // convert x into W
    double W = std::sqrt(jpsiMass2/x);
    // fill TGraph
    gr->SetPoint(gr->GetN(),W,sig);
  }
  ifs.close();
}
//_____________________________________________________________________________
// read model from Phys.Lett.B 817 (2021) 136306, b-BK-A and b-BK--GG
void readBK(TGraph *grA, TGraph *grGG)
{
  int nlines = 100;
  
  // tmp storage of input data
  double x = 0;
  double sig = 0;
  double W = 0;

  // reading b-BK-A file
  ifstream ifsA("predictions/Pb_A_jpsi_Q2_0.05_x_W_cs.txt");
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifsA >> x;
    ifsA >> W;
    ifsA >> sig;
    // fill TGraph
    grA->SetPoint(grA->GetN(),W,sig);
  }
  ifsA.close();
  
  // reading b-BK-GG file
  ifstream ifsGG("predictions/Pb_GG_jpsi_Q2_0.05_x_W_cs.txt");
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifsGG >> x;
    ifsGG >> W;
    ifsGG >> sig;
    // fill TGraph
    grGG->SetPoint(grGG->GetN(),W,sig);
  }
  ifsGG.close();
}
//_____________________________________________________________________________
// read model from STARlight Impulse Approximation
void readImpulseApproximation(TGraph *gr)
{
  // tmp storage of input data
  double y = 0;
  double W = 0;
  double flux = 0;
  double sig = 0;
  double dsig = 0;
  double W2 = 0;
  double flux2 = 0;
  double sig2 = 0;
  double dsig2 = 0;
  double boh = 0;

  // reading input file
  ifstream ifs("predictions/ImpulseApproximation.txt");
  int nlines = 820;
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> y;
    ifs >> W;
    ifs >> flux;
    ifs >> sig;
    ifs >> dsig;
    ifs >> W2;
    ifs >> flux2;
    ifs >> sig2;
    ifs >> dsig2;
    ifs >> boh;
    gr->SetPoint(gr->GetN(),W,sig*1000);
  }
  ifs.close();
}
//_____________________________________________________________________________
// read model STARlight 
void readSTARlight(TGraph *gr)
{
  // tmp storage of input data
  double y = 0;
  double W = 0;
  double flux = 0;
  double sig = 0;
  double dsig = 0;
  double W2 = 0;
  double flux2 = 0;
  double sig2 = 0;
  double dsig2 = 0;
  double boh = 0;

  // reading input file
  ifstream ifs("predictions/STARlight-xnxn.txt");
  int nlines = 820;
  for (int i=0;i<nlines;i++) {
    // read values in this line
    ifs >> y;
    ifs >> W;
    ifs >> flux;
    ifs >> sig;
    ifs >> dsig;
    ifs >> W2;
    ifs >> flux2;
    ifs >> sig2;
    ifs >> dsig2;
    ifs >> boh;
    gr->SetPoint(gr->GetN(),W,sig);
  }
  ifs.close();
}
//_____________________________________________________________________________
void SetStyle(TH1D* h, Color_t color, Style_t style, Width_t width = 2)
{
  h->SetLineColor(color);
  h->SetLineStyle(style);
  h->SetLineWidth(width);
}
//_____________________________________________________________________________
void SetStyle(TGraph* g, Color_t color, Style_t style, Width_t width = 2)
{
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
}
//_____________________________________________________________________________
vector<double> getSquareSum(double inputA[], double inputB[], int size) 
{
  vector<double> vec;
  for (int i = 0; i < size; i++)
  {
    vec.push_back( sqrt( pow(inputA[i],2)+pow(inputB[i],2) ) );
  }
  return vec;
}
//_____________________________________________________________________________
TGraph *getAvergaGraph(TGraph *inputA, TGraph *inputB)
{
  // create the difference graph
  vector<double> v;
  v.insert( v.end(), inputA->GetX(), (inputA->GetX() + inputA->GetN()) );
  v.insert( v.end(), inputB->GetX(), (inputB->GetX() + inputB->GetN()) );
  sort( v.begin(), v.end() );
  v.erase( unique( v.begin(), v.end() ), v.end() );
  int n = v.size();
  TGraph *average = new TGraph(n);
  for (Int_t i = 0; i < n; i++)
    average->SetPoint( i, v[i], (inputA->Eval(v[i]) + inputB->Eval(v[i]))/2 );
  return (TGraph*) average;
}
//_____________________________________________________________________________
void drawGammaPb()
{
  // configure layout of the plot
  //--- ranges
  double gxmin = 12., gymin = 3, gxmax = 1200., gymax = 1700;

  //--- style, title...
  TCanvas* c1 = new TCanvas("c1","c1",1000,700);
  TH1F* frame = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
  
  float fSiz = 0.045;
  int iFont = 42;
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetTitleOffset(1.25);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetXaxis()->SetTitleSize(fSiz);
  frame->GetYaxis()->SetTitleSize(fSiz);
  frame->GetXaxis()->SetLabelSize(fSiz);
  frame->GetYaxis()->SetLabelSize(fSiz);
  frame->GetXaxis()->SetTitleFont(iFont);
  frame->GetYaxis()->SetTitleFont(iFont);
  frame->GetXaxis()->SetLabelFont(iFont);
  frame->GetYaxis()->SetLabelFont(iFont);
  frame->SetTitleSize(fSiz);       
  frame->SetLabelSize(fSiz);
  frame->SetTitleSize(fSiz, "Y");  
  frame->SetLabelSize(fSiz, "Y");
  frame->SetTitle(";W_{#gammaPb,n} (GeV);#sigma(#gammaPb) (#mub)");
  
  gPad->SetMargin(0.10,0.05,0.12,0.10);
  gPad->SetLogy();
  gPad->SetLogx();
  
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(2.3);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetTextSize(fSiz);
  gStyle->SetTextFont(iFont);

  TGaxis::SetMaxDigits(3);

  //--- upper horizontal axis with Bjorken-x
  double bxmin = pow((mjpsi/gxmax),2.);
  double bxmax = pow((mjpsi/gxmin),2.);
  TF1 *fbx = new TF1("fbx","pow(([0]/x),2.)", bxmin, bxmax);
  fbx->SetParameter(0, mjpsi);

  TGaxis *axis = new TGaxis(gxmax, gymax, gxmin, gymax, "fbx", 510, "+G");
  axis->SetTextFont(iFont);
  axis->SetLabelFont(iFont);
  axis->SetTitleSize(fSiz); 
  axis->SetLabelSize(fSiz);
  axis->SetLabelOffset(-0.035);
  axis->SetTitle("                                                                                              Bjorken-#it{x}");
  axis->Draw("same");


  // Theory prediction 
  //--- read
  TGraph *predHS = new TGraph(); // prediction from the hot-spot model
  readHS(predHS);
  TGraph *predBKA = new TGraph(); // prediction from the b-BK-A equation
  TGraph *predBKGG = new TGraph(); // prediction from the b-BK-GG equation  
  readBK(predBKA,predBKGG);
  TGraph *predIA = new TGraph(); // prediction from the IA model
  readImpulseApproximation(predIA);
  TGraph *predSL = new TGraph(); // prediction from the STARlight model
  readSTARlight(predSL);
  TFile *f = new TFile("predictions/gkz.root"); // prediction from the EPS09 model
  TGraph *ltaW = (TGraph*) f->Get("cs_LTA_weak");
  TGraph *ltaS = (TGraph*) f->Get("cs_LTA_strong");
  TGraph *epsW = (TGraph*) f->Get("cs_EPS09_weak");
  TGraph *epsC = (TGraph*) f->Get("cs_EPS09_central");
  TGraph *epsS = (TGraph*) f->Get("cs_EPS09_strong");

  TGraph *ltaAverage = getAvergaGraph(ltaW, ltaS);
  
  //--- style and plot
  SetStyle(predHS                 ,kRed-4 ,2);
  predHS->Draw("csame");
  SetStyle(predBKA                ,kRed+2,8);
  predBKA->Draw("csame");
  // SetStyle(predBKGG            ,kCyan ,7);
  // predBKGG->Draw("csame");
  SetStyle(predIA                 ,kAzure+1 ,7);
  predIA->Draw("csame");
  SetStyle(predSL                 ,kBlue-2  ,5);
  predSL->Draw("csame");
  epsC->Scale(1000,"y");     
  SetStyle(epsC                  ,kGreen+2 ,1);
  epsC->Draw("csame");
  ltaAverage->Scale(1000,"y");
  SetStyle(ltaAverage            ,kTeal+2 ,4);
  ltaAverage->Draw("csame");


  // ALICE
  // Input ALICE data
  int nP = 9;
  double ALICE_W[9]    = { 19.12, 813.05, 24.55, 633.21, 31.53, 493.14, 97.11, 160.1, 124.69  };
  double ALICE_sig[9]  = { 8.80, 60.23, 13.86, 47.80, 16.78, 46.00, 20.21, 27.03, 24.10       };
  double ALICE_unc[9]  = { 0.30, 19.94, 0.23, 6.50, 0.59, 6.32, 5.14, 7.40, 0.70              };
  double ALICE_cor[9]  = { 0.67, 8.21, 1.10, 10.76, 1.31, 5.3,3.13, 4.96, 1.36                };
  double ALICE_mig[9]  = { 0.12, 15.22, 0.19, 8.73, 0.35, 6.72, 7.53, 10.94, 0.23             };
  double ALICE_flu[9]  = { 0.06, 8.81,  0.11, 5.16, 0.25, 4.22, 3.81, 5.44, 0.15              };

  // Plotting
  vector<double> rgbRaw = {0,0,0};
  vector<double> rgb = {0,0,0};
  //--- define colors
  int cData = TColor::GetFreeColorIndex();
  rgbRaw = {165, 0, 0};
  for (int i = 0; i < 3; i++) rgb[i] = rgbRaw[i]/255;
  TColor *colorData = new TColor(cData, rgb[0], rgb[1], rgb[2]);
  cData = kBlack;

  int cError = TColor::GetFreeColorIndex();
  rgbRaw = {255, 183, 183};
  for (int i = 0; i < 3; i++) rgb[i] = rgbRaw[i]/255;
  TColor *colorError = new TColor(cError, rgb[0], rgb[1], rgb[2]);
  cError = kGray;

  //--- plot data
  double ySize[nP];
  for (int i = 0; i < nP; i++) ySize[i] = 0.07*ALICE_W[i];

  TGraphErrors *aliceDataErr2  = new TGraphErrors(nP, ALICE_W, ALICE_sig, ySize, ALICE_flu);
  aliceDataErr2->SetFillColor(cError);
  aliceDataErr2->Draw("2");

  TGraphErrors *aliceDataErr1  = new TGraphErrors(nP, ALICE_W, ALICE_sig, ySize, &getSquareSum(ALICE_cor,ALICE_mig,nP)[0]);
  aliceDataErr1->SetFillColor(cData);
  aliceDataErr1->SetFillStyle(0);
  aliceDataErr1->SetLineWidth(2);
  aliceDataErr1->Draw("2");

  TGraphErrors *aliceData = new TGraphErrors(nP, ALICE_W, ALICE_sig, NULL, ALICE_unc);
  aliceData->SetMarkerColor(cData);
  aliceData->SetMarkerStyle(kFullCircle);
  aliceData->SetMarkerSize(1.2);
  aliceData->SetLineColor(cData);
  aliceData->SetLineWidth(2);
  aliceData->Draw("zpsame");


  //legend 
  double xl = 0.14, dxl = 0.2;
  double yl = 0.58, dyl = 0.28;
  TLegend* leg1 = new TLegend(xl, yl, xl+dxl, yl+dyl);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(fSiz-0.005);
  leg1->AddEntry(aliceData,"ALICE","EP");
  leg1->AddEntry(predIA,"Impulse approximation","l");
  leg1->AddEntry(predSL,"STARlight","l");
  leg1->AddEntry(epsC,"EPS09 LO","l");
  leg1->AddEntry(ltaAverage,"LTA","l");
  leg1->AddEntry(predHS,"GG-HS","l"); //(PRC97 (2018) 024901)
  leg1->AddEntry(predBKA,"b-BK-A","l"); // (PLB817(2021)136306)
  // leg1->AddEntry(predBKGG,"b-BK-GG (PLB817(2021)136306)","l");  
  leg1->Draw("same");

  TLatex* latex5 = new TLatex();
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.49,0.82,"ALICE Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"); //  #gamma+Pb#rightarrow J/#psi+Pb,

  c1->SaveAs("gammaPb.eps");
}//main
