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
// read model from STARlight Impulse Approximation
void readImpulseApproximation(TGraph *gr, TString path)
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
  ifstream ifs(path.Data());
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
    if (i % 2 != 0 ) continue; // Skip half of the poins because the data is too dense for TGraph to draw smoothly
    gr->SetPoint(gr->GetN(),y,boh);
  }
  ifs.close();
}
//_____________________________________________________________________________
// read model STARlight 
void readSTARlight(TGraph *gr, TString path)
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
  ifstream ifs(path.Data());
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
    if (i % 2 != 0 ) continue; // Skip half of the poins because the data is too dense for TGraph to draw smoothly
    gr->SetPoint(gr->GetN(),y,boh);
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
TGraph *removeLastValue(TGraph *inputA)
{
  TGraph *res = new TGraph();
  int counter = inputA->GetMaxSize() - 1;
  for (int i = 0; i < counter; i++) 
  {
    res->SetPoint(res->GetN(),inputA->GetPointX(i),inputA->GetPointY(i));
  }
  return (TGraph*) res;
}
//_____________________________________________________________________________
void drawUpcCs(Int_t selectionFlag = 1)
{
  // configure layout of the plot
  //--- ranges
  double gxmin = -4.5, gymin = 0., gxmax = 1.2, gymax = 5.;

  if (selectionFlag == 1) gymax = 11.;
  if (selectionFlag == 2) gymax = 3.5;
  if (selectionFlag == 3) gymax = 1.6;
  if (selectionFlag == 4) gymax = 1.25;

  //--- style, title...
  TCanvas* c1 = new TCanvas("c1","c1",1000,700);
  TH1F* frame = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
  
  float fSiz = 0.045;
  int iFont = 42;
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetXaxis()->SetTitleSize(fSiz);
  frame->GetYaxis()->SetTitleSize(fSiz);
  frame->GetXaxis()->SetLabelSize(fSiz);
  frame->GetYaxis()->SetLabelSize(fSiz);
  frame->GetXaxis()->SetTitleFont(iFont);
  frame->GetYaxis()->SetTitleFont(iFont);
  frame->GetXaxis()->SetLabelFont(iFont);
  frame->GetYaxis()->SetLabelFont(iFont);
  frame->GetYaxis()->SetDecimals(1); 
  frame->SetTitleSize(fSiz);       
  frame->SetLabelSize(fSiz);
  frame->SetTitleSize(fSiz, "Y");  
  frame->SetLabelSize(fSiz, "Y");
  frame->SetTitle(";y;d#sigma/dy (mb)");
  
  gPad->SetMargin(0.10,0.03,0.12,0.03);

  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(2.3);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetTextSize(fSiz);
  gStyle->SetTextFont(iFont);
  gStyle->SetJoinLinePS(1); // Set the PDF line join 

  TGaxis::SetMaxDigits(3);


  // Theory prediction 
  //--- read
  TFile *fileGKZ = new TFile("predictions/predictions_with_neutrons.root");
  TFile *filePrg = new TFile("predictions/NeutronPredictions.root");

  //--- 0N0N
  // EPS 
  TGraph *eps_0N0N_c = (TGraph*) fileGKZ->Get("geps00_0");
  // LTA 
  TGraph *lta_0N0N_h = (TGraph*) fileGKZ->Get("glta00_h");
  TGraph *lta_0N0N_l = (TGraph*) fileGKZ->Get("glta00_l");
  TGraph *ltaAverage_0N0N = getAvergaGraph(lta_0N0N_h, lta_0N0N_l);
  // IA 
  TGraph *predIA_0N0N = new TGraph(); 
  readImpulseApproximation(predIA_0N0N,"predictions/IA-onon.txt");
  // STARlight
  TGraph *predSL_0N0N = new TGraph(); 
  readSTARlight(predSL_0N0N,"predictions/SL-0N0N.txt");
  // Hot Spot 
  TGraph *predHS_0N0N = (TGraph*) filePrg->Get("gGGHS_0n0n");
  // b-BK-A 
  TGraph *predBKA_0N0N = (TGraph*) filePrg->Get("gBBK_0n0n");

  //--- 0NXN
  // EPS 
  TGraph *eps_0NXN_source = (TGraph*) fileGKZ->Get("geps0X_0");
  TGraph *eps_0NXN_c = removeLastValue(eps_0NXN_source);
  // LTA 
  TGraph *lta_0NXN_h = (TGraph*) fileGKZ->Get("glta0X_h");
  TGraph *lta_0NXN_l = (TGraph*) fileGKZ->Get("glta0X_l");
  TGraph *ltaAverage_0NXN = getAvergaGraph(lta_0NXN_h, lta_0NXN_l);
  // IA 
  TGraph *predIA_0NXN = new TGraph(); 
  readImpulseApproximation(predIA_0NXN,"predictions/IA-xnon+onxn.txt");
  // STARlight
  TGraph *predSL_0NXN = new TGraph(); 
  readSTARlight(predSL_0NXN,"predictions/SL-XN0N+0NXN.txt");
  // Hot Spot 
  TGraph *predHS_0NXN = (TGraph*) filePrg->Get("gGGHS_Xn0n");
  // b-BK-A 
  TGraph *predBKA_0NXN = (TGraph*) filePrg->Get("gBBK_Xn0n");

  //--- XNXN
  // EPS 
  TGraph *eps_XNXN_source = (TGraph*) fileGKZ->Get("gepsXX_0");
  TGraph *eps_XNXN_c = removeLastValue(eps_XNXN_source);
  // LTA 
  TGraph *lta_XNXN_h = (TGraph*) fileGKZ->Get("gltaXX_h");
  TGraph *lta_XNXN_l = (TGraph*) fileGKZ->Get("gltaXX_l");
  TGraph *ltaAverage_XNXN = getAvergaGraph(lta_XNXN_h, lta_XNXN_l);
  // IA 
  TGraph *predIA_XNXN = new TGraph(); 
  readImpulseApproximation(predIA_XNXN,"predictions/IA-xnxn.txt");
  // STARlight
  TGraph *predSL_XNXN = new TGraph(); 
  readSTARlight(predSL_XNXN,"predictions/SL-XNXN.txt");
  // Hot Spot 
  TGraph *predHS_XNXN = (TGraph*) filePrg->Get("gGGHS_XnXn");
  // b-BK-A 
  TGraph *predBKA_XNXN = (TGraph*) filePrg->Get("gBBK_XnXn");

  //--- style and plot
  Int_t c0N0N = TColor::GetFreeColorIndex();
  TColor *color0N0N = new TColor(c0N0N, 0.8, 0.8784, 0.8);

  SetStyle(eps_0N0N_c                       ,kGreen+2 ,1);
  SetStyle(ltaAverage_0N0N                  ,kTeal+2  ,4);
  SetStyle(predIA_0N0N                      ,kAzure+1 ,7);
  SetStyle(predSL_0N0N                      ,kBlue-2  ,5);
  SetStyle(predHS_0N0N                      ,kRed-4   ,2);
  SetStyle(predBKA_0N0N                     ,kRed+2   ,8);

  SetStyle(eps_0NXN_c                       ,kGreen+2 ,1);
  SetStyle(ltaAverage_0NXN                  ,kTeal+2  ,4);
  SetStyle(predIA_0NXN                      ,kAzure+1 ,7);
  SetStyle(predSL_0NXN                      ,kBlue-2  ,5);
  SetStyle(predHS_0NXN                      ,kRed-4   ,2);
  SetStyle(predBKA_0NXN                     ,kRed+2   ,8);

  SetStyle(eps_XNXN_c                       ,kGreen+2 ,1);
  SetStyle(ltaAverage_XNXN                  ,kTeal+2  ,4);
  SetStyle(predIA_XNXN                      ,kAzure+1 ,7);
  SetStyle(predSL_XNXN                      ,kBlue-2  ,5);
  SetStyle(predHS_XNXN                      ,kRed-4   ,2);
  SetStyle(predBKA_XNXN                     ,kRed+2   ,8);

  if (selectionFlag == 1) {
    eps_0N0N_c->Draw("csame");
    ltaAverage_0N0N->Draw("csame");
    predSL_0N0N->Draw("csame");
    predIA_0N0N->Draw("csame");
    predHS_0N0N->Draw("csame");
    predBKA_0N0N->Draw("csame");
  }
  if (selectionFlag == 2) {  
    eps_0NXN_c->Draw("csame");
    ltaAverage_0NXN->Draw("csame");
    predSL_0NXN->Draw("csame");
    predIA_0NXN->Draw("csame");
    predHS_0NXN->Draw("csame");
    predBKA_0NXN->Draw("csame");
  }
  if (selectionFlag == 3) {
    eps_0NXN_c->Scale(0.5,"y");
    eps_0NXN_c->Draw("csame");
    ltaAverage_0NXN->Scale(0.5,"y"); 
    ltaAverage_0NXN->Draw("csame");
    predSL_0NXN->Scale(0.5,"y"); 
    predSL_0NXN->Draw("csame");
    predIA_0NXN->Scale(0.5,"y");
    predIA_0NXN->Draw("csame");
    predHS_0NXN->Scale(0.5,"y"); 
    predHS_0NXN->Draw("csame");
    predBKA_0NXN->Scale(0.5,"y"); 
    predBKA_0NXN->Draw("csame");
  }
  if (selectionFlag == 4) {
    eps_XNXN_c->Draw("csame");
    ltaAverage_XNXN->Draw("csame");
    predSL_XNXN->Scale(0.001,"y");
    predSL_XNXN->Draw("csame");
    predIA_XNXN->Scale(0.001,"y");
    predIA_XNXN->Draw("csame");
    predHS_XNXN->Draw("csame");
    predBKA_XNXN->Draw("csame");
  }

  // Input ALICE data
  int nP;
  if (selectionFlag == 1) nP = 6;
  if (selectionFlag == 2) nP = 3;
  if (selectionFlag == 3) nP = 3;
  if (selectionFlag == 4) nP = 6;
  
  double ALICE_y_0N0N[6]      = { 0.50, 0.00, -0.50, -2.75, -3.25, -3.75 };
  double ALICE_y_0NXNXN0N[3]  = { 0.50, 0.00, -0.50                      };
  double ALICE_y_XN0N[3]      = {                    -2.75, -3.25, -3.75 };
  double ALICE_y_XNXN[6]      = { 0.50, 0.00, -0.50, -2.75, -3.25, -3.75 };

  double ySize_0N0N[6]      = { 0.30, 0.20, 0.30, 0.25, 0.25, 0.25 };
  double ySize_0NXNXN0N[3]  = { 0.30, 0.20, 0.30                   };
  double ySize_XN0N[3]      = {                   0.25, 0.25, 0.25 };
  double ySize_XNXN[6]      = { 0.30, 0.20, 0.30, 0.25, 0.25, 0.25 };

  double ALICE_sig_0N0N[6]  = { 2.880, 3.110, 2.880, 2.657, 2.316, 1.586 };
  double ALICE_sta_0N0N[6]  = { 0.070, 0.090, 0.070, 0.076, 0.033, 0.049 };
  double ALICE_cor_0N0N[6]  = { 0.043, 0.047, 0.043, 0.011, 0.010, 0.010 };
  double ALICE_unc_0N0N[6]  = { 0.151, 0.163, 0.151, 0.198, 0.173, 0.118 };
  double ALICE_mig_0N0N[6]  = { 0.058, 0.062, 0.058, 0.027, 0.023, 0.016 };

  double ALICE_sig_0NXNXN0N[3]  = { 0.810, 0.730, 0.810 };
  double ALICE_sta_0NXNXN0N[3]  = { 0.040, 0.050, 0.040 };
  double ALICE_cor_0NXNXN0N[3]  = { 0.017, 0.015, 0.017 };
  double ALICE_unc_0NXNXN0N[3]  = { 0.050, 0.045, 0.050 };
  double ALICE_mig_0NXNXN0N[3]  = { 0.041, 0.036, 0.041 };

  double ALICE_sig_XN0N[3]  = { 0.236, 0.168, 0.099 };
  double ALICE_sta_XN0N[3]  = { 0.021, 0.010, 0.008 };
  double ALICE_cor_XN0N[3]  = { 0.003, 0.002, 0.001 };
  double ALICE_unc_XN0N[3]  = { 0.018, 0.013, 0.007 };
  double ALICE_mig_XN0N[3]  = { 0.006, 0.004, 0.003 };

  double ALICE_sig_XNXN[6]  = { 0.300, 0.250, 0.300, 0.268, 0.169, 0.083 };
  double ALICE_sta_XNXN[6]  = { 0.030, 0.030, 0.030, 0.025, 0.012, 0.012 };
  double ALICE_cor_XNXN[6]  = { 0.006, 0.005, 0.006, 0.005, 0.006, 0.002 };
  double ALICE_unc_XNXN[6]  = { 0.019, 0.016, 0.019, 0.026, 0.016, 0.008 };
  double ALICE_mig_XNXN[6]  = { 0.030, 0.025, 0.030, 0.027, 0.017, 0.008 };

  // ALICE data plotting
  vector<double> rgbRaw = {0,0,0};
  vector<double> rgb = {0,0,0};
  //--- define colors
  int cData = TColor::GetFreeColorIndex();
  // rgbRaw = {165, 0, 0};
  // for (int i = 0; i < 3; i++) rgb[i] = rgbRaw[i]/255;
  // TColor *colorData = new TColor(cData, rgb[0], rgb[1], rgb[2]);
  cData = kBlack;

  //--- plot data
  TGraphErrors *aliceDataErr1;
  if (selectionFlag == 1) aliceDataErr1  = new TGraphErrors(nP, ALICE_y_0N0N, ALICE_sig_0N0N, ySize_0N0N, &getSquareSum(ALICE_cor_0N0N,ALICE_mig_0N0N,nP)[0]);
  if (selectionFlag == 2) aliceDataErr1  = new TGraphErrors(nP, ALICE_y_0NXNXN0N, ALICE_sig_0NXNXN0N, ySize_0NXNXN0N, &getSquareSum(ALICE_cor_0NXNXN0N,ALICE_mig_0NXNXN0N,nP)[0]);
  if (selectionFlag == 3) aliceDataErr1  = new TGraphErrors(nP, ALICE_y_XN0N, ALICE_sig_XN0N, ySize_XN0N, &getSquareSum(ALICE_cor_XN0N,ALICE_mig_XN0N,nP)[0]);
  if (selectionFlag == 4) aliceDataErr1  = new TGraphErrors(nP, ALICE_y_XNXN, ALICE_sig_XNXN, ySize_XNXN, &getSquareSum(ALICE_cor_XNXN,ALICE_mig_XNXN,nP)[0]);
  aliceDataErr1->SetFillColor(cData);
  aliceDataErr1->SetFillStyle(0);
  aliceDataErr1->SetLineWidth(2);
  aliceDataErr1->Draw("2");

  TGraphErrors *aliceData;
  if (selectionFlag == 1) aliceData = new TGraphErrors(nP, ALICE_y_0N0N, ALICE_sig_0N0N, NULL, &getSquareSum(ALICE_sta_0N0N,ALICE_unc_0N0N,nP)[0]);
  if (selectionFlag == 2) aliceData = new TGraphErrors(nP, ALICE_y_0NXNXN0N, ALICE_sig_0NXNXN0N, NULL, &getSquareSum(ALICE_sta_0NXNXN0N,ALICE_unc_0NXNXN0N,nP)[0]);
  if (selectionFlag == 3) aliceData = new TGraphErrors(nP, ALICE_y_XN0N, ALICE_sig_XN0N, NULL, &getSquareSum(ALICE_sta_XN0N,ALICE_unc_XN0N,nP)[0]);
  if (selectionFlag == 4) aliceData = new TGraphErrors(nP, ALICE_y_XNXN, ALICE_sig_XNXN, NULL, &getSquareSum(ALICE_sta_XNXN,ALICE_unc_XNXN,nP)[0]);
  aliceData->SetMarkerColor(cData);
  aliceData->SetMarkerStyle(kFullCircle);
  aliceData->SetMarkerSize(1.2);
  aliceData->SetLineColor(cData);
  aliceData->SetLineWidth(2);
  aliceData->Draw("zpsame");


  //legend for experimental data
  Double_t xl = 0.15, dxl = 0.2;
  Double_t yl = 0.60, dyl = 0.25;
  if ( selectionFlag != 0 ) yl = 0.62, dyl = 0.25;
  TLegend* leg = new TLegend(xl, yl, xl+dxl, yl+dyl);
  leg->SetFillStyle(0);
  leg->SetTextSize(fSiz-0.005);
  if ( selectionFlag == 0 || selectionFlag == 1 ) {
    leg->AddEntry(aliceData,"ALICE 0n0n", "EP");
    leg->AddEntry(predIA_0N0N,"Impulse approximation", "l");
    leg->AddEntry(predSL_0N0N,"STARlight", "l");
    leg->AddEntry(eps_0N0N_c,"EPS09 LO", "l");
    leg->AddEntry(ltaAverage_0N0N,"LTA", "l");
    leg->AddEntry(predHS_0N0N,"GG-HS", "l");
    leg->AddEntry(predBKA_0N0N,"b-BK-A", "l"); 
  }
  if ( selectionFlag == 0 || selectionFlag == 2 ) {
    leg->AddEntry(aliceData,"ALICE 0nXn+Xn0n", "EP");
    leg->AddEntry(predIA_0NXN,"Impulse approximation", "l");
    leg->AddEntry(predSL_0NXN,"STARlight", "l");
    leg->AddEntry(eps_0NXN_c,"EPS09 LO", "l");
    leg->AddEntry(ltaAverage_0NXN,"LTA", "l");
    leg->AddEntry(predHS_0NXN,"GG-HS", "l");
    leg->AddEntry(predBKA_0NXN,"b-BK-A", "l");
  }
  if ( selectionFlag == 0 || selectionFlag == 3 ) {
    leg->AddEntry(aliceData,"ALICE Xn0n", "EP");
    leg->AddEntry(predIA_0NXN,"Impulse approximation", "l");
    leg->AddEntry(predSL_0NXN,"STARlight", "l");
    leg->AddEntry(eps_0NXN_c,"EPS09 LO", "l");
    leg->AddEntry(ltaAverage_0NXN,"LTA", "l");
    leg->AddEntry(predHS_0NXN,"GG-HS", "l");
    leg->AddEntry(predBKA_0NXN,"b-BK-A", "l");
  }
  if ( selectionFlag == 0 || selectionFlag == 4 ) {
    leg->AddEntry(aliceData,"ALICE XnXn", "EP");
    leg->AddEntry(predIA_XNXN,"Impulse approximation", "l");
    leg->AddEntry(predSL_XNXN,"STARlight", "l");
    leg->AddEntry(eps_XNXN_c,"EPS09 LO", "l");
    leg->AddEntry(ltaAverage_XNXN,"LTA", "l");
    leg->AddEntry(predHS_XNXN,"GG-HS", "l");
    leg->AddEntry(predBKA_XNXN,"b-BK-A", "l");    
  }
  leg->Draw();

  TLatex* latex5 = new TLatex();
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.17,0.90,"ALICE Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");

  c1->SaveAs(Form("upcCs_%d.eps", selectionFlag));
}
