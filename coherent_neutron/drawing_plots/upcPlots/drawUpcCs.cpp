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
void SetStyle(TH1D* h, Color_t color, Style_t style, Width_t width = 3)
{
  h->SetLineColor(color);
  h->SetLineStyle(style);
  h->SetLineWidth(width);
}
//_____________________________________________________________________________
void SetStyle(TGraph* g, Color_t color, Style_t style, Width_t width = 3)
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
// read IA model by GKZ
TGraph *getIAgkz()
{
  ifstream f("predictions/gkzIA.d");
  const int n=41;
  double w[n];
  double x[n];
  double cs[n];
  double y[n];

  for (int i=0;i<n;i++){
    f >> w[i] >> x[i] >> cs[i];
    y[i] = -std::log(w[i]*w[i]/(mjpsi*5020));
    // cout << w[i] << ",  " << y[i] << ",  " << cs[i] << endl;
  }
  TGraph* g = new TGraph(n,w,cs);  

  return (TGraph*) g;
}
//_____________________________________________________________________________
// get n00n flux
TGraph *getn00nFlux(int neutronClass)
{

  double n00n_y[33]          = { -4.00, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50, -2.25, -2.00, -1.75, -1.50, -1.25, -1.00, -0.75, -0.50, -0.25, 0.00,  0.25,  0.50,  0.75,  1.00,  1.25,  1.50,  1.75,  2.00,  2.25,  2.50,  2.75,  3.00,  3.25,  3.50,  3.75,  4.00 };
  double n00n_Flux[33]       = { 210.736000,  203.064000,  195.280000,  187.505000,  179.787000,  171.913000,  164.258000,  156.438000,  148.696000,  140.944000,  133.103000,  125.418000,  117.568000,  109.870000,  102.099000,  94.316600, 86.639400, 78.816800, 71.250700, 63.593600, 56.068400, 48.692900, 41.374600, 34.510300, 27.827800, 21.742100, 16.222800, 11.398100, 7.565360,  4.511350,  2.502040,  1.199790,  0.490969 };
  double n00n_Frac_0n0n[33]  = { 0.883772,  0.879343,  0.874520,  0.869308,  0.863675,  0.857466,  0.850805,  0.843400,  0.835302,  0.826347,  0.816351,  0.805350,  0.792914,  0.779125,  0.763551,  0.745999,  0.726345,  0.703957,  0.679173,  0.651144,  0.620075,  0.585867,  0.548189,  0.508406,  0.465739,  0.422223,  0.378088,  0.334363,  0.293320,  0.254114,  0.219903,  0.189250,  0.163063 };
  double n00n_Frac_0nXn[33]  = { 0.0861876, 0.0894742, 0.0930528, 0.0969181, 0.1010920, 0.1056940, 0.1106280, 0.1161080, 0.1220960, 0.1287100, 0.1360820, 0.1441750, 0.1532980, 0.1633710, 0.1746900, 0.1873590, 0.2014130, 0.2172530, 0.2345020, 0.2536680, 0.2743960, 0.2965240, 0.3200280, 0.3435360, 0.3673230, 0.3895610, 0.4098290, 0.4272380, 0.4403070, 0.4496250, 0.4538860, 0.4542070, 0.4508290 };
  double n00n_Frac_XnXn[33]  = { 0.0300406, 0.0311824, 0.0324271, 0.0337744, 0.0352326, 0.0368405, 0.0385679, 0.0404919, 0.0426017, 0.0449426, 0.0475674, 0.0504756, 0.0537880, 0.0575044, 0.0617595, 0.0666425, 0.0722426, 0.0787896, 0.0863253, 0.0951879, 0.1055290, 0.1176100, 0.1317830, 0.1480580, 0.1669380, 0.1882160, 0.2120830, 0.2383990, 0.2663740, 0.2962610, 0.3262110, 0.3565430, 0.3861080 };


  double flx[33];
  for (int i = 0; i < 33; ++i)
  {
    if (neutronClass == 0) flx[i] = n00n_Flux[i]*n00n_Frac_0n0n[i];
    if (neutronClass == 1) flx[i] = n00n_Flux[i]*n00n_Frac_0nXn[i];
    if (neutronClass == 2) flx[i] = n00n_Flux[i]*n00n_Frac_XnXn[i];
  }

  TGraph* g = new TGraph(33,n00n_y,flx);  

  return (TGraph*) g;
}
//_____________________________________________________________________________
// compute W from y
double getWfromY(double y)
{
  double w = sqrt(mjpsi*5020*exp(-y));
  // y[i] = -std::log(pow(w[i],2)/(mjpsi*5020));
  return (double) w;
}
//_____________________________________________________________________________
// get n00n flux
TGraph *getUpcCsIaN00n(int neutronClass)
{
  TGraph *predIAgkz = getIAgkz(); // IA from GKZ
  TGraph *n00nFlux = getn00nFlux(neutronClass);

  double y[23] = { -4.10, -4.00, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50, -2.25, -2.00, -1.75, -1.50, -1.25, -1.00, -0.75, -0.50, -0.25, 0.00,  0.25,  0.50,  0.75,  1.00, 1.20 };
  double cs[23];

  for (int i = 0; i < 23; ++i)
  {
    cs[i] = (n00nFlux->Eval(y[i])) * (predIAgkz->Eval(getWfromY(-y[i])))/1000 + (n00nFlux->Eval(-y[i])) * (predIAgkz->Eval(getWfromY(y[i])))/1000;
    // cout << y[i]  << ",  " << getWfromY(-y[i]) << ",  " << n00nFlux->Eval(y[i]) << ",  " << (predIAgkz->Eval(getWfromY(-y[i])))/1000  << ",  " << (n00nFlux->Eval(y[i])) * (predIAgkz->Eval(getWfromY(-y[i])))/1000 << ",  " << getWfromY(y[i]) << ",  " << n00nFlux->Eval(-y[i]) << ",  " << (predIAgkz->Eval(getWfromY(y[i])))/1000  << ",  " << (n00nFlux->Eval(-y[i])) * (predIAgkz->Eval(getWfromY(y[i])))/1000  << ",  " << cs[i] << endl;
  }

  TGraph* g = new TGraph(23,y,cs);  

  return (TGraph*) g;
}
//_____________________________________________________________________________
void drawUpc(Int_t selectionFlag = 4)
{
  // configure layout of the plot
  //--- ranges
  double gxmin = -4.5, gymin = 0., gxmax = 1.2, gymax = 5.;

  if (selectionFlag == 1) gymax = 8.5;
  if (selectionFlag == 2) gymax = 2.5;
  if (selectionFlag == 3) gymax = 1.3;
  if (selectionFlag == 4) gymax = 1.05;

  //--- style, title...
  TCanvas* c1 = new TCanvas("c1","c1",1000,700);
  TH1F* frame = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
  
  float fSiz = 0.055;
  int iFont = 42;
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleOffset(0.9);
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
  frame->SetTitle(";#it{y};d#sigma/d#it{y} (mb)");
  
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
  // TGraph *predIA_0N0N_SL = new TGraph(); 
  // readImpulseApproximation(predIA_0N0N_SL,"predictions/IA-onon.txt");
  TGraph *predIA_0N0N = getUpcCsIaN00n(0);
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
  // TGraph *predIA_0NXN = new TGraph(); 
  // readImpulseApproximation(predIA_0NXN,"predictions/IA-xnon+onxn.txt");
  TGraph *predIA_0NXN = getUpcCsIaN00n(1);
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
  // TGraph *predIA_XNXN = new TGraph(); 
  // readImpulseApproximation(predIA_XNXN,"predictions/IA-xnxn.txt");
  TGraph *predIA_XNXN = getUpcCsIaN00n(2);
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
  // SetStyle(predIA_0N0N_SL                      ,kAzure+2 ,9);
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
    // predIA_0N0N_SL->Draw("csame");
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
    // predIA_XNXN->Scale(0.001,"y");
    predIA_XNXN->Draw("csame");
    predHS_XNXN->Draw("csame");
    predBKA_XNXN->Draw("csame");
  }

  // Input ALICE data
  int nP;
  if (selectionFlag == 1) nP = 5;
  if (selectionFlag == 2) nP = 2;
  if (selectionFlag == 3) nP = 3;
  if (selectionFlag == 4) nP = 5;
  
  double ALICE_y_0N0N[5]      = { 0.00, -0.50, -2.75, -3.25, -3.75 };
  double ALICE_y_0NXNXN0N[2]  = { 0.00, -0.50                      };
  double ALICE_y_XN0N[3]      = {              -2.75, -3.25, -3.75 };
  double ALICE_y_XNXN[5]      = { 0.00, -0.50, -2.75, -3.25, -3.75 };

  double ySize_0N0N[5]      = { 0.20, 0.30, 0.25, 0.25, 0.25 };
  double ySize_0NXNXN0N[2]  = { 0.20, 0.30                   };
  double ySize_XN0N[3]      = {             0.25, 0.25, 0.25 };
  double ySize_XNXN[5]      = { 0.20, 0.30, 0.25, 0.25, 0.25 };

  double ALICE_sig_0N0N[5]  = { 3.130, 2.900, 2.668, 2.322, 1.590 };
  double ALICE_sta_0N0N[5]  = { 0.090, 0.070, 0.076, 0.033, 0.049 };
  double ALICE_cor_0N0N[5]  = { 0.047, 0.044, 0.011, 0.010, 0.010 };
  double ALICE_unc_0N0N[5]  = { 0.164, 0.152, 0.199, 0.173, 0.119 };
  double ALICE_mig_0N0N[5]  = { 0.122, 0.104, 0.009, 0.005, 0.003 };

  double ALICE_sig_0NXNXN0N[2]  = { 0.730, 0.800 };
  double ALICE_sta_0NXNXN0N[2]  = { 0.050, 0.040 };
  double ALICE_cor_0NXNXN0N[2]  = { 0.015, 0.017 };
  double ALICE_unc_0NXNXN0N[2]  = { 0.045, 0.050 };
  double ALICE_mig_0NXNXN0N[2]  = { 0.025, 0.025 };

  double ALICE_sig_XN0N[3]  = { 0.242, 0.172, 0.101 };
  double ALICE_sta_XN0N[3]  = { 0.021, 0.010, 0.009 };
  double ALICE_cor_XN0N[3]  = { 0.003, 0.002, 0.001 };
  double ALICE_unc_XN0N[3]  = { 0.018, 0.013, 0.008 };
  double ALICE_mig_XN0N[3]  = { 0.009, 0.006, 0.003 };

  double ALICE_sig_XNXN[5]  = { 0.250, 0.300, 0.256, 0.161, 0.079 };
  double ALICE_sta_XNXN[5]  = { 0.024, 0.029, 0.024, 0.011, 0.011 };
  double ALICE_cor_XNXN[5]  = { 0.005, 0.006, 0.005, 0.005, 0.002 };
  double ALICE_unc_XNXN[5]  = { 0.016, 0.019, 0.024, 0.015, 0.008 };
  double ALICE_mig_XNXN[5]  = { 0.002, 0.003, 0.009, 0.006, 0.003 };

  // Mirrored data point
  double ALICE_y_0N0N_mirror[1]      = { 0.50 };
  double ALICE_y_0NXNXN0N_mirror[1]  = { 0.50 };
  double ALICE_y_XNXN_mirror[1]      = { 0.50 };

  double ySize_0N0N_mirror[1]      = { 0.30 };
  double ySize_0NXNXN0N_mirror[1]  = { 0.30 };
  double ySize_XNXN_mirror[1]      = { 0.30 };

  double ALICE_sig_0N0N_mirror[1]  = { 2.900 };
  double ALICE_sta_0N0N_mirror[1]  = { 0.070 };
  double ALICE_cor_0N0N_mirror[1]  = { 0.044 };
  double ALICE_unc_0N0N_mirror[1]  = { 0.152 };
  double ALICE_mig_0N0N_mirror[1]  = { 0.104 };

  double ALICE_sig_0NXNXN0N_mirror[1]  = { 0.800 };
  double ALICE_sta_0NXNXN0N_mirror[1]  = { 0.040 };
  double ALICE_cor_0NXNXN0N_mirror[1]  = { 0.017 };
  double ALICE_unc_0NXNXN0N_mirror[1]  = { 0.050 };
  double ALICE_mig_0NXNXN0N_mirror[1]  = { 0.025 };

  double ALICE_sig_XNXN_mirror[1]  = { 0.300 };
  double ALICE_sta_XNXN_mirror[1]  = { 0.029 };
  double ALICE_cor_XNXN_mirror[1]  = { 0.006 };
  double ALICE_unc_XNXN_mirror[1]  = { 0.019 };
  double ALICE_mig_XNXN_mirror[1]  = { 0.003 };

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
  TGraphErrors *aliceDataErr1_mirror;
  if (selectionFlag == 1) 
  {
    aliceDataErr1  = new TGraphErrors(nP, ALICE_y_0N0N, ALICE_sig_0N0N, ySize_0N0N, &getSquareSum(ALICE_cor_0N0N,ALICE_mig_0N0N,nP)[0]);
    aliceDataErr1_mirror  = new TGraphErrors(1, ALICE_y_0N0N_mirror, ALICE_sig_0N0N_mirror, ySize_0N0N_mirror, &getSquareSum(ALICE_cor_0N0N_mirror,ALICE_mig_0N0N_mirror,1)[0]);
  }
  if (selectionFlag == 2) 
  {
    aliceDataErr1  = new TGraphErrors(nP, ALICE_y_0NXNXN0N, ALICE_sig_0NXNXN0N, ySize_0NXNXN0N, &getSquareSum(ALICE_cor_0NXNXN0N,ALICE_mig_0NXNXN0N,nP)[0]);
    aliceDataErr1_mirror  = new TGraphErrors(1, ALICE_y_0NXNXN0N_mirror, ALICE_sig_0NXNXN0N_mirror, ySize_0NXNXN0N_mirror, &getSquareSum(ALICE_cor_0NXNXN0N_mirror,ALICE_mig_0NXNXN0N_mirror,1)[0]);
  }
  if (selectionFlag == 3) aliceDataErr1  = new TGraphErrors(nP, ALICE_y_XN0N, ALICE_sig_XN0N, ySize_XN0N, &getSquareSum(ALICE_cor_XN0N,ALICE_mig_XN0N,nP)[0]);
  if (selectionFlag == 4) 
  {
    aliceDataErr1  = new TGraphErrors(nP, ALICE_y_XNXN, ALICE_sig_XNXN, ySize_XNXN, &getSquareSum(ALICE_cor_XNXN,ALICE_mig_XNXN,nP)[0]);
    aliceDataErr1_mirror  = new TGraphErrors(1, ALICE_y_XNXN_mirror, ALICE_sig_XNXN_mirror, ySize_XNXN_mirror, &getSquareSum(ALICE_cor_XNXN_mirror,ALICE_mig_XNXN_mirror,1)[0]);
  }
  aliceDataErr1->SetFillColor(cData);
  aliceDataErr1->SetFillStyle(0);
  aliceDataErr1->SetLineWidth(2);
  aliceDataErr1->Draw("2");
  if (selectionFlag != 3)
  {
    aliceDataErr1_mirror->SetFillColor(cData);
    aliceDataErr1_mirror->SetFillStyle(0);
    aliceDataErr1_mirror->SetLineWidth(2);
    aliceDataErr1_mirror->Draw("2");
  }


  TGraphErrors *aliceData;
  TGraphErrors *aliceData_mirror;
  if (selectionFlag == 1) 
  {
    aliceData = new TGraphErrors(nP, ALICE_y_0N0N, ALICE_sig_0N0N, NULL, &getSquareSum(ALICE_sta_0N0N,ALICE_unc_0N0N,nP)[0]);
    aliceData_mirror = new TGraphErrors(1, ALICE_y_0N0N_mirror, ALICE_sig_0N0N_mirror, NULL, &getSquareSum(ALICE_sta_0N0N_mirror,ALICE_unc_0N0N_mirror,1)[0]);
  }
  if (selectionFlag == 2) 
  {
    aliceData = new TGraphErrors(nP, ALICE_y_0NXNXN0N, ALICE_sig_0NXNXN0N, NULL, &getSquareSum(ALICE_sta_0NXNXN0N,ALICE_unc_0NXNXN0N,nP)[0]);
    aliceData_mirror = new TGraphErrors(1, ALICE_y_0NXNXN0N_mirror, ALICE_sig_0NXNXN0N_mirror, NULL, &getSquareSum(ALICE_sta_0NXNXN0N_mirror,ALICE_unc_0NXNXN0N_mirror,1)[0]);
  }
  if (selectionFlag == 3) aliceData = new TGraphErrors(nP, ALICE_y_XN0N, ALICE_sig_XN0N, NULL, &getSquareSum(ALICE_sta_XN0N,ALICE_unc_XN0N,nP)[0]);
  if (selectionFlag == 4) 
  {
    aliceData = new TGraphErrors(nP, ALICE_y_XNXN, ALICE_sig_XNXN, NULL, &getSquareSum(ALICE_sta_XNXN,ALICE_unc_XNXN,nP)[0]);
    aliceData_mirror = new TGraphErrors(1, ALICE_y_XNXN_mirror, ALICE_sig_XNXN_mirror, NULL, &getSquareSum(ALICE_sta_XNXN_mirror,ALICE_unc_XNXN_mirror,1)[0]);
  }
  aliceData->SetMarkerColor(cData);
  aliceData->SetMarkerStyle(kFullCircle);
  aliceData->SetMarkerSize(1.2);
  aliceData->SetLineColor(cData);
  aliceData->SetLineWidth(3);
  aliceData->Draw("zpsame");
  if (selectionFlag != 3)
  {
    aliceData_mirror->SetMarkerColor(cData);
    aliceData_mirror->SetMarkerStyle(kOpenCircle);
    aliceData_mirror->SetMarkerSize(1.2);
    aliceData_mirror->SetLineColor(cData);
    aliceData_mirror->SetLineWidth(3);
    aliceData_mirror->Draw("zpsame");
  }

  //legend for experimental data
  Double_t xl = 0.13, dxl = 0.2;
  Double_t yl = 0.60, dyl = 0.25;
  if ( selectionFlag != 0 ) yl = 0.54, dyl = 0.33;
  TLegend* leg = new TLegend(xl, yl, xl+dxl, yl+dyl);
  leg->SetFillStyle(0);
  leg->SetTextSize(fSiz-0.005);
  if ( selectionFlag == 0 || selectionFlag == 1 ) {
    leg->AddEntry(aliceData,"ALICE 0n0n", "EP");
    leg->AddEntry(predIA_0N0N,"Impulse approximation", "l");
    // leg->AddEntry(predIA_0N0N_SL,"Impulse approximation SL", "l");
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
  // latex5->SetTextSize(fSiz-0.005);
  latex5->DrawLatex(0.15,0.90,"ALICE Pb#kern[0.2]{#minus}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");

  c1->SaveAs(Form("upcCs_%d.eps", selectionFlag));
  c1->SaveAs(Form("upcCs_%d.pdf", selectionFlag));
}


void drawUpcCs()
{
  for (int i = 1; i < 5; ++i)
  {
    drawUpc(i);
  }
}