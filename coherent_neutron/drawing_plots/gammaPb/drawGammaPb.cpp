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
    y[i] = -std::log(pow(w[i],2)/(mjpsi*5020));
    // cout << w[i] << ",  " << y[i] << ",  " << cs[i] << endl;
  }
  TGraph* g = new TGraph(n,w,cs);  

  return (TGraph*) g;
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
  double gxmin = 12., gymin = 3, gxmax = 1200., gymax = 1300;

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
  frame->SetTitle(";#it{W}_{#gammaPb,n} (GeV);#sigma(#gammaPb) (#mub)");
  
  gPad->SetMargin(0.10,0.05,0.12,0.10);
  gPad->SetLogy();
  gPad->SetTicky(1);
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
  // TGraph *predIAsl = new TGraph(); // prediction from the IA model from STARlight
  // readImpulseApproximation(predIAsl);
  TGraph *predIA = getIAgkz(); // prediction from the IA model by GKZ
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


  // ALICE data
  int nP = 9;
  double ALICE_W[9]    = { 19.12, 813.05, 24.55, 633.21, 31.53, 493.14, 97.11, 160.1, 124.69  };
  double ALICE_sig[9]  = { 8.84, 57.32, 13.89, 46.58, 16.89, 44.68, 21.73, 25.00, 24.15       };
  double ALICE_unc[9]  = { 0.30, 20.77, 0.23,  6.61,  0.59,  6.38,  5.12,  7.33,  0.69        };
  double ALICE_cor[9]  = { 0.68, 7.57,  1.08,  5.73,  1.32,  5.15,  3.12,  4.88,  1.37        };
  double ALICE_mig[9]  = { 0.02, 6.41,  0.05,  3.77,  0.11,  2.73,  4.32,  5.43,  0.50        };
  double ALICE_flu[9]  = { 0.04, 6.56,  0.08,  3.63,  0.18,  2.97,  2.73,  3.91,  0.06        };

  // CMS data
  double CMS_W[6]      = { 40.5, 45.9, 52.0, 299.0, 339.0, 384.0 };
  double CMS_sig[6]    = { 20.22, 22.44, 22.07, 33.90, 29.50, 32.80 };
  double CMS_sta[6]    = { 0.38, 0.43, 0.77, 2.9, 1.8, 2.1 };
  double CMS_teo[6]    = { 0.59, 0.76, 1.14, 5.5, 4.8, 5.3 };
  double CMS_sys[6]    = { 1.2, 1.39, 1.86, 3.2, 3.4, 3.1 };

  // Data from Guillermo paper Phys.Rev.C 96 (2017) 1, 015203: https://inspirehep.net/literature/1491190
  double guillermo_W[3]      = {    18.2,    92.4,    469.5  };
  double guillermo_sig[3]    = {     5.2,    17.9,     38.1  };
  double guillermo_sta_h[3]  = {     1.0,     0.0,     15.0  };
  double guillermo_sta_l[3]  = {     1.0,     0.0,     15.0  };
  double guillermo_sys_h[3]  = {     1.0,     2.6,      9.9  };
  double guillermo_sys_l[3]  = {     1.0,     1.8,     11.3  };
  double guillermo_stasys_h[3];
  double guillermo_stasys_l[3];
    for(int i = 0; i < 3; i++){
    guillermo_stasys_h[i] = TMath::Sqrt(guillermo_sta_h[i]*guillermo_sta_h[i]+guillermo_sys_h[i]*guillermo_sys_h[i]);
    guillermo_stasys_l[i] = TMath::Sqrt(guillermo_sta_l[i]*guillermo_sta_l[i]+guillermo_sys_l[i]*guillermo_sys_l[i]);
  }

  // Data from Evgeny paper  Phys.Lett.B 726 (2013) 290-295: https://inspirehep.net/literature/1232397
  double evgeny_W[2]         = {   19.6,    92.4   };
  double evgeny_sig[2]       = {    6.1,    17.6   };
  double evgeny_stasys_h[2]  = {    1.8,     2.7   };
  double evgeny_stasys_l[2]  = {    2.0,     2.0   };

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

  //--- plot ALICE data
  double ySize[nP];
  double CMSySize[6];
  for (int i = 0; i < nP; i++) ySize[i] = 0.07*ALICE_W[i];
  for (int i = 0; i < 6; i++) CMSySize[i] = 0.05*CMS_W[i];

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

  //--- plot CMS data
  TGraphErrors *CMSDataErr2  = new TGraphErrors(6, CMS_W, CMS_sig, CMSySize, CMS_sys);
  CMSDataErr2->SetFillColor(cError);
  // CMSDataErr2->Draw("2");

  TGraphErrors *CMSDataErr1  = new TGraphErrors(6, CMS_W, CMS_sig, CMSySize, CMS_teo);
  CMSDataErr1->SetFillColor(cData);
  CMSDataErr1->SetFillStyle(0);
  CMSDataErr1->SetLineWidth(2);
  // CMSDataErr1->Draw("2");

  TGraphErrors *cmsData = new TGraphErrors(6, CMS_W, CMS_sig, NULL, CMS_sta);
  cmsData->SetMarkerColor(kRed+1);
  cmsData->SetMarkerStyle(kOpenCircle);
  cmsData->SetMarkerSize(1.0);
  cmsData->SetLineColor(kRed+1);
  cmsData->SetLineWidth(2);
  // cmsData->Draw("zpsame");


  //--- plot data from Guillermo paper Phys.Rev.C 96 (2017) 1, 015203: https://inspirehep.net/literature/1491190
  TGraphAsymmErrors* guillermoData = new TGraphAsymmErrors(3,guillermo_W, guillermo_sig,NULL,NULL,guillermo_stasys_l,guillermo_stasys_h);
  guillermoData->SetMarkerColor(kBlue-1);
  guillermoData->SetMarkerStyle(kOpenTriangleUp);
  guillermoData->SetMarkerSize(1.2);
  guillermoData->SetLineColor(kBlue-1);
  guillermoData->SetLineWidth(2);
  guillermoData->Draw("zpsame");

  //--- plot data from Evgeny paper  Phys.Lett.B 726 (2013) 290-295: https://inspirehep.net/literature/1232397
  TGraphAsymmErrors* evgenyData = new TGraphAsymmErrors(2,evgeny_W, evgeny_sig,NULL,NULL,evgeny_stasys_l,evgeny_stasys_h);
  evgenyData->SetMarkerColor(kBlue-3);
  evgenyData->SetMarkerStyle(kOpenSquare);
  evgenyData->SetMarkerSize(1.2);
  evgenyData->SetLineColor(kBlue-3);
  evgenyData->SetLineWidth(2);
  evgenyData->Draw("zpsame");

  //legend 
  double xl = 0.14, dxl = 0.2;
  double yl = 0.58, dyl = 0.28;
  // yl = 0.44, dyl = 0.42; // without CMS
  // yl = 0.652596, dyl = 0.2074; //without models
  yl = 0.393333, dyl = 0.4666666;
  TLegend* leg1 = new TLegend(xl, yl, xl+dxl, yl+dyl);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(fSiz-0.010);
  leg1->AddEntry(aliceData,"ALICE, #scale[0.9]{Pb#kern[0.15]{#minus}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}","EP");
  // leg1->AddEntry(cmsData,"CMS, #scale[0.9]{Pb#kern[0.15]{#minus}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV} #scale[0.8]{(arXiv:2303.16984})","EP");
  leg1->AddEntry(evgenyData,"Guzey et al., #scale[0.9]{using ALICE Pb#kern[0.15]{#minus}Pb #sqrt{#it{s}_{NN}} = 2.76 TeV} #scale[0.8]{(PLB 726 (2013) 290-295)} ","EP");
  leg1->AddEntry(guillermoData,"Contreras, #scale[0.9]{using ALICE Pb#kern[0.15]{#minus}Pb #sqrt{#it{s}_{NN}} = 2.76 TeV} #scale[0.8]{(PRC 96 (2017) 015203)}","EP");
  leg1->AddEntry(predIA,"Impulse approximation","l");
  leg1->AddEntry(predSL,"STARlight","l");
  leg1->AddEntry(epsC,"EPS09 LO","l");
  leg1->AddEntry(ltaAverage,"LTA","l");
  leg1->AddEntry(predHS,"GG-HS","l"); //(PRC97 (2018) 024901)
  leg1->AddEntry(predBKA,"b-BK-A","l"); // (PLB817(2021)136306)
  // leg1->AddEntry(predBKGG,"b-BK-GG (PLB817(2021)136306)","l");  
  leg1->Draw("same");

  // TLatex* latex5 = new TLatex();
  // latex5->SetTextAlign(11);
  // latex5->SetNDC();
  // latex5->DrawLatex(0.49,0.82,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"); //  #gamma+Pb#rightarrow J/#psi+Pb,

  c1->SaveAs("gammaPb.eps");
  c1->SaveAs("gammaPb.pdf");
}//main
