// #include "Riostream.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TVectorD.h"


// https://root-forum.cern.ch/t/tdecompsvd-not-working/22313
//_____________________________________
Double_t FittingPhotonuclear(Int_t element, Int_t mode = 0, Int_t rap = 0)
{
  //==========================
  // Solving linear equations
  // with SVD technique
  //==========================
  //
  // We have an equation like
  // b = Ax
  // where b are the measured
  // cross sections,
  // A is the flux matrix
  // x is our unknown vector
  // of photonuclear sigmas
  // for positive and
  // negative rapidity


  // 4 neutron emission classes
  // 2 unknown photonuclear xsecs
  TMatrixD Amatrix(4,2);
  TMatrixD Amatrix0(4,2);
  TMatrixD Amatrix1(4,2);
  TMatrixD Amatrix2(4,2);
  Amatrix.Zero();
  Amatrix0.Zero();
  Amatrix1.Zero();
  Amatrix2.Zero();


  // Fluxes from Michal
  /*
  y = 0 0n0n:  Low = 61.537 High = 61.537
  y = 0 0nXn:  Low = 16.1841 High = 16.1841
  y = 0 XnXn:  Low = 5.03313 High = 5.03313
  y = 0.25 0n0n:  Low = 68.7589 High = 54.2806
  y = 0.25 0nXn:  Low = 16.3828 High = 15.8818
  y = 0.25 XnXn:  Low = 5.05204 High = 4.99478
  y = 0.575 0n0n:  Low = 78.3588 High = 45.2674
  y = 0.575 0nXn:  Low = 16.5604 High = 15.3431
  y = 0.575 XnXn:  Low = 5.06777 High = 4.92338
  y = 2.75 0n0n:  Low = 142.948 High = 3.66191
  y = 2.75 0nXn:  Low = 16.7132 High = 4.28895
  y = 2.75 XnXn:  Low = 5.04786 High = 2.04812
  y = 3.25 0n0n:  Low = 157.789 High = 1.07192
  y = 3.25 0nXn:  Low = 16.6849 High = 1.70239
  y = 3.25 XnXn:  Low = 5.03916 High = 0.949336
  y = 3.75 0n0n:  Low = 172.551 High = 0.200768
  y = 3.75 0nXn:  Low = 16.6506 High = 0.421108
  y = 3.75 XnXn:  Low = 5.02952 High = 0.275626
  */


  // -4 < y < -2.5
  Amatrix(0,0) = 157.789;
  Amatrix(0,1) = 1.07192;
  Amatrix(1,0) = 16.6849*0.5;
  Amatrix(1,1) = 1.70239*0.5;
  Amatrix(2,0) = 16.6849*0.5;
  Amatrix(2,1) = 1.70239*0.5;
  Amatrix(3,0) = 5.03916;
  Amatrix(3,1) = 0.949336;
  // -4 < y < -3.5
  Amatrix0(0,0) = 172.551;
  Amatrix0(0,1) = 0.200768;
  Amatrix0(1,0) = 16.6506*0.5;
  Amatrix0(1,1) = 0.421108*0.5;
  Amatrix0(2,0) = 16.6506*0.5;
  Amatrix0(2,1) = 0.421108*0.5;
  Amatrix0(3,0) = 5.02952;
  Amatrix0(3,1) = 0.275626;
  // -3.5 < y < -3.0
  Amatrix1(0,0) = 157.789;
  Amatrix1(0,1) = 1.07192;
  Amatrix1(1,0) = 16.6849*0.5;
  Amatrix1(1,1) = 1.70239*0.5;
  Amatrix1(2,0) = 16.6849*0.5;
  Amatrix1(2,1) = 1.70239*0.5;
  Amatrix1(3,0) = 5.03916;
  Amatrix1(3,1) = 0.949336;
  // -3.0 < y < -2.5
  Amatrix2(0,0) = 142.948;
  Amatrix2(0,1) = 3.66191;
  Amatrix2(1,0) = 16.7132*0.5;
  Amatrix2(1,1) = 4.28895*0.5;
  Amatrix2(2,0) = 16.7132*0.5;
  Amatrix2(2,1) = 4.28895*0.5;
  Amatrix2(3,0) = 5.04786;
  Amatrix2(3,1) = 2.04812;


  TDecompSVD svd(Amatrix,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd0(Amatrix0,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd1(Amatrix1,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  TDecompSVD svd2(Amatrix2,0.00001); //Decomp my A matrix with a tolerance of 1e-5
  Bool_t ok;
  Bool_t ok0;
  Bool_t ok1;
  Bool_t ok2;


  //====================
  // Corrected sigmas
  // -4 < y < -2.5
  //====================
  // CorrSigma  stat.
  // 2.2390    0.0307
  // 0.1184    0.0041
  // 0.1667    0.0074
  // 0.1695    0.0088
  TVectorD xsec4to25(4);
  xsec4to25.Zero();
  xsec4to25(0) = 2.2390;
  xsec4to25(1) = 0.1184;
  xsec4to25(2) = 0.1667;
  xsec4to25(3) = 0.1695;

  if( mode == 1 ){
    xsec4to25(0) = 2.2390 +TMath::Sqrt( 0.0307*0.0307 +  0.0855*0.0855*2.2390*2.2390 );
    xsec4to25(1) = 0.1184 +TMath::Sqrt( 0.0041*0.0041 +  0.2005*0.2005*0.1184*0.1184 );
    xsec4to25(2) = 0.1667 +TMath::Sqrt( 0.0074*0.0074 +  0.1085*0.1085*0.1667*0.1667 );
    xsec4to25(3) = 0.1695 +TMath::Sqrt( 0.0088*0.0088 +  0.0924*0.0924*0.1695*0.1695 );
  } else if (mode == 2){
    xsec4to25(0) = 2.2390 -TMath::Sqrt( 0.0307*0.0307 +  0.0855*0.0855*2.2390*2.2390 );
    xsec4to25(1) = 0.1184 -TMath::Sqrt( 0.0041*0.0041 +  0.2005*0.2005*0.1184*0.1184 );
    xsec4to25(2) = 0.1667 -TMath::Sqrt( 0.0074*0.0074 +  0.1085*0.1085*0.1667*0.1667 );
    xsec4to25(3) = 0.1695 -TMath::Sqrt( 0.0088*0.0088 +  0.0924*0.0924*0.1695*0.1695 );

  }









  TVectorD xsec4to35(4);
  xsec4to35.Zero();
  xsec4to35(0) = 1.5860;
  xsec4to35(1) = 0.0430;
  xsec4to35(2) = 0.0993;
  xsec4to35(3) = 0.0830;

  if( mode == 1 ){
    xsec4to35(0) = 1.5860 +TMath::Sqrt( 0.0307*0.0307 +  0.0690*0.0690*1.5860*1.5860 );
    xsec4to35(1) = 0.0430 +TMath::Sqrt( 0.0489*0.0489 +  0.3872*0.3872*0.0430*0.0430 );
    xsec4to35(2) = 0.0993 +TMath::Sqrt( 0.0084*0.0084 +  0.1092*0.1092*0.0993*0.0993 );
    xsec4to35(3) = 0.0830 +TMath::Sqrt( 0.0119*0.0119 +  0.0761*0.0761*0.0830*0.0830 );
  } else if (mode == 2){
    xsec4to35(0) = 1.5860 -TMath::Sqrt( 0.0307*0.0307 +  0.0690*0.0690*1.5860*1.5860 );
    xsec4to35(1) = 0.0430 -TMath::Sqrt( 0.0489*0.0489 +  0.3872*0.3872*0.0430*0.0430 );
    xsec4to35(2) = 0.0993 -TMath::Sqrt( 0.0084*0.0084 +  0.1092*0.1092*0.0993*0.0993 );
    xsec4to35(3) = 0.0830 -TMath::Sqrt( 0.0119*0.0119 +  0.0761*0.0761*0.0830*0.0830 );

  }







  TVectorD xsec35to3(4);
  xsec35to3.Zero();
  xsec35to3(0) = 2.3150;
  xsec35to3(1) = 0.1092;
  xsec35to3(2) = 0.1681;
  xsec35to3(3) = 0.1693;

  if( mode == 1 ){
    xsec35to3(0) = 2.3150 +TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.3150*2.3150 );
    xsec35to3(1) = 0.1092 +TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.1092*0.1092 );
    xsec35to3(2) = 0.1681 +TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.1681*0.1681 );
    xsec35to3(3) = 0.1693 +TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.1693*0.1693 );
  } else if (mode == 2){
    xsec35to3(0) = 2.3150 -TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.3150*2.3150 );
    xsec35to3(1) = 0.1092 -TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.1092*0.1092 );
    xsec35to3(2) = 0.1681 -TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.1681*0.1681 );
    xsec35to3(3) = 0.1693 -TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.1693*0.1693 );

  }










  TVectorD xsec3to25(4);
  xsec3to25.Zero();
  xsec3to25(0) = 2.6570;
  xsec3to25(1) = 0.2302;
  xsec3to25(2) = 0.2359;
  xsec3to25(3) = 0.2681;

  if( mode == 1 ){
    xsec3to25(0) = 2.6570 +TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.6570*2.6570 );
    xsec3to25(1) = 0.2302 +TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.2302*0.2302 );
    xsec3to25(2) = 0.2359 +TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.2359*0.2359 );
    xsec3to25(3) = 0.2681 +TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.2681*0.2681 );
  } else if (mode == 2){
    xsec3to25(0) = 2.6570 -TMath::Sqrt( 0.0307*0.0307 +  0.0679*0.0679*2.6570*2.6570 );
    xsec3to25(1) = 0.2302 -TMath::Sqrt( 0.0489*0.0489 +  0.2082*0.2082*0.2302*0.2302 );
    xsec3to25(2) = 0.2359 -TMath::Sqrt( 0.0084*0.0084 +  0.0975*0.0975*0.2359*0.2359 );
    xsec3to25(3) = 0.2681 -TMath::Sqrt( 0.0119*0.0119 +  0.0782*0.0782*0.2681*0.2681 );

  }



  const TVectorD photonuclear4to25 = svd.Solve(xsec4to25,ok);
  const TVectorD photonuclear4to35 = svd0.Solve(xsec4to35,ok0);
  const TVectorD photonuclear35to3 = svd1.Solve(xsec35to3,ok1);
  const TVectorD photonuclear3to25 = svd2.Solve(xsec3to25,ok2);
  cout << "Was it fine? " << ok << ok0 << ok1 << ok2 << endl;
  photonuclear4to25.Print();
  photonuclear4to35.Print();
  photonuclear35to3.Print();
  photonuclear3to25.Print();

  Double_t Result = -999.;
  if( rap == 0){
    if (element == 0) {
      Result = photonuclear4to25(0);
    } else {
      Result = photonuclear4to25(1);
    }
  } else if (rap == 1) {
    if (element == 0) {
      Result = photonuclear4to35(0);
    } else {
      Result = photonuclear4to35(1);
    }
  } else if (rap == 2) {
    if (element == 0) {
      Result = photonuclear35to3(0);
    } else {
      Result = photonuclear35to3(1);
    }
  } else if (rap == 3) {
    if (element == 0) {
      Result = photonuclear3to25(0);
    } else {
      Result = photonuclear3to25(1);
    }
  } else {
    Result = -999.;
  }
  return Result;
}
//_____________________________________________________________________________
void Computing(){
  Double_t HighSolution  = FittingPhotonuclear(0,0,0);
  Double_t HighSolutionU = FittingPhotonuclear(0,1,0)-FittingPhotonuclear(0,0,0);
  Double_t HighSolutionD = FittingPhotonuclear(0,2,0)-FittingPhotonuclear(0,0,0);
  Double_t LowSolution   = FittingPhotonuclear(1,0,0);
  Double_t LowSolutionU  = FittingPhotonuclear(1,1,0)-FittingPhotonuclear(1,0,0);
  Double_t LowSolutionD  = FittingPhotonuclear(1,2,0)-FittingPhotonuclear(1,0,0);

  // cout << "Element 0 (-4 < y < -2.5)= " << HighSolution << " + " << HighSolutionU << " - " << HighSolutionD << endl;
  // cout << "Element 1 (-4 < y < -2.5)= " << LowSolution  << " + " << LowSolutionU  << " - " << LowSolutionD  << endl;



  Double_t HighSolution0  = FittingPhotonuclear(0,0,1);
  Double_t HighSolutionU0 = FittingPhotonuclear(0,1,1)-FittingPhotonuclear(0,0,1);
  Double_t HighSolutionD0 = FittingPhotonuclear(0,2,1)-FittingPhotonuclear(0,0,1);
  Double_t LowSolution0   = FittingPhotonuclear(1,0,1);
  Double_t LowSolutionU0  = FittingPhotonuclear(1,1,1)-FittingPhotonuclear(1,0,1);
  Double_t LowSolutionD0  = FittingPhotonuclear(1,2,1)-FittingPhotonuclear(1,0,1);

  // cout << "Element 0 (-4 < y < -2.5)= " << HighSolution << " + " << HighSolutionU << " - " << HighSolutionD << endl;
  // cout << "Element 1 (-4 < y < -2.5)= " << LowSolution  << " + " << LowSolutionU  << " - " << LowSolutionD  << endl;
  //
  // cout << "Element 0 (-4 < y < -3.5)= " << HighSolution0 << " + " << HighSolutionU0 << " - " << HighSolutionD0 << endl;
  // cout << "Element 1 (-4 < y < -3.5)= " << LowSolution0  << " + " << LowSolutionU0  << " - " << LowSolutionD0  << endl;



  Double_t HighSolution1  = FittingPhotonuclear(0,0,2);
  Double_t HighSolutionU1 = FittingPhotonuclear(0,1,2)-FittingPhotonuclear(0,0,2);
  Double_t HighSolutionD1 = FittingPhotonuclear(0,2,2)-FittingPhotonuclear(0,0,2);
  Double_t LowSolution1   = FittingPhotonuclear(1,0,2);
  Double_t LowSolutionU1  = FittingPhotonuclear(1,1,2)-FittingPhotonuclear(1,0,2);
  Double_t LowSolutionD1  = FittingPhotonuclear(1,2,2)-FittingPhotonuclear(1,0,2);

  // cout << "Element 0 (-4 < y < -2.5)= " << HighSolution << " + " << HighSolutionU << " - " << HighSolutionD << endl;
  // cout << "Element 1 (-4 < y < -2.5)= " << LowSolution  << " + " << LowSolutionU  << " - " << LowSolutionD  << endl;
  //
  // cout << "Element 0 (-4 < y < -3.5)= " << HighSolution0 << " + " << HighSolutionU0 << " - " << HighSolutionD0 << endl;
  // cout << "Element 1 (-4 < y < -3.5)= " << LowSolution0  << " + " << LowSolutionU0  << " - " << LowSolutionD0  << endl;
  //
  // cout << "Element 0 (-3.5 < y < -3)= " << HighSolution1 << " + " << HighSolutionU1 << " - " << HighSolutionD1 << endl;
  // cout << "Element 1 (-3.5 < y < -3)= " << LowSolution1  << " + " << LowSolutionU1  << " - " << LowSolutionD1  << endl;



  Double_t HighSolution2  = FittingPhotonuclear(0,0,3);
  Double_t HighSolutionU2 = FittingPhotonuclear(0,1,3)-FittingPhotonuclear(0,0,3);
  Double_t HighSolutionD2 = FittingPhotonuclear(0,2,3)-FittingPhotonuclear(0,0,3);
  Double_t LowSolution2   = FittingPhotonuclear(1,0,3);
  Double_t LowSolutionU2  = FittingPhotonuclear(1,1,3)-FittingPhotonuclear(1,0,3);
  Double_t LowSolutionD2  = FittingPhotonuclear(1,2,3)-FittingPhotonuclear(1,0,3);

  cout << "Element 0 (-4 < y < -2.5)= " << HighSolution << " + " << HighSolutionU << " - " << HighSolutionD << endl;
  cout << "Element 1 (-4 < y < -2.5)= " << LowSolution  << " + " << LowSolutionU  << " - " << LowSolutionD  << endl;

  cout << "Element 0 (-4 < y < -3.5)= " << HighSolution0 << " + " << HighSolutionU0 << " - " << HighSolutionD0 << endl;
  cout << "Element 1 (-4 < y < -3.5)= " << LowSolution0  << " + " << LowSolutionU0  << " - " << LowSolutionD0  << endl;

  cout << "Element 0 (-3.5 < y < -3)= " << HighSolution1 << " + " << HighSolutionU1 << " - " << HighSolutionD1 << endl;
  cout << "Element 1 (-3.5 < y < -3)= " << LowSolution1  << " + " << LowSolutionU1  << " - " << LowSolutionD1  << endl;

  cout << "Element 0 (-3 < y < -2.5)= " << HighSolution2 << " + " << HighSolutionU2 << " - " << HighSolutionD2 << endl;
  cout << "Element 1 (-3 < y < -2.5)= " << LowSolution2  << " + " << LowSolutionU2  << " - " << LowSolutionD2  << endl;












  TH1F* histo = new TH1F("histo", "histo", 1000, -0.5, 999.5);
  Double_t rapidity[8] = {-3.25, -3.75, -3.25, -2.75, 3.25, 3.75, 3.25, 2.75};
  // Double_t rapidity[8] = {3.25, 3.75, 3.25, 2.75, -3.25, -3.75, -3.25, -2.75};
  Double_t Wgp[8];
  for (Int_t i = 0; i < 8; i++) {
    Wgp[i] = TMath::Sqrt(2.*2760*3.1*TMath::Exp(rapidity[i]));
  }
  Double_t y[8]   = {HighSolution, HighSolution0, HighSolution1, HighSolution2, LowSolution, LowSolution0, LowSolution1, LowSolution2};
  Double_t ey[8]  = {HighSolutionU, HighSolutionU0, HighSolutionU1, HighSolutionU2, LowSolutionU, LowSolutionU0, LowSolutionU1, LowSolutionU2};
  for (Int_t i = 0; i < 8; i++) {
    Int_t bin = histo->GetXaxis()->FindBin(Wgp[i]);
    cout << "bin = " << bin << endl;
    histo->SetBinContent( histo->GetXaxis()->FindBin(Wgp[i]),  y[i]);
    histo->SetBinError(   histo->GetXaxis()->FindBin(Wgp[i]), ey[i]);
  }
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTopMargin(0.14);
  // gPad->SetTickx(1);
  // gPad->SetTicky(1);
  // gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
  histo->SetTitle("");
  histo->GetXaxis()->SetTitleOffset(1.25);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->GetXaxis()->SetTitleSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.045);
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetYaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelFont(42);

  histo->GetXaxis()->SetTitle("W_{#gammaPb} [GeV]");
  histo->GetYaxis()->SetTitle("#sigma_{#gammaPb} [mb]");
  histo->GetYaxis()->SetRangeUser(0.002,0.2);
  // histo->GetYaxis()->SetRangeUser(-0.1,0.2);
  histo->GetXaxis()->SetRangeUser(10.,1000.);
  histo->SetLineWidth(5);
  histo->SetLineColor(2);
  histo->Draw("same");
  Double_t mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  Double_t bxmin = TMath::Power((mjpsi/1000.),2.);
  Double_t bxmax = TMath::Power((mjpsi/10.),2.);
  TF1 *fbx = new TF1("fbx","TMath::Power(([0]/x),2.)", bxmin, bxmax);
  fbx->SetParameter(0, mjpsi);

  TGaxis *axis = new TGaxis(1000., 0.2, 10., 0.2, "fbx", 510, "+G");
  axis->SetTextFont(42);
  axis->SetLabelFont(42);
  Double_t siz = 0.045;
  axis->SetTitleSize(siz); axis->SetLabelSize(siz);
  axis->SetLabelOffset(-0.035);
  axis->SetTitleOffset();
  axis->SetTitle("Bjorken-#it{x}");
  axis->Draw("same");

  TLatex *bxtit = new TLatex();
  bxtit->SetTextFont(42);
  bxtit->SetTextSize(siz);
  bxtit->SetTextAlign(31);
  bxtit->DrawLatex(1000., 0.2+450, "Bjorken-#it{x}");
  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  // latex5->DrawLatex(0.31,0.94,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // latex5->DrawLatex(0.31,0.7,"This thesis");
  // latex5->DrawLatex(0.31,0.8,"Coherent J/#psi");
  latex5->DrawLatex(0.31,0.78,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.7,"This thesis");
  latex5->DrawLatex(0.31,0.62,"Coherent J/#psi");

  gPad->SaveAs("pngResults/solutions.pdf", "recreate");














  // ==================
  // Suppression factor
  // ==================
  TH1F* SuppressionFactor = new TH1F("SuppressionFactor", "SuppressionFactor", 1000, -0.5, 999.5);
  // Double_t rapidity[8] = {-3.25, -3.75, -3.25, -2.75, 3.25, 3.75, 3.25, 2.75};
  // // Double_t rapidity[8] = {3.25, 3.75, 3.25, 2.75, -3.25, -3.75, -3.25, -2.75};
  // Double_t Wgp[8];
  // for (Int_t i = 0; i < 8; i++) {
  //   Wgp[i] = TMath::Sqrt(2.*2760*3.1*TMath::Exp(rapidity[i]));
  // }
  // Double_t y[8]   = {HighSolution, HighSolution0, HighSolution1, HighSolution2, LowSolution, LowSolution0, LowSolution1, LowSolution2};
  // Double_t ey[8]  = {HighSolutionU, HighSolutionU0, HighSolutionU1, HighSolutionU2, LowSolutionU, LowSolutionU0, LowSolutionU1, LowSolutionU2};
  // ========================================
  // Wgp(1) (GeV) k*(dn_1/dk)  sigma_2(gam+A)
  // ----------------------------------------
  // 2.1082E+01   1.9869E+02   1.5119E-02
  // 2.7060E+01   1.8317E+02   1.9706E-02
  // 3.4050E+01   1.6889E+02   2.3864E-02
  // 5.1928E+02   1.0599E+01   1.4670E-01
  // 6.6677E+02   4.0781E+00   1.7259E-01
  // 8.5615E+02   1.0018E+00   2.0305E-01
  Double_t IA[8]  = { 1.9706E-02, 1.5119E-02, 1.9706E-02, 2.3864E-02, 1.7259E-01, 2.0305E-01, 1.7259E-01, 1.4670E-01 };
  for (Int_t i = 0; i < 8; i++) {
    Int_t bin = histo->GetXaxis()->FindBin(Wgp[i]);
    cout << "bin = " << bin << endl;
    SuppressionFactor->SetBinContent( SuppressionFactor->GetXaxis()->FindBin(Wgp[i]),  TMath::Sqrt(y[i]/IA[i]));
    // SuppressionFactor->SetBinError(   SuppressionFactor->GetXaxis()->FindBin(Wgp[i]),  TMath::Sqrt(ey[i]/IA[i]));
    SuppressionFactor->SetBinError(   SuppressionFactor->GetXaxis()->FindBin(Wgp[i]),  ey[i]/(2.*TMath::Sqrt(y[i]*IA[i])) );
  }
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTopMargin(0.14);
  // gPad->SetTickx(1);
  // gPad->SetTicky(1);
  // gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
  SuppressionFactor->SetTitle("");
  SuppressionFactor->GetXaxis()->SetTitleOffset(1.25);
  SuppressionFactor->GetYaxis()->SetTitleOffset(1.25);
  SuppressionFactor->GetXaxis()->SetTitleSize(0.045);
  SuppressionFactor->GetYaxis()->SetTitleSize(0.045);
  SuppressionFactor->GetXaxis()->SetLabelSize(0.045);
  SuppressionFactor->GetYaxis()->SetLabelSize(0.045);
  SuppressionFactor->GetXaxis()->SetTitleFont(42);
  SuppressionFactor->GetYaxis()->SetTitleFont(42);
  SuppressionFactor->GetXaxis()->SetLabelFont(42);
  SuppressionFactor->GetYaxis()->SetLabelFont(42);

  SuppressionFactor->GetXaxis()->SetTitle("W_{#gammaPb} [GeV]");
  SuppressionFactor->GetYaxis()->SetTitle("S_{Pb} [a.u.]");
  SuppressionFactor->GetYaxis()->SetRangeUser(-0.3,1.5);
  SuppressionFactor->GetXaxis()->SetRangeUser(10.,1000.);
  SuppressionFactor->SetLineWidth(5);
  SuppressionFactor->SetLineColor(2);
  SuppressionFactor->Draw("same");
  // Double_t mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  // Double_t bxmin = TMath::Power((mjpsi/1000.),2.);
  // Double_t bxmax = TMath::Power((mjpsi/10.),2.);
  // TF1 *fbx = new TF1("fbx","TMath::Power(([0]/x),2.)", bxmin, bxmax);
  // fbx->SetParameter(0, mjpsi);

  TGaxis *axis2 = new TGaxis(1000., 1.5, 10., 1.5, "fbx", 510, "+G");
  axis2->SetTextFont(42);
  axis2->SetLabelFont(42);
  // Double_t siz = 0.045;
  axis2->SetTitleSize(siz); axis->SetLabelSize(siz);
  axis2->SetLabelOffset(-0.035);
  axis2->SetTitleOffset();
  axis2->SetTitle("Bjorken-#it{x}");
  axis2->Draw("same");

  TLatex *bxtit2 = new TLatex();
  bxtit2->SetTextFont(42);
  bxtit2->SetTextSize(siz);
  bxtit2->SetTextAlign(31);
  bxtit2->DrawLatex(1000., 1.5+450, "Bjorken-#it{x}");
  TLatex* latex6 = new TLatex();
  latex6->SetTextSize(0.045);
  latex6->SetTextFont(42);
  latex6->SetTextAlign(11);
  latex6->SetNDC();
  // latex5->DrawLatex(0.31,0.94,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // latex5->DrawLatex(0.31,0.7,"This thesis");
  // latex5->DrawLatex(0.31,0.8,"Coherent J/#psi");
  latex6->DrawLatex(0.4,0.78,"LHC18qr, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex6->DrawLatex(0.4,0.7,"This thesis");
  latex6->DrawLatex(0.4,0.62,"Coherent J/#psi");

  gPad->SaveAs("pngResults/SuppressionFactor.pdf", "recreate");

}
