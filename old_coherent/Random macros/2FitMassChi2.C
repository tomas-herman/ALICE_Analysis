
//
// program to fit the invariant mass distribution using Chi2 minimization
// and standard DQ formulas
//

// external headers
#include "CandRoot.h"


// jpsi mass
const Double_t mJPsiPDG = 3.096916;

// define cuts
#include "Cuts.h"

/////////////////////////////////////////////////////////////
// function from AliAnalysisMuMuJpsiResult::FitFunctionNA60New
/////////////////////////////////////////////////////////////

Double_t FitFunctionNA60New(Double_t *x,Double_t *par)
{
  // New Formulation of NA60 : 11 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = p1Left
  // par[4] = p2Left
  // par[5] = p3Left
  // par[6] = p1Right
  // par[7] = p2Right
  // par[8] = p3Right
  // par[9] = alphaLeft
  // par[10] = alphaRight

  
  const Double_t t = (x[0]-par[1])/par[2];
  
  Double_t sigmaRatio(0.);
  if( t < par[9] ) sigmaRatio = ( 1.0 + TMath::Power( par[3]*(par[9]-t), par[4]-par[5]*TMath::Sqrt(par[9] - t) ) );
  else if( t >= par[9] && t < par[10] ) sigmaRatio = 1;
  else if( t >= par[10] ) sigmaRatio = ( 1.0 + TMath::Power( par[6]*(t-par[10]), par[7]-par[8]*TMath::Sqrt(t - par[10]) ) );
  
  return par[0]*TMath::Exp( -(1/2.)*TMath::Power(t/sigmaRatio,2.));
  
}


/////////////////////////////////////////////////////////////
// function from AliAnalysisMuMuJpsiResult::FitFunctionSignalCrystalBallExtended
/////////////////////////////////////////////////////////////

Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par)
{
  // Extended Crystal Ball : 7 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n
  // par[5] = alpha'
  // par[6] = n'
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }
  
  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  
  if (t >= absAlpha2) //right tail
  {
    
    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }
  
  return 0. ;
}

/////////////////////////////////////////////////////////////
// function for gg
/////////////////////////////////////////////////////////////
Double_t FitFunctionExp(Double_t *x,Double_t *par)
{
    Double_t m = x[0];
    return par[0]*TMath::Exp(-par[1]*m);
}

/////////////////////////////////////////////////////////////
// function for gg, j/psi, psip
/////////////////////////////////////////////////////////////
Double_t FitFunctionAll(Double_t *x,Double_t *par)
{
  // par[0-1] exp
  // par[2-9] cb jpsi
  // par2 cb psip
  Double_t m = x[0];
  Double_t gg = par[0]*TMath::Exp(-par[1]*m);
  Double_t jpsi = FitFunctionSignalCrystalBallExtended(x,&par[2]);
  Double_t par2[7] = {
    par[9],
    3.68609+(par[3]-3.096916)/3.096916*3.68609,
     par[4]/3.096916*3.68609,
    par[5],
    par[6],
    par[7],
    par[8]
  };
  Double_t psip = FitFunctionSignalCrystalBallExtended(x,par2);
  
  return gg+jpsi+psip;
}

/////////////////////////////////////////////////////////////
// read data into histo
/////////////////////////////////////////////////////////////
void ReadData(char *InFile, TH1* histo, Double_t y_min, Double_t y_max)
{
  Double_t pt,y,m;
  ifstream ifs;
  ifs.open(InFile);
  while (!ifs.eof()) {
    ifs>>pt;
    ifs>>y;
    ifs>>m;
    if (y_min>y) continue;
    if (y_max<y) continue;
    if (PtMax<pt) continue;
    histo->Fill(m);
  }
  ifs.close();
}

/////////////////////////////////////////////////////////////
// fills histo with appropriate data
/////////////////////////////////////////////////////////////

void FillHisto(TH1* Mass, Int_t OptFile, Double_t y_min, Double_t y_max)
{
    // read data into histo
  if (OptFile == 1) {
    ReadData("TXT/LHC16r_std.txt",Mass, y_min, y_max);
  }
  if (OptFile == 11) {
    ReadData("TXT/LHC16s_std.txt",Mass, y_min, y_max);
  }
  if (OptFile == 20) {
    ReadData("TXT/LHC16r_std_MCgp.txt",Mass, y_min, y_max);
  }
  if (OptFile == 21) {
    ReadData("TXT/LHC16s_std_MCgp.txt",Mass, y_min, y_max);
  }
}

TF1* SetFitModel(Int_t OptFile, Int_t OptFit, TH1* Mass)
{
  TF1 *FitModel = NULL;
  if (OptFit == 1) {
    FitModel = new TF1("FitModel",FitFunctionSignalCrystalBallExtended, mMinFitMass, mMaxFitMass,7);
    FitModel->SetParNames("N","#mu","#sigma","#alpha","n","#alpha^{'}","n^{'}");
    FitModel->SetParameters(Mass->GetEntries(),mJPsiPDG,0.08,0.01,5,2,5);
    FitModel->SetParLimits(0,1.0,2.0*Mass->GetEntries());
    FitModel->SetParLimits(1, 3.0,3.2);
    FitModel->SetParLimits(2, 0.05,0.15);
    FitModel->SetParLimits(3, 0.0,2.0);
    FitModel->SetParLimits(4, 0.0,20.0);    
    FitModel->SetParLimits(5, 0.0,5.0);
    FitModel->SetParLimits(6, 0.0,20.0);    
  }
  if (OptFit == 2) {
    FitModel = new TF1("FitModel",FitFunctionNA60New, mMinFitMass, mMaxFitMass,11);
    FitModel->SetParNames("N","#mu","#sigma","p1L","p2L","p3L","p1R","p2R","p3R","#alpha_{L}","#alpha_{R}");
    FitModel->SetParameters(Mass->GetEntries(),mJPsiPDG,0.08,0.1,5,2,5,20,30,0.1,0.1);
    FitModel->SetParLimits(0,1.0,2.0*Mass->GetEntries());
    FitModel->SetParLimits(1, 3.0,3.2);
    FitModel->SetParLimits(2, 0.05,0.15);
    FitModel->SetParLimits(3, 0,100);
    FitModel->SetParLimits(4,0,100);
    FitModel->SetParLimits(5, 0,100);
    FitModel->SetParLimits(6, 0,100);
    FitModel->SetParLimits(7,0,100);
    FitModel->SetParLimits(8, 0,100);
    FitModel->SetParLimits(9, 0.0,2.0);
    FitModel->SetParLimits(10, 0.0,2.0);    
  }
  if (OptFit == 3) {
    FitModel = new TF1("FitModel",FitFunctionAll, mMinFitMass, mMaxFitMass,10);
    FitModel->SetParNames("N_{exp}","#lambda","N","#mu","#sigma","#alpha","n","#alpha^{'}","n^{'}",
			  "N_{#psi^{'}}");
    FitModel->SetParameters(0.5*Mass->GetEntries(),1.0,
			    Mass->GetEntries(),mJPsiPDG,0.08,0.01,5,2,5,
			    0.1*Mass->GetEntries());
    FitModel->SetParLimits(0,1.0,2.0*Mass->GetEntries());
    FitModel->SetParLimits(1, 0,10);
    FitModel->SetParLimits(2,1.0,2.0*Mass->GetEntries());    
    FitModel->SetParLimits(3, 3.0,3.2);
    FitModel->SetParLimits(4, 0.05,0.15);
    FitModel->SetParLimits(5, 0.0,2.0);
    FitModel->SetParLimits(6, 0.0,20.0);    
    FitModel->SetParLimits(7, 0.0,5.0);
    FitModel->SetParLimits(8, 0.0,20.0);
    FitModel->SetParLimits(9,0.0,2.0*Mass->GetEntries());
    if (OptFile == 1) {
      FitModel->FixParameter(5,1.235);
      FitModel->FixParameter(6,2.659);
      FitModel->FixParameter(7,2.62);
      FitModel->FixParameter(8,3.05);
    } else if (OptFile == 11) {
      FitModel->FixParameter(5,1.11);
      FitModel->FixParameter(6,4.588);
      FitModel->FixParameter(7,2.707);
      FitModel->FixParameter(8,3.339);
    }
  }
  return FitModel;
}

/////////////////////////////////////////////////////////////
// main function
/////////////////////////////////////////////////////////////

void FitMassChi2(Int_t OptFile, Int_t OptFit, Double_t y_min, Double_t y_max)
{
  // define histogram to be fitted
  Int_t nBins = (mMaxFitMass-mMinFitMass)*20; // 20 bins per GeV
  TH1D *Mass = new TH1D("Mass","Mass",nBins,mMinFitMass,mMaxFitMass);
  Mass->SetTitle(";m_{#mu^{+}#mu^{-}} GeV/#it{c}^{2};Dimuon candidates/(50 MeV/#it{c^{2}})");

  // read data into histo
  FillHisto(Mass,OptFile, y_min, y_max);
  Mass->SetBinErrorOption(TH1::kPoisson);
  // set fit model
  TF1 *FitModel = SetFitModel(OptFile,OptFit,Mass);

  // draw and fit
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);    
  Mass->SetMarkerStyle(20);
  Mass->SetMarkerColor(kBlack);  Mass->SetLineColor(kBlack);
  Mass->Draw("e");
  TFitResultPtr r = Mass->Fit("FitModel","RlS","",1.5,4.5);
  Double_t chi2 = FitModel->GetChisquare();
  chi2/=((Double_t)FitModel->GetNDF());
  if (OptFit == 3) {
    // get error matrix
    TMatrixDSym cov = r->GetCovarianceMatrix();
    Double_t *fullmat;
    fullmat = cov.GetMatrixArray();
    // --> set matrices for gg, j/psi and for psip
    Double_t gg_cov[4];
    Double_t jpsi_cov[49];
    Double_t psip_cov[49];
    for(Int_t i=0;i<2;i++){
      for(Int_t j=0;j<2;j++) gg_cov[2*i+j] = fullmat[i*10+j];
    }
    for(Int_t i=0;i<7;i++){
      for(Int_t j=0;j<7;j++) psip_cov[7*i+j] = jpsi_cov[7*i+j] = fullmat[22+i*10+j];
    }
    // change first line 
    for(Int_t j=1;j<7;j++) psip_cov[j] = fullmat[9*10+2+j];
    psip_cov[0] = fullmat[99];
    // change first row
    for(Int_t i=1;i<7;i++) psip_cov[7*i] = fullmat[29+i*10];

    // function for the j/psi
    TF1 *jpsi_fit = new TF1("jpsi_fit",FitFunctionSignalCrystalBallExtended, mMinFitMass, mMaxFitMass,7);
    Double_t jpsi_par[7];
    for(Int_t i=0;i<7;i++) {
      jpsi_fit->SetParameter(i,FitModel->GetParameter(i+2));
      jpsi_par[i] = FitModel->GetParameter(i+2);
    }
    jpsi_fit->SetLineColor(kBlue);
    jpsi_fit->Draw("same");
    Double_t njpsi = jpsi_fit->Integral(mMinFitMass, mMaxFitMass)/Mass->GetBinWidth(1);
    Double_t njpsi_err =
      jpsi_fit->IntegralError(mMinFitMass, mMaxFitMass,jpsi_par,jpsi_cov)/Mass->GetBinWidth(1);
    cout << " N_jpsi = " << njpsi << " +/- " << njpsi_err << endl;
    cout << " (2.8,3.2) N_jpsi = " << jpsi_fit->Integral(2.8, 3.2)/Mass->GetBinWidth(1)
	 << " +/- " << jpsi_fit->IntegralError(2.8, 3.2,jpsi_par,jpsi_cov)/Mass->GetBinWidth(1) << endl;
    
    // function for the psi prime
    TF1 *psip_fit = new TF1("psip_fit",FitFunctionSignalCrystalBallExtended, mMinFitMass, mMaxFitMass,7);
    Double_t psip_par[7];
    for(Int_t i=0;i<7;i++) {
      psip_fit->SetParameter(i,FitModel->GetParameter(i+2));
      psip_par[i] = FitModel->GetParameter(i+2);
    }
    psip_fit->SetParameter(0,FitModel->GetParameter(9));
    psip_par[0] = FitModel->GetParameter(9);
    psip_fit->SetParameter(1,3.68609+(FitModel->GetParameter(3)-3.096916)/3.096916*3.68609);
    psip_par[1] = 3.68609+(FitModel->GetParameter(3)-3.096916)/3.096916*3.68609;
    psip_fit->SetParameter(2,FitModel->GetParameter(4)/3.096916*3.68609);
    psip_par[2] = FitModel->GetParameter(4)/3.096916*3.68609;
    psip_fit->SetLineColor(kCyan);
    psip_fit->Draw("same");
    Double_t npsip = psip_fit->Integral(mMinFitMass, mMaxFitMass)/Mass->GetBinWidth(1);
    Double_t npsip_err =
      psip_fit->IntegralError(mMinFitMass, mMaxFitMass,psip_par,psip_cov)/Mass->GetBinWidth(1);
    cout << " N_psip = " << npsip << " +/- " << npsip_err << endl;
    
    // function for the gg
    TF1 *gg_fit = new TF1("gg_fit",FitFunctionExp, mMinFitMass, mMaxFitMass,2);
    Double_t gg_par[2];
    for(Int_t i=0;i<2;i++) {
      gg_fit->SetParameter(i,FitModel->GetParameter(i));
      gg_par[i] = FitModel->GetParameter(i);
    }
    Double_t ngg = gg_fit->Integral(mMinFitMass, mMaxFitMass)/Mass->GetBinWidth(1);
    Double_t ngg_err =
      gg_fit->IntegralError(mMinFitMass, mMaxFitMass,gg_par,gg_cov)/Mass->GetBinWidth(1);
					       
    cout << " N_gg = " << ngg << " +/- " << ngg_err << endl;
    cout << " (2.8,3.2) N_gg = " << gg_fit->Integral(2.8, 3.2)/Mass->GetBinWidth(1)
	<< " +/- " << gg_fit->IntegralError(2.8, 3.2, gg_par,gg_cov)/Mass->GetBinWidth(1) << endl;

    // print results in plot
    Double_t ymax = 1.2*Mass->GetMaximum();
    char text3[100];
    sprintf(text3,"#chi^{2} = %3.2f",chi2);
    TLatex *l0 = new TLatex(1.55,ymax,text3);
    l0->Draw();
   sprintf(text3,"N_{J/#psi}= %3.0f #pm %2.0f",njpsi,njpsi_err);
    TLatex *l1 = new TLatex(2.55,ymax,text3);
    l1->Draw();
    sprintf(text3,"N_{#psi(2S)}= %3.0f #pm %2.0f",npsip,npsip_err);
    TLatex *l2 = new TLatex(3.55,ymax,text3);
    l2->Draw();

  gPad->SetLogy();

  }
}


/*
fitting from 2.7 to 3.35, CB gives chi2= 1126/123
   1  N            5.38840e+04   5.21155e+01   4.46390e-08  -2.70282e-02
   2  #mu          3.10957e+00   9.93014e-05   2.03010e-07  -4.94607e-03
   3  #sigma       7.16021e-02   8.31791e-05  -1.64981e-06  -7.36400e-03
   4  #alpha       9.22823e-01   4.10942e-03  -1.68928e-06  -1.94667e-03
   5  n            9.09453e+00   2.53389e-01   5.31124e-06   5.01778e-04
   6  #alpha^{'}   2.35165e+00   3.85075e-02  -9.71681e-05   1.40846e-03
   7  n^{'}        1.03635e+01   3.35085e+00   2.04929e-03  -2.61955e-03

NA chi2= 1000/119
   1  N            5.32024e+04   5.04656e+01   9.53736e-06  -1.56643e+04
   2  #mu          3.11234e+00   1.31231e-04  -8.89718e-04   1.14707e+03
   3  #sigma       7.06413e-02   8.74137e-05   1.24411e-03  -1.64451e+03
   4  p1L          1.48456e-01   1.55868e-03   5.26846e-04  -6.79863e+04
   5  p2L          3.38600e+00   1.08120e-01  -5.70566e-03   9.00221e+03
   6  p3L          7.97959e-01   3.75671e-02  -3.04638e-03  -9.12370e+03
   7  p1R          2.55344e+01   3.58798e-02  -1.02746e-03   0.00000e+00
   8  p2R          1.57850e+01   7.01281e-03   1.33727e-04   7.58787e-02
   9  p3R          3.53374e+01   7.46450e-03  -1.16611e-04  -9.93906e+03
  10  #alpha_{L}   1.26187e+00   7.44457e-02   1.02224e-01   6.09460e+02
  11  #alpha_{R}   1.68104e-02   7.17141e-05   2.55838e-04   0.00000e+00
*/
