//
// program to test SVD for coherent production in neutron classes
//

//-------------------------------------------------------
// Header files
//-------------------------------------------------------

// c++ headers
#include <iostream>
#include <fstream>

// root headers
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"

//-------------------------------------------------------
// Globals
//-------------------------------------------------------

const int nXS = 2; // two cross sections at each rapidity

//-------------------------------------------------------
// Functions
//-------------------------------------------------------

//-------------------------------------------------------
// set some reasonable uncertainties in mb
void getUpcErr(const int optRap, double *upcErr)
{
  if (optRap==1) {
    upcErr[0] = 0.114;
    upcErr[1] = 0.013;
    upcErr[2] = 0.014;
  } else   if (optRap==2) {
    upcErr[0] = 0.160;
    upcErr[1] = 0.018;
    upcErr[2] = 0.018;
  } else   if (optRap==3) {
    upcErr[0] = 0.183;
    upcErr[1] = 0.024;
    upcErr[2] = 0.024;
  }

}


//-------------------------------------------------------
// compute the upc cross sections for each neutron class
void getUpcXS( const int n, double *gXS, TMatrixD &fluxM, double *upcXS)
{
  for(int i=0;i<n;i++) {
    upcXS[i] = fluxM(i,0)*gXS[0]+fluxM(i,1)*gXS[1];
  }
}

//-------------------------------------------------------
// set the flux matrix
// fluxes from Simone who got them from Starlight
void setFlux(const int optRap,TMatrixD &fluxM)
{
  if (optRap == 1) {
      fluxM(0,0) = 172.551;
      fluxM(0,1) = 0.200768;
      fluxM(1,0) = 16.6506*0.5;
      fluxM(1,1) = 0.421108*0.5;
      fluxM(2,0) = 5.02952;
      fluxM(2,1) = 0.275626;
  } else if (optRap == 2) {
      fluxM(0,0) = 157.789;
      fluxM(0,1) = 1.07192;
      fluxM(1,0) = 16.6849*0.5;
      fluxM(1,1) = 1.70239*0.5;
      fluxM(2,0) = 5.03916;
      fluxM(2,1) = 0.94933;
  } else if (optRap == 3) {
      fluxM(0,0) = 142.948;
      fluxM(0,1) = 3.66191;
      fluxM(1,0) = 16.7132*0.5;
      fluxM(1,1) = 4.28895*0.5;
      fluxM(2,0) = 5.04786;
      fluxM(2,1) = 2.04812;
  }
}


//-------------------------------------------------------
// set number of neutron classes
int setClasses(const int optRap)
{
  if (optRap>0 && optRap<4) return 3;
  return 4;
}

//-------------------------------------------------------
// set photonuclear cross section according to the hs model
// in mb
void setPredictedXS(const int optRap, double *xs)
{
  if (optRap==1) { //y=3.8
    xs[0] = 0.00904619;
    xs[1] = 0.0665287;
  } else  if (optRap==2) { // y=3.2
    xs[0] = 0.0107257;
    xs[1] = 0.0586511;
  } else if (optRap==3) { // y=2.8
    xs[0] = 0.0121004;
    xs[1] = 0.0536615;
  }

}

//-------------------------------------------------------
// set rapidity range
void setRap(const int optRap, double *rap)
{
  if (optRap==1) {
    rap[0] = -4;
    rap[1] = -3.5;
  } else if (optRap==2) {
    rap[0] = -3.5;
    rap[1] = -3.;
  } else if (optRap==3) {
    rap[0] = -3;
    rap[1] = -2.5;
  }

}

//-------------------------------------------------------
// entry point
void testSVD(int optRap)
// optRap: chose the rapidity bin to be studied
{
  // control printing
  bool doPrint = kTRUE;

  // check input
  if (optRap<1 || optRap>3) {
    if (doPrint) cout << " Not a valid option. Bye." << endl;
    return;
  }
  
  //-----------------
  // initialisation
  //-----------------

  // for clarity leave some space in the output
  cout << endl;
  cout << endl;

  // set rapidity range
  double rap[2];
  setRap(optRap, rap);
  if (doPrint) {
    cout << " Rapidity range: (" << rap[0] << "," << rap[1] << ")" << endl;
  }
  
  // set input photonuclear cross sections
  double gXS[2];
  setPredictedXS(optRap,gXS);
  if (doPrint) {
    cout << endl;
    cout << " Input cross sections are: " << gXS[0] << " and " << gXS[1] << " in mb" << endl;
  }
      
  // set number of classes to be used
  const int nClasses = setClasses(optRap);
  if (doPrint) {
    //    cout << " Using " << nClasses << " neutron classes "  << endl;
  }
  

  // set the flux matrix
  TMatrixD fluxM(nClasses, nXS);
  setFlux(optRap,fluxM);
  if (doPrint) {
    cout << endl;
    cout << " Flux matrix (0n0n,Xn0n,XnXn)x(low W, high W): " << endl;
    for(int i=0;i<nClasses;i++) cout << "   " << fluxM(i,0) << " " << fluxM(i,1) << endl;
  }

  // compute the upc predictions
  double upcXS[nClasses];
  getUpcXS(nClasses,gXS,fluxM,upcXS);
  if (doPrint) {
    cout << endl;
    cout << " Predicted UPC cross sections are: " << endl;
    cout << "   "; 
    for(int i=0;i<nClasses;i++) cout << upcXS[i] << " ";
    cout << "mb"<< endl;
  }

  // compute the upc predictions
  double upcErr[nClasses];
  getUpcErr(optRap,upcErr);
  if (doPrint) {
    cout << endl;
    cout << " Uncertainties to be used are: " << endl;
    cout << "   "; 
    for(int i=0;i<nClasses;i++) cout << upcErr[i] << " (" << (upcErr[i]/upcXS[i]) << " %) ";
    cout << endl;
  }

  //-----------------
  // analysis part
  //-----------------

  // for clarity leave some space in the output
  cout << endl;
  cout << endl;
  
  // number of tries
  const int nTries = 2500;

  // space to store the results of the 'measured' cross sections
  TH1D *lowH = new TH1D("lowH",Form("y in (%0.2f,%0.2f), xs = %0.5f mb; #sigma(#gammaPb) (mb); Entries",
				    rap[0],rap[1],gXS[0]),100,0.25*gXS[0],1.75*gXS[0]);
  TH1D *highH = new TH1D("highH",Form("y in (%0.2f,%0.2f), xs = %0.5f mb; #sigma(#gammaPb) (mb); Entries",
				      rap[0],rap[1],gXS[1]),100,0.1*gXS[1],1.9*gXS[1]);
  
  // decompose the flux matrix
  double tol = 0.00001;
  TDecompSVD svd(fluxM,tol); //Decomp my A matrix with a tolerance of 1e-5

  // perform the tries
  gRandom->SetSeed(0); // set a new seed
  TVectorD mUpcXs(nClasses); // to store the 'measured' upc cross section
  for(int i=0;i<nTries;i++) {
    // randomise the predicted upc cross sections
    bool isNegative = kFALSE;
    for(int j=0;j<nClasses;j++) {
      mUpcXs(j) = gRandom->Gaus(upcXS[j],upcErr[j]);
      if (mUpcXs(j)<0) isNegative=kTRUE;
    }
    if (isNegative) continue;
    
    // get the 'measured' photonuclear cross sections
    bool isOK = kFALSE;
    const TVectorD svdXS = svd.Solve(mUpcXs,isOK);
    if (isOK) { // store the results
      lowH->Fill(svdXS(0));
      highH->Fill(svdXS(1));      
    }
  }

  //-----------------
  // show results
  //-----------------
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  TCanvas *c = new TCanvas("c","c",1200,600);
  c->Divide(2,1);
  c->cd(1);
  lowH->Fit("gaus");
  lowH->Draw();
  c->cd(2);
  highH->Fit("gaus");
  highH->Draw();
  
}

