//
// program to fit data on coherent j/psi photoproduction in PbPb
// in forward neutron classes to extract the gPb cross section
//
// It uses the prescription in
// Eur.Phys.J. C63 (2009) 625-678, section 9, eq (31)
// https://inspirehep.net/record/817368?ln=en
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
#include <TMinuit.h>

// -------------------------------------------------------------------------------
//  data (global variables)
// -------------------------------------------------------------------------------
#include "inputData.h"

// -------------------------------------------------------------------------------
//  other global variables
// -------------------------------------------------------------------------------
double totalUncorErr[nData];
double totalCorErr[nData];

int nPhysPar = 2; // number of physics (as opposed to nuisance) parameters
double gChi2 = 0; // to keep track of the chi2
double gChi2Point[nData];
int gFitType = 0; // to keep track of the current type of fit
bool gPrintMessages = kFALSE; // print or not messages
bool gMidBin = kFALSE;

// -------------------------------------------------------------------------------
//  define chi2 without correlated uncertainties
// -------------------------------------------------------------------------------

void fcnChi2ModelSimple(int &npar, double *gin, double &f, double *par, int iflag)

{
  // ch2 = sum_i [m_i-mu_i]^2/D
  // D = (d_i_stat*m_i)^2+(d_i,unc*m_i)^2
  // mu_i is the measurement
  // d_i = relative uncertainty (either stat or uncorr)

  // get chi2
  double chi2 = 0;
  for (int i=0; i<nData;i++) {
    double mu_i = crossSection[i];
    double m_i = 0;
    if (gMidBin) {
      m_i = 2.0*fluxL[i]*par[0];
    } else {
      m_i = fluxL[i]*par[0] + fluxH[i]*par[1];
    }
    double d_stat_i =  statErr[i]/crossSection[i];
    double d_unc_i = totalUncorErr[i];
    double D = (d_stat_i*d_stat_i*m_i*m_i)+((d_unc_i*m_i)*(d_unc_i*m_i));
    gChi2Point[i] = ((m_i-mu_i)*(m_i-mu_i)/D);
    chi2 += gChi2Point[i];
    if ((m_i-mu_i)<0) gChi2Point[i] *= -1;
  }
  f = chi2;
  gChi2 = f;
}


// -------------------------------------------------------------------------------
//  define chi2 for the full model
// -------------------------------------------------------------------------------

void fcnChi2ModelFull(int &npar, double *gin, double &f, double *par, int iflag)

{
    // ch2 = sum_i [m_i-mu_i-Sij]^2/D + Sbj
    // Sij = sum_j g_ij*m_i*b_j
    // Sbj = sum_j (b_j)^2
    // D = d_i_stat^2*mu_i*(m_i-Sij)+(d_i,unc*m_i)^2
    // mu_i is the measurement
    // g_ij relative normalization uncertainty at point i from source j
    // d_i = relative uncertainty (either stat or uncorr)

    // set bj
    int nCorErr = nGloCorErr;
    if (gFitType == 2) nCorErr = nGloCorErr+nLocCorErr;
    double bj[nCorErr];
    for(int j=0;j<nCorErr;j++) bj[j] = par[nPhysPar+j];

    // set Sbj
    double Sbj = 0;
    for(int j=0;j<nCorErr;j++) Sbj += (bj[j]*bj[j]);

    // get chi2
    double chi2 = 0;
    for (int i=0; i<nData;i++) {
        double mu_i = crossSection[i];
	double m_i = 0;
	if (gMidBin) {
	  m_i = 2.0*fluxL[i]*par[0];
	} else {
	  m_i = fluxL[i]*par[0] + fluxH[i]*par[1];
	}
        double Sij = 0;
        for(int j=0;j<nGloCorErr;j++) Sij += m_i*bj[j]*gloCorErr[j][i];
        if (gFitType == 2) {
            for(int j=nGloCorErr;j<nCorErr;j++) Sij += m_i*bj[j]*locCorErr[j-nGloCorErr][i];
        }
        double d_stat_i =  statErr[i]/crossSection[i];
        double d_unc_i = totalUncorErr[i];
        double D = (d_stat_i*d_stat_i*mu_i*(m_i-Sij))+((d_unc_i*m_i)*(d_unc_i*m_i));
        chi2 += ((m_i-mu_i-Sij)*(m_i-mu_i-Sij)/D);
    }
    chi2 += Sbj;
    f = chi2;
    gChi2 = f;
}

// -------------------------------------------------------------------------------
//  do the fit
// -------------------------------------------------------------------------------

void doFit(double *fitResult)
// -------------------------------------------------------------------------------
//  gFitType =
//   0 => only uncorrelated uncertainties
//   1 => uncorrelated and globally correlated uncertainties
//   2 => uncorrelated, globally and locally correlated uncertainties
// -------------------------------------------------------------------------------
{
    ///////////////////////////////////////////////////////////////////////////////
    // set up minuit
    ///////////////////////////////////////////////////////////////////////////////

    // initialize minuit with a the number of parameters
    int nCorErr = nGloCorErr;
    if (gFitType == 2) nCorErr = nGloCorErr+nLocCorErr;
    int nPar = nPhysPar+nCorErr;
    if (gFitType == 0) nPar = nPhysPar;
    TMinuit *myMinuit = new TMinuit(nPar);
      
    // set the function with the minimization process
    if (gFitType == 0) myMinuit->SetFCN(fcnChi2ModelSimple);
    else myMinuit->SetFCN(fcnChi2ModelFull);
    
    // define parameters
    // args = parNo, name, initVal, initErr, lowerLim, upperLim
    if (gMidBin) {
      myMinuit->DefineParameter(0,"#sigma_{L}",0.05,0.001,0.0,1.);
    } else {
      myMinuit->DefineParameter(0,"#sigma_{L}",0.05,0.001,0.0,1.);
      myMinuit->DefineParameter(1,"#sigma_{H}",0.05,0.005,0.0,1.);
    }
    if(gFitType>0) {
        for(int j=0;j<nCorErr;j++) {
            myMinuit->DefineParameter(nPhysPar+j,Form("b_{%d}",j),0.001,0.0001,-1,1.);
        }
    } // end gFitType>0
    
    // do the fit with migrad
    myMinuit->SetMaxIterations(500);
    myMinuit->Migrad();

    ///////////////////////////////////////////////////////////////////////////////
    // get results
    ///////////////////////////////////////////////////////////////////////////////
    double Cov[nPar*nPar];
    myMinuit->mnemat(Cov,nPar);
    double xsL = 0;  // photo-nuclear cross section for low energy photon
    double xsLerr = 0;  // its error
    myMinuit->GetParameter(0,xsL,xsLerr);
    double xsH = 0;  // photo-nuclear cross section for high energy photon
    double xsHerr = 0;  //its error
    double covLL = 0;
    double covHH = 0;
    double covLH = 0;
    double rho = 0;
    if (gMidBin) {
      xsH = xsL;
      xsHerr = xsLerr;
    } else {
      myMinuit->GetParameter(1,xsH,xsHerr);
      covLL = Cov[0];
      covHH = Cov[nPar+1];
      covLH = Cov[1];
      rho = covLH/TMath::Sqrt(covLL*covHH);
    }
    // prepare result (use mub for xs)
    fitResult[0] = xsL*1000.;
    fitResult[1] = xsLerr*1000.;
    fitResult[2] = xsH*1000.;
    fitResult[3] = xsHerr*1000.;
    fitResult[4] = rho;
    fitResult[5] = gChi2;
}

// -------------------------------------------------------------------------------
// Print input data for mid rap
// -------------------------------------------------------------------------------
void printInputData()
{
  // header
  cout << endl << " rapidity bin $(" << yMin <<", "<<yMax<<")$" << endl;
  cout << "\\hline" << endl;
  if (yMin>2) // cout << " Source & 0n0n & 0nXn & Xn0n & XnXn \\\\" << endl;
    cout << " Source & 0n0n & Xn0n & XnXn \\\\" << endl;
  else  cout << " Source & 0n0n & 0nXn + Xn0n & XnXn \\\\" << endl;
  cout << "\\hline" << endl;
  // flux
  cout << " Flux low energy ";
  for(int j=0;j<nData;j++) cout << Form("& %4.2f ",fluxL[j]);
  cout <<  " \\\\ "<< endl;
  cout << " Flux high energy ";
  for(int j=0;j<nData;j++) cout << Form("& %4.2f ",fluxH[j]);
  cout <<  " \\\\ "<< endl;
  cout << "\\hline" << endl;
  // uncorrelated
  cout << " Cross section (mb)  ";
  for(int j=0;j<nData;j++) cout << Form("& %4.2f ",crossSection[j]);
  cout <<  " \\\\ "<< endl;
  cout << " Stat. (mb)  ";
  for(int j=0;j<nData;j++) cout << Form("& %4.2f ",statErr[j]);
  cout <<  " \\\\ "<< endl;
  for(int i=0;i<nUncorErr;i++) {
    cout << Form(" %s ",uncorErrNames[i]);
    for(int j=0;j<nData;j++) cout << Form("& %4.2f ",100*uncorErr[i][j]);
    cout <<  " \\\\ "<< endl;
  }
  cout << "\\hline" << endl;
  // global
  for(int i=0;i<nGloCorErr;i++) {
    cout << Form(" %s ",gloCorErrNames[i]);
    for(int j=0;j<nData;j++) cout << Form("& %4.2f ",100*gloCorErr[i][j]);
    cout <<  " \\\\ "<< endl;
  }
  cout << "\\hline" << endl;
  // local
  for(int i=0;i<nLocCorErr;i++) {
    cout << Form(" %s ",locCorErrNames[i]);
    for(int j=0;j<nData;j++) cout << Form("& %4.2f ",100*locCorErr[i][j]);
    cout <<  " \\\\ "<< endl;
  }
  cout << "\\hline" << endl;
  // migration
  for(int i=0;i<nMigErr;i++) {
    cout << Form(" %s ",migErrNames[i]);
    for(int j=0;j<nData;j++) cout << Form("& %4.2f ",100*migErr[i][j]);
    cout <<  " \\\\ "<< endl;
  }
  cout << "\\hline" << endl;

}

// -------------------------------------------------------------------------------
// Print header for table of results
// -------------------------------------------------------------------------------
void printFitHeader()
{
    cout << endl << " rapidity bin $(" << yMin <<", "<<yMax<<")$" << endl;
    cout << "\\hline" << endl;
    cout << " Fit & low energy ($\\mu$b) & high energy ($\\mu$b) & correlation & $\\chi^2$ \\\\" << endl;
    cout << "\\hline" << endl;
}


// -------------------------------------------------------------------------------
// Print result of the fit
// -------------------------------------------------------------------------------
void printFitResult(const char *fitType, double *fitResult)
{
    // latex format
  cout << fitType << " & "
       << Form("%4.2f $\\pm$ %4.2f", fitResult[0], fitResult[1]) << " & "
       << Form("%4.2f $\\pm$ %4.2f", fitResult[2], fitResult[3]) << " & "
       << Form("%4.2f", fitResult[4]) << " & "
       << Form("%4.2f", fitResult[5]) << " \\\\ "
       << endl;
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
// Entry point
// -------------------------------------------------------------------------------

void fitgPb()
{
    ///////////////////////////////////////////////////////////////////////////////
    // adapt the input
    ///////////////////////////////////////////////////////////////////////////////
    // get the total uncorrelated error
    for(int j=0;j<nData;j++) {
        double unc = 0;
        for(int i=0;i<nUncorErr;i++) unc+=TMath::Power(uncorErr[i][j],2);
        totalUncorErr[j]=TMath::Sqrt(unc);
    }

    // set up the fluxes
    setFluxes(nData,fluxL,fluxL_noon_mn,fluxH,fluxH_noon_mn);

    // check if mid rapidity bin
    if (yMin<0) {
      gMidBin = kTRUE;
      nPhysPar = 1;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // print input cross section
    ///////////////////////////////////////////////////////////////////////////////
    bool printInputXS = true;
    if (printInputXS) {
      const char *nClass[nData] = {"0n0n", "Xn0n", "XnXn"};
      // get total correlated error
      for(int j=0;j<nData;j++) {
        double cor = 0;
	// skip lumi
        for(int i=1;i<nGloCorErr;i++) cor+=TMath::Power(gloCorErr[i][j],2);
	for(int i=0;i<nLocCorErr;i++) cor+=TMath::Power(locCorErr[i][j],2);
        totalCorErr[j]=TMath::Sqrt(cor);
      }
      // print in "latex" format
      for(int i=0;i<nData;i++) {
	cout << nClass[i] 
	     << Form(" $%4.3f \\pm %4.3f  \\pm %4.3f  \\pm %4.3f \\pm %4.3f$",
		     crossSection[i], statErr[i],
		     (crossSection[i]*totalUncorErr[i]),
		     (crossSection[i]*totalCorErr[i]),
		     (crossSection[i]*migErr[0][i])) 
	     << endl;
      }
    }

    ///////////////////////////////////////////////////////////////////////////////
    // perform the first three fits
    ///////////////////////////////////////////////////////////////////////////////

    // Fit A
    gFitType = 0; // only uncorrelated uncertainties
    double fitResultA[6]; // to store the output
    doFit(fitResultA);

    // Fit B
    gFitType = 1; // uncorrelated and globally correlated uncertainties
    double fitResultB[6]; // to store the output
    doFit(fitResultB);

    // Fit C
    gFitType = 2; // uncorrelated, globally and locally correlated uncertainties.
    double fitResultC[6]; // to store the output
    doFit(fitResultC);

    ///////////////////////////////////////////////////////////////////////////////
    // uncertainties from the flux
    ///////////////////////////////////////////////////////////////////////////////
    // Fit Fup
    setFluxes(nData,fluxL,fluxL_noon_up,fluxH,fluxH_noon_up);
    gFitType = 2; // uncorrelated, globally and locally correlated uncertainties.
    double fitResultFup[6]; // to store the output
    doFit(fitResultFup);
     // Fit Fdn
    setFluxes(nData,fluxL,fluxL_noon_dn,fluxH,fluxH_noon_dn);
    gFitType = 2; // uncorrelated, globally and locally correlated uncertainties.
    double fitResultFdn[6]; // to store the output
    doFit(fitResultFdn);

    // reset  the fluxes
    setFluxes(nData,fluxL,fluxL_noon_mn,fluxH,fluxH_noon_mn);
    
    ///////////////////////////////////////////////////////////////////////////////
    // modify the input and perform the rest of the fits
    ///////////////////////////////////////////////////////////////////////////////

    // save original data
    double crossSectionOri[nData];
    double statErrOri[nData];
    for(int i=0;i<nData;i++) {
        crossSectionOri[i]=crossSection[i];
        statErrOri[i]=statErr[i];
    }

    // loop over all anti-correlated sources
    double sign = 1;
    double fitResultM[6][nMigErr];
    double fitResultTmp[6];
    for(int j=0;j<nMigErr;j++) {    
      // modify data 
      for(int i=0;i<nData;i++) {
        crossSection[i] = crossSectionOri[i]*(1.0+sign*migErr[j][i]);
        statErr[i] = statErrOri[i]*(1.0+sign*migErr[j][i]);
      }
      gFitType = 2; // uncorrelated, globally and locally correlated uncertainties.
      doFit(fitResultTmp);
      for(int k=0;k<6;k++) fitResultM[k][j]=fitResultTmp[k];
    } // end anti-correlates sys

    // reset data
    for(int i=0;i<nData;i++) {
        crossSection[i]=crossSectionOri[i];
        statErr[i]=statErrOri[i];
    }

    ///////////////////////////////////////////////////////////////////////////////
    // print input data
    ///////////////////////////////////////////////////////////////////////////////
    printInputData();
    
    ///////////////////////////////////////////////////////////////////////////////
    // print fit results
    ///////////////////////////////////////////////////////////////////////////////

    printFitHeader();
    printFitResult("A",fitResultA);
    printFitResult("B",fitResultB);
    printFitResult("C",fitResultC);
    printFitResult("Fup",fitResultFup);
    printFitResult("Fdn",fitResultFdn);
    for(int j=0;j<nMigErr;j++) {    
      for(int k=0;k<6;k++) fitResultTmp[k]=fitResultM[k][j];
      printFitResult(Form("M"),fitResultTmp);
    }
    cout << "\\hline" << endl;

    ///////////////////////////////////////////////////////////////////////////////
    // obtain final results results
    ///////////////////////////////////////////////////////////////////////////////

    // kinematics
    double W[2] = {0,0};
    getW(yMin,yMax,W);
    
    // uncertainties
    double uU[2] = {0,0}; // uncorrelated
    double cU[2] = {0,0}; // correlated   
    double fU[2] = {0,0}; // flux fractions
    double mU[2] = {0,0};

    // uncertainties: uncorrelated
    //   get unc. in percentage of fit A and apply that to fitC
    for(int i=0;i<2;i++) {
      double frac = fitResultA[1+2*i]/fitResultA[2*i];
      uU[i] = fitResultC[2*i]*frac;
    }
    
    // uncertainties: correlated
    //   subtract corr unc from fitC
    for(int i=0;i<2;i++) {
      cU[i] = TMath::Sqrt((fitResultC[1+2*i]*fitResultC[1+2*i])
			  -(uU[i]*uU[i]));
    }

    // uncertainties: flux
    //   compare fitC with fitF and take the largest variation
    for(int i=0;i<2;i++) {
      double up = TMath::Abs(fitResultC[2*i]-fitResultFup[2*i]);
      double dn = TMath::Abs(fitResultC[2*i]-fitResultFdn[2*i]);
      fU[i] = TMath::Max(up,dn)/TMath::Sqrt(2);
    }
    
    // migrations
    for(int i=0;i<2;i++) {
      for(int j=0;j<nMigErr;j++) {
	double diff = TMath::Abs(fitResultM[2*i][j]-fitResultC[2*i])/TMath::Sqrt(2);
	if (diff>mU[i]) mU[i] = diff;
      }
    }

    
    ///////////////////////////////////////////////////////////////////////////////
    // print final results results
    ///////////////////////////////////////////////////////////////////////////////

    
    cout << endl << " rapidity bin $(" << yMin <<", "<<yMax<<")$" << endl;
    cout << "\\hline" << endl;
    cout << " $\\WgPb$ (GeV) & $\\sigma$ ($\\mu$b) & unc. ($\\mu$b) "
	 << " & corr. ($\\mu$b) & mig. ($\\mu$b)  & flux frac. ($\\mu$b)  \\\\" << endl;
    cout << "\\hline" << endl;
    for(int i=0;i<2;i++) {    
      cout << Form("%4.2f &",W[i]);
      cout << Form("%4.2f &",fitResultC[2*i]);
      cout << Form("%4.2f &",uU[i]);            
      cout << Form("%4.2f &",cU[i]);            
      cout << Form("%4.2f &",mU[i]);            
      cout << Form("%4.2f  \\\\",fU[i]);            
      cout << endl;
    }
    cout << "\\hline" << endl;
     
}
