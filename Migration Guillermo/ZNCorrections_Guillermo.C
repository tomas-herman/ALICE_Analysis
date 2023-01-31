
//
// this programs computes the migrations in ZN classes due to pile-up and
// efficiency in ZN
//
// The equation to be solved is \vec{x}=(\matrix{A})^-1\vector{b}
// here \vector{b} are the cross sections before ZN corrections,
// \vector{x} are the cross sections after the ZN corrections
// the \matrix{A} contains the corrections for each term
//
// indices 1,2, 3 -> 0n0n, Xn0n+0nXn, XnXn

// c++ headers
#include <iostream>
#include <fstream>
using namespace std;

// root headers
#include <TMatrixD.h>

Double_t xs00_dh;
Double_t xsx0_dh;  
Double_t xsxx_dh; 

/////////////////////////////////////////////////////////////////////////////

Double_t Get_11(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  return (1.0-pA*(1.0-pC) - pC*(1.0-pA) - pA*pC);
}

/////////////////////////////////////////////////////////////////////////////

Double_t Get_12(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  Double_t c_x0 = 0.5*(1.0-eA)*(1.0-pA)*(1.0-pC);
  Double_t c_0x = 0.5*(1.0-eC)*(1.0-pA)*(1.0-pC);  
  return (c_x0+c_0x);
}

/////////////////////////////////////////////////////////////////////////////

Double_t Get_13(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  Double_t c_xx = (1.0-eA)*(1.0-eC)*(1.0-pA)*(1.0-pC);
  return (c_xx);
}

/////////////////////////////////////////////////////////////////////////////

Double_t Get_21(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  Double_t c_0x = pC*(1-pA);
  Double_t c_x0 = pA*(1-pC);  
  return (c_x0+c_0x);
}

/////////////////////////////////////////////////////////////////////////////

Double_t Get_22(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  //  Double_t c_x0 = 0.5*(1.0-(1.0-eA)*(1.0-pA)*((1.0-pC)+pC)
  Double_t c_x0 = 0.5*(1.0-(1.0-eA)*(1.0-pA)*(1.0-pC)
		       -(eA*pC+(1.0-eA)*pA*pC));
  //  Double_t c_0x = 0.5*(1.0-(1.0-eC)*(1.0-pC)*((1.0-pA)+pA)
  Double_t c_0x = 0.5*(1.0-(1.0-eC)*(1.0-pC)*(1.0-pA)
		       -(eC*pA+(1.0-eC)*pA*pC));  
  return (c_x0+c_0x);
}

/////////////////////////////////////////////////////////////////////////////

Double_t Get_23(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  Double_t c_x0 = (1.0-eC)*(1.0-pC)*(eA+(1.0-eA)*pA);
  Double_t c_0x = (1.0-eA)*(1.0-pA)*(eC+(1.0-eC)*pC);  
  return (c_x0+c_0x);
}

/////////////////////////////////////////////////////////////////////////////

Double_t Get_31(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  return (pA*pC);
}

/////////////////////////////////////////////////////////////////////////////

Double_t Get_32(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  Double_t c_x0 = 0.5*(eC*pA+(1-eC)*pA*pC);
  Double_t c_0x = 0.5*(eA*pC+(1-eA)*pA*pC);  
  return (c_x0+c_0x);
}



/////////////////////////////////////////////////////////////////////////////

Double_t Get_33(Double_t eA, Double_t eC, Double_t pA, Double_t pC)
{
  Double_t c1 = (1.0-eA)*(1.0-pA)*(eC+(1-eC)*pC);
  Double_t c2 = (1.0-eC)*(1.0-pC)*(eA+(1-eA)*pA);
  Double_t c3 = (1.0-eA)*(1.0-eC)*(1.0-pA)*(1.0-pC);    
  return (1.0-(c1+c2+c3));
}

/////////////////////////////////////////////////////////////////////////////

void ZNCorrections(Double_t eA, // efficiency ZNA
		   Double_t eC, // efficiency ZNC
		   Double_t pA, // pile-up ZNA
		   Double_t pC, // pile-up ZNC
		   Double_t XS_00, // cross section in 0n0n
		   Double_t XS_X0, // cross section in Xn0n+0nXn
		   Double_t XS_XX, // cross section in XnXn		   
		   Double_t *XS_corr)
{
  // print out input
  cout << " efficiency: ZNA = " << eA << ", ZNC = " << eC << endl;
  cout << " pile-up: ZNA = " << pA << ", ZNC = " << pC << endl;
  cout << " XS before corrections: 0n0n = " << XS_00  << endl;
  cout << " XS before corrections: 0nXn+Xn0n = " << XS_X0  << endl;
  cout << " XS before corrections: XnXn = " << XS_XX  << endl << endl;
  
  // compute the different coefficients
  Double_t a11 = Get_11(eA,eC,pA,pC);
  Double_t a12 = Get_12(eA,eC,pA,pC);
  Double_t a13 = Get_13(eA,eC,pA,pC);  
  Double_t a21 = Get_21(eA,eC,pA,pC);
  Double_t a22 = Get_22(eA,eC,pA,pC);
  Double_t a23 = Get_23(eA,eC,pA,pC);  
  Double_t a31 = Get_31(eA,eC,pA,pC);
  Double_t a32 = Get_32(eA,eC,pA,pC);
  Double_t a33 = Get_33(eA,eC,pA,pC);

  // check
  Double_t d00 = a11 + a21 + a31;
  Double_t dxn = a22 + a12 + a32;
  Double_t dxx = a33 + a13 + a23;
  cout << " Intermediate check: " << d00 << ", " << dxn << ", "<< dxx << ", "<< endl;

  // fill matrix (row, col)
  TMatrixD A(3,3);
  A(0,0)=a11;A(0,1)=a12;A(0,2)=a13;
  A(1,0)=a21;A(1,1)=a22;A(1,2)=a23;
  A(2,0)=a31;A(2,1)=a32;A(2,2)=a33;
  // Define inverse
  TMatrixD Ainv(A);
  Ainv.Invert();
  // test inverse
  cout << endl;
  (A*Ainv).Print();
  cout << endl;

  // define inhomogeneous term
  TMatrixD b(3,1);
  b(0,0) = XS_00; b(1,0) = XS_X0; b(2,0) = XS_XX;

  // solve the equation
  TMatrixD x(Ainv*b);

  // final check
  cout << " sum input " << (XS_00+XS_X0+XS_XX) << " sum output " << (x(0,0)+x(1,0)+x(2,0)) << endl;

  // final printout
  cout << endl << endl << " Corrected cross sections: " << endl;
  cout << " xs_00 = " << x(0,0) << " xs_X0 = " << x(1,0) << " xs_XX = " << x(2,0) << endl;
  cout << " Percentual change (after/before) " << endl;
  cout << "0n0n => " << x(0,0)/XS_00 <<" Xn0n+0nXn => "
       << x(1,0)/XS_X0 << " XnXn => " << x(2,0)/XS_XX << endl;
  cout << " Percentual change (after/dh) " << endl;
  cout << "0n0n => " << x(0,0)/xs00_dh <<" Xn0n+0nXn => "
       << x(1,0)/xsx0_dh << " XnXn => " << x(2,0)/xsxx_dh << endl;

  // store found cross sections
  XS_corr[0] = x(0,0);
  XS_corr[1] = x(1,0);
  XS_corr[2] = x(2,0);  
}

void GetSys(Int_t opt)
{
  // cross sections from David
  Double_t xs00=417.6; // |y|<0.2
  Double_t xsx0=77.5;  
  Double_t xsxx=13.3; 
  if (opt==1) {xs00=420.1;xsx0=76.0;xsxx=13.3;} // 0.2<|y|<0.45
  if (opt==2) {xs00=423.0;xsx0=77.6;xsxx=13.9;}  // 0.45<|y|<0.08
    
  Double_t pA = 0.031;
  Double_t pAunc = 0.003;
  Double_t pC = 0.031;
  Double_t pCunc = 0.003;
  Double_t eA = 0.93;
  Double_t eAunc = 0.01;
  Double_t eC = 0.93;
  Double_t eCunc = 0.01;
  // reserve space
  Double_t *xs_std = new Double_t[3];
  xs_std[0]=xs_std[1]=xs_std[2]=0.0;

  Double_t *xs_pu_up = new Double_t[3];
  xs_pu_up[0]=xs_pu_up[1]=xs_pu_up[2]=0.0;
 
  Double_t *xs_pu_dn = new Double_t[3];
  xs_pu_dn[0]=xs_pu_dn[1]=xs_pu_dn[2]=0.0;

  Double_t *xs_eff_up = new Double_t[3];
  xs_eff_up[0]=xs_eff_up[1]=xs_eff_up[2]=0.0;
 
  Double_t *xs_eff_dn = new Double_t[3];
  xs_eff_dn[0]=xs_eff_dn[1]=xs_eff_dn[2]=0.0;

  // compute systematics
  ZNCorrections(eA,eC,pA,pC,xs00,xsx0,xsxx,xs_std);
  ZNCorrections(eA+eAunc,eC+eCunc,pA,pC,xs00,xsx0,xsxx,xs_eff_up);
  ZNCorrections(eA-eAunc,eC-eCunc,pA,pC,xs00,xsx0,xsxx,xs_eff_dn);  
  ZNCorrections(eA,eC,pA+pAunc,pC+pCunc,xs00,xsx0,xsxx,xs_pu_up);
  ZNCorrections(eA,eC,pA-pAunc,pC-pCunc,xs00,xsx0,xsxx,xs_pu_dn);  

  // print out
  cout << endl << " =======================std/sys======================= "  << endl << endl;
  cout << " eff_up: " << (xs_std[0]/xs_eff_up[0]) << ",  "
       << (xs_std[1]/xs_eff_up[1]) << ",  "
       << (xs_std[2]/xs_eff_up[2]) << endl;
  cout << " eff_dn: " << (xs_std[0]/xs_eff_dn[0]) << ",  "
       << (xs_std[1]/xs_eff_dn[1]) << ",  "
       << (xs_std[2]/xs_eff_dn[2]) << endl;
  cout << " pu_up: " << (xs_std[0]/xs_pu_up[0]) << ",  "
       << (xs_std[1]/xs_pu_up[1]) << ",  "
       << (xs_std[2]/xs_pu_up[2]) << endl;
  cout << " pu_dn: " << (xs_std[0]/xs_pu_dn[0]) << ",  "
       << (xs_std[1]/xs_pu_dn[1]) << ",  "
       << (xs_std[2]/xs_pu_dn[2]) << endl;
  cout << endl << " ===================================================== " << endl;  

  // clean up
  delete [] xs_std;
  delete [] xs_pu_up;
  delete [] xs_pu_dn;
  delete [] xs_eff_up;
  delete [] xs_eff_dn;  

}

void GetSysNew(Int_t opt,Int_t unc)
{
  // cross sections from David
  xs00_dh=417.6; // |y|<0.2
  xsx0_dh=77.5;  
  xsxx_dh=13.3; 
  if (opt==1) {xs00_dh=420.1;xsx0_dh=76.0;xsxx_dh=13.3;} // 0.2<|y|<0.45
  if (opt==2) {xs00_dh=423.0;xsx0_dh=77.6;xsxx_dh=13.9;}  // 0.45<|y|<0.08

  // new: correct cross sections for EMD and signals in AD/V0
  // Double_t cf_sEDM = 1.0-0.29;
  // Double_t cf_mEDM = 1.0-0.48;
  Double_t cf_sEDM = 1.0-0.26;
  Double_t cf_mEDM = 1.0-0.43;
  Double_t cfe_sEDM = 0.04;  
  Double_t cfe_mEDM = 0.05;
  
  Double_t xs00=xs00_dh; // |y|<0.2
  Double_t xsx0=xsx0_dh;  
  Double_t xsxx=xsxx_dh;
  Double_t allB = xs00+xsx0+xsxx;  // before
  if (unc==0) {
    xsx0 = xsx0/cf_sEDM; 
    xsxx = xsxx/cf_mEDM;
  } else if  (unc==1) {
    xsx0 = xsx0/(cf_sEDM+cfe_sEDM); 
    xsxx = xsxx/(cf_mEDM); // 
  } else if  (unc==2) {
    xsx0 = xsx0/(cf_sEDM-cfe_sEDM); 
    xsxx = xsxx/(cf_mEDM);
  } else if  (unc==3) {
    xsx0 = xsx0/(cf_sEDM); 
    xsxx = xsxx/(cf_mEDM+cfe_mEDM); 
  } else if   (unc==4) {
    xsx0 = xsx0/(cf_sEDM); 
    xsxx = xsxx/(cf_mEDM-cfe_mEDM); 
  }

  Double_t allA = xs00+xsx0+xsxx; // after
  // new pile-up 0.031 -> 0.026 for ZNA and ZNC
  Double_t pA = 0.026;
  Double_t pAunc = 0.003;
  Double_t pC = 0.026;
  Double_t pCunc = 0.003;
  Double_t eA = 0.93;
  Double_t eAunc = 0.01;
  Double_t eC = 0.93;
  Double_t eCunc = 0.01;
  // reserve space
  Double_t *xs_std = new Double_t[3];
  xs_std[0]=xs_std[1]=xs_std[2]=0.0;

  Double_t *xs_pu_up = new Double_t[3];
  xs_pu_up[0]=xs_pu_up[1]=xs_pu_up[2]=0.0;
 
  Double_t *xs_pu_dn = new Double_t[3];
  xs_pu_dn[0]=xs_pu_dn[1]=xs_pu_dn[2]=0.0;

  Double_t *xs_eff_up = new Double_t[3];
  xs_eff_up[0]=xs_eff_up[1]=xs_eff_up[2]=0.0;
 
  Double_t *xs_eff_dn = new Double_t[3];
  xs_eff_dn[0]=xs_eff_dn[1]=xs_eff_dn[2]=0.0;

  // compute systematics
  ZNCorrections(eA,eC,pA,pC,xs00,xsx0,xsxx,xs_std);
  /*
  ZNCorrections(eA+eAunc,eC+eCunc,pA,pC,xs00,xsx0,xsxx,xs_eff_up);
  ZNCorrections(eA-eAunc,eC-eCunc,pA,pC,xs00,xsx0,xsxx,xs_eff_dn);  
  ZNCorrections(eA,eC,pA+pAunc,pC+pCunc,xs00,xsx0,xsxx,xs_pu_up);
  ZNCorrections(eA,eC,pA-pAunc,pC-pCunc,xs00,xsx0,xsxx,xs_pu_dn);  

  // print out
  cout << endl << " =======================std/sys======================= "  << endl << endl;
  cout << " eff_up: " << (xs_std[0]/xs_eff_up[0]) << ",  "
       << (xs_std[1]/xs_eff_up[1]) << ",  "
       << (xs_std[2]/xs_eff_up[2]) << endl;
  cout << " eff_dn: " << (xs_std[0]/xs_eff_dn[0]) << ",  "
       << (xs_std[1]/xs_eff_dn[1]) << ",  "
       << (xs_std[2]/xs_eff_dn[2]) << endl;
  cout << " pu_up: " << (xs_std[0]/xs_pu_up[0]) << ",  "
       << (xs_std[1]/xs_pu_up[1]) << ",  "
       << (xs_std[2]/xs_pu_up[2]) << endl;
  cout << " pu_dn: " << (xs_std[0]/xs_pu_dn[0]) << ",  "
       << (xs_std[1]/xs_pu_dn[1]) << ",  "
       << (xs_std[2]/xs_pu_dn[2]) << endl;
  cout << endl << " ===================================================== " << endl;  
  */
  cout << " allA = " << allA << " allA/allB = " << (allA/allB) << endl;

  // clean up
  delete [] xs_std;
  delete [] xs_pu_up;
  delete [] xs_pu_dn;
  delete [] xs_eff_up;
  delete [] xs_eff_dn;  

}

void GetXeXe(Int_t opt,Int_t inc)
{
  if (inc==0) {
    /*
  xs00_dh=1868.06; //±20 (incoherent subtraction)
  xsx0_dh=131.35; // ±4 (incoherent subtraction)
  xsxx_dh=10.83; // ±0.5 (incoherent subtraction)
    */
  xs00_dh=1715.55; 
  xsx0_dh=151.23; 
  xsxx_dh=16.12;    
  } 
  if (inc==1) {
    /*
  xs00_dh=1888.06; //±20 (incoherent subtraction)
  xsx0_dh=131.35; // ±4 (incoherent subtraction)
  xsxx_dh=10.83; // ±0.5 (incoherent subtraction)
    */
    xs00_dh=1715.55; 
  xsx0_dh=139.2; 
  xsxx_dh=14.84;    
  }   if (inc==2) {
    /*
  xs00_dh=1848.06; //±20 (incoherent subtraction)
  xsx0_dh=131.35; // ±4 (incoherent subtraction)
  xsxx_dh=10.83; // ±0.5 (incoherent subtraction)
    */
    xs00_dh=1715.55; 
  xsx0_dh=163.2; 
  xsxx_dh=17.4;    
  } 
  if (inc==3) {
    /*
  xs00_dh=1868.06; //±20 (incoherent subtraction)
  xsx0_dh=135.35; // ±4 (incoherent subtraction)
  xsxx_dh=11.33; // ±0.5 (incoherent subtraction)
    */
    xs00_dh=1690.00; 
  xsx0_dh=151.2; 
  xsxx_dh=16.1;    
    
  } 
  if (inc==4) {
    /*
  xs00_dh=1868.06; //±20 (incoherent subtraction)
  xsx0_dh=127.35; // ±4 (incoherent subtraction)
  xsxx_dh=10.33; // ±0.5 (incoherent subtraction)
    */
    xs00_dh=1741.13; 
  xsx0_dh=151.2; 
  xsxx_dh=16.12;    
  } 


  
  Double_t xs00=xs00_dh; // |y|<0.2
  Double_t xsx0=xsx0_dh;  
  Double_t xsxx=xsxx_dh;
  Double_t allB = xs00+xsx0+xsxx;  // before

  Double_t allA = xs00+xsx0+xsxx; // after
  Double_t pA = 0.0047;
  Double_t pAunc = 0.0002;
  Double_t pC = 0.0044;
  Double_t pCunc = 0.0002;
  Double_t eA = 0.91;
  Double_t eAunc = 0.02;
  Double_t eC = 0.92;
  Double_t eCunc = 0.02;
  // reserve space
  Double_t *xs_std = new Double_t[3];
  xs_std[0]=xs_std[1]=xs_std[2]=0.0;

  Double_t *xs_pu_up = new Double_t[3];
  xs_pu_up[0]=xs_pu_up[1]=xs_pu_up[2]=0.0;
 
  Double_t *xs_pu_dn = new Double_t[3];
  xs_pu_dn[0]=xs_pu_dn[1]=xs_pu_dn[2]=0.0;

  Double_t *xs_eff_up = new Double_t[3];
  xs_eff_up[0]=xs_eff_up[1]=xs_eff_up[2]=0.0;
 
  Double_t *xs_eff_dn = new Double_t[3];
  xs_eff_dn[0]=xs_eff_dn[1]=xs_eff_dn[2]=0.0;

  // compute systematics
  if (opt == 0) ZNCorrections(eA,eC,pA,pC,xs00,xsx0,xsxx,xs_std);
  
  if (opt == 1)ZNCorrections(eA+eAunc,eC+eCunc,pA,pC,xs00,xsx0,xsxx,xs_std);
  if (opt == 2)ZNCorrections(eA-eAunc,eC-eCunc,pA,pC,xs00,xsx0,xsxx,xs_std);  
  if (opt == 3)ZNCorrections(eA,eC,pA+pAunc,pC+pCunc,xs00,xsx0,xsxx,xs_std);
  if (opt == 4)ZNCorrections(eA,eC,pA-pAunc,pC-pCunc,xs00,xsx0,xsxx,xs_std);  
  /*
  // print out
  cout << endl << " =======================std/sys======================= "  << endl << endl;
  cout << " eff_up: " << (xs_std[0]/xs_eff_up[0]) << ",  "
       << (xs_std[1]/xs_eff_up[1]) << ",  "
       << (xs_std[2]/xs_eff_up[2]) << endl;
  cout << " eff_dn: " << (xs_std[0]/xs_eff_dn[0]) << ",  "
       << (xs_std[1]/xs_eff_dn[1]) << ",  "
       << (xs_std[2]/xs_eff_dn[2]) << endl;
  cout << " pu_up: " << (xs_std[0]/xs_pu_up[0]) << ",  "
       << (xs_std[1]/xs_pu_up[1]) << ",  "
       << (xs_std[2]/xs_pu_up[2]) << endl;
  cout << " pu_dn: " << (xs_std[0]/xs_pu_dn[0]) << ",  "
       << (xs_std[1]/xs_pu_dn[1]) << ",  "
       << (xs_std[2]/xs_pu_dn[2]) << endl;
  cout << endl << " ===================================================== " << endl;  
  cout << " allA = " << allA << " allA/allB = " << (allA/allB) << endl;
  */
  cout << " xs_std = " << xs_std[0] << ", "
       << xs_std[1] << ", " << xs_std[2] << ", "<<endl;
  Double_t sum_std = xs_std[0]+xs_std[1]+xs_std[2];
  Double_t f0 = xs_std[0]/sum_std;
  Double_t f1 = xs_std[1]/sum_std;
  Double_t f2 = xs_std[2]/sum_std;

  cout << sum_std << " " << f0 << " " << f1 << " " << f2 << endl;
  // clean up
  delete [] xs_std;
  delete [] xs_pu_up;
  delete [] xs_pu_dn;
  delete [] xs_eff_up;
  delete [] xs_eff_dn;  

}
