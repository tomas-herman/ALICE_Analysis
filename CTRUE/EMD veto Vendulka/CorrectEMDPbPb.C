// C++ headers
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>      // std::setprecision


void CorrectEMDPbPb(Int_t opt_n =2, Int_t opt_pu=0,Int_t opt_emd=0)
{
  // input data
  Double_t nXn0n_V0A=107;     //V0A Decision == 1
  Double_t nXn0n_noV0A=1154;  //V0A Decision == 0

  Double_t nXnXn_V0A=111;     //V0A Decision == 1
  Double_t nXnXn_noV0A=802;   //V0A Decision == 0
  
 if (opt_n == 1){nXn0n_V0A+=13; nXn0n_noV0A+=43; nXnXn_V0A+=43; nXnXn_noV0A+=35;} 
 if (opt_n == 2){nXn0n_V0A-=13; nXn0n_noV0A-=43; nXnXn_V0A-=43; nXnXn_noV0A-=35;}


  // Pile-up inefficiency (different pair) factors for V0A: 
  Double_t puV0A = 0.0385;

  if(opt_pu == 1) {puV0A+=0.0015;}
  if(opt_pu == 2) {puV0A-=0.0015;}

  if(opt_pu == 1) {cout <<"The case for upper limit of pile-up factor" << endl;}
  if(opt_pu == 2) {cout <<"The case for lower limit of pile-up factor" << endl;}

  // EMD inefficiency (same pair) factors for V0A: 
  Double_t Xn0n_vetoV0A = 0.38;
  Double_t XnXn_vetoV0A = 0.35;

  if (opt_emd == 1) {Xn0n_vetoV0A+=0.03;XnXn_vetoV0A+=0.05;}
  if (opt_emd == 2) {Xn0n_vetoV0A-=0.03;XnXn_vetoV0A-=0.05;}
  
  if(opt_emd == 1) {cout <<"The case for upper limit of EMD veto factor" << endl;}
  if(opt_emd == 2) {cout <<"The case for lower limit of EMD veto factor" << endl;}
  
  // compute and print
  Double_t F_Xn0n = nXn0n_noV0A/(1-puV0A) / ( nXn0n_V0A/(1+puV0A) / (1-Xn0n_vetoV0A) + nXn0n_noV0A/(1-puV0A) );
  Double_t F_XnXn = nXnXn_noV0A/(1-puV0A) / ( nXnXn_V0A/(1+puV0A) / (1-XnXn_vetoV0A) + nXnXn_noV0A/(1-puV0A) );
  
  cout << "F_Xn0n: " << F_Xn0n << endl;
  cout << "F_XnXn: " << F_XnXn << endl;

}




