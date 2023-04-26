
// bin definition
const double yMin = -0.2;
const double yMax = 0.2;

// fluxes in this bin
const int nData = 3; // 0n0n, Xn0n, XnXn
double fluxL[nData];
double fluxH[nData];
// noon
double fluxH_noon_mn[nData] = {86.603*0.725891,86.603*0.201711,86.603*0.0723976};
double fluxL_noon_mn[nData] = {86.603*0.725891,86.603*0.201711,86.603*0.0723976};
double fluxH_noon_up[nData] = {86.603*0.71808,86.603*0.204728,86.603*0.0771925};
double fluxL_noon_up[nData] = {86.603*0.71808,86.603*0.204728,86.603*0.0771925};
double fluxH_noon_dn[nData] = {86.603*0.734033,86.603*0.198374,86.603*0.0675925};
double fluxL_noon_dn[nData] = {86.603*0.734033,86.603*0.198374,86.603*0.0675925};


// SL
double fluxL_sl[nData] = {62.374, 17.565, 6.4383};
double fluxH_sl[nData] = {62.374, 17.565, 6.4383};

// cross section
double crossSection[nData] = { // mb
    3.13, 0.73, 0.25
};

// uncertainties
double statErr[nData] = { 0.09,0.05,0.024 }; // mb

const int nUncorErr = 2;
const double uncorErr[nUncorErr][nData] = {
  0.015, 0.015, 0.015, //signal extraction
  0.001, 0.015, 0.013 //incoherent fraction
};
const char *uncorErrNames[nUncorErr] = {
  "Signal extraction","Incoherent fraction"
};

const int nGloCorErr = 3;
const double gloCorErr[nGloCorErr][nData] = {
  0.02, 0.02, 0.02, // flux
  0.025, 0.025, 0.025, //lumi vdM
  0.005, 0.005, 0.005 //BR
};
const char *gloCorErrNames[nGloCorErr] = {
  "Photon flux","Luminosity","Branching ratio"
};

const int nLocCorErr = 8;
const double locCorErr[nLocCorErr][nData] = {
  0.001, 0.008, 0.006, // Incoherent b-slope
  0.015, 0.015, 0.015, // trigger live time
  0.028, 0.028, 0.028, // its-tpc matching
  0.007, 0.007, 0.007, // TOF trig
  0.01, 0.01, 0.01, // SPD trigger
  0.006, 0.006, 0.006, // fD
  0.0, 0.032, 0.035, // EMD pile-up
  0.03, 0.03, 0.03 // pile-up
};
const char *locCorErrNames[nLocCorErr] = {
  "Coherent shape","Trigger live time","ITS-TPC matching","TOF trigger",
  "SPD trigger","Feed down","EMD pile-up","Pile-up"
};

/*
const int nMigErr = 6;
const double migErr[nMigErr][nData] = {
  -0.013,0.005, 0.009, // eff upup
  0.013,-0.004,-0.009, // eff dndn
  -0.039,0.034,0.004, // pu upup
  0.038,-0.033,-0.005, // pu dndn
  0.003, -0.003, -0.001, // pu updn
  -0.004, 0.003, 0.0 // pu dnup  
};


const char *migErrNames[nMigErr] = {
  "eff upup", "eff dndn", "pu upup", "pu dndn", "pu updn", "pu dnup"
};

*/

const int nMigErr = 1;
const double migErr[nMigErr][nData] = {
  -0.039,0.034, 0.009
};


const char *migErrNames[nMigErr] = {
   "Migrations"
};

