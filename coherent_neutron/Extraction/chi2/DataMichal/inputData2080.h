
// bin definition
const double yMin = 0.2;
const double yMax = 0.80;

// fluxes in this bin
const int nData = 3; // 0n0n, Xn0n, XnXn
double fluxL[nData];
double fluxH[nData];
// noon
double fluxH_noon_mn[nData] = {71.2488*0.678053,71.2488*0.235154,71.2488*0.0867926};
double fluxL_noon_mn[nData] = {102.088*0.762836,102.088*0.175178,102.088*0.0619866};
double fluxH_noon_up[nData] = {71.2488*0.669012,71.2488*0.238461,71.2488*0.092527};
double fluxL_noon_up[nData] = {102.088*0.756026,102.088*0.177878,102.088*0.0660961};
double fluxH_noon_dn[nData] = {71.2488*0.687487,71.2488*0.231468,71.2488*0.0810441};
double fluxL_noon_dn[nData] = {102.088*0.769928,102.088*0.172203,102.088*0.0578688};

// starlight
double fluxL_sl[nData] = {77.709, 17.962, 6.4902};
double fluxH_sl[nData] = {47.873, 16.789, 6.3229};

// cross section
double crossSection[nData] = { // mb
    2.90, 0.8, 0.30
};

// uncertainties
double statErr[nData] = { 0.07,0.04,0.029 }; // mb

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
  -0.014,0.004, 0.01, // eff upup
  0.015,-0.004,-0.011, // eff dndn
  -0.036,0.031,0.005, // pu upup
  0.035,-0.03,-0.005, // pu dndn
  0.003, -0.002, -0.002, // pu updn
  -0.003, 0.003, 0.0 // pu dnup  
};


const char *migErrNames[nMigErr] = {
   "eff upup", "eff dndn", "pu upup", "pu dndn", "pu updn", "pu dnup"
};
*/

const int nMigErr = 1;
const double migErr[nMigErr][nData] = {
  -0.036,0.031, 0.011
};


const char *migErrNames[nMigErr] = {
   "Migrations"
};




