
// bin definition
const double yMin = 3.5;
const double yMax = 4.0;

// fluxes in this bin
const int nData = 3; // 0n0n, Xn0n, XnXn
double fluxL[nData];
double fluxH[nData];
// noon
double fluxH_noon_mn[nData] = {1.34335*0.190331,0.5*1.34335*0.453302,1.34335*0.356367};
double fluxL_noon_mn[nData] = {203.026*0.87924,0.5*203.026*0.0895507,203.026*0.0312093};
double fluxH_noon_up[nData] = {1.34335*0.176846,0.5*1.34335*0.446474,1.34335*0.376679};
double fluxL_noon_up[nData] = {203.026*0.875745,0.5*203.026*0.0909746,203.026*0.0332802};
double fluxH_noon_dn[nData] = {1.34335*0.205025,0.5*1.34335*0.459334,1.34335*0.335641};
double fluxL_noon_dn[nData] = {203.026*0.882878,0.5*203.026*0.0879878,203.026*0.0291345};

// starlight
double fluxL_sl[nData] = {178.28,  0.5*18.312, 6.5199};
double fluxH_sl[nData] = {0.2149,  0.5*0.5266, 0.425};

// cross section
double totalCrossSection = (1.615+1.938)/2.0; // mb
double crossSection[nData] = { // mb
    1.5900,  0.1008 , 0.0791
};

// uncertainties
double statErr[nData] = { 0.0490,  0.0085, 0.0114}; // mb

const int nUncorErr = 2;
const double uncorErr[nUncorErr][nData] = {
  0.0054,  0.0053, 0.009, //signal extraction
  0.0037,  0.0053, 0.0223 //incoherent fraction
};
const char *uncorErrNames[nUncorErr] = {
  "Signal extraction","Incoherent fraction"
};

const int nGloCorErr = 3;
const double gloCorErr[nGloCorErr][nData] = {
  0.02, 0.02, 0.02, // flux
  0.025,  0.025, 0.025, //lumi vdM
  0.006,  0.006, 0.006 //BR
};
const char *gloCorErrNames[nGloCorErr] = {
  "Photon flux","Luminosity","Branching ratio"
};

const int nLocCorErr = 7;
const double locCorErr[nLocCorErr][nData] = {
    0.03,  0.03, 0.03, // trk
    0.062,  0.062, 0.062, // trig
    0.01,  0.01, 0.01, // matching
    0.001,  0.001, 0.001, // pt-shape
    0.007,  0.007, 0.007, // fD
    0.0,  0.0114, 0.0595, // EMD pile-up
    0.002,  0.002, 0.002 // pile-up
};
const char *locCorErrNames[nLocCorErr] = {
  "Tracking","Trigger","Matching",
  "Coherent shape","Feed down","EMD pile-up","Pile-up"
};

/*
const int nMigErr = 6;
const double migErr[nMigErr][nData] = {
  0.0012,-0.008,-0.0012, // pu upup
  -0.0001,-0.0094,0.0005, // pu updn
  0.0001,0.0094,-0.0006, // pu dnup
  0.0016,-0.0041,-0.0355, // eff upup
  0.0006,-0.0322,0.0004, // eff updn
  -0.007, 0.0329,0.0003 // eff dnup
};


const char *migErrNames[nMigErr] = {
  "pu upup", "pu updn", "pu dnup", "eff upup", "eff updn", "eff dnup"
};
*/

const int nMigErr = 1;
const double migErr[nMigErr][nData] = {
  -0.0016,0.0329, 0.0355
};


const char *migErrNames[nMigErr] = {
   "Migrations"
};

