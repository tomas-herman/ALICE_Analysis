
// bin definition
const double yMin = 3.0;
const double yMax = 3.5;

// fluxes in this bin
const int nData = 3; // 0n0n, Xn0n, XnXn
double fluxL[nData];
double fluxH[nData];
// noon
double fluxH_noon_mn[nData] = {4.78854*0.255502,0.5*4.78854*0.448296,4.78854*0.296202};
double fluxL_noon_mn[nData] = {187.521*0.869203,0.5*187.521*0.096995,187.521*0.033802};
double fluxH_noon_up[nData] = {4.78854*0.241193,0.5*4.78854*0.44493,4.78854*0.313877};
double fluxL_noon_up[nData] = {187.521*0.865418,0.5*187.521*0.0985375,187.521*0.0360449};
double fluxH_noon_dn[nData] = {4.78854*0.270912,0.5*4.78854*0.450825,4.78854*0.278263};
double fluxL_noon_dn[nData] = {187.521*0.873143,0.5*187.521*0.095302,187.521*0.0315548};

// starlight
double fluxL_sl[nData] = {162.75,  0.5*18.313, 6.5209};
double fluxH_sl[nData] = {1.0973,  0.5*1.9936, 1.3546};

// cross section
double totalCrossSection = (2.377+2.831)/2.0; // mb
double crossSection[nData] = { // mb
    2.3220,  0.1716 , 0.1613
};

// uncertainties
double statErr[nData] = { 0.0334,  0.0099, 0.0114}; // mb

const int nUncorErr = 2;
const double uncorErr[nUncorErr][nData] = {
  0.0010, 0.0063, 0.0067, //signal extraction
  0.0041, 0.0090, 0.0332 //incoherent fraction
};
const char *uncorErrNames[nUncorErr] = {
  "Signal extraction","Incoherent fraction"
};

const int nGloCorErr = 3;
const double gloCorErr[nGloCorErr][nData] = {
  0.02,  0.02, 0.02, // flux
  0.025,  0.025, 0.025, //lumi vdM
  0.006, 0.006, 0.006 //BR
};
const char *gloCorErrNames[nGloCorErr] = {
  "Photon flux","Luminosity","Branching ratio"
};

const int nLocCorErr = 7;
const double locCorErr[nLocCorErr][nData] = {
    0.03, 0.03, 0.03, // trk
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
  0.0012,-0.0068,-0.0011, // pu upup
  -0.0001,-0.0082,0.0003, // pu updn
  0.0001,0.0082,-0.0003, // pu dnup
  0.0022,-0.0013,-0.0355, // eff upup
  0.0004,-0.0351,0.0004, // eff updn
  -0.006, 0.0357,0.0003 // eff dnup
};


const char *migErrNames[nMigErr] = {
  "pu upup", "pu updn", "pu dnup", "eff upup", "eff updn", "eff dnup"
};
*/
const int nMigErr = 1;
const double migErr[nMigErr][nData] = {
  -0.0022,0.0357, 0.0355
};


const char *migErrNames[nMigErr] = {
   "Migrations"
};

