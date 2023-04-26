
// bin definition
const double yMin = 2.5;
const double yMax = 3.0;

// fluxes in this bin
const int nData = 3; // 0n0n, Xn0n, XnXn
double fluxL[nData];
double fluxH[nData];
// noon
double fluxH_noon_mn[nData] = {11.6585*0.335105,0.5*11.6585*0.426104,11.6585*0.238791};
double fluxL_noon_mn[nData] = {171.994*0.857367,0.5*171.994*0.105767,171.994*0.0368668};
double fluxH_noon_up[nData] = {11.6585*0.320662,0.5*11.6585*0.425762,11.6585*0.253577};
double fluxL_noon_up[nData] = {171.994*0.853239,0.5*171.994*0.107448,171.994*0.0393131};
double fluxH_noon_dn[nData] = {11.6585*0.350503,0.5*11.6585*0.425647,11.6585*0.22385};
double fluxL_noon_dn[nData] = {171.994*0.861663,0.5*171.994*0.103921,171.994*0.0344159};

// starlight
double fluxL_sl[nData] = {147.23, 0.5*18.311, 6.5225};
double fluxH_sl[nData] = {3.6998, 0.5*4.8253, 2.7674};

// cross section
double totalCrossSection = (3.018+3.531)/2.0; // mb
double crossSection[nData] = { // mb
    2.6680,  0.2418 , 0.2555
};

// uncertainties
double statErr[nData] = { 0.0763,  0.0213, 0.0237}; // mb

const int nUncorErr = 2;
const double uncorErr[nUncorErr][nData] = {
  0.0020,  0.0126, 0.0081, //signal extraction
  0.0038,  0.0061, 0.0159 //incoherent fraction
};
const char *uncorErrNames[nUncorErr] = {
  "Signal extraction","Incoherent fraction"
};

const int nGloCorErr = 3;
const double gloCorErr[nGloCorErr][nData] = {
  0.02,  0.02, 0.02, // flux
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
  0.0012,-0.0054,-0.0012, // pu upup
  -0.0001,-0.0068,0.0001, // pu updn
  0.0001,0.0068,-0.0002, // pu dnup
  0.0032,0.0007,-0.0355, // eff upup
  0.0001,-0.0372,0.0004, // eff updn
  -0.002, 0.0379,0.0003 // eff dnup
};


const char *migErrNames[nMigErr] = {
  "pu upup", "pu updn", "pu dnup", "eff upup", "eff updn", "eff dnup"
};
*/

const int nMigErr = 1;
const double migErr[nMigErr][nData] = {
  -0.0032,0.0379, 0.0355
};


const char *migErrNames[nMigErr] = {
   "Migrations"
};

