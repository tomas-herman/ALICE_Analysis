////////////////////////////////////////
// Define global variabels
////////////////////////////////////////

//--Jpsi mass windows
const Double_t m_min_jpsi = 2.85; 
const Double_t m_max_jpsi = 3.35; 
//--Psi' mass windows
const Double_t m_min_psi2 = 3.5; 
const Double_t m_max_psi2 = 3.9; 
//--Mass range
const Double_t m_range_min = 2; 
const Double_t m_range_max = 6; 
//--Pt range
const Double_t pt_range_min = 0.; 
const Double_t pt_range_max = 2.7; 
//--Pt cut for mass fit
const Double_t pt_cut = 0.25; 

//Roofit variabels
RooRealVar fMCMuMuM("fMCMuMuM","fMCMuMuM",0);
RooRealVar fMCMuMuPt("fMCMuMuPt","fMCMuMuPt",0);
RooRealVar fMCMuMuY("fMCMuMuY","fMCMuMuY",0);
RooRealVar fMCRunNum("fMCRunNum","fMCRunNum",0);

RooRealVar fMuMuM("fMuMuM","fMuMuM",0);
RooRealVar fMuMuPt("fMuMuPt","fMuMuPt",0);
RooRealVar fMuMuY("fMuMuY","fMuMuY",0);
RooRealVar fMuMuPhi("fMuMuPhi","fMuMuPhi",0);
RooRealVar fRunNum("fRunNum","fRunNum",0);
RooRealVar fADADecision("fADADecision","fADADecision",0);
RooRealVar fADCDecision("fADCDecision","fADCDecision",0);
RooRealVar fV0ADecision("fV0ADecision","fV0ADecision",0);
RooRealVar fV0CDecision("fV0CDecision","fV0CDecision",0);
RooRealVar fV0CFiredCells("fV0CFiredCells","fV0CFiredCells",0);
RooRealVar fIsZNAFired("fIsZNAFired","fIsZNAFired",0);
RooRealVar fIsZNCFired("fIsZNCFired","fIsZNCFired",0);
RooRealVar fV0AOfflineTrigger("fV0AOfflineTrigger","fV0AOfflineTrigger", 0);
RooRealVar fV0COfflineTrigger("fV0COfflineTrigger","fV0COfflineTrigger", 0);
RooRealVar fTracklets("fTracklets","fTracklets",0);
RooRealVar fMuPt1("fMuPt1","fMuPt1",0); 
RooRealVar fMuPt2("fMuPt2","fMuPt2",0);
RooRealVar fMuEta1("fMuEta1","fMuEta1",0); 
RooRealVar fMuEta2("fMuEta2","fMuEta2",0);
RooRealVar fMuPhi1("fMuPhi1","fMuPhi1",0); 
RooRealVar fMuPhi2("fMuPhi2","fMuPhi2",0);

//--Rapidity bins
const Int_t n_y_bins = 10;
// -------------------------------0------1-----2------3------4------5------6------7------8------9----
Float_t y_bins_min[n_y_bins] = {-4.0, -4.00, -3.75, -3.50, -3.25, -3.00, -2.75, -4.00, -3.50, -3.00};
Float_t y_bins_max[n_y_bins] = {-2.5, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50, -3.50, -3.00, -2.50};

// --------------------------------------------------------------------------------
// Define cuts
// --For fits
char cuts_pt_mc[500];
char cuts_m_mc[500];
char cuts_m_mc_PtAll[500];

char cuts_pt_data[500];
char cuts_m_data[500];
char cuts_m_data_PtAll[500];

// --For efficiency
char cuts_m_mc_Rec_PtCut[500];
char cuts_m_mc_Rec_PtAll[500];

char cuts_m_mc_Gen[200];