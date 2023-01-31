// --------------------------------
// tree variables
// --------------------------------
TTree *fAnaTree; //! analysis tree
Int_t fRunNum;
UShort_t fBC;
Int_t fTracklets;
Int_t fCtrue;
Int_t fC1zed;	

UInt_t fL0inputs;
UInt_t fL1inputs;	

Double_t fZem1Energy;
Double_t fZem2Energy; 	

Double_t fZNCEnergy; 
Double_t fZNAEnergy;
Double_t fZPCEnergy; 
Double_t fZPAEnergy;
Double_t fZNATDC[4];
Double_t fZNCTDC[4];
Double_t fZPATDC[4];
Double_t fZPCTDC[4];
Int_t fV0ADecision; 
Int_t fV0CDecision;
Int_t fADADecision; 
Int_t fADCDecision;

TBits *fIR1Map;
TBits *fIR2Map;
