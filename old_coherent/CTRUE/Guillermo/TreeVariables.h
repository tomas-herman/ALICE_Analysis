// --------------------------------
// tree variables
// --------------------------------
TTree *fAnaTree; //! analysis tree
Int_t fRunNum;
UShort_t fBC;
Int_t fTracklets;
Int_t fGoodTracks;
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

// --------------------------------
// setup tree
// --------------------------------

void SetupTree(TFile *fIn)
{
  fAnaTree = (TTree*) fIn->Get("CtrueTask/fAnaTree");
  fAnaTree ->SetBranchAddress("fRunNum", &fRunNum);
  fAnaTree ->SetBranchAddress("fBCrossNum", &fBC);  
  fAnaTree ->SetBranchAddress("fTracklets", &fTracklets);	
  fAnaTree ->SetBranchAddress("fCtrue", &fCtrue);
  fAnaTree ->SetBranchAddress("fC1zed", &fC1zed);  
  fAnaTree ->SetBranchAddress("fL0inputs", &fL0inputs);
  fAnaTree ->SetBranchAddress("fL1inputs", &fL1inputs);
  fAnaTree ->SetBranchAddress("fZem1Energy", &fZem1Energy);
  fAnaTree ->SetBranchAddress("fZem2Energy", &fZem2Energy);  
  fAnaTree ->SetBranchAddress("fZNCEnergy", &fZNCEnergy);  
  fAnaTree ->SetBranchAddress("fZNAEnergy", &fZNAEnergy);
  fAnaTree ->SetBranchAddress("fZPCEnergy", &fZPCEnergy);
  fAnaTree ->SetBranchAddress("fZPAEnergy", &fZPAEnergy);  
  fAnaTree ->SetBranchAddress("fZNATDC", &fZNATDC[0]);
  fAnaTree ->SetBranchAddress("fZNCTDC", &fZNCTDC[0]);  
  fAnaTree ->SetBranchAddress("fZPATDC", &fZPATDC[0]);
  fAnaTree ->SetBranchAddress("fZPCTDC", &fZPCTDC[0]);  
  fAnaTree ->SetBranchAddress("fV0ADecision", &fV0ADecision);
  fAnaTree ->SetBranchAddress("fV0CDecision", &fV0CDecision);  
  fAnaTree ->SetBranchAddress("fADADecision", &fADADecision);
  fAnaTree ->SetBranchAddress("fADCDecision", &fADCDecision);  
  fAnaTree ->SetBranchAddress("fIR1Map", &fIR1Map);
  fAnaTree ->SetBranchAddress("fIR2Map", &fIR2Map);   

}
