//
//  global variables with the definition of leaves in the input chain
//

// --- for reconstructed data
 int           fRunNum;
 uint          fL0inputs;
 int           fTracklets;
 float         fZNCEnergy;
 float         fZNAEnergy;
 float         fZNATDC[4];
 float         fZNCTDC[4];
 int           fV0ADecision;
 int           fV0CDecision;
 int           fV0AFiredCells;
 int           fV0CFiredCells;
 bool          fV0AOfflineTrigger[32];
 bool          fV0COfflineTrigger[32];
 int           fADADecision;
 int           fADCDecision;
 int           fIsZNAFired;
 int           fIsZNCFired;
 float         fMuMuPt;
 float         fMuMuPhi;
 float         fMuMuY;
 float         fMuMuM;
 float         fMuPt1;
 float         fMuPt2;
 float         fMuEta1;
 float         fMuEta2;
 float         fMuPhi1;
 float         fMuPhi2;
 float         fGenMuMuPt;
 float         fGenMuMuPhi;
 float         fGenMuMuY;
 float         fGenMuMuM;
 float         fGenMuPt1;
 float         fGenMuPt2;
 float         fGenMuEta1;
 float         fGenMuEta2;
 float         fGenMuPhi1;
 float         fGenMuPhi2;
 int           fCMUP6Decision;
 int           fCMUP10Decision;
 int           fCMUP11Decision;

// --- for generated data
 int           fMCRunNum;
 float         fMCMuMuPt;
 float         fMCMuMuPhi;
 float         fMCMuMuY;
 float         fMCMuMuM;

//-----------------------------------------------
// set up the addresses for reconstructed data
void setChainBranches(TChain *chain)
{
   if (!chain) return;

   chain->SetBranchAddress("fRunNum", &fRunNum);
   chain->SetBranchAddress("fL0inputs", &fL0inputs);
   chain->SetBranchAddress("fTracklets", &fTracklets);
   chain->SetBranchAddress("fZNCEnergy", &fZNCEnergy);
   chain->SetBranchAddress("fZNAEnergy", &fZNAEnergy);
   chain->SetBranchAddress("fZNATDC", fZNATDC);
   chain->SetBranchAddress("fZNCTDC", fZNCTDC);
   chain->SetBranchAddress("fV0ADecision", &fV0ADecision);
   chain->SetBranchAddress("fV0CDecision", &fV0CDecision);
   chain->SetBranchAddress("fV0AFiredCells", &fV0AFiredCells);
   chain->SetBranchAddress("fV0CFiredCells", &fV0CFiredCells);
   chain->SetBranchAddress("fV0AOfflineTrigger", fV0AOfflineTrigger);
   chain->SetBranchAddress("fV0COfflineTrigger", fV0COfflineTrigger);
   chain->SetBranchAddress("fADADecision", &fADADecision);
   chain->SetBranchAddress("fADCDecision", &fADCDecision);
   chain->SetBranchAddress("fIsZNAFired", &fIsZNAFired);
   chain->SetBranchAddress("fIsZNCFired", &fIsZNCFired);
   chain->SetBranchAddress("fMuMuPt", &fMuMuPt);
   chain->SetBranchAddress("fMuMuPhi", &fMuMuPhi);
   chain->SetBranchAddress("fMuMuY", &fMuMuY);
   chain->SetBranchAddress("fMuMuM", &fMuMuM);
   chain->SetBranchAddress("fMuPt1", &fMuPt1);
   chain->SetBranchAddress("fMuPt2", &fMuPt2);
   chain->SetBranchAddress("fMuEta1", &fMuEta1);
   chain->SetBranchAddress("fMuEta2", &fMuEta2);
   chain->SetBranchAddress("fMuPhi1", &fMuPhi1);
   chain->SetBranchAddress("fMuPhi2", &fMuPhi2);
   chain->SetBranchAddress("fGenMuMuPt", &fGenMuMuPt);
   chain->SetBranchAddress("fGenMuMuPhi", &fGenMuMuPhi);
   chain->SetBranchAddress("fGenMuMuY", &fGenMuMuY);
   chain->SetBranchAddress("fGenMuMuM", &fGenMuMuM);
   chain->SetBranchAddress("fGenMuPt1", &fGenMuPt1);
   chain->SetBranchAddress("fGenMuPt2", &fGenMuPt2);
   chain->SetBranchAddress("fGenMuEta1", &fGenMuEta1);
   chain->SetBranchAddress("fGenMuEta2", &fGenMuEta2);
   chain->SetBranchAddress("fGenMuPhi1", &fGenMuPhi1);
   chain->SetBranchAddress("fGenMuPhi2", &fGenMuPhi2);
   chain->SetBranchAddress("fCMUP6Decision", &fCMUP6Decision);
   chain->SetBranchAddress("fCMUP10Decision", &fCMUP10Decision);
   chain->SetBranchAddress("fCMUP11Decision", &fCMUP11Decision);
}

//-----------------------------------------------
// set up the addresses for generated data
void setChainBranchesGenerated(TChain *chain)
{
   if (!chain) return;

   chain->SetBranchAddress("fMCRunNum", &fMCRunNum);
   chain->SetBranchAddress("fMCMuMuPt", &fMCMuMuPt);
   chain->SetBranchAddress("fMCMuMuPhi", &fMCMuMuPhi);
   chain->SetBranchAddress("fMCMuMuY", &fMCMuMuY);
   chain->SetBranchAddress("fMCMuMuM", &fMCMuMuM);

}

//-----------------------------------------------
// get the right file and add it to the chain
void addFileToChain(int iData, TChain *chain)
{
  // real data
  if (iData == 0) { // LHC18q and LHC18r
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18q_AnalysisResults.root");
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18r_AnalysisResults.root");    
  } else if (iData == 1) { // LHC18q 
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18q_AnalysisResults.root");
  } else if (iData == 2) { // LHC18q 
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/594_20220123-1648-Data/LHC18r_AnalysisResults.root");
  }

  // mc
  else if (iData == 10) { // LHC18l7 CohJpsiToMu
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/CohJpsiToMu_AnalysisResults.root");
  } else if (iData == 11) { // LHC18l7 IncohJpsiToMu
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/IncohJpsiToMu_AnalysisResults.root");
  } else if (iData == 12) { // LHC18l7 CohPsi2sToMu
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/CohPsi2sToMu_AnalysisResults.root");
  } else if (iData == 13) { // LHC18l7 TwoGammaToMuMedium
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/TwoGammaToMuMedium_AnalysisResults.root");
  } else if (iData == 14) { // LHC18l7 IncohPsi2sToMu
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/IncohPsi2sToMu_AnalysisResults.root");
  } else if (iData == 15) { // LHC18l7 CohPsi2sToMuPi
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/CohPsi2sToMuPi_AnalysisResults.root");
  } else if (iData == 16) { // LHC18l7 IncohPsi2sToMuPi
    chain->AddFile("/mnt/Data/Data_processing/Processed_data/TrainResults/595_20220123-1654-LHC18l7/IncohPsi2sToMuPi_AnalysisResults.root");
  }

  // unkown case
  else {
    cout << "Data set " << iData << " not known. Bye!" << endl;
    exit(-1);
  }
}
