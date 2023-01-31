//
// this program creates a tree to be used by roofit
// as input takes files produced by the lego train
// and produces an output with a flag with different selections
//

// -----------------------------------------------------------------
// all headers are defined here

// c++ headers
#include <iostream>
#include <fstream>

// root headers
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>

// my headers
#include "chainVariables.h"
#include "selection.h"


//-----------------------------------------------
// tag the event as 0n0n, 0nXn, Xn0n or XnXn
// MC has iData 10 or larger
int getEventClass(int iData)
{
  if (iData>9) return -1; // MC: not defined
  if ((fIsZNAFired == 0) &&  (fIsZNCFired == 0)) return 0; // 0n0n
  if ((fIsZNAFired == 0) &&  (fIsZNCFired == 1)) return 1; // 0nxn
  if ((fIsZNAFired == 1) &&  (fIsZNCFired == 0)) return 2; // xn0n     
  if ((fIsZNAFired == 1) &&  (fIsZNCFired == 1)) return 3; // xnxn
  // we should not reach this point
  cout << " getEventClass: we should not reach this point ..."
       << " fIsZNAFired = " << fIsZNAFired
       << " fIsZNCFired = " << fIsZNCFired
       << " iData = " << iData
       << endl;
 exit(-1);
}

// -----------------------------------------------------------------
// auxiliary function to relate a muon to a v0c cell
Int_t getV0CIdx(Double_t eta, Double_t phi)
{
  Int_t i_eta = -1;
  Int_t idx = -1;
  if (eta>-3.7 && eta<-3.2) i_eta = 0;
  else if (eta>-3.2 && eta<-2.7) i_eta = 1;
  else if (eta>-2.7 && eta<-2.2) i_eta = 2;  
  
  if (i_eta> -1) {
    Int_t i_phi = (Int_t) ((4.0*phi/TMath::Pi()));
    idx = i_eta*8+i_phi; 
  }

  return idx;
}

// -----------------------------------------------------------------
// Get number of muons with an active V0C cell
int getMatchedV0C()
{
  int n=0;

  // get the v0c cell indices
  Int_t idx1 = getV0CIdx(fMuEta1,fMuPhi1);
  Int_t idx2 = getV0CIdx(fMuEta2,fMuPhi2);
  // check if they have signal
  if (idx1>-1 && fV0COfflineTrigger[idx1]) n++;
  if (idx2>-1 && fV0COfflineTrigger[idx2]) n++;  

  if (idx1>-1 && idx2>-1 && idx1==idx2) {
    return (100+idx1);
  }
  // return number of matched tracks
  return n;
}

//-----------------------------------------------
// produce a tree to be used in RooFit
void produceFitTree(int iData = 16, int iSel = 1)
// iData = specify input data
// iSel = specify the selection
  
{
  // open input file and get the info
  TChain *chain = new TChain("NanoMUONCMUP6/fRecTree");
  addFileToChain(iData, chain);
  cout << " Data found. Chain has " << chain->GetEntries() << " entries " << endl;
  // set up the branches
  setChainBranches(chain);

  // set up the output file and tree
  TFile *fOut = new TFile(Form("/mnt/Data/Data_processing/Processed_data/TrainResults/fitTrees/fitTree_Data_%d_%d.root",iData,iSel),"recreate");
  TTree* fitTree = new TTree("fitTree", "fitTree");
  double mass; // GeV/c^2
  fitTree ->Branch("mass", &mass, "mass/D");
  double rap;
  fitTree ->Branch("rap", &rap, "rap/D");
  double pt; // GeV/c
  fitTree ->Branch("pt", &pt, "pt/D");
  int znclass; // negative for MC, 0=>0n0n, 1=>0nXn, 2=>Xn0n, 3=>XnXn
  fitTree ->Branch("znclass", &znclass, "znclass/I");
  
  // fill the output tree
  // --> loop over all entries
  int nSelected = 0;
  int nEntries = chain->GetEntries();
  for (int iEntry = 0; iEntry<nEntries; iEntry++) {
    // get the entry
    chain->GetEntry(iEntry);

    // apply trigger for real events
    if (iData<10 && fCMUP6Decision != 1) continue;
    
    // set the event class
    // MC has iData 10 or larger
    znclass = getEventClass(iData);

    // get number of v0 cells
    // gOfflineCellsV0A and gOfflineCellsV0C are global
    gOfflineCellsV0A = 0;
    gOfflineCellsV0C = 0;
    for (int i=0;i<32;i++) {
      if (fV0AOfflineTrigger[i]) gOfflineCellsV0A ++;
      if (fV0COfflineTrigger[i]) gOfflineCellsV0C ++;    
    }

    // get number of matched muon tracks in this event
    gMatchedCells = getMatchedV0C();
    
    // select the event
    bool isGood = SelectEvent(iSel, znclass);
    if (!isGood) continue;
    nSelected++;
    
    // compute the new variables
    mass = fMuMuM;
    rap = fMuMuY;
    pt = fMuMuPt;

    // fill the new tree
    fitTree->Fill();
  } // end loop over all entries

  // print out
  cout << "  " << nSelected << " entries were selected " << endl;

  // write the tree to the output file
  fOut->cd();
  fitTree->Write();
  fOut->Close();
 
}
