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
// produce a tree to be used in RooFit
void produceFitTree(int iData = 0, int iSel = 1)
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
  TString fullpath = Form("/mnt/Data/Data_processing/Processed_data/TrainResults/fitTrees/fitTree_Data_%d_%d.root",iData,iSel);
  gSystem->mkdir(fullpath, kTRUE);
  TFile *fOut = new TFile(fullpath,"recreate");
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
