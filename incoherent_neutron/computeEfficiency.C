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
#include <filesystem>
#include <string>


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
void computeOneEfficiency(int iData, int iSel, 
                          double minY, double maxY )
// iData = specify input data
// iSel = specify the selection
{
  // open input file and get the info
  TChain *chainReco = new TChain("NanoMUONCMUP6/fRecTree");
  addFileToChain(iData, chainReco);
  cout << " Data found. Chain has " << chainReco->GetEntries() << " entries " << endl;
  // set up the branches
  setChainBranches(chainReco);

  // open input file and get the info
  TChain *chainGen = new TChain("NanoMUONCMUP6/fGenTree");
  addFileToChain(iData, chainGen);
  cout << " Data found. Chain has " << chainGen->GetEntries() << " entries " << endl;
  // set up the branches
  setChainBranchesGenerated(chainGen);

  // Set binning
  int nBins = 7;
  double xBins [8] = {0.0, 0.3, 0.5, 0.7, 0.9, 1.2, 1.5, 3.0};

  //-------------------
  // --> loop over all Reconstructed entries
  TH1D recoHist("recoHist","recoHist",nBins,xBins);
  TH1D recoHistIntegrated("recoHistIntegrated","recoHistIntegrated",1,0.3,1.5);
  TH1D recoHistAll("recoHistAll","recoHistAll",1,0.0,3.0);
  int nEntries = chainReco->GetEntries();
  int znclass; // negative for MC, 0=>0n0n, 1=>0nXn, 2=>Xn0n, 3=>XnXn
  cout << "Entering the reconstructed events loop with " << nEntries << " events." << endl;
  for (int iEntry = 0; iEntry<nEntries; iEntry++) {
    // get the entry
    chainReco->GetEntry(iEntry);

    // apply trigger for real events
    if (iData<10 && fCMUP6Decision != 1) continue;
    
    // get number of v0 cells
    // gOfflineCellsV0A and gOfflineCellsV0C are global
    gOfflineCellsV0A = 0;
    gOfflineCellsV0C = 0;
    for (int i=0;i<32;i++) {
      if (fV0AOfflineTrigger[i]) gOfflineCellsV0A ++;
      if (fV0COfflineTrigger[i]) gOfflineCellsV0C ++;    
    }

    // set the event class
    // MC has iData 10 or larger
    znclass = getEventClass(iData);

    // get number of matched muon tracks in this event
    gMatchedCells = getMatchedV0C();
    
    // select the event
    bool isGood = SelectEvent(iSel, znclass);
    if (!isGood) continue;
    if (fMuMuY < minY || fMuMuY > maxY ) continue;
    recoHist.Fill(fMuMuPt);
    recoHistIntegrated.Fill(fMuMuPt);
    recoHistAll.Fill(fMuMuPt);
  } // end loop over all entries

  //-------------------
  // --> loop over all generated entries
  TH1D genHist("genHist","genHist",nBins,xBins);
  TH1D genHistIntegrated("genHistIntegrated","genHistIntegrated",1,0.3,1.5);
  TH1D genHistAll("genHistAll","genHistAll",1,0.0,3.0);
  nEntries = chainGen->GetEntries();
  cout << "Entering the generated events loop with " << nEntries << " events." << endl;
  for (int iEntry = 0; iEntry<nEntries; iEntry++) {
    // get the entry
    chainGen->GetEntry(iEntry);

    // select the event
    if (fMCMuMuY < minY || fMCMuMuY > maxY ) continue;
    genHist.Fill(fMCMuMuPt);
    genHistIntegrated.Fill(fMCMuMuPt);
    genHistAll.Fill(fMCMuMuPt);
  } // end loop over all entries

  //-------------------
  // Set overflow reco bin to 0 to avoid tefficeincy error due to more reco than gen events
  recoHist.SetBinContent(nBins+1,0);
  recoHistIntegrated.SetBinContent(2,0);
  recoHistAll.SetBinContent(2,0);

  //-------------------
  // Create efficiency histogram
  TEfficiency *effHist = new TEfficiency(recoHist,genHist);
  effHist->SetName("effHist");
  effHist->SetTitle("effHist;p_{T} (GeV/#it{c}); N_{rec}/N_{gen}");

  TEfficiency *effHistIntegrated = new TEfficiency(recoHistIntegrated,genHistIntegrated);
  effHistIntegrated->SetName("effHistIntegrated");
  effHistIntegrated->SetTitle("effHistIntegrated;p_{T} (GeV/#it{c}); N_{rec}/N_{gen}");

  TEfficiency *effHistAll = new TEfficiency(recoHistAll,genHistAll);
  effHistAll->SetName("effHistAll");
  effHistAll->SetTitle("effHistAll;p_{T} (GeV/#it{c}); N_{rec}/N_{gen}");

  cout << "Finished computing efficiency in rapidity range ("  << minY << "," << maxY << ")" <<endl;
  for (int i = 0; i < nBins; ++i)
  {
    cout << "Efficiency in pt range (" << xBins[i] << "," << xBins[i+1] << ") is: " << effHist->GetEfficiency(i+1) << endl;
  }
  cout << "Efficiency in pt range (" << 0.3 << "," << 1.5 << ") is: " << effHistIntegrated->GetEfficiency(1) << endl;
  cout << "Efficiency in pt range (" << 0.0 << "," << 3.0 << ") is: " << effHistAll->GetEfficiency(1) << endl;

  TString fullpath = Form("efficiency/Selection_%i/eff_data_%i_%.2f_y_%.2f.root",iSel,iData,abs(minY),abs(maxY));
  gSystem->mkdir(fullpath, kTRUE);
  TFile *fOut = new TFile(fullpath,"recreate");

  // write the histograms to the output file 
  fOut->cd();
  recoHist.Write();
  genHist.Write();
  recoHistIntegrated.Write();
  genHistIntegrated.Write();
  recoHistAll.Write();
  genHistAll.Write();
  effHist->Write();
  effHistIntegrated->Write();
  effHistAll->Write();
  fOut->Close();

}

void computeEfficiency()
{

  int iData = 11;
  int iSel = 1;
  double minRap[3] = {-4.0, -4.00, -3.25 };
  double maxRap[3] = {-2.5, -3.25, -2.50 };

  for (int i = 0; i < 3; ++i)
  {
    computeOneEfficiency(iData, iSel, minRap[i], maxRap[i]);
  }

}
