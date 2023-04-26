



/ *********************************** Reconstruction efficiency ****************************************
  Printf("Calculating reconstruction efficiency...");

  TEfficiency *efficiencyRec  = new TEfficiency("efficiencyOMU ","efficiencyOMU ",1,0.,1.);

  for (Int_t iEn=0;iEn<tInput->GetEntries();iEn++){
    tInput->GetEntry(iEn);
    // particle generated in wanted rapidity/pt2 interval
    if (!NanoGenSelected()) continue;
    // track was reconstructed
    if (fPt==0.) {
      efficiencyRec->Fill(kFALSE,0.5);
    }
    else {
      efficiencyRec->Fill(kTRUE,0.5);
    }
  }

  effRec = efficiencyRec->GetEfficiency(1);
  errUpStatEffRec = efficiencyRec->GetEfficiencyErrorUp(1);
  errLowStatEffRec = efficiencyRec->GetEfficiencyErrorLow(1);