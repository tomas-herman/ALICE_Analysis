
#include "GoodRuns.h"

void GetSeenTriggers(Int_t opt)
{
  TFile *fIn = NULL;
  if (opt == 0) fIn = new TFile("../AnalysisResultsUD_pPb_AOD_78_20161215-2003.root");
  else fIn = new TFile("../AnalysisResultsUD_pPb_AOD_83_20170501-2053.root");

  TDirectoryFile *dirfil = (TDirectoryFile*) fIn->Get("UpcFilter");
  TList *list = (TList *) dirfil->FindObjectAny("HistList");
  TH2I *h = (TH2I*) list->At(1);

  Int_t nRuns;
  if (opt == 0) nRuns = nGoodRunsLHC16r;
  else nRuns = nGoodRunsLHC16s;

  Int_t iy;
  if (opt == 0) iy = 26+1; // cmup14-b
  else iy = 32+1; // cmup23

  for(Int_t i=0;i<nRuns;i++) {
    Int_t r;
    if (opt ==0) r= GoodRunsLHC16r[i];
    else r = GoodRunsLHC16s[i];
    Int_t i0 = h->GetXaxis()->GetBinCenter(1);
    Int_t ix = r-i0+1;
    //cout << i0 << " " << r << " " << ix << " " <<  h->GetXaxis()->GetBinCenter(ix) << endl;
    /*
    for (Int_t iz = 0; iz<h->GetNbinsY(); iz++) {
      Int_t c = h->GetBinContent(ix,iz+1);
      if (c>0)    cout << r << " "<< (iz+1) << " " << c << endl;
    }
    */
    cout  << h->GetBinContent(ix,iy) << "," << endl;
  }
	
}
