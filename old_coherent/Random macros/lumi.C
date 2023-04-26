#include <TH1.h>
#include "GoodRuns.h"
#include "SeenTriggers.h"

void lumi(Int_t opt = 0){
  
  TString className;
  Int_t nRuns;
  if (opt == 0) {
    className.Append("CMUP14-B-NOPF-MUFAST");
    nRuns = nGoodRunsLHC16r;
  }  else {
    className.Append("CMUP23-B-NOPF-MUFAST");
    nRuns = nGoodRunsLHC16s;
  }

  TChain* t = new TChain("trending");
  t->AddFile("trending.root");
  TObjArray* classes = new TObjArray();
  Double_t  class_lumi[100] = {0};
  Double_t  class_ds[100] = {0};
  ULong64_t  class_l2a[100] = {0};    
  Int_t run;
  Double_t mu = 0;
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("classes",&classes);
  t->SetBranchAddress("class_lumi",&class_lumi);
  t->SetBranchAddress("class_ds",&class_ds);
  t->SetBranchAddress("class_l2a",&class_l2a);  
  t->BuildIndex("run");
  
  TH1D* hLumi = new TH1D("hLumi","",nRuns,0,nRuns);
  TH1D* hLumiS = new TH1D("hLumiS","",nRuns,0,nRuns);
  TH1D* hLumiSPu = new TH1D("hLumiSPu","",nRuns,0,nRuns);

  TH1D* hScale = new TH1D("hScale","",nRuns,0,nRuns);
  TH1D* hTrigger_ds = new TH1D("hTrigger_ds","hTrigger_ds",nRuns,0,nRuns);  
  TH1D* hMu = new TH1D("hMu","#mu(INEL)",nRuns,0,nRuns);
  TH1D* hPileup = new TH1D("hPileup","Pile-up",nRuns,0,nRuns);  

  Double_t lumi_tot=0;
  Double_t lumiS_tot=0;
  Double_t lumiSPu_tot=0;    
  for (Int_t i=0;i<nRuns;i++){
    Int_t r;
    Int_t s;
    if (opt ==0) {
      r = GoodRunsLHC16r[i];
      s = SeenLHC16r[i];
    } else {
      r = GoodRunsLHC16s[i];
      s = SeenLHC16s[i];
    }
    char* srun = Form("%i",r);
    t->GetEntryWithIndex(r);
    //  if (r==265589) classes->ls(); // run in LHC16s    
    //    if (r==266405) classes->ls(); // run in LHC16s
    AliTriggerClass* cl = (AliTriggerClass*) classes->FindObject(className.Data());
    if (!cl) cout << " Run = " << r << " idx " << i << " missing" <<endl;
    if (!cl) continue;
    Int_t iclass = classes->IndexOf(cl);
    
    Double_t l2a = (Double_t) class_l2a[iclass];
    Double_t scale;
    if (l2a < 1) cout << " Run = " << r << " idx " << i << " has no events" <<endl;
    if (l2a < 1) continue;
    Double_t scale = ((Double_t) s)/l2a;
    hScale->Fill(srun,scale);
    // cout << r << " " << l2a << " " << s << " " << scale << endl;
    hLumi->Fill(srun,class_lumi[iclass]);
    hLumiS->Fill(srun,scale*class_lumi[iclass]);    
    hTrigger_ds->Fill(srun,class_ds[iclass]);
    hMu->Fill(srun,mu);
    Double_t pu = 1.0-TMath::Exp(-mu);
    Double_t Cpu = 1-pu;
    hPileup->Fill(srun,pu);
    hLumiSPu->Fill(srun,scale*Cpu*class_lumi[iclass]);
    lumi_tot+=class_lumi[iclass];
    lumiS_tot+=(scale*class_lumi[iclass]);
    lumiSPu_tot+=(scale*Cpu*class_lumi[iclass]);    
    /*
    cout << " run = " << r << " ds = " << class_ds[iclass] << " mu = " << mu
	 << " lumi = " << class_lumi[iclass]
	 << "  correction " << scale
	 << endl;
    */
    cout << (scale*Cpu*class_lumi[iclass]/1000) << "," << endl; // lumi in 1/nb
    // cout << i << " " << (scale*Cpu*class_lumi[iclass]/1000) << "," << endl; // lumi in 1/nb
  }

  cout << " Integrated lumi = " << lumi_tot << " 1/mub"<< endl;
  cout << " Seen integrated lumi = " << lumiS_tot << " 1/mub"<< endl;
  cout << " Pile-up corrected seen integrated lumi = "
       << lumiSPu_tot << " 1/mub"<< endl;

  // plots
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TFile* f = NULL;
  if (opt == 0) f = new TFile("lumi_CMUP14.root","recreate");
  else f = new TFile("lumi_CMUP23.root","recreate");
  hLumi->Write();
  hLumiS->Write();
  hScale->Write();
  hLumiSPu->Write();      
  hMu->Write();
  hTrigger_ds->Write();
  hPileup->Write();
  f->Close();
  
  TCanvas *cPU = new TCanvas("cPU","cPU",1200,450);
  hPileup->SetXTitle("");
  hPileup->SetYTitle("pile-up probability");
  hPileup->Draw();

  TCanvas *cMu = new TCanvas("cMu","cMu",1200,450);
  hMu->SetXTitle("");
  hMu->SetYTitle("#mu");
  hMu->Draw();

  TCanvas *cLumi = new TCanvas("cLumi","cLumi",1200,450);
  hLumi->SetXTitle("");
  hLumi->SetYTitle("Luminosity (1/#mub)");
  hLumi->Draw("");
  hLumiS->SetMarkerColor(kRed);
  hLumiS->SetMarkerStyle(20);  
  hLumiS->Draw("same,p");  
  hLumiSPu->SetMarkerColor(kBlue);
  hLumiSPu->SetMarkerStyle(24);  
  hLumiSPu->Draw("same,p");  
  TLegend* legLumi = new TLegend(0.6,0.7,0.9,0.9);
  legLumi->SetFillColor(kWhite);
  legLumi->AddEntry(hLumi,Form("Luminosity: %.3f 1/#mub",lumi_tot),"l");
  legLumi->AddEntry(hLumiS,Form("Analysed luminosity: %.3f 1/#mub",lumiS_tot),"p");
  legLumi->AddEntry(hLumiSPu,Form("Pileup corrected analysed luminosity: %.3f 1/#mub",lumiSPu_tot),"p");
  legLumi->Draw();
  TCanvas *cDS = new TCanvas("cDS","cDS",1200,450);
  hTrigger_ds->SetXTitle("");
  hTrigger_ds->SetYTitle("Downscale factor");
  gPad->SetLogy();
  hTrigger_ds->Draw();
  
  TCanvas *cS = new TCanvas("cS","cS",1200,450);
  hScale->SetXTitle("");
  hScale->SetYTitle("Scale factor");
  //  gPad->SetLogy();
  hScale->Draw();
  
}
