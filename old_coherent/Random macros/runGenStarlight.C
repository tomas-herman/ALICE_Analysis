#include "/home/ekryshen/alice/AliDPG/MC/GeneratorConfig.C"
#include "/home/ekryshen/alice/AliRoot/STARLIGHT/AliStarLight/AliGenStarLight.h"

TString processConfig = "kCohJpsiToMu";
//TString processConfig = "kIncohJpsiToMu";
//TString processConfig = "kTwoGammaToMuLow";
//TString processConfig = "kIncohPsi2sToMu";
TString comment = "";
Int_t seedConfig = 1;

TString systemConfig = "Pb-Pb";
Double_t yminConfig = -7.0;
Double_t ymaxConfig = +7.0; // !!! ymaxConfig must be equal to -yminConfig 
Double_t energyConfig = 5020;
Int_t nYbins = (ymaxConfig-yminConfig)/0.1; // binWidth = 0.1;
Double_t cs = 38.789; // mb

void runGenStarlight(Int_t nEvents=10000000){
  
  AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"recreate");
  rl->MakeTree("E");
  rl->MakeStack();
  AliStack* stack = rl->Stack();
  gAlice->SetRunLoader(rl);

  AliGenStarLight* gener = (AliGenStarLight*) GeneratorStarlight();
  gener->SetStack(stack);
  gener->Init();

  TH1D* hY  = new TH1D(
      Form("hY_%s_%s_%.0f"                                 ,processConfig.Data(),systemConfig.Data(),energyConfig),
      Form("%s, %s@%.3f TeV, #sigma=%.3f ub;y; #sigma (mb)",processConfig.Data(),systemConfig.Data(),energyConfig/1000.,cs),
      nYbins,-yminConfig,ymaxConfig
    );
  
  TH2D* hYM = new TH2D(
      Form("hYM_%s_%s_%.0f"                                ,processConfig.Data(),systemConfig.Data(),energyConfig),
      Form("%s, %s@%.3f TeV, #sigma=%.3f ub;y; #sigma (mb)",processConfig.Data(),systemConfig.Data(),energyConfig/1000.,cs),
      nYbins,-yminConfig,ymaxConfig, 
      146,0.4,15.0
    );
  
  TLorentzVector p;
  for (Int_t iev=0; iev<nEvents; iev++){
    if (iev%1000==0) printf("Event: %i\n",iev);
    stack->Reset();
    gener->Generate();
    TLorentzVector pAll;
    for (Int_t i=0;i<stack->GetNtrack();i++){
      TParticle* part = stack->Particle(i);
      part->Momentum(p);
      pAll+=p;
    }
    hY->Fill(pAll.Rapidity());
    hYM->Fill(pAll.Rapidity(),pAll.M());
  }
  new TCanvas;
  hY->Draw();

  hY->Scale(cs/hY->Integral());

  new TCanvas;
  hYM->Draw("colz");
  hYM->Scale(cs/hYM->Integral());
  
  TFile* f = new TFile("cs.root","update");
  hY->Write("",TObject::kOverwrite);
  hYM->Write("",TObject::kOverwrite);
  f->Close();
}
