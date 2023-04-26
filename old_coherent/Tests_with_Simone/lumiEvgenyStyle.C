#include "TChain.h"
#include "TH1D.h"
#include "AliTriggerClass.h"
#include "TFile.h"
#include "TH2D.h"
// #include "../task/runlist.txt"
#include "TCanvas.h"
#include "TLegend.h"

// Int_t nRuns = 361;
Int_t nRuns = 97; // 128 97
Int_t runList[] = {
    /*244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151, 245152,
    245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347, 245353,
    245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504, 245505,
    245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700, 245705,
    245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793, 245829,
    245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003, 246012,
    246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113, 246115,
    246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220, 246222,
    246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428, 246431,
    246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750, 246751,
    246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805, 246806,
    246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855, 246859,
    246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948, 246949,
    246980, 246982, 246984, 246989, 246991, 246994,*/
    /*295585, 295586, 295587, 295588, 295589, 295612, 295615, 295665, 295666, 295667,
    295668, 295671, 295673, 295675, 295676, 295677, 295714, 295716, 295717, 295718,
    295719, 295723, 295725, 295753, 295754, 295755, 295758, 295759, 295762, 295763,
    295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829,
    295831, 295854, 295855, 295856, 295859, 295860, 295861, 295863, 295881, 295908,
    295909, 295910, 295913, 295936, 295937, 295941, 295942, 295943, 295945, 295947,
    296061, 296062, 296063, 296065, 296066, 296068, 296123, 296128, 296132, 296133,
    296134, 296135, 296142, 296143, 296191, 296192, 296194, 296195, 296196, 296197,
    296198, 296241, 296242, 296243, 296244, 296246, 296247, 296269, 296270, 296273,
    296279, 296280, 296303, 296304, 296307, 296309, 296312, 296377, 296378, 296379,
    296380, 296381, 296383, 296414, 296419, 296420, 296423, 296424, 296433, 296472,
    296509, 296510, 296511, 296514, 296516, 296547, 296548, 296549, 296550, 296551,
    296552, 296553, 296615, 296616, 296618, 296619, 296622, 296623 */
    296690, 296691, 296694, 296749, 296750, 296781, 296784, 296785, 296786, 296787,
    296791, 296793, 296794, 296799, 296836, 296838, 296839, 296848, 296849, 296850,
    296851, 296852, 296890, 296894, 296899, 296900, 296903, 296930, 296931, 296932,
    296934, 296935, 296938, 296941, 296966, 296967, 296968, 296969, 296971, 296975,
    296976, 296979, 297029, 297031, 297035, 297085, 297117, 297118, 297119, 297123,
    297124, 297128, 297129, 297132, 297133, 297193, 297194, 297196, 297218, 297219,
    297221, 297222, 297278, 297310, 297312, 297315, 297317, 297363, 297366, 297367,
    297372, 297379, 297380, 297405, 297408, 297413, 297414, 297415, 297441, 297442,
    297446, 297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540,
    297541, 297542, 297544, 297558, 297588, 297590, 297595 
};


// void lumiEvgenyStyle(TString className = "CMUL7-B-NOPF-MUFAST"){
void lumiEvgenyStyle(TString className = "CMUP6-B-NOPF-MUFAST"){
//void lumiEvgenyStyle(TString className = "CCUP9-B-NOPF-CENTNOTRD"){
  TH1D* hTrigRecorded = new TH1D("hTrigRecorded","",nRuns,0,nRuns);
  TH1D* hTrigAnalysed = new TH1D("hTrigAnalysed","",nRuns,0,nRuns);
  TH1D* hLumiRecorded = new TH1D("hLumiRecorded","",nRuns,0,nRuns);
  TH1D* hLumiAnalysed = new TH1D("hLumiAnalysed","",nRuns,0,nRuns);
  Int_t run;
  TObjArray* classes = new TObjArray();
  ULong64_t class_l0b[100];
  ULong64_t class_l2a[100];
  Double_t  class_lumi[100];
  Double_t mu = 0;
  TChain* t = new TChain("trending");
  // t->AddFile("trending2015.root");
  // t->AddFile("trending_fixed.root");
  // t->AddFile("trending2018.root"); //- old
  t->AddFile("Lumi/trending_merged_PbPb2018.root"); //- old
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("classes",&classes);
  t->SetBranchAddress("class_lumi",&class_lumi);
  t->SetBranchAddress("class_l2a",&class_l2a);
  t->SetBranchAddress("class_l0b",&class_l0b);
  t->BuildIndex("run");

  AliTriggerClass* cl     = 0x0;
  AliTriggerClass* clRef  = 0x0;
  AliTriggerClass* clLT   = 0x0;
  TFile* f = new TFile("AnalysisResultsLHC18qr_ptCoherent.root");
  // TFile* f = new TFile("AnalysisResultsLHC15o14082019.root");
  TList* list = (TList*) f->Get("MyTask/ADcheck_"); // provide your list here
  TH2F* h2 = (TH2F*) list->FindObject("fTriggersVsRunH");


  for (Int_t i=0;i<nRuns;i++){
    Int_t r = runList[i];
    char* srun = Form("%i",r);
    t->GetEntryWithIndex(r);
    printf("%i ",run);
    className = "CMUP6-B-NOPF-MUFAST";
    // if (run<=245542) className = "CMUP10-B-NOPF-MUFAST";
    // AliTriggerClass* cl     = (AliTriggerClass*) classes->FindObject(className.Data());
    // AliTriggerClass* clRef  = (AliTriggerClass*) classes->FindObject("C0V0M-B-NOPF-CENTNOTRD");
    // AliTriggerClass* clLT   = (AliTriggerClass*) classes->FindObject("CMUL7-B-NOPF-MUFAST");
    cl     = (AliTriggerClass*) classes->FindObject(className.Data());
    clRef  = (AliTriggerClass*) classes->FindObject("C0V0M-B-NOPF-CENTNOTRD");
    clLT   = (AliTriggerClass*) classes->FindObject("CMUL7-B-NOPF-MUFAST");
    Int_t icl        = classes->IndexOf(cl);
    Int_t iclassRef  = classes->IndexOf(clRef);
    Int_t iclassLT   = classes->IndexOf(clLT);
    Double_t l0bRef = class_l0b[iclassRef];
    Double_t l0bLT  = class_l0b[iclassLT];
    Double_t l2aLT  = class_l2a[iclassLT];
    Double_t l2a    = class_l2a[icl];
    Double_t lumi   = class_lumi[iclassLT];
    Double_t eff = 0.615;
    if (295585 <=run && run<= 295589) eff = 0.536;
    if (295612 <=run && run<= 295615) eff = 0.532;
    if (295665 <=run && run<= 296198) eff = 0.518;
    if (296241 <=run && run<= 297595) eff = 0.514;
    Double_t sigma = 7.66*eff;
    if (295671 <=run && run<= 295715) l0bRef/=0.135;
    if (295717 <=run && run<= 295720) l0bRef/=0.135;

    Double_t lumiSeen = l0bRef/sigma*mu*eff/(1-exp(-mu*eff));
    Double_t lt = l2aLT/l0bLT;
    Double_t lumiMy = lumiSeen*lt/1000000.;
    printf("%f %f %f\n",lumiMy,lumi,lumiMy/lumi);
//    if (!cl) continue;
//    Int_t iclass = classes->IndexOf(cl);
//    if (iclass<0) continue;
//    Double_t l0b = class_l0b[iclass];
//    Double_t l2a = class_l2a[iclass];
//
//    Double_t lumiMy = l0bRef/4.*l2a/l0b;
////    if ((run>=295671 && run<=295715) || (run>=295717 && run<=295720)) class_lumi[iclass]/=0.135;
//    Double_t lumi = class_lumi[iclass];
//    printf("%f\n",lumiMy/1000000./lumi);
    hTrigRecorded->Fill(srun,l2a);
    hLumiRecorded->Fill(srun,lumiMy);


//     TFile* f = new TFile("AnalysisResultsLHC18qr15o30082019.root");
//     TList* list = (TList*) f->Get("MyTask/MyOutputContainer"); // provide your list here
//     // TH2D* fTriggersVsRun = (TH2D*) list->FindObject("fTriggersVsRun");
//     // TH2F* fTriggersVsRun = (TH2F*) list->FindObject("fTriggersVsRunH");
//
//
//     // TFile* f = new TFile(Form("../data18/AnalysisResults.%i.root",r));
// //    TFile* f = new TFile(Form("../data18/AnalysisResults.%i.root",r));
// //    TFile* f = new TFile(Form("/eos/user/e/ekryshen/trees/pb2018/mup/task/merged/AnalysisResults.000%i.root",r));
//     // TList* list = (TList*) f->Get("UpcTree/histos");
//     TH2F* h2 = (TH2F*) list->FindObject("fTriggersVsRunH");
//     Int_t nTrigAnalysed =  h2->GetBinContent(run<=245542 ? 5 : 3,h2->GetYaxis()->FindFixBin(run+0.5));
//     hTrigAnalysed->Fill(srun,nTrigAnalysed);
// //    printf("l2a=%0.f nTrigAnalysed=%i\n",l2a,nTrigAnalysed);
//     hLumiAnalysed->Fill(srun,lumi>1e-10 ? Double_t(nTrigAnalysed)/l2a*lumiMy : 0);



    Int_t class_bin = 3;
    // check for consistency with fTriggersVsRun contents
    if (className.Contains("CMUP6-B-NOPF-MUFAST")) class_bin = 3;
    // if (className.Contains("CMUP11-B-NOPF-MUFAST")) class_bin = 1;
    Int_t run_bin       = h2->GetYaxis()->FindFixBin(run+0.5);
    cout << "run bin = " << run_bin << endl;
    Int_t nTrigAnalysed = h2->GetBinContent(class_bin,h2->GetYaxis()->FindFixBin(run+0.5));
    cout << "ana bin = " << nTrigAnalysed << endl;
    Int_t nTrigRecorded = l2a;
    // Double_t lumiAnalysed = nTrigRecorded>=1 ? Double_t(nTrigAnalysed)/nTrigRecorded*lumiRecorded : 0;
    // hTrigRecorded->Fill(srun,l2a);
    // hTrigAnalysed->Fill(srun,nTrigAnalysed);
    // hLumiRecorded->Fill(srun,lumiRecorded/1000000.); // b -> ub
    // hLumiAnalysed->Fill(srun,lumiAnalysed/1000000.); // b -> ub
    // Int_t nTrigAnalysed =  h2->GetBinContent(run<=245542 ? 5 : 3,h2->GetYaxis()->FindFixBin(run+0.5));
    hTrigAnalysed->Fill(srun,nTrigAnalysed);
//    printf("l2a=%0.f nTrigAnalysed=%i\n",l2a,nTrigAnalysed);
    hLumiAnalysed->Fill(srun,lumi>1e-10 ? Double_t(nTrigAnalysed)/l2a*lumiMy : 0);

  }
  new TCanvas;
  hTrigRecorded->SetLineColor(kBlue);
  hTrigAnalysed->SetLineColor(kRed);
  hTrigRecorded->Draw("hist");
  hTrigAnalysed->Draw("hist same");

  new TCanvas;
  hLumiRecorded->SetLineColor(kBlue);
  hLumiAnalysed->SetLineColor(kRed);
  hLumiRecorded->Draw("hist");
  hLumiAnalysed->Draw("hist same");

  TLegend* l = new TLegend(0.5,0.7,0.9,0.9);
  l->AddEntry(hLumiRecorded,Form("Recorded = %.3f /mb",hLumiRecorded->Integral()));
  l->AddEntry(hLumiAnalysed,Form("Analysed = %.3f /mb",hLumiAnalysed->Integral()));
  l->Draw();

  printf("%f\n",hLumiAnalysed->Integral());
  TFile* f2 = new TFile("lumi18.root","recreate");
  hLumiRecorded->Write();
  hLumiAnalysed->Write();
  f2->Close();

  for( Int_t iBins = 1; iBins < hLumiAnalysed->GetNbinsX(); iBins++) {

  // else if ( fRunNum == 246809 ) { fLumiPerRun = 2.47068;  }

  cout << "else if ( fRunNum == "  << hLumiAnalysed->GetXaxis()->GetBinLabel(iBins);
  cout << " ) { fLumiPerRun = " << hLumiAnalysed->GetBinContent(iBins) << "; }" << endl;
}

}
