// Function to load generated or reconstructed data and create datasets
//-----------------------------------------------------------------------------------

void LoadData(char const *path, TString tree_opt, TString folder, RooDataSet *dataset)
{
  // --------------------------------------------------------------------------------
  // Load data
  TFile *f=new TFile(path); 
  TTree *t;
  if(tree_opt.Contains("Gen")){ // Import generated data
    f->GetObject(folder+"/fGenTree",t); 
    RooDataSet *metadata = new RooDataSet("metadata", "metadata", RooArgSet(fMCMuMuM,fMCMuMuY,fMCMuMuPt,fMCRunNum), Import(*t));  // Read data from Tree
    dataset->append(*metadata);
    if (metadata) {delete metadata;}
  } else if(tree_opt.Contains("DataRec")){ // Import reconstructed data
    f->GetObject(folder+"/fRecTree",t);
    RooArgSet Arguments(fMuMuM,fMuMuY,fMuMuPt,fRunNum,fADADecision,fADCDecision,fV0ADecision,fV0CDecision,fV0CFiredCells);
    Arguments.add( RooArgSet(fIsZNAFired,fIsZNCFired,fV0AOfflineTrigger,fV0COfflineTrigger,fTracklets,fMuMuPhi) ); // You can only add 9 arguments at a time
    Arguments.add( RooArgSet(fMuPt1,fMuPt2,fMuEta1,fMuEta2,fMuPhi1,fMuPhi2) );
    RooDataSet *metadata = new RooDataSet("metadata", "metadata", Arguments, Import(*t));  // Read data from Tree
    dataset->append(*metadata); 
    if (metadata) {delete metadata;}
  } else if(tree_opt.Contains("MCRec")){ // Import reconstructed data
    f->GetObject(folder+"/fRecTree",t);
    RooArgSet Arguments(fMuMuM,fMuMuY,fMuMuPt,fRunNum,fADADecision,fADCDecision,fV0ADecision,fV0CDecision,fV0CFiredCells);;
    RooDataSet *metadata = new RooDataSet("metadata", "metadata", Arguments, Import(*t));  // Read data from Tree
    Arguments.add( RooArgSet(fIsZNAFired,fIsZNCFired,fV0AOfflineTrigger,fV0COfflineTrigger,fTracklets,fMuMuPhi) ); // You can only add 9 arguments at a time
    Arguments.add( RooArgSet(fMuPt1,fMuPt2,fMuEta1,fMuEta2,fMuPhi1,fMuPhi2) );
    dataset->append(*metadata); 
    if (metadata) {delete metadata;}
  } else {
      cout << "Option for loading data not valid! Bye." << endl;
      gROOT->ProcessLine(".q");
  }
  if (f) {delete f;}
}