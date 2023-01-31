int nPtBins;
vector<double> ptBinBoundariesVec;

// -----------------------------------------------------------------
// get the input tree for a given data input and selection
TTree *getTree(int iData, int iSel)
{
  // get the file
  TFile *dataFile = new TFile(Form("/mnt/Data/Data_processing/Processed_data/TrainResults/fitTrees/fitTree_Data_%d_%d.root",iData,iSel));
  if (!dataFile) {
    cout << " Data file not found " << endl;
    return NULL;
  }
  
  // return the tree
  return ((TTree *) dataFile->Get("fitTree"));
}

//-----------------------------------------------
// Get the mass fit results
double getMassFitResults(int iSel, int znSelection, TString valOrErr,TString paramName,double minPt,double maxPt,double minRap,double maxRap)
{
  TFile *inFile = TFile::Open(Form("massFitResults/Selection_%i/Data_ZNclass%d_%.2f_y_%.2f_%.2f_pt_%.2f.root",iSel,znSelection,abs(minRap),abs(maxRap),minPt,maxPt));

  TVector *ptrParam;
  if (paramName.Contains("nCbParam")) { ptrParam = (TVector*)inFile->Get("nCbParam"); }
  else if (paramName.Contains("nBkgdParam")) { ptrParam = (TVector*)inFile->Get("nBkgdParam"); }
  else {
    cout << "Parameter " << paramName.Data() << " not found! Exiting." << endl;
    throw;
  }

  TVector Param = *ptrParam;
  double par;
  // return the histo
  if (valOrErr.Contains("val")){ par = Param[0]; }
  else if (valOrErr.Contains("err")){ par = Param[1]; }
  else {
    cout << "Select val or err. " << valOrErr.Data() << " is not known! Exiting." << endl;
    throw;
  }
  return (par);  
}

//-----------------------------------------------
// Get the pt fit results
double getPtFitResults(int iSel, int znSelection, TString valOrErr,TString paramName,double minMass,double maxMass,double minRap,double maxRap)
{
  TFile *inFile = TFile::Open(Form("ptFitResults/Selection_%i/Data_ZNclass%d_%.2f_y_%.2f_%.2f_mass_%.2f.root",iSel,znSelection,abs(minRap),abs(maxRap),minMass,maxMass));

  TVector *ptrParam;
  if (paramName.Contains("fiParam")) { ptrParam = (TVector*)inFile->Get("fiParam"); }
  else {
    cout << "Parameter " << paramName.Data() << " not found! Exiting." << endl;
    throw;
  }

  TVector Param = *ptrParam;
  double par;
  // return the histo
  if (valOrErr.Contains("val")){ par = Param[0]; }
  else if (valOrErr.Contains("err")){ par = Param[1]; }
  else {
    cout << "Select val or err. " << valOrErr.Data() << " is not known! Exiting." << endl;
    throw;
  }
  return (par);  
}

// -----------------------------------------------------------------
// set bin ranges and #bins for pt template fits
void setPtBinning(int znSelection)
{
  // Set the bin sizes and ranges
  vector<double> BinsWidths;
  vector<double> BinsRanges;


  if (znSelection==0) {
    BinsWidths = {0.02, 0.025, 0.10, 0.20, 0.25, 0.50}; // GeV 
    BinsRanges = {0.00,  0.20, 0.40, 0.60, 1.00, 2.00, 3.00}; // GeV
  }
  if (znSelection==1) {
    BinsWidths = {0.04, 0.05, 0.10, 0.50}; // GeV 
    BinsRanges = {0.00, 0.20, 0.60, 2.00, 3.00}; // GeV
  }
  if (znSelection==2) {
    BinsWidths = {0.03, 0.05, 0.20, 0.50, 1.50}; // GeV 
    BinsRanges = {0.00, 0.15, 0.40, 1.00, 1.50, 3.00}; // GeV
  }
  if (znSelection==3) {
    BinsWidths = {0.05, 0.10, 0.20, 0.40, 0.50}; // GeV 
    BinsRanges = {0.00, 0.20, 0.40, 0.60, 1.00, 3.00}; // GeV
  }

  int nBinsTypes = BinsWidths.size();

  // Fill the vector containing the calculated boundaries
  double ptVal = 0.0;
  ptBinBoundariesVec.push_back(ptVal);
  for (int iBinType = 0; iBinType < nBinsTypes; iBinType++) {
    double nBinsThisType = (BinsRanges[iBinType+1] - BinsRanges[iBinType]) / BinsWidths[iBinType];
    double iBin = 0.0;
    while (iBin < nBinsThisType) {
      ptVal += BinsWidths[iBinType];
      ptBinBoundariesVec.push_back(ptVal);
      iBin++;
    }
  }
  // set the number of bins
  nPtBins = ptBinBoundariesVec.size() - 1;
}