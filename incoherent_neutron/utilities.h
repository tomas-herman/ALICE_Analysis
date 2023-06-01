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
// Get the yields frm gamma mass fit results
double getYieldFitResults(int iSel, int znSelection, TString valOrErr,TString paramName, int bin)
{
  TFile *inFile = TFile::Open(Form("zdcClassYields/Selection_%i/Yields_ZNclass%d.root",iSel,znSelection));

  TVector *paramVec;
  if (paramName.Contains(Form("nCbYield_%i",bin))) { paramVec = (TVector*)inFile->Get(Form("nCbYield_%i",bin)); }
  else if (paramName.Contains(Form("nBkgdYield_%i",bin))) { paramVec = (TVector*)inFile->Get(Form("nBkgdYield_%i",bin)); }
  else {
    cout << "Parameter " << paramName.Data() << " not found! Exiting." << endl;
    throw;
  }

  double par;
  // return the value
  if (valOrErr.Contains("val")){ par = paramVec->operator[](0); }
  else if (valOrErr.Contains("err")){ par = paramVec->operator[](1); }
  else {
    cout << "Select val or err. " << valOrErr.Data() << " is not known! Exiting." << endl;
    throw;
  }
  return (par);  
}

//-----------------------------------------------
// Get the mass fit results
double getMassFitResults(int iSel, int znSelection, TString valOrErr,TString paramName,double minPt,double maxPt,double minRap,double maxRap)
{
  TFile *inFile = TFile::Open(Form("massFitResults/Selection_%i/Data_ZNclass%d_%.2f_y_%.2f_%.2f_pt_%.2f.root",iSel,znSelection,abs(minRap),abs(maxRap),minPt,maxPt));

  TVector *paramVec;
  if (paramName.Contains("nCbParam")) { paramVec = (TVector*)inFile->Get("nCbParam"); }
  else if (paramName.Contains("nBkgdParam")) { paramVec = (TVector*)inFile->Get("nBkgdParam"); }
  else {
    cout << "Parameter " << paramName.Data() << " not found! Exiting." << endl;
    throw;
  }

  double par;
  // return the value
  if (valOrErr.Contains("val")){ par = paramVec->operator[](0); }
  else if (valOrErr.Contains("err")){ par = paramVec->operator[](1); }
  else {
    cout << "Select val or err. " << valOrErr.Data() << " is not known! Exiting." << endl;
    throw;
  }
  return (par);  
}

//-----------------------------------------------
// Get a MC parameter for mass fit
double getMcMassParam(int iData,TString valOrErr,TString paramName,double minPt,double maxPt,double minRap,double maxRap)
{
  TFile *inFile = TFile::Open(Form("massMcParameters/mcParam_data_%i_%.2f_y_%.2f_%.2f_pt_%.2f.root",iData,abs(minRap),abs(maxRap),minPt,maxPt));

  TVector *paramVec;
  if (paramName.Contains("m0Param")) { paramVec = (TVector*)inFile->Get("m0Param"); }
  else if (paramName.Contains("sigmaLParam")) { paramVec = (TVector*)inFile->Get("sigmaLParam"); }
  else if (paramName.Contains("sigmaRParam")) { paramVec = (TVector*)inFile->Get("sigmaRParam"); }
  else if (paramName.Contains("nLParam")) { paramVec = (TVector*)inFile->Get("nLParam"); }
  else if (paramName.Contains("nRParam")) { paramVec = (TVector*)inFile->Get("nRParam"); }
  else if (paramName.Contains("alphaLParam")) { paramVec = (TVector*)inFile->Get("alphaLParam"); }
  else if (paramName.Contains("alphaRParam")) { paramVec = (TVector*)inFile->Get("alphaRParam"); }

  else if (paramName.Contains("lambdaParam")) { paramVec = (TVector*)inFile->Get("lambdaParam"); }
  else if (paramName.Contains("a2Param")) { paramVec = (TVector*)inFile->Get("a2Param"); }
  else if (paramName.Contains("a3Param")) { paramVec = (TVector*)inFile->Get("a3Param"); }
  else if (paramName.Contains("a4Param")) { paramVec = (TVector*)inFile->Get("a4Param"); }

  else {
    cout << "Parameter " << paramName.Data() << " not found! Exiting." << endl;
    throw;
  }
  
  double par;
  // return the value
  if (valOrErr.Contains("val")){ par = paramVec->operator[](0); }
  else if (valOrErr.Contains("err")){ par = paramVec->operator[](1); }
  else {
    cout << "Select val or err. " << valOrErr.Data() << " is not known! Exiting." << endl;
    throw;
  }
  return (par);  
}

//-----------------------------------------------
// Get the efficiency value
double getEfficiency(int iData, int iSel, TString valOrErr, TString histName, double minRap, double maxRap, int bin)
{
  TFile *inFile = TFile::Open(Form("efficiency/Selection_%i/eff_data_%i_%.2f_y_%.2f.root",iSel,iData,abs(minRap),abs(maxRap)));

  if (!inFile) {
    cout << " File not found: " << Form("efficiency/Selection_%i/eff_data_%i_%.2f_y_%.2f.root",iSel,iData,abs(minRap),abs(maxRap)) << endl;
    throw;
  }

  TEfficiency *h = nullptr;
  double eff;

  // This is the binning used in the efficeincy computation
  // int nBins = 7;
  // double xBins [8] = {0.0, 0.3, 0.5, 0.7, 0.9, 1.2, 1.5, 3.0};

  if (histName.Contains("effHistIntegrated")) // Efficiency in pt range (0.3,1.5)
  {
    h = inFile->Get<TEfficiency>("effHistIntegrated"); 
    if (valOrErr.Contains("val")){ eff = h->GetEfficiency(1); }
  }
  else if (histName.Contains("effHistAll")) // Efficiency in pt range (0.0,3.0)
  {
    h = inFile->Get<TEfficiency>("effHistAll");
    if (valOrErr.Contains("val")){ eff = h->GetEfficiency(1); }
  }
  else
  {
    h = inFile->Get<TEfficiency>("effHist"); // Efficiency in given pt bin
    if (valOrErr.Contains("val")){ eff = h->GetEfficiency(bin+1); }
  }
  
  // return the value
  return (eff);
}

//-----------------------------------------------
// Get the pt fit results
double getPtFitResults(int iSel, int znSelection, TString valOrErr,TString paramName,double minMass,double maxMass,double minRap,double maxRap)
{
  TFile *inFile = TFile::Open(Form("ptFitResults/Selection_%i/Data_ZNclass%d_%.2f_y_%.2f_%.2f_mass_%.2f.root",iSel,znSelection,abs(minRap),abs(maxRap),minMass,maxMass));

  TVector *paramVec;
  if (paramName.Contains("fiParam")) { paramVec = (TVector*)inFile->Get("fiParam"); }
  else {
    cout << "Parameter " << paramName.Data() << " not found! Exiting." << endl;
    throw;
  }

  double par;
  // return the value
  if (valOrErr.Contains("val")){ par = paramVec->operator[](0); }
  else if (valOrErr.Contains("err")){ par = paramVec->operator[](1); }
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

  if (znSelection==5) {
    BinsWidths = {0.02, 0.025, 0.10, 0.20, 0.25, 0.50}; // GeV 
    BinsRanges = {0.00,  0.20, 0.40, 0.60, 1.00, 2.00, 3.00}; // GeV
  }
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