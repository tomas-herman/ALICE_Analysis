//
// this headers contains the global constants used in the
// selection of events as well as the selection functions
//

// -----------------------------------------------------------------
// global constants

// general kinematic selection to ve applied to all events
const double minMuMuM = 2.0; // GeV/c2
const double maxMuMuM = 6.0; // GeV/c2 
const double maxMuMuPt = 5.0; // GeV/c
const double beamEnergy = 2510; // GeV

// auxiliary variables
int gOfflineCellsV0A = 0;
int gOfflineCellsV0C = 0;
int gMatchedCells = 0;

// -----------------------------------------------------------------
// selection functions

//-----------------------------------------------
// general kinematics to be applied to all events
bool generalKinematics()
{
  if (fMuMuM < minMuMuM) return false;
  if (fMuMuM > maxMuMuM) return false;
  if (fMuMuPt > maxMuMuPt) return false;
  return true;
}

//-----------------------------------------------
// New ZDC class criteria 
//-----------------------------------------------

//-----------------------------------------------
// checks if the event passes the XnXn selection
bool passXnXnSelection(int extraTracklets, int extraCellsV0C, float maxNeutronsZN)
{
  // tracklets
  if (fTracklets>extraTracklets) return false;
  // ad
  if (!(fADADecision == 0 || fADADecision == 1)) return false;
  if (!(fADCDecision == 0 || fADCDecision == 1)) return false;
  // v0
  if (fV0ADecision != 0) return false;
  if (gOfflineCellsV0C > (gMatchedCells+extraCellsV0C) ) return false;
  // zn
  // if (((fZNAEnergy/beamEnergy)+(fZNCEnergy/beamEnergy)) > maxNeutronsZN) return false;
  // event passes selection
  return true;
}

//-----------------------------------------------
// checks if the event passes the Xn0n selection
bool passXn0nSelection(int extraTracklets, int extraCellsV0C, float maxNeutronsZN)
{
  // tracklets
  if (fTracklets>extraTracklets) return false;
  // ad
  if (!(fADADecision == 0 || fADADecision == 1)) return false;
  // if (fADCDecision != 0) return false;
  // v0
  if (fV0ADecision != 0) return false;
  if (gOfflineCellsV0C > (gMatchedCells+extraCellsV0C) ) return false;
  // zn
  // if ((fZNAEnergy/beamEnergy) > maxNeutronsZN) return false;
  // event passes selection
  return true;
}

//-----------------------------------------------
// checks if the event passes the 0nXn selection
bool pass0nXnSelection(int extraTracklets, int extraCellsV0C, float maxNeutronsZN)
{
  // tracklets
  if (fTracklets>extraTracklets) return false;
  // ad
  // if (fADADecision != 0) return false;
  if (!(fADCDecision == 0 || fADCDecision == 1)) return false;
  // v0
  if (fV0ADecision != 0) return false;
  if (gOfflineCellsV0C > (gMatchedCells+extraCellsV0C) ) return false;
  // zn
  // if ((fZNCEnergy/beamEnergy) > maxNeutronsZN) return false;
  // event passes selection
  return true;
}

//-----------------------------------------------
// checks if the event passes the 0n0n selection
bool pass0n0nSelection(int extraTracklets, int extraCellsV0C)
{
  // tracklets
  if (fTracklets>extraTracklets) return false;
  // ad
  // if (fADADecision != 0) return false;
  // if (fADCDecision != 0) return false;
  // v0
  if (fV0ADecision != 0) return false;
  if (gOfflineCellsV0C > (gMatchedCells+extraCellsV0C) ) return false;
  // event passes the selection
  return true;
}

//-----------------------------------------------
// Original ZDC criteria used in the Analysis
//-----------------------------------------------

// checks if the event passes the selection
bool passOriginSelection(int extraTracklets, int extraCellsV0C, float maxNeutronsZN)
{
  // tracklets
  if (fTracklets>extraTracklets) return false;
  // ad
  // if (!(fADADecision == 0 || fADADecision == 1)) return false;
  // if (!(fADCDecision == 0 || fADCDecision == 1)) return false;
  // v0
  if (fV0ADecision != 0) return false;
  if (fV0CDecision != 0 && fV0CDecision != 1) return false;
  if (fV0CFiredCells >= 3) return false;
  // if (gOfflineCellsV0C > (gMatchedCells+extraCellsV0C) ) return false;
  // zn
  // if (((fZNAEnergy/beamEnergy)+(fZNCEnergy/beamEnergy)) > maxNeutronsZN) return false;
  // event passes selection
  return true;
}

//-----------------------------------------------
// perform different event selections
// it depends on the ZN class
bool SelectEvent(int iSel, int iClassZN)
{
  // to have variations in the selection
  int extraTracklets = 500;
  int extraCellsV0C = 0;
  float maxNeutronsZN = 100;
  
  // general kinematics applied to all events
  if (!generalKinematics()) return false;
  
  // select
  bool isSel = false;
  if (iSel == 0) { // standard selection
    if ((iClassZN<0) || (iClassZN==0)) { // mc or 0n0n
      isSel = pass0n0nSelection(extraTracklets,extraCellsV0C);
    } else if (iClassZN==1) { // 0nxn
      isSel = pass0nXnSelection(extraTracklets,extraCellsV0C,maxNeutronsZN);
    } else if (iClassZN==2) { // xn0n
      isSel = passXn0nSelection(extraTracklets,extraCellsV0C,maxNeutronsZN);
    } else if (iClassZN==3) { // xnxn
      isSel = passXnXnSelection(extraTracklets,extraCellsV0C,maxNeutronsZN);
    }
  } // end of: iSel == 0

  if (iSel == 1) { // Simone/original selection, same for all classes
    extraTracklets = 999999999;
    isSel = passOriginSelection(extraTracklets,extraCellsV0C, maxNeutronsZN);
  } // end of: iSel == 1

  if (iSel == 2) { // Simone/original selection with a tracklet cut 2, same for all classes
    extraTracklets = 2;
    isSel = passOriginSelection(extraTracklets,extraCellsV0C, maxNeutronsZN);
  } // end of: iSel == 2

  if (iSel == 3) { // Simone/original selection with tracklet cut 0, same for all classes
    extraTracklets = 0;
    isSel = passOriginSelection(extraTracklets,extraCellsV0C, maxNeutronsZN);
  } // end of: iSel == 3
  
  if (iSel == 4) { // Simone/original selection with tracklet cut 1, same for all classes
    extraTracklets = 1;
    isSel = passOriginSelection(extraTracklets,extraCellsV0C, maxNeutronsZN);
  } // end of: iSel == 4

  // end of selections
  return isSel;
}

