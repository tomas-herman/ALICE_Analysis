// Function to set cuts
//-----------------------------------------------------------------------------------

void SetCuts(Double_t y_min, Double_t y_max, TString ZDC_class, TString ADveto)
{
  // For fits ------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------
  // With AD vetos
  if(ADveto.Contains("ADvetoOn")){
    sprintf(cuts_pt_mc,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                     "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                     "fADADecision == 0 && fADCDecision == 0"
                     ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
    sprintf(cuts_m_mc,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                     "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                     "fADADecision == 0 && fADCDecision == 0"
                     ,y_min, y_max, pt_cut, m_range_min, m_range_max);
    sprintf(cuts_m_mc_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                     "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                     "fADADecision == 0 && fADCDecision == 0"
                     ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    if(ZDC_class.Contains("all")){ 
      // ---ZDC all classes
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired > -1 && fIsZNCFired > -1 && "
                           "fADADecision == 0 && fADCDecision == 0"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired > -1 && fIsZNCFired > -1 && "
                          "fADADecision == 0 && fADCDecision == 0"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired > -1 && fIsZNCFired > -1 && "
                                "fADADecision == 0 && fADCDecision == 0"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("0N0N")){ 
      // ---0N0N
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 0 && fIsZNCFired == 0 && "
                           "fADADecision == 0 && fADCDecision == 0"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 0 && fIsZNCFired == 0 && "
                          "fADADecision == 0 && fADCDecision == 0"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 0 && fIsZNCFired == 0 && "
                                "fADADecision == 0 && fADCDecision == 0"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("0NXN")){ 
      // ---0NXN
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 0 && fIsZNCFired == 1 && "
                           "fADADecision == 0 && fADCDecision == 0"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 0 && fIsZNCFired == 1 && "
                          "fADADecision == 0 && fADCDecision == 0"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 0 && fIsZNCFired == 1 && "
                                "fADADecision == 0 && fADCDecision == 0"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("XN0N")){ 
      // ---XN0N
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 1 && fIsZNCFired == 0 && "
                           "fADADecision == 0 && fADCDecision == 0"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 1 && fIsZNCFired == 0 && "
                          "fADADecision == 0 && fADCDecision == 0"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 1 && fIsZNCFired == 0 && "
                                "fADADecision == 0 && fADCDecision == 0"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("XNXN")){ 
      // ---XNXN
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 1 && fIsZNCFired == 1 && "
                           "fADADecision == 0 && fADCDecision == 0"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 1 && fIsZNCFired == 1 && "
                          "fADADecision == 0 && fADCDecision == 0"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 1 && fIsZNCFired == 1 && "
                                "fADADecision == 0 && fADCDecision == 0"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);                                                                
    } else {
          cout << "Option for ZDC class not valid! Bye." << endl;
          gROOT->ProcessLine(".q");
    }
  } 
  // -------------------------
  // Without AD vetos
  if(ADveto.Contains("ADvetoOff")){
    sprintf(cuts_pt_mc,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                     "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 "
                     ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
    sprintf(cuts_m_mc,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                     "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 "
                     ,y_min, y_max, pt_cut, m_range_min, m_range_max);
    sprintf(cuts_m_mc_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                     "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 "
                     ,y_min, y_max, pt_range_max, m_range_min, m_range_max);    
    if(ZDC_class.Contains("all")){ 
      // ---ZDC all classes
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired > -1 && fIsZNCFired > -1 "
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired > -1 && fIsZNCFired > -1 "
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired > -1 && fIsZNCFired > -1 "
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("0N0N")){ 
      // ---0N0N
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 0 && fIsZNCFired == 0"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 0 && fIsZNCFired == 0"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 0 && fIsZNCFired == 0"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("0NXN")){ 
      // ---0NXN
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 0 && fIsZNCFired == 1"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 0 && fIsZNCFired == 1"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 0 && fIsZNCFired == 1"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("XN0N")){ 
      // ---XN0N
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 1 && fIsZNCFired == 0"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 1 && fIsZNCFired == 0"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 1 && fIsZNCFired == 0"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
    } else if(ZDC_class.Contains("XNXN")){ 
      
      // ---XNXN
      sprintf(cuts_pt_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                           "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                           "fIsZNAFired == 1 && fIsZNCFired == 1"
                           ,y_min, y_max, pt_range_max, m_min_jpsi, m_max_jpsi);
      sprintf(cuts_m_data,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                          "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                          "fIsZNAFired == 1 && fIsZNCFired == 1"
                          ,y_min, y_max, pt_cut, m_range_min, m_range_max);
      sprintf(cuts_m_data_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && "
                                "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                                "fIsZNAFired == 1 && fIsZNCFired == 1"
                                ,y_min, y_max, pt_range_max, m_range_min, m_range_max);                                                                
    } else {
          cout << "Option for ZDC class not valid! Bye." << endl;
          gROOT->ProcessLine(".q");
    }
  }

  // For efficiency ------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------
  // With AD vetos
  if(ADveto.Contains("ADvetoOn")){
    sprintf(cuts_m_mc_Rec_PtCut,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                       "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                       "fADADecision == 0 && fADCDecision == 0"
                       ,y_min, y_max, pt_cut, m_range_min, m_range_max);
    sprintf(cuts_m_mc_Rec_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                       "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 && "
                       "fADADecision == 0 && fADCDecision == 0"
                       ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
  } 
  // -------------------------
  // Without AD vetos
  if(ADveto.Contains("ADvetoOff")){
    sprintf(cuts_m_mc_Rec_PtCut,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                       "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 "
                       ,y_min, y_max, pt_cut, m_range_min, m_range_max);
    sprintf(cuts_m_mc_Rec_PtAll,"fMuMuY>%f && fMuMuY<%f && fMuMuPt<%f && fMuMuM>%f && fMuMuM< %f && " 
                       "fV0ADecision == 0 && (fV0CDecision == 0 || fV0CDecision == 1) && fV0CFiredCells < 3 "
                       ,y_min, y_max, pt_range_max, m_range_min, m_range_max);
  }

  sprintf(cuts_m_mc_Gen,"fMCMuMuY>%f && fMCMuMuY<%f",y_min, y_max);

}