// Function to load values for efficiencies
//-----------------------------------------------------------------------------------

void LoadEff(Double_t y_min, Double_t y_max, 
             double &Eff_CohJpsiToMu_PtCut, double &Eff_CohPsi2sToMu_PtCut, double &Eff_CohPsi2sToMuPi_PtCut,
             double &Eff_CohJpsiToMu_PtAll, double &Eff_CohPsi2sToMu_PtAll, double &Eff_CohPsi2sToMuPi_PtAll,
             TString ADveto, TString trigger)
{
  char FileName[120];
  sprintf(FileName,"/home/tomas/cernbox/work/ALICE_Analysis/Efficiency/%s/%s_%.2f_%.2f.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));

  // Open file for reading efficiency
  fstream EffFile;
  EffFile.open(FileName);
  // Define variabels for reading 
  string line;
  string eff_value;
  size_t pos;

  while (std::getline(EffFile, line)) {
    if (line.find("=") != string::npos){
      pos = line.find("=");
      eff_value = line.substr(pos + 2, 6);
      if (line.find("Eff_CohJpsiToMu_PtCut ") != string::npos){
        Eff_CohJpsiToMu_PtCut = atof(eff_value.c_str());
      } else if (line.find("Eff_CohJpsiToMu_PtAll ") != string::npos){
        Eff_CohJpsiToMu_PtAll = atof(eff_value.c_str());
      } else if (line.find("Eff_CohPsi2sToMu_PtCut ") != string::npos){
        Eff_CohPsi2sToMu_PtCut = atof(eff_value.c_str());
      } else if (line.find("Eff_CohPsi2sToMu_PtAll ") != string::npos){
        Eff_CohPsi2sToMu_PtAll = atof(eff_value.c_str());  
      } else if (line.find("Eff_CohPsi2sToMuPi_PtCut ") != string::npos){
        Eff_CohPsi2sToMuPi_PtCut = atof(eff_value.c_str());
      } else if (line.find("Eff_CohPsi2sToMuPi_PtAll ") != string::npos){
        Eff_CohPsi2sToMuPi_PtAll = atof(eff_value.c_str());  
      }  
    } 
  }

  EffFile.close();
}