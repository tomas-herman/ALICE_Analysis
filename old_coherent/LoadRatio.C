// Function to load values for efficiencies
//-----------------------------------------------------------------------------------

void LoadRatio(Double_t y_min, Double_t y_max, TString ZDC_class, TString trigger, TString ADveto, double &ratio, double &ratio_err)
{
  char FileName[120];
  sprintf(FileName,"Ratio/%s/%s_%.2f_%.2f_%s.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max),ZDC_class.Data());

  // Open file for reading efficiency
  fstream RatioFile;
  RatioFile.open(FileName);
  // Define variabels for reading 
  string line;
  string ratio_value;
  size_t pos;

  while (std::getline(RatioFile, line)) {
    if (line.find("=") != string::npos){
      pos = line.find("=");
      ratio_value = line.substr(pos + 2, 10);
      if (line.find("ratio ") != string::npos){
        ratio = atof(ratio_value.c_str());
      } else if (line.find("ratio_err ") != string::npos){
        ratio_err = atof(ratio_value.c_str());
      }  
    } 
  }
  RatioFile.close();
}