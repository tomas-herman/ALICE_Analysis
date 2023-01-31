// Function to load values for efficiencies
//-----------------------------------------------------------------------------------

void LoadPtParam(Double_t y_min, Double_t y_max, TString ADveto, TString ZDC_class, TString trigger, 
                double &b_val, double &b_err,
                double &n_val, double &n_err)
{ 
  char FileName[120];

  sprintf(FileName,"Fit_param/%s/%s_%s_%.2f_%.2f_pt_param.txt",trigger.Data(),ADveto.Data(),ZDC_class.Data(),abs(y_min),abs(y_max));

  // Open file for reading efficiency
  fstream File;
  File.open(FileName);
  // Define variabels for reading 
  string line;
  string value;
  size_t pos;

  while (std::getline(File, line)) {
    if (line.find("=") != string::npos){
      pos = line.find("=");
      value = line.substr(pos + 2, 10);
      if (line.find("b_val") != string::npos){
        b_val = atof(value.c_str());
      } 
      if (line.find("b_err") != string::npos){
        b_err = atof(value.c_str());
      }
      if (line.find("n_val") != string::npos){
        n_val = atof(value.c_str());
      } 
      if (line.find("n_err ") != string::npos){
        n_err = atof(value.c_str());
      }
    } 
  }
  File.close();
}