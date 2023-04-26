// Function to load values for efficiencies
//-----------------------------------------------------------------------------------

void LoadMCmassParam(Double_t y_min, Double_t y_max, TString ADveto, TString trigger, TString Pt_str, 
                    double &lambda_val, double &lambda_err,
                    double &a2_val, double &a2_err,
                    double &a3_val, double &a3_err,
                    double &a4_val, double &a4_err,
                    double &n_cb_val, double &n_cb_err,
                    double &alpha_val, double &alpha_err,
                    double &sigmaMC_val, double &sigmaMC_err,
                    double &sigma2MC_val, double &sigma2MC_err)
{ 
  char FileName[120];

  if (Pt_str.Contains("PtAll")) sprintf(FileName,"Fit_param/%s/%s_%.2f_%.2f_MC_m_param_PtAll.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));
  if (Pt_str.Contains("PtCut")) sprintf(FileName,"Fit_param/%s/%s_%.2f_%.2f_MC_m_param_PtCut.txt",trigger.Data(),ADveto.Data(),abs(y_min),abs(y_max));

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
      if (line.find("lambda_val") != string::npos){
        lambda_val = atof(value.c_str());
      }
      if (line.find("lambda_err") != string::npos){
        lambda_err = atof(value.c_str());
      }
      if (line.find("a2_val") != string::npos){
        a2_val = atof(value.c_str());
      } 
      if (line.find("a2_err") != string::npos){
        a2_err = atof(value.c_str());
      }
      if (line.find("a3_val") != string::npos){
        a3_val = atof(value.c_str());
      } 
      if (line.find("a3_err") != string::npos){
        a3_err = atof(value.c_str());
      }
      if (line.find("a4_val") != string::npos){
        a4_val = atof(value.c_str());
      } 
      if (line.find("a4_err") != string::npos){
        a4_err = atof(value.c_str());
      }
      if (line.find("n_cb_val") != string::npos){
        n_cb_val = atof(value.c_str());
      } 
      if (line.find("n_cb_err") != string::npos){
        n_cb_err = atof(value.c_str());
      }
      if (line.find("alpha_val") != string::npos){
        alpha_val = atof(value.c_str());
      } 
      if (line.find("alpha_err") != string::npos){
        alpha_err = atof(value.c_str());
      }
      if (line.find("sigmaMC_val") != string::npos){
        sigmaMC_val = atof(value.c_str());
      } 
      if (line.find("sigmaMC_err") != string::npos){
        sigmaMC_err = atof(value.c_str());
      }
      if (line.find("sigma2MC_val") != string::npos){
        sigma2MC_val = atof(value.c_str());
      } 
      if (line.find("sigma2MC_err") != string::npos){
        sigma2MC_err = atof(value.c_str());
      }
    } 
  }
  File.close();
}