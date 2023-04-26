// Macro to load good runs into a char
//-----------------------------------------------------------------------------------

void LoadGoodRuns(char *goodruns, int size, TString period)
{
  // opt = 0; 2018 q || opt = 90; 2018 q,r MC - LHC18l7
  // opt = 1; 2018 r 
  // opt = 2; 2015 o || opt = 91; 2015 o MC - LHC16b2

  //-----------------------------------------------------------------------------------
  // Check if the size of the char is large enough to contain all good runs

  // ROOT::v5::TFormula::SetMaxima(RunsCharSize);
  
  if(sizeof("fMCRunNum == 295585 || _")*Vector_GoodRuns15o.size()*1.1>size){

    cout << "########################################################" << endl;
    cout << "The size of the char for list of good runs is too small." << endl;
    cout << "Size = " << sizeof("fRunNum == 295585 || ")*Vector_GoodRuns15o.size()*1.1 << endl;
    cout << "########################################################" << endl;
    cout << endl;

    gROOT->ProcessLine(".q");
  }

  //-----------------------------------------------------------------------------------
  // Create lists of good runs

  char GoodRuns18q[size];
  char GoodRuns18r[size];
  char GoodRuns15o[size];
  char GoodRunsAll[size];
  char GoodRunsMCAll[size];

  sprintf(GoodRuns18q," ");
  sprintf(GoodRuns18r," ");
  sprintf(GoodRuns15o," ");
  sprintf(GoodRunsAll," ");
  sprintf(GoodRunsMCAll," ");


  // Fill lists of good runs from the header file
  char OneGoodRun[sizeof("fMCRunNum == 295585 || _")];

  for(auto i : Vector_GoodRuns18q){
    if(i == Vector_GoodRuns18q.back()){
      sprintf(OneGoodRun, "fRunNum == %d ", i);
      strcat(GoodRuns18q, OneGoodRun);
    } else {
      sprintf(OneGoodRun, "fRunNum == %d || ", i);
      strcat(GoodRuns18q, OneGoodRun);
    }
  }

  for(auto i : Vector_GoodRuns18r){
    if(i == Vector_GoodRuns18r.back()){
      sprintf(OneGoodRun, "fRunNum == %d ", i);
      strcat(GoodRuns18r, OneGoodRun);
    } else {
      sprintf(OneGoodRun, "fRunNum == %d || ", i);
      strcat(GoodRuns18r, OneGoodRun);
    }
  }

  for(auto i : Vector_GoodRuns15o){
    if(i == Vector_GoodRuns15o.back()){
      sprintf(OneGoodRun, "fRunNum == %d ", i);
      strcat(GoodRuns15o, OneGoodRun);
    } else {
      sprintf(OneGoodRun, "fRunNum == %d || ", i);
      strcat(GoodRuns15o, OneGoodRun);
    }
  }

  for(auto i : Vector_GoodRunsAll){
    if(i == Vector_GoodRunsAll.back()){
      sprintf(OneGoodRun, "fRunNum == %d ", i);
      strcat(GoodRunsAll, OneGoodRun);
    } else {
      sprintf(OneGoodRun, "fRunNum == %d || ", i);
      strcat(GoodRunsAll, OneGoodRun);
    }    
  }

  for(auto i : Vector_GoodRunsAll){
    if(i == Vector_GoodRunsAll.back()){
      sprintf(OneGoodRun, "fMCRunNum == %d ", i);
      strcat(GoodRunsMCAll, OneGoodRun);
    } else {
      sprintf(OneGoodRun, "fMCRunNum == %d || ", i);
      strcat(GoodRunsMCAll, OneGoodRun);
    }  
  }

  //-----------------------------------------------------------------------------------
  // Copy the required list to the output string

  if(period.Contains("18q")){
    for(int i=0; i < size; ++i){
      goodruns[i] = GoodRuns18q[i];
    }
  } else if(period.Contains("18r")){
    for(int i=0; i < size; ++i){
      goodruns[i] = GoodRuns18r[i];
    }
  } else if(period.Contains("15o")){
    for(int i=0; i < size; ++i){
      goodruns[i] = GoodRuns15o[i];
    }
  } else if(period.Contains("DataAll")){
    for(int i=0; i < size; ++i){
      goodruns[i] = GoodRunsAll[i];
    }
  } else if(period.Contains("MCAll")){
    for(int i=0; i < size; ++i){
      goodruns[i] = GoodRunsMCAll[i];
    }    
  } else {
      cout << "Option for the period of loading good runs not valid! Bey." << endl;
      gROOT->ProcessLine(".q");
  }

}
   
