// #ifndef __CINT__
// #include "RooGlobalFunc.h"
// #endif

// c++ headers
#include <iostream>
#include <fstream>

// ROOT headers
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TLatex.h"

// RooFit headers
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
using namespace RooFit;

// my headers
// #include "GoodRuns.h"



// --------------------------------------------------------------------------------
//Constants and selection parameters
const Double_t mMin0 = 2.85; //2.85
const Double_t mMax0 = 3.35; //3.35
const Double_t mMin = 3.5; //2.85
const Double_t mMax = 3.9; //3.35
const Double_t mrangeMin = 2.; //2.85
const Double_t mrangeMax = 6.; //3.35
const Double_t ptMin = 0.; 
const Double_t ptMax = 3.; //3.
const Double_t ptcut = 3; 

const Double_t n_bins = 120; //200, 150
const Double_t n_binsM = 80; //150

const Double_t scale = 1; //4,2,1,0.5,0.25,0.1,0.05
Double_t scalesMC[7] = {1, 0.0385, 0.1415, 0.2663, 0.2975, 0.2001, 0.0561};
const Double_t scaleMC = scalesMC[0]; //0.0385, 0.1415, 0.2663, 0.2975, 0.2001, 0.0561

const Int_t data_opt = 1; // 1 = data, 2 = MC 

const Int_t nYBins = 10;
Float_t gYMin[nYBins] = {-4.0, -4.00, -3.75, -3.50, -3.25, -3.00, -2.75, -4.00, -3.50, -3.00};
Float_t gYMax[nYBins] = {-2.5, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50, -3.50, -3.00, -2.50};
// Float_t gYMin[nYBins] = {-4.0, -4.00, 000, -3.50, 000, -3.00, 000};
// Float_t gYMax[nYBins] = {-2.5, -3.50, 000, -3.00, 000, -2.50, 000};
const Int_t YBin = 0;


void DoFit(Double_t y_min, Double_t y_max)
{

RooBinning bin(n_bins,ptMin,ptMax);
RooBinning binM(n_binsM,mrangeMin,mrangeMax);
// --------------------------------------------------------------------------------
  // Create roofit variables

  RooRealVar jRecM("jRecM","jRecM",2,15) ;
  RooRealVar jRecPt("jRecPt","jRecPt",0,5) ;
  RooRealVar jRecY("jRecY","jRecY",-4,-2.5) ;

  RooRealVar fMuMuM("fMuMuM","fMuMuM",0,15) ;
  RooRealVar fMuMuPt("fMuMuPt","fMuMuPt",0,5) ;
  RooRealVar fMuMuY("fMuMuY","fMuMuY",-4,-2.5) ;
  RooRealVar fRunNum("fRunNum","fRunNum",240000,300000) ;
  RooRealVar fV0ADecision("fV0ADecision","fV0ADecision",0,5) ;
  RooRealVar fADADecision("fADADecision","fADADecision",0,5) ;
  RooRealVar fADCDecision("fADCDecision","fADCDecision",0,5) ;
  RooRealVar fV0CFiredCells("fV0CFiredCells","fV0CFiredCells",0,25) ;

 
  jRecPt.setBinning(bin);
  jRecM.setBinning(binM);

// --------------------------------------------------------------------------------
// Cuts

char strReduce[120];
sprintf(strReduce,"jRecY>%f && jRecY<%f && jRecPt<%f && jRecM>%f && jRecM< %f",y_min,y_max,ptMax, mMin0, mMax0);

char strReduceM[120];
sprintf(strReduceM,"jRecY>%f && jRecY<%f && jRecPt<%f && jRecM>%f && jRecM< %f",y_min,y_max,ptcut, mrangeMin, mrangeMax);

char strReduceData[500];
sprintf(strReduceData,"jRecY>%f && jRecY<%f && jRecPt<%f && jRecM>%f && jRecM< %f &&"
                  "fV0ADecision == 0 && fV0CFiredCells < 3 && fADADecision == 0 && fADCDecision == 0",y_min,y_max,ptMax, mMin0, mMax0);

char strReduceMData[500];
sprintf(strReduceMData,"jRecY>%f && jRecY<%f && jRecPt<%f && jRecM>%f && jRecM< %f &&"
                   "fV0ADecision == 0 && fV0CFiredCells < 3 && fADADecision == 0 && fADCDecision == 0",y_min,y_max,ptcut, mrangeMin, mrangeMax);

char GoodRuns15o[3000];
sprintf(GoodRuns15o,
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d",
    244918, 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
    245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
    245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
    245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
    245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
    245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
    246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
    246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
    246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
    246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
    246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
    246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
    246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
    246949, 246980, 246982, 246984, 246989, 246991, 246994
    );

char GoodRuns18q[3000];
sprintf(GoodRuns18q,
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d",
    295585, 295586, 295587, 295588, 295589, 295612, 295615, 295665, 295666, 295667,
    295668, 295671, 295673, 295675, 295676, 295677, 295714, 295716, 295717, 295718,
    295719, 295723, 295725, 295753, 295754, 295755, 295758, 295759, 295762, 295763,
    295786, 295788, 295791, 295816, 295818, 295819, 295822, 295825, 295826, 295829,
    295831, 295854, 295855, 295856, 295859, 295860, 295861, 295863, 295881, 295908,
    295909, 295910, 295913, 295936, 295937, 295941, 295942, 295943, 295945, 295947,
    296061, 296062, 296063, 296065, 296066, 296068, 296123, 296128, 296132, 296133,
    296134, 296135, 296142, 296143, 296191, 296192, 296194, 296195, 296196, 296197,
    296198, 296241, 296242, 296243, 296244, 296246, 296247, 296269, 296270, 296273,
    296279, 296280, 296303, 296304, 296307, 296309, 296312, 296377, 296378, 296379, 
    296380, 296381, 296383, 296414, 296419, 296420, 296423, 296424, 296433, 296472, 
    296509, 296510, 296511, 296514, 296516, 296547, 296548, 296549, 296550, 296551, 
    296552, 296553, 296615, 296616, 296618, 296619, 296622, 296623
    );


char GoodRuns18r[3000];
sprintf(GoodRuns18r,  
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d ||"
   "fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d || fRunNum == %d",
    296690, 296691, 296694, 296749, 296750, 296781, 296784, 296785, 296786, 296787,
    296791, 296793, 296794, 296799, 296836, 296838, 296839, 296848, 296849, 296850,
    296851, 296852, 296890, 296894, 296899, 296900, 296903, 296930, 296931, 296932,
    296934, 296935, 296938, 296941, 296966, 296967, 296968, 296969, 296971, 296975,
    296976, 296979, 297029, 297031, 297035, 297085, 297117, 297118, 297119, 297123, 
    297124, 297128, 297129, 297132, 297133, 297193, 297194, 297196, 297218, 297219, 
    297221, 297222, 297278, 297310, 297312, 297315, 297317, 297363, 297366, 297367, 
    297372, 297379, 297380, 297405, 297408, 297413, 297414, 297415, 297441, 297442,
    297446, 297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540,
    297541, 297542, 297544, 297558, 297588, 297590, 297595
    );


// --------------------------------------------------------------------------------
  // import Coherent jpsi MC for PDF

 TFile *f1=new TFile("/mnt/data/Processed_data/Old/LHC16b2a-h_Train337_child_1.root"); 
 TTree *t1;
 f1->GetObject("jRecTree",t1); 

 RooDataSet *dataINCoh = new RooDataSet("dataINCoh", "dataINCoh", RooArgSet(jRecM,jRecY,jRecPt), Import(*t1));  // Read data from Tree

 RooAbsData *dataCoh = dataINCoh->reduce(Cut(strReduce),EventRange((t1->GetEntries())/2,(t1->GetEntries())));
 RooDataSet *MCdataCoh = (RooDataSet*) dataINCoh->reduce(EventRange(0,scale*8300));


 RooDataHist cohDH("cohDH","cohDH",jRecPt,*dataCoh,1);


 // --------------------------------------------------------------------------------
  // import Incoherent jpsi MC for PDF

 TFile *f2=new TFile("/mnt/data/Processed_data/Old/LHC16b2a-h_Train337_child_2.root"); 
 TTree *t2;
 f2->GetObject("jRecTree",t2); 

 RooDataSet *dataINIncoh = new RooDataSet("dataINIncoh", "dataINIncoh", RooArgSet(jRecM,jRecY,jRecPt), Import(*t2));  // Read data from Tree

 RooAbsData *dataIncoh = dataINIncoh->reduce(Cut(strReduce),EventRange((t2->GetEntries())/2,(t2->GetEntries())));
 RooDataSet *MCdataIncoh = (RooDataSet*) dataINIncoh->reduce(EventRange(0,scale*970));

 RooDataHist incohDH("incohDH","incohDH",jRecPt,*dataIncoh,1);


 // --------------------------------------------------------------------------------
  // import Coherent Psi2S MC for PDF

 TFile *f3=new TFile("/mnt/data/Processed_data/Old/LHC16b2a-h_Train337_child_3.root"); 
 TTree *t3;
 f3->GetObject("jRecTree",t3); 

 RooDataSet *dataINDirectCoh2 = new RooDataSet("dataINDirectCoh2", "dataINDirectCoh2", RooArgSet(jRecM,jRecY,jRecPt), Import(*t3));  // Read data from Tree

 RooAbsData *dataDirectCoh2 = dataINDirectCoh2->reduce(Cut(strReduce),EventRange((t3->GetEntries())/2,(t3->GetEntries())));
 RooDataSet *MCdataDirectCoh2 = (RooDataSet*) dataINDirectCoh2->reduce(EventRange(0,scale*970));

 RooDataHist DirectCoh2DH("DirectCoh2DH","DirectCoh2DH",jRecPt,*dataDirectCoh2,1);


 // --------------------------------------------------------------------------------
  // import Incoherent Psi2S MC for PDF

 TFile *f4=new TFile("/mnt/data/Processed_data/Old/LHC16b2a-h_Train337_child_4.root"); 
 TTree *t4;
 f4->GetObject("jRecTree",t4); 

 RooDataSet *dataINDirectIncoh2 = new RooDataSet("dataINDirectIncoh2", "dataINDirectIncoh2", RooArgSet(jRecM,jRecY,jRecPt), Import(*t4));  // Read data from Tree

 RooAbsData *dataDirectIncoh2 = dataINDirectIncoh2->reduce(Cut(strReduce),EventRange((t4->GetEntries())/2,(t4->GetEntries())));
 RooDataSet *MCdataDirectIncoh2 = (RooDataSet*) dataINDirectIncoh2->reduce(EventRange(0,scale*970));

 RooDataHist DirectIncoh2DH("DirectIncoh2DH","DirectIncoh2DH",jRecPt,*dataDirectIncoh2,1);


 // --------------------------------------------------------------------------------
  // import Coherent psi2s to jpsi MC for PDF

 TFile *f5=new TFile("/mnt/data/Processed_data/Old/LHC16b2a-h_Train337_child_5.root"); 
 TTree *t5;
 f5->GetObject("jRecTree",t5); 

 RooDataSet *dataINCoh2 = new RooDataSet("dataINCoh2", "dataINCoh2", RooArgSet(jRecM,jRecY,jRecPt), Import(*t5));  // Read data from Tree

 RooAbsData *dataCoh2 = dataINCoh2->reduce(Cut(strReduce),EventRange((t5->GetEntries())/2,(t5->GetEntries())));
 RooDataSet *MCdataCoh2 = (RooDataSet*) dataINCoh2->reduce(EventRange(0,scale*790));

 RooDataHist cohDH2("cohDH2","cohDH2",jRecPt,*dataCoh2,1);


 // --------------------------------------------------------------------------------
  // import Incoherent psi2s to jpsi MC for PDF

 TFile *f6=new TFile("/mnt/data/Processed_data/Old/LHC16b2a-h_Train337_child_6.root"); 
 TTree *t6;
 f6->GetObject("jRecTree",t6); 

 RooDataSet *dataINIncoh2 = new RooDataSet("dataINIncoh2", "dataINIncoh2", RooArgSet(jRecM,jRecY,jRecPt), Import(*t6));  // Read data from Tree

 RooAbsData *dataIncoh2 = dataINIncoh2->reduce(Cut(strReduce),EventRange((t6->GetEntries())/2,(t6->GetEntries())));
 RooDataSet *MCdataIncoh2 = (RooDataSet*) dataINIncoh2->reduce(EventRange(0,0));

 RooDataHist incohDH2("incohDH2","incohDH2",jRecPt,*dataIncoh2,1);


// --------------------------------------------------------------------------------
  // import gg from 2 // and 4 GeV MC for PDF - not using the 4 GeV

 TFile *f7=new TFile("/mnt/data/Processed_data/Old/LHC16b2a-h_Train337_child_7.root"); 
 TTree *t7;
 f7->GetObject("jRecTree",t7); 

 RooDataSet *dataINGam2GeV = new RooDataSet("dataINGam2GeV", "dataINGam2GeV", RooArgSet(jRecM,jRecY,jRecPt), Import(*t7));  // Read data from Tree

 // TFile *f8=new TFile("/home/alidock/Data_processing/Processed_data/LHC16b2a-h_Train337_child_8.root"); 
 // TTree *t8;
 // f8->GetObject("jRecTree",t8); 

 // RooDataSet *dataINGam4GeV = new RooDataSet("dataINGam4GeV", "dataINGam4GeV", RooArgSet(jRecM,jRecY,jRecPt), Import(*t8));  // Read data from Tree

 // dataINGam2GeV->append(*dataINGam4GeV); //Merging 2 GeV and 4GeV gamma gamma distributions leads to 4GeV peak in gg background

 RooAbsData *dataGam2GeV = dataINGam2GeV->reduce(Cut(strReduce),EventRange((t7->GetEntries())/2,(t7->GetEntries())));
 RooDataSet *MCdataGam2GeV = (RooDataSet*) dataINGam2GeV->reduce(EventRange(0,scale*7160));

 RooDataHist ggDH("ggDH","ggDH",jRecPt,*dataGam2GeV,1);  
    

// --------------------------------------------------------------------------------   
  // create templates
  // --> template coherent jpsi signal
  RooHistPdf cohPdf("cohPdf", "cohPdf", jRecPt,cohDH,0); 
  // --> template incoherent jspi signal
  RooHistPdf incohPdf("incohPdf", "incohPdf", jRecPt,incohDH,0); 
  // --> template Coherent psi2S signal
  RooHistPdf DirectcohPdf2("DirectcohPdf2", "DirectcohPdf2", jRecPt,DirectCoh2DH,0); 
  // --> template incoherent Psi2S signal
  RooHistPdf DirectincohPdf2("DirectincohPdf2", "DirectincohPdf2", jRecPt,DirectIncoh2DH,0); 
  // --> template coherent psi2s to jpsi signal
  RooHistPdf cohPdf2("cohPdf2", "cohPdf2", jRecPt,cohDH2,0); 
    // --> template incoherent psi2s to jspi signal
  RooHistPdf incohPdf2("incohPdf2", "incohPdf2", jRecPt,incohDH2,0);   
  // --> template gg
  RooHistPdf ggPdf("ggPdf", "ggPdf", jRecPt,ggDH,0);

 // --------------------------------------------------------------------------------
  // Create Incoherent + disocitation PDF

  RooRealVar b("b","b",1.79, 0, 1.8); // 1, -10,50)
  RooRealVar n("n","n",3.58, 0,8);   // 4, 0.01,50)
  b.setConstant(kTRUE);
  n.setConstant(kTRUE); 
  
  RooGenericPdf incoh_disPdf("incoh_disPdf","jRecPt*pow((1+pow(jRecPt,2)*b/n),-n)",RooArgSet(jRecPt, b, n)); 
  RooDataSet *MCINincoh_dis = incoh_disPdf.generate(jRecPt,scale*1926);
  RooDataSet *MCincoh_dis = (RooDataSet*) MCINincoh_dis->reduce("jRecPt<3");

 // --------------------------------------------------------------------------------
  // Create CB1 PDF

  RooRealVar mMC("mMC","massMC",3.09867,3.05,3.15) ;
  RooRealVar sigMC("#sigma_{MC}","resolution",0.0854761,0.01,0.2) ;
  RooRealVar n_cbMC("n_{cb,MC}","n_cbMC",93.1194,0,120) ;
  RooRealVar alphaMC("#alpha_{MC}","alphaMC",1.12159,0,10) ;

  mMC.setConstant(kTRUE);
  sigMC.setConstant(kTRUE);
  n_cbMC.setConstant(kTRUE);
  alphaMC.setConstant(kTRUE);
  
  RooCBShape cb1MC("cb1MC","Crystal Ball PDF",jRecM,mMC,sigMC,alphaMC,n_cbMC) ;
  RooDataSet *CB1MC = cb1MC.generate(jRecM,scaleMC*20998);
  
 // --------------------------------------------------------------------------------
  // Create CB2 PDF

  RooFormulaVar m2MC("m2MC","mMC+3.68609-3.096916",RooArgList(mMC));
  RooFormulaVar sig2MC("sig2MC","#sigma_{MC}*1.09",RooArgList(sigMC));
  
  RooCBShape cb2MC("cb2","Crystal Ball PDF",jRecM,m2MC,sig2MC,alphaMC,n_cbMC) ;
  RooDataSet *CB2MC = cb2MC.generate(jRecM,scaleMC*512);
  
 // --------------------------------------------------------------------------------
  // Create BKGD PDF

  RooRealVar lambdaMC("#lambda_{MC}","exponent",-0.934756,-5.,0.);
  RooRealVar a1MC("a_{1,MC}","parameter a1",0.648354,0,2);
  RooRealVar a2MC("a_{2,MC}","parameter a2",0.900351,0,2);
  RooRealVar a3MC("a_{3,MC}","parameter a3",0.239378,0,2);

  lambdaMC.setConstant(kTRUE);
  a1MC.setConstant(kTRUE);
  a2MC.setConstant(kTRUE);
  a3MC.setConstant(kTRUE);

  cout << lambdaMC.getVal() << endl;
  cout << a1MC.getVal() << endl;
  cout << a2MC.getVal() << endl;
  cout << a3MC.getVal() << endl;
  
  RooGenericPdf BkgdMC("BkgdMC","jRecM>4? exp(jRecM*#lambda_{MC}) : exp(jRecM*#lambda_{MC})*(1+a_{1,MC}*pow((jRecM-4),2)+a_{2,MC}*pow((jRecM-4),3)+a_{3,MC}*pow((jRecM-4),4))",RooArgSet(jRecM,lambdaMC,a1MC,a2MC,a3MC));
  RooDataSet *bgMC = BkgdMC.generate(jRecM,scaleMC*29008);
  
// --------------------------------------------------------------------------------
  // Import real data
 TFile *f18q=new TFile("/mnt/data/Processed_data/LHC18q_AnalysisResults_pass3_AOD_501.root"); 
 TTree *t18q;

 f18q->GetObject("NanoMUON/fRecTree",t18q);
 
 RooDataSet *dataINRaw18q = new RooDataSet("dataIN", "dataIN", RooArgSet(fMuMuM,fMuMuY,fMuMuPt, fRunNum, fADADecision, fADCDecision, fV0ADecision, fV0CFiredCells), Import(*t18q));  // Read data from Tree
 dataINRaw18q->changeObservableName("fMuMuM","jRecM");
 dataINRaw18q->changeObservableName("fMuMuY","jRecY");
 dataINRaw18q->changeObservableName("fMuMuPt","jRecPt");
 RooDataSet *dataIN18q = (RooDataSet*) dataINRaw18q->reduce(Cut(GoodRuns18q));
// --------------------------------------------------------------------------------
 // TFile *f18r=new TFile("/mnt/data/Processed_data/LHC18r_AnalysisResults_pass3_AOD_502.root"); 
 // TTree *t18r;

 // f18r->GetObject("NanoMUON/fRecTree",t18r);

 // RooDataSet *dataINRaw18r = new RooDataSet("dataIN", "dataIN", RooArgSet(fMuMuM,fMuMuY,fMuMuPt, fRunNum, fADADecision, fADCDecision, fV0ADecision, fV0CFiredCells), Import(*t18r));  // Read data from Tree
 // dataINRaw18r->changeObservableName("fMuMuM","jRecM");
 // dataINRaw18r->changeObservableName("fMuMuY","jRecY");
 // dataINRaw18r->changeObservableName("fMuMuPt","jRecPt");
 // RooDataSet *dataIN18r = (RooDataSet*) dataINRaw18r->reduce(Cut(GoodRuns18r));
 // --------------------------------------------------------------------------------
 // TFile *f15o=new TFile("/mnt/data/Processed_data/LHC15o_AnalysisResults.root"); 
 // TTree *t15o;

 // f15o->GetObject("NanoMUON/fAnaTree",t15o);
 

 // RooDataSet *dataINRaw15o = new RooDataSet("dataIN", "dataIN", RooArgSet(fMuMuM,fMuMuY,fMuMuPt, fRunNum, fADADecision, fADCDecision, fV0ADecision, fV0CFiredCells), Import(*t15o));  // Read data from Tree
 // dataINRaw15o->changeObservableName("fMuMuM","jRecM");
 // dataINRaw15o->changeObservableName("fMuMuY","jRecY");
 // dataINRaw15o->changeObservableName("fMuMuPt","jRecPt");
 // RooDataSet *dataIN15o = (RooDataSet*) dataINRaw15o->reduce(Cut(GoodRuns15o));

 //  // --------------------------------------------------------------------------------
 // TFile *f18l1=new TFile("/mnt/data/Processed_data/MCtestAnalysisResults.root"); 
 // TTree *t18l1MC;
 // TTree *t18l1ana;

 // f18l1->GetObject("NanoMUON/fMCTree",t18l1MC);
 // f18l1->GetObject("NanoMUON/fAnaTree",t18l1ana);
 

 // RooDataSet *dataINRaw18l1MC = new RooDataSet("dataIN", "dataIN", RooArgSet(fMuMuM,fMuMuY,fMuMuPt, fRunNum), Import(*t18l1MC));  // Read data from Tree
 // RooDataSet *dataINRaw18l1ana = new RooDataSet("dataIN", "dataIN", RooArgSet(fMuMuM,fMuMuY,fMuMuPt, fRunNum), Import(*t18l1MC));  // Read data from Tree
 // dataINRaw15o->changeObservableName("fMuMuM","jRecM");
 // dataINRaw15o->changeObservableName("fMuMuY","jRecY");
 // dataINRaw15o->changeObservableName("fMuMuPt","jRecPt");
 // RooDataSet *dataIN15o = (RooDataSet*) dataINRaw15o->reduce(Cut(GoodRuns15o));


 // Create MC data set
 RooDataSet* MCdataIN = new RooDataSet("MCdataIN", "MCdataIN", RooArgSet(jRecM,jRecY,jRecPt));  // Create MC data 
 RooDataSet* MCdataINM = new RooDataSet("MCdataINM", "MCdataINM", RooArgSet(jRecM,jRecY,jRecPt));  // Create MC data M
 RooDataSet* MCdataINM_pure_sig = new RooDataSet("MCdataINM_pure_sig", "MCdataINM_pure_sig", RooArgSet(jRecM,jRecY,jRecPt));  // Create MC data M pure signal
 
 MCdataIN->append(*MCdataCoh);
 MCdataIN->append(*MCdataIncoh);
 MCdataIN->append(*MCdataCoh2);
 // // MCdataIN->append(*MCdataIncoh2); 
 MCdataIN->append(*MCincoh_dis);
 MCdataIN->append(*MCdataGam2GeV);
// ------------------------------------
 MCdataINM->append(*MCdataCoh);
 MCdataINM->append(*MCdataIncoh);
 MCdataINM->append(*MCdataCoh2);
 // MCdataIN->append(*MCdataIncoh2); 
 MCdataINM->append(*MCdataGam2GeV);
 // MCdataINM->append(*MCincoh_dis);
// -------------------------------------
 MCdataINM_pure_sig->append(*MCdataCoh);
 MCdataINM_pure_sig->append(*MCdataIncoh);
 MCdataINM_pure_sig->append(*MCdataCoh2);
 
 // dataIN18q->append(*dataIN18r);
 // dataIN18q->append(*dataIN15o);

 CB1MC->append(*CB2MC);
 CB1MC->append(*bgMC);

 RooAbsData* data;
 RooAbsData* dataM;

  // Get data and do cuts  -- select between real and MC data
if (data_opt == 1) {
  data = dataIN18q->reduce(Cut(strReduceData));
  dataM = dataIN18q->reduce(Cut(strReduceMData));
}
// -------------------------------------------------------
else if (data_opt == 2) {
  data = MCdataIN->reduce(strReduce);
  dataM = MCdataINM->reduce(strReduceM);
  // dataM = CB1MC;
}
// -------------------------------------------------------

  RooAbsData* dataM_pure_sig = MCdataINM_pure_sig->reduce(strReduceM);
 



// // *********************************** Reconstruction efficiency ****************************************
//   Printf("Calculating reconstruction efficiency...");

//   TEfficiency *efficiencyRec  = new TEfficiency("efficiencyOMU ","efficiencyOMU ",1,0.,1.);

//   for (Int_t iEn=0;iEn<tInput->GetEntries();iEn++){
//     t->GetEntry(iEn);
//     // particle generated in wanted rapidity/pt2 interval
//     if (!NanoGenSelected()) continue;
//     // track was reconstructed
//     if (fMuMuPt==0.) {
//       efficiencyRec->Fill(kFALSE,0.5);
//     }
//     else {
//       efficiencyRec->Fill(kTRUE,0.5);
//     }
//   }

//   effRec = efficiencyRec->GetEfficiency(1);
//   errUpStatEffRec = efficiencyRec->GetEfficiencyErrorUp(1);
//   errLowStatEffRec = efficiencyRec->GetEfficiencyErrorLow(1);





  // #####################################################################################
  // Invariant mass fit
  // #####################################################################################
  // Extract parameters from Evgeny's file
  TFile* fin = new TFile("EfitM.root");
  TH1D* hM_1c = (TH1D*) fin->Get(Form("1c_%i_%i_1",YBin,0));
  // TH1D* hM_2c = (TH1D*) fin->Get(Form("2c_%i_%i_1",iy,ipt));
  TH1D* hM_gl = (TH1D*) fin->Get(Form("gl_%i_%i_1",YBin,0));


  TF1 *fBG;
  // TF1 *fCB1;
  fBG  = (TF1*) hM_gl->GetFunction("f2")->Clone("fBG");
  // fCB1 = (TF1*) hM_1c->GetFunction("fCB")->Clone("fCB1");
  // fCB2 = (TF1*) hM_2c->GetFunction("fCB")->Clone("fCB2");
  // Double_t nE = fCB1->GetParameter(0);
  Double_t bE = fBG->GetParameter(1);
  // Double_t m0 = fBG->GetParameter(2);
  Double_t a2E = fBG->GetParameter(3);
  Double_t a3E = fBG->GetParameter(4);
  Double_t a4E = fBG->GetParameter(5);

  // cout << "a2: " << a2E << endl;
  // cout << "a3: " << a3E << endl;
  // cout << "a4: " << a4E << endl;


  Int_t nEvts = dataM->numEntries();
  RooRealVar m0("m","mass",3.096,3.05,3.15) ;
  RooRealVar n_cb("n_{cb}","n_cb",4,0,50) ;
  RooRealVar sig("#sigma","resolution",0.092,0.01,0.2) ;
  RooRealVar alpha("#alpha","alpha",0.97,0,20) ;
  RooFormulaVar m02("m02","m+3.68609-3.096916",RooArgList(m0));
  RooFormulaVar sig2("sig2","#sigma*1.09",RooArgList(sig));

  // RooRealVar cb("cb","cb_new",2,5) ;
  // RooFormulaVar cb("cb","(((jRecM-m0)/sig)>=-alpha)? exp(-0.5*pow((jRecM-m0)/sig,2)) : exp(0.5*pow(alpha,2)+alpha*(jRecM-m0)/sig)",RooArgList(jRecM,m0,sig,alpha));     
  // RooFormulaVar cb("cb","((jRecM-m0)/sig)>=-alpha? exp(-0.5*pow((jRecM-m0)/sig,2)) : exp(-0.5*pow((jRecM-m0)/sig,2))",RooArgList(jRecM,m0,sig,alpha));     

  RooCBShape cb("cb","Crystal Ball PDF",jRecM,m0,sig,alpha,n_cb) ;
  RooCBShape cb2("cb2","Crystal Ball PDF",jRecM,m02,sig2,alpha,n_cb) ;


  RooRealVar Ncb("Ncb","Number of CB pure",0.2*nEvts,0,nEvts);
  RooAddPdf cb_pure("cb_pure","Crystal Ball pure PDF", RooArgList(cb), RooArgList(Ncb));

  // RooFitResult* rM_pure_sig = cb_pure.fitTo(*dataM_pure_sig,Extended(kTRUE),Range(2,6),Save());

  // n_cb.setConstant(kTRUE);
  
  // TCanvas *cM_pure = new TCanvas("cM_pure","cM_pure",900,600);
  // cM_pure->cd(1);
  // RooPlot* frameM_pure = jRecM.frame(Title("Mass fit")) ; 
  // dataM_pure_sig->plotOn(frameM_pure,Name("dataM_pure_sig"),Binning(binM),MarkerStyle(20),MarkerSize(1.));
  // cb_pure.plotOn(frameM_pure,Name("cb_pure"),LineColor(kBlue)) ;
  // frameM_pure->Draw();

  //  Double_t chi2M_pure = frameM_pure->chiSquare("cb_pure","dataM_pure_sig",rM_pure_sig->floatParsFinal().getSize());

  //  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  //  leg->SetFillStyle(0);
  //  leg->SetBorderSize(0);
  //  leg->SetHeader("Crystal Ball - Pure Signal","C"); 
  //  leg->AddEntry((TObject*)0,Form("Rapidity: %.1f<y<%.1f",y_min,y_max),"");   
  //  leg->AddEntry("dataM_pure_sig","Data", "P");
  //  leg->AddEntry("cb_pure","Crystal Ball","L");
  //  leg->AddEntry((TObject*)0,Form("#chi^{2}/NDF: %.3f",chi2M_pure),"");
  //  // leg->AddEntry((TObject*)0,Form("%s = %.0f+/-%.0f ",FitModel->GetParName(5),FitModel->GetParameter(5),FitModel->GetParError(5)),"");
  //  leg->AddEntry((TObject*)0,Form("N_{J/#psi} = %3.0f #pm %3.0f",Ncb.getVal(),Ncb.getError()),"");
  //  leg->AddEntry((TObject*)0,Form("# Entries = %.0d",dataM_pure_sig->numEntries()),"");
  //  leg->SetTextSize(0.027);
  //  leg->Draw();

  // create a PDF to describe the background
  RooRealVar lambda("#lambda","exponent",-bE,-5.,0.);
  // lambda.setConstant(kTRUE);

  RooRealVar a1("a_{1}","parameter a1",a2E,0.0001,2);
  RooRealVar a2("a_{2}","parameter a2",a3E,0.0001,2);
  RooRealVar a3("a_{3}","parameter a3",a4E,0.0001,2);
  a1.setConstant(kTRUE);
  a2.setConstant(kTRUE);
  a3.setConstant(kTRUE);

  // RooFormulaVar Bkgd("Bkgd","(jRecM<=4)*exp(jRecM*lambda)*(1+a1*pow((jRecM-4),2)+a2*pow((jRecM-4),3)+a3*pow((jRecM-4),4))+(jRecM>4)*exp(jRecM*lambda)",RooArgList(jRecM,lambda,a1,a2,a3));
  // RooFormulaVar Bkgd("Bkgd","exp(jRecM*lambda)",RooArgList(jRecM,lambda));
  // RooFormulaVar Bkgd("Bkgd","exp(jRecM*lambda)*(1+a1*pow((jRecM-4),2)+a2*pow((jRecM-4),3)+a3*pow((jRecM-4),4))",RooArgList(jRecM,lambda,a1,a2,a3));
  RooGenericPdf Bkgd("Bkgd","jRecM>4? exp(jRecM*#lambda) : exp(jRecM*#lambda)*(1+a_{1}*pow((jRecM-4),2)+a_{2}*pow((jRecM-4),3)+a_{3}*pow((jRecM-4),4))",RooArgSet(jRecM,lambda,a1,a2,a3));

  
  RooRealVar nCB("N_{J/#psi}","Number of CB events",0.5*nEvts,0,nEvts);
  RooRealVar nCB2("N_{#psi'}","Number of CB events",0.1*nEvts,0,nEvts);
  RooRealVar nBG("N_{bg}","Number of BG events",0.6*nEvts,0,nEvts);
  RooAddPdf CbBG("CbBG","Crystal Ball plus Background PDF", RooArgList(Bkgd,cb,cb2), RooArgList(nBG,nCB,nCB2));

  RooFitResult* rM = CbBG.fitTo(*dataM,Extended(kTRUE),Range(2,6),Save());

  // compute background in range
  Double_t nBkgd[2];
  Double_t nBkgd0[2];
  Double_t nCBM[2];
  Double_t nCBM2[2];
  jRecM.setRange("JPsiMassRange",mMin,mMax);
  jRecM.setRange("JPsiMassRange0",mMin0,mMax0);
  RooAbsReal *iBG = Bkgd.createIntegral(jRecM,NormSet(jRecM),Range("JPsiMassRange"));
  RooAbsReal *iBG0 = Bkgd.createIntegral(jRecM,NormSet(jRecM),Range("JPsiMassRange0"));
  RooAbsReal *iCB = cb.createIntegral(jRecM,NormSet(jRecM),Range(2.,6.));
  RooAbsReal *iCB2 = cb2.createIntegral(jRecM,NormSet(jRecM),Range(2.,6.));

  nCBM[0] = iCB->getVal()*nCB.getVal();
  nCBM[1] =  iCB->getVal()*nCB.getError();
  nCBM2[0] = iCB2->getVal()*nCB2.getVal();
  nCBM2[1] =  iCB2->getVal()*nCB2.getError();
  nBkgd[0] =  iBG->getVal()*nBG.getVal();
  nBkgd[1] =  iBG->getVal()*nBG.getError();
  nBkgd0[0] =  iBG0->getVal()*nBG.getVal();
  nBkgd0[1] =  iBG0->getVal()*nBG.getError();  

  Double_t ggN_min=(nBkgd[0]);//-3*nBkgd[1];
  Double_t ggN_max=(nBkgd[0]);//+3*nBkgd[1];


  // cout << "Number of coh entries: " << (MCdataCoh->reduce(strReduce))->numEntries() << endl;
  // cout << "Number of incoh entries: " << (MCdataIncoh->reduce(strReduce))->numEntries() << endl;
  // cout << "Number of coh2 entries: " << (MCdataCoh2->reduce(strReduce))->numEntries() << endl;  
  // cout << "Number of gg entries: " << (MCdataGam2GeV->reduce(strReduce))->numEntries() << endl;
  // cout << "Number of disoc entries: " << (MCincoh_dis->reduce("jRecPt<3"))->numEntries() << endl;

// #####################################################################################
// Transverse momentum fit
// #####################################################################################

// --------------------------------------------------------------------------------
  // create the model
  RooRealVar cohN("J#psi_{coh}","number of coherent jpsi signal",0.75*data->numEntries(),0,data->numEntries());   
  RooRealVar incohN("J#psi_{incoh}","number of incoherent jpsi signal",0.1*data->numEntries(),0*data->numEntries(),data->numEntries()); //800 min to force 
  RooRealVar DirectcohN2("Direct #psi(2S)_{coh}","number of direct coherent psi2s signal",0.75*data->numEntries(),0,data->numEntries());   
  RooRealVar DirectincohN2("Direct #psi(2S)_{incoh}","number of direct incoherent psi2s signal",0.1*data->numEntries(),0*data->numEntries(),data->numEntries()); 
  RooFormulaVar cohN2("#psi'_{coh}","J#psi_{coh}*0.085",RooArgList(cohN));
  RooFormulaVar incohN2("#psi'_{incoh}","J#psi_{incoh}*0.085",RooArgList(incohN));
  // RooRealVar cohN2("#psi(2S)_{coh}","number of coherent psi2s signal",0.1*data->numEntries(),0,0.5*data->numEntries()); //0.1*data->numEntries(),0,0.5*data->numEntries() 
  // RooRealVar incohN2("#psi(2S)_{incoh}","number of incoherent psi2s signal",0.01*data->numEntries(),0,0.5*data->numEntries()); //0.01*data->numEntries(),0,0.5*data->numEntries()
  RooRealVar ggN("#gamma#gamma","number of gg", nBkgd0[0],nBkgd0[0],nBkgd0[0]); //0.3*data->numEntries(),0,data->numEntries()
  RooRealVar incoh_disN("N_{disoc}","number of incoherent jpsi disociated signal",0.1*data->numEntries(),0,data->numEntries());
  ggN.setConstant(kTRUE);

  // Create the model as the sum of the templates
  // --- For J/Psi
    RooAddPdf Model("Model","Sum of background and coherent signal",RooArgList(cohPdf,incohPdf, cohPdf2, incohPdf2, ggPdf, incoh_disPdf), RooArgList(cohN,incohN, cohN2, incohN2, ggN, incoh_disN));
  // --- For Psi(2S)
    // RooAddPdf Model("Model","Sum of background and coherent signal",RooArgList(DirectcohPdf2, DirectincohPdf2, ggPdf, incoh_disPdf), RooArgList(DirectcohN2, DirectincohN2, ggN, incoh_disN));
  // RooAddPdf Model("Model","Sum of background and coherent signal",RooArgList(cohPdf,incohPdf, cohPdf2,  ggPdf, incoh_disPdf), RooArgList(cohN,incohN, cohN2, ggN, incoh_disN));
  // RooAddPdf Model("Model","Sum of background and coherent signal",RooArgList(cohPdf,incohPdf, cohPdf2, ggPdf), RooArgList(cohN,incohN, cohN2, ggN));

// --------------------------------------------------------------------------------
  // perform fit
  RooFitResult* r = Model.fitTo(*data,Extended(kTRUE),Save(),Range(0,ptMax));
 

 // -------------------------------------------------------------------------------- 
  // plot results
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  // gStyle->SetPadTickX(1);
  // gStyle->SetPadTickY(1);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetFrameLineWidth(1);

  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTextSize(0.04);


  TCanvas *cPt = new TCanvas("cPt","cPt",1600,800);
  cPt->Divide(2);  

  //-----------------------------------------------------------------------------------
// Draw Correlation Matrix
   cPt->cd(2); 
   TPad *cPt_2 = new TPad("cPt_2", "cPt_2",0.001,0.001,0.999,0.999);
   
   cPt_2->SetRightMargin(0.15);
   cPt_2->SetLeftMargin(0.15);
   cPt_2->Draw();
   cPt_2->cd();

   TH2* hcorr = r->correlationHist();
   hcorr->SetMarkerSize(1.6);
   hcorr->GetXaxis()->SetLabelSize(0.045);
   hcorr->GetYaxis()->SetLabelSize(0.045);
   hcorr->Draw("zcol,text");

   //-----------------------------------------------------------------------------------
// Draw Histogram  
  cPt->cd(1);
  TPad *cPt_1 = new TPad("cPt_1", "cPt_1",0.001,0.001,0.999,0.999);
  
  cPt_1->SetLogy();
  cPt_1->SetLeftMargin(0.14);
  cPt_1->SetRightMargin(0.01);
  cPt_1->SetBottomMargin(0.12);
  cPt_1->Draw(); 
  cPt_1->cd();

  RooPlot* frame = jRecPt.frame(Title("Pt fit")) ;  
  
  data->plotOn(frame,Name("data"),Binning(bin),MarkerStyle(20),MarkerSize(0.9));
  Model.plotOn(frame,Name("Model"),LineColor(kBlack), LineWidth(2)) ;
  // Model.paramOn(frame,data);
  // ------Psi(2S) 
  // Model.plotOn(frame,Name("DirectcohPdf2"), Components(DirectcohPdf2), LineColor(kBlue+1), LineWidth(2));
  // Model.plotOn(frame,Name("DirectincohPdf2"),Components(DirectincohPdf2), LineColor(kRed+1), LineWidth(2));
  // ------J/Psi
  Model.plotOn(frame,Name("cohPdf"), Components(cohPdf), LineColor(kBlue+1), LineWidth(2));
  Model.plotOn(frame,Name("incohPdf"),Components(incohPdf), LineColor(kRed+1), LineWidth(2));
  Model.plotOn(frame,Name("cohPdf2"),Components(cohPdf2), LineColor(kCyan+1), LineWidth(2));
  Model.plotOn(frame,Name("incohPdf2"),Components(incohPdf2), LineColor(kOrange), LineWidth(2));
  
  Model.plotOn(frame,Name("incoh_disPdf"),Components(incoh_disPdf), LineColor(kMagenta+1), LineWidth(2));
  Model.plotOn(frame,Name("ggPdf"),Components(ggPdf), LineColor(kGreen+2), LineWidth(2));

  frame->SetAxisRange(0,3,"X");
  frame->SetAxisRange(1,110000,"Y");
  frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  frame->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c} ", (ptMax-ptMin)/n_bins*1000));
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->Draw();


  Double_t ncoh[2];
  Double_t ncoh2[2];
  Double_t nincoh[2];
  Double_t nincoh2[2];
  Double_t ndisoc[2];  
  jRecPt.setRange("JPsiptRange",0,ptcut);
  
  RooAbsReal *icoh = cohPdf.createIntegral(jRecPt,NormSet(jRecPt),Range("JPsiptRange"));
  RooAbsReal *icoh2 = cohPdf2.createIntegral(jRecPt,NormSet(jRecPt),Range("JPsiptRange"));
  RooAbsReal *iincoh = incohPdf.createIntegral(jRecPt,NormSet(jRecPt),Range("JPsiptRange"));
  RooAbsReal *iincoh2 = incohPdf2.createIntegral(jRecPt,NormSet(jRecPt),Range("JPsiptRange"));
  RooAbsReal *idisoc = incoh_disPdf.createIntegral(jRecPt,NormSet(jRecPt),Range("JPsiptRange"));

  ncoh[0] = icoh->getVal()*cohN.getVal();
  ncoh[1] = icoh->getVal()*cohN.getError();
  ncoh2[0] = icoh2->getVal()*cohN2.getVal();

  nincoh[0] = iincoh->getVal()*incohN.getVal();
  nincoh[1] = iincoh->getVal()*incohN.getError();
  nincoh2[0] = iincoh2->getVal()*incohN2.getVal();

  ndisoc[0] = idisoc->getVal()*incoh_disN.getVal();
  ndisoc[1] = idisoc->getVal()*incoh_disN.getError();

//----------------------------------------------------------------------
// Calculate chi2
  Double_t chi2 = frame->chiSquare("Model","data",r->floatParsFinal().getSize()); 

  TLatex * text1 = new TLatex (0.1,0.45*frame->GetMaximum(),"This thesis");
  text1->Draw();
  if (data_opt == 1) {
	  TLatex * text2 = new TLatex (0.85,0.45*frame->GetMaximum(),"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
	  text2->Draw();
	}	  
  if (data_opt == 2) {
  	  TLatex * text2 = new TLatex (0.85,0.45*frame->GetMaximum(),"#bf{ALICE MC, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
	  text2->Draw();
  }
  TLatex * text3 = new TLatex (0.85,0.22*frame->GetMaximum(),Form("#bf{%.1f < y < %.1f, %.2f < #it{m_{#mu^{+}#mu^{-}}} < %.2f}",y_min,y_max, mMin0, mMax0));
  text3->Draw();
  TLatex * text4 = new TLatex (0.1,0.11*frame->GetMaximum(),Form("#bf{f_{I} = %.1f%% #pm %.1f%%}",100*(nincoh[0]+nincoh2[0]+ndisoc[0])/(ncoh[0]+ncoh2[0]), 100*sqrt( (nincoh[1]/nincoh[0])*(nincoh[1]/nincoh[0])+(ndisoc[1]/ndisoc[0])*(ndisoc[1]/ndisoc[0])+(ncoh[1]/ncoh[0])*(ncoh[1]/ncoh[0])  ) * (nincoh[0]+nincoh2[0]+ndisoc[0])/(ncoh[0]+ncoh2[0])));
  text4->Draw();



  TLegend *leg1 = new TLegend(0.53,0.40,0.95,0.79);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.04);
  leg1->AddEntry("data","Data", "P");
  leg1->AddEntry("Model","Model","L");
  // ------Psi(2S)
  // leg1->AddEntry("DirectcohPdf2","Coherent #psi(2S)", "L");
  // leg1->AddEntry("DirectincohPdf2","Incoherent #psi(2S)", "L");
  // ------J/Psi
  leg1->AddEntry("cohPdf","Coherent J/#psi", "L");
  leg1->AddEntry("incohPdf","Incoherent J/#psi", "L");
  leg1->AddEntry("cohPdf2","Coh. #psi' #rightarrow  J/#psi", "L");
  leg1->AddEntry("incohPdf2","Incoh. #psi' #rightarrow  J/#psi", "L");
  
  leg1->AddEntry("incoh_disPdf","Nucleon disoc.", "L");
  leg1->AddEntry("ggPdf","#gamma#gamma #rightarrow #mu#mu", "L");
  leg1->AddEntry((TObject*)0,Form("#chi^{2}/NDF: %.3f",chi2),"");
  // TLegendEntry* l1 = leg1->AddEntry((TObject*)0,Form("This thesis"),""); 
  // l1->SetTextColor(922);

  leg1->Draw();

  // save as pdf
  char ptFit[120];
  sprintf(ptFit,"Results/OLD_pass3_Pt_fit_%.2f_y_%.2f.pdf",y_min,y_max);
  // sprintf(ptFit,"Results/diploma/Pt_fit_%.2f_y_%.2f.pdf",y_min,y_max);
  cPt->Print(ptFit);

  // --------------------------------------------------------------------------------


TCanvas *cM = new TCanvas("cM","cM",1600,800);
cM->Divide(2);  

//-----------------------------------------------------------------------------------
// Draw Correlation Matrix
 cM->cd(2); 
 TPad *cM_2 = new TPad("cM_2", "cM_2",0.001,0.001,0.999,0.999);
 
 cM_2->SetRightMargin(0.15);
 cM_2->SetLeftMargin(0.15);
 cM_2->Draw();
 cM_2->cd();

 TH2* hcorrM = rM->correlationHist();
 hcorrM->SetMarkerSize(1.2);
 hcorrM->GetXaxis()->SetLabelSize(0.045);
 hcorrM->GetYaxis()->SetLabelSize(0.045);
 hcorrM->Draw("zcol,text");

//-----------------------------------------------------------------------------------
// Draw Histogram  
cM->cd(1);
TPad *cM_1 = new TPad("cM_1", "cM_1",0.001,0.001,0.999,0.999);
cM_1->SetLeftMargin(0.15);
cM_1->SetRightMargin(0.01);
cM_1->SetBottomMargin(0.12);
cM_1->Draw(); 
cM_1->cd();

RooPlot* frameM = jRecM.frame(Title("Mass fit")) ; 
dataM->plotOn(frameM,Name("dataM"),Binning(binM),MarkerStyle(20),MarkerSize(0.5));
CbBG.plotOn(frameM,Name("CbBG"),LineColor(kBlack), LineWidth(1)) ;
CbBG.plotOn(frameM,Name("cb"), Components(cb), LineColor(kRed+1),LineWidth(1));
CbBG.plotOn(frameM,Name("cb2"), Components(cb2), LineColor(kGreen+2), LineWidth(1));
CbBG.plotOn(frameM,Name("Bkgd"), Components(Bkgd), LineColor(kBlue+1), LineWidth(1));
frameM->GetXaxis()->SetTitle("#it{m_{#mu^{+}#mu^{-}}} (GeV/#it{c}^{2})");
frameM->GetYaxis()->SetTitle(Form("Counts per  %.0f MeV/#it{c^{2}} ", (mrangeMax-mrangeMin)/n_binsM*1000));
frameM->GetYaxis()->SetTitleOffset(1.6);
frameM->SetAxisRange(0,nCBM[0]/3,"Y");
frameM->Draw();

Double_t chi2M = frameM->chiSquare("CbBG","dataM",rM->floatParsFinal().getSize());


TLatex * text1M = new TLatex (2.2,0.9*frameM->GetMaximum(),"This thesis");
text1M->Draw();
if (data_opt == 1) {
TLatex * text2M = new TLatex (3.5,0.9*frameM->GetMaximum(),"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
text2M->Draw();
}
if (data_opt == 2) {
TLatex * text2M = new TLatex (3.35,0.9*frameM->GetMaximum(),"#bf{ALICE MC, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
text2M->Draw();
}	
TLatex * text3M = new TLatex (3.5,0.84*frameM->GetMaximum(),Form("#bf{%.2f < y < %.2f, #it{p}_{T} < %.2f}",y_min,y_max, ptcut));
text3M->Draw();
TLatex * text4M = new TLatex (2.2,0.79*frameM->GetMaximum(),Form("#bf{f_{D} = %.1f%% #pm %.1f%%}",100*nCBM2[0]/nCBM[0]*0.05961*0.12000/(0.00800*0.15800-0.61400*0.07200*0.05961*nCBM2[0]/nCBM[0])*0.61400*0.07200/0.12000, 100*sqrt( (nCBM2[1]/nCBM2[0])*(nCBM2[1]/nCBM2[0])+(nCBM[1]/nCBM[0])*(nCBM[1]/nCBM[0]) ) * nCBM2[0]/nCBM[0]*0.05961*0.12000/(0.00800*0.15800-0.61400*0.07200*0.05961*nCBM2[0]/nCBM[0])*0.61400*0.07200/0.12000));
text4M->Draw();




TLegend *leg2 = new TLegend(0.41,0.33,0.85,0.75);
leg2->SetFillStyle(0);
leg2->SetBorderSize(0);
leg2->SetTextSize(0.038);
// leg2->AddEntry("dataM","Data", "P");
// leg2->AddEntry("CbBG","Model","L");
// leg2->AddEntry("cb","Signal J/#psi","L");
// leg2->AddEntry("cb2","Signal #psi(2S)","L");
// leg2->AddEntry("Bkgd","Background","L");
leg2->AddEntry((TObject*)0,Form("#chi^{2}/NDF: %.3f",chi2M),"");
// leg->AddEntry((TObject*)0,Form("%s = %.0f+/-%.0f ",FitModel->GetParName(5),FitModel->GetParameter(5),FitModel->GetParError(5)),"");
leg2->AddEntry((TObject*)0,Form("N_{J/#psi} = %3.0f #pm %3.0f",nCBM[0],nCBM[1]),"");
leg2->AddEntry((TObject*)0,Form("N_{#psi'} = %3.0f #pm %3.0f",nCBM2[0],nCBM2[1]),"");
leg2->AddEntry((TObject*)0,Form("N_{#psi'}/N_{J/#psi} = %3.4f #pm %3.4f",nCBM2[0]/nCBM[0],sqrt( ( (nCBM2[1]/nCBM[0])*(nCBM2[1]/nCBM[0]) ) +
                                                                                    ( (nCBM2[0]/nCBM[0]*nCBM[1]/nCBM[0])*(nCBM2[0]/nCBM[0]*nCBM[1]/nCBM[0]) ) -
                                                                                    ( (nCBM2[1]/nCBM[0])*(nCBM2[0]/nCBM[0]*nCBM[1]/nCBM[0]) ) )),"");
leg2->AddEntry((TObject*)0,Form("N_{bg}(2,6) = %3.0f #pm %3.0f", nBG.getVal(),nBG.getError()),"");
leg2->AddEntry((TObject*)0,Form("N_{bg}(%3.2f,%3.2f) = %3.0f #pm %3.0f", mMin0, mMax0, nBkgd0[0],nBkgd0[1]),"");
leg2->AddEntry((TObject*)0,Form("N_{bg}(%3.2f,%3.2f) = %3.0f #pm %3.0f", mMin, mMax, nBkgd[0],nBkgd[1]),"");

leg2->AddEntry((TObject*)0,Form("# Entries = %.0i",dataM->numEntries()),"");
leg2->Draw();

// save as pdf
char massFit[120];
sprintf(massFit,"Results/OLD_pass3_Mass_fit_%.2f_y_%.2f.pdf",y_min,y_max);
// sprintf(massFit,"Results/diploma/Mass_fit_%.2f_y_%.2f.pdf",y_min,y_max);
cM->Print(massFit);


// RooFormulaVar fd("f_{D}","J#psi_{incoh}/J#psi_{coh}*0.05961*0.12000/(0.00800*0.15800-0.61400*0.07200*0.05961*J#psi_{incoh}/J#psi_{coh})*0.61400*0.07200/0.12000",RooArgList(incohN, cohN));


  cout << "###################################################################" << endl;
  cout << "###################################################################" << endl;
  cout << "Number of all entries: " << data->numEntries() << endl;
  cout << "Number of coh entries: " << (MCdataCoh->reduce(strReduce))->numEntries() << endl;
  cout << "Number of incoh entries: " << (MCdataIncoh->reduce(strReduce))->numEntries() << endl;
  cout << "Number of coh2 entries: " << (MCdataCoh2->reduce(strReduce))->numEntries() << endl;  
  cout << "Number of gg entries: " << (MCdataGam2GeV->reduce(strReduce))->numEntries() << endl;
  cout << "Number of disoc entries: " << (MCincoh_dis->reduce("jRecPt<3"))->numEntries() << endl;
  cout << "cohN: " << (Int_t) cohN.getVal() << " \\pm " << (Int_t) cohN.getError() << endl;  
  cout << "incohN: " << (Int_t) incohN.getVal() << " \\pm " << (Int_t) incohN.getError() << endl;
  cout << "cohN2: " << (Int_t) cohN2.getVal() << " \\pm "   << endl;  
  cout << "incohN2: " << (Int_t) incohN2.getVal() << " \\pm " << endl;
  cout << "ggN: " << (Int_t) ggN.getVal() << " \\pm " << (Int_t) ggN.getError() << endl;
  cout << "incoh_disN: " << (Int_t) incoh_disN.getVal() << " \\pm " << (Int_t) incoh_disN.getError() << endl;
  cout << "b: " << (Double_t) b.getVal() << " \\pm " << (Double_t) b.getError() << endl;
  cout << "n: " << (Double_t) n.getVal() << " \\pm " << (Double_t) n.getError() << endl;
  cout << "Chi2: " << chi2 << endl ;
  cout << "---------------------------------------------" << endl;
  cout << "ggN_range: " << ggN_min << " - " << ggN_max << " , " << (nBkgd[0]) << " +/- " << nBkgd[1] << endl;
  cout << "J/Psi signal from invariant mass fit " << nCBM[0] << " +/- " << nCBM[1] << endl;
  cout << "---------------------------------------------" << endl;
  cout << "#alpha: " << (Double_t) alpha.getVal() << " \\pm " << (Double_t) alpha.getError() << endl;   
  cout << "n_cb: " << (Double_t) n_cb.getVal() << " \\pm " << (Double_t) n_cb.getError() << endl; 
  cout << "m: " << (Double_t) m0.getVal() << " \\pm " << (Double_t) m0.getError() << endl;
  cout << "#sigma: " << (Double_t) sig.getVal() << " \\pm " << (Double_t) sig.getError() << endl;
  cout << "nCB2: " << (Int_t) nCB2.getVal() << " \\pm " << (Int_t) nCB2.getError() << endl;
  cout << "nCB: " << (Int_t) nCB.getVal() << " \\pm " << (Int_t) nCB.getError() << endl;
  cout << "nBG: " << (Int_t) nBG.getVal() << " \\pm " << (Int_t) nBG.getError() << endl;
  cout << "#lambda: " << (Double_t) lambda.getVal() << " \\pm " << (Double_t) lambda.getError() << endl; 
  cout << "a_2: " << (Double_t) a1.getVal() << " \\pm " << (Double_t) a1.getError() << endl; 
  cout << "a_3: " << (Double_t) a2.getVal() << " \\pm " << (Double_t) a2.getError() << endl; 
  cout << "a_4: " << (Double_t) a3.getVal() << " \\pm " << (Double_t) a3.getError() << endl;  
  cout << "###################################################################" << endl;
  cout << "###################################################################" << endl;
  cout << "f_D: " << nCBM2[0]/nCBM[0]*0.05961*0.12000/(0.00800*0.15800-0.61400*0.07200*0.05961*nCBM2[0]/nCBM[0])*0.61400*0.07200/0.12000 << " \\pm " << sqrt( (nCBM2[1]/nCBM2[0])*(nCBM2[1]/nCBM2[0])+(nCBM[1]/nCBM[0])*(nCBM[1]/nCBM[0]) ) * nCBM2[0]/nCBM[0]*0.05961*0.12000/(0.00800*0.15800-0.61400*0.07200*0.05961*nCBM2[0]/nCBM[0])*0.61400*0.07200/0.12000<< endl;  
  cout << "f_I: " << (nincoh[0]+nincoh2[0]+ndisoc[0])/(ncoh[0]+ncoh2[0]) << " \\pm " << sqrt( (nincoh[1]/nincoh[0])*(nincoh[1]/nincoh[0])+(ndisoc[1]/ndisoc[0])*(ndisoc[1]/ndisoc[0])+(ncoh[1]/ncoh[0])*(ncoh[1]/ncoh[0])  ) * (nincoh[0]+nincoh2[0]+ndisoc[0])/(ncoh[0]+ncoh2[0]) << endl;  



	 
}

void OldAnalysisJPsi(){

 DoFit(gYMin[YBin], gYMax[YBin]);

}
