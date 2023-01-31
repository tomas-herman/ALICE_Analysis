/*
Runs for forward and central analyses.
Also runs with TRD.
*/


// runs to be used for the forward analyses

const Int_t nGoodRunsLHC16r_fwd = 57;
Int_t GoodRunsLHC16r_fwd[nGoodRunsLHC16r_fwd] = {
  265589, 265594, 265596, 265607, 265669, 265691, 265694, 
  265697, 265698, 265700, 265701, 265709, 265713, 265714, 265740,
  265741, 265742, 265744, 265746, 265754, 265756, 265785, 265787,
  265788, 265789, 265792, 265795, 265797, 265840, 266022, 266023,
  266025, 266034, 266074, 266076, 266081, 266084, 266085, 266086,
  266117, 266187, 266189, 266190, 266193, 266196, 266197, 266208,
  266234, 266235, 266296, 266299, 266300, 266304, 266305, 266312,
  266316, 266318};

const Int_t nGoodRunsLHC16s_fwd = 76;
Int_t GoodRunsLHC16s_fwd[nGoodRunsLHC16s_fwd] = {
    266439, 266441, 266472, 266479,
  266480, 266487, 266514, 266516, 266518, 266520, 266522, 266523,
  266525, 266533, 266534, 266539, 266543, 266549, 266584, 266587,
  266588, 266591, 266593, 266595, 266613, 266614, 266618,
  266621, 266630, 266657, 266658, 266659, 266665, 266668, 266669,
  266674, 266676, 266702, 266703, 266706, 266708, 266775, 266776,
  266800, 266805, 266807, 266857, 266878, 266880, 266882, 266883,
  266885, 266886, 266912, 266915, 266940, 266942, 266943, 266944,
  266988, 266993, 266994, 266997, 266998, 267020, 267022, 267062,
  267063, 267067, 267070, 267072, 267077, 267109, 267110, 267130,
  267131};

// runs to be used for the mid and semi-forward analyses
const Int_t nGoodRunsLHC16r_mid = 45;
Int_t GoodRunsLHC16r_mid[nGoodRunsLHC16r_mid] = {
  266318, 266317, 266316, 266305, 266304, 266300, 266299, 266296,
  266208, 266197, 266196, 266193, 266190, 266189, 266187, 266117,
  266086, 266085, 266084, 266083, 266081, 266076, 266074, 266034,
  265797, 265795, 265789, 265788, 265756, 265754, 265746, 265744,
  265742, 265741, 265714, 265713, 265709, 265705, 265701, 265700,
  265698, 265697, 265607, 265596, 265594};

const Int_t nGoodRunsLHC16s_mid = 24;
Int_t GoodRunsLHC16s_mid[nGoodRunsLHC16s_mid] = {
  267110, 267081, 267077, 267072, 267070, 267030, 266998, 266997,
  266994, 266993, 266944, 266886, 266885, 266883, 266882, 266881,
  266549, 266543, 266492, 266487, 266441, 266439, 266438, 266437};

// run with TRD 
const Int_t nTRDLHC16rRun = 45;
Int_t TRDLHC16rRun[nTRDLHC16rRun] = {
  265589, 265594, 265596, 265607, 265669,  265697, 265698, 265700,
  265701, 265709, 265713, 265714, 265741, 265742, 265744, 265746, 265754,
  265756, 265788, 265789, 265792, 265795, 265797, 266034, 266074, 266076,
  266081, 266084, 266085, 266086, 266117, 266187, 266189, 266190, 266193,
  266196, 266197, 266208, 266296, 266299, 266300, 266304, 266305, 266316,
  266318};

const Int_t nTRDLHC16sRun = 19;
Int_t TRDLHC16sRun[nTRDLHC16sRun] = {
  266439, 266441, 266487, 266543, 266549, 266630,
  266882, 266883, 266885, 266886, 266944, 266993, 266994, 266997, 266998,
  267070, 267072, 267077, 267110};


const Int_t FirstRun = 265115;
const Int_t LastRun = 267166;
const Int_t nRuns = LastRun-FirstRun+1;
Int_t GoodRuns[nRuns];
Int_t GoodRunsIdx[nRuns];
Int_t TRDRunsIdx[nRuns];

void SetTRDRunsIdx(Int_t per)
{
  for(Int_t i=0;i<nRuns;i++) TRDRunsIdx[i] = -1;
  if (per == 0) {
    for(Int_t i=0;i<nTRDLHC16rRun;i++)
      TRDRunsIdx[TRDLHC16rRun[i]-FirstRun] = i;
  } else if (per==1) {
    for(Int_t i=0;i<nTRDLHC16sRun;i++)
    TRDRunsIdx[TRDLHC16sRun[i]-FirstRun] = i;
  }
}

void SetGoodRunsIdx_fwd(Int_t per)
{
  for(Int_t i=0;i<nRuns;i++) GoodRunsIdx[i] = -1;
  if (per == 0) {
  for(Int_t i=0;i<nGoodRunsLHC16r_fwd;i++)
    GoodRunsIdx[GoodRunsLHC16r_fwd[i]-FirstRun] = i;
  } else if (per==1) {
  for(Int_t i=0;i<nGoodRunsLHC16s_fwd;i++)
    GoodRunsIdx[GoodRunsLHC16s_fwd[i]-FirstRun] = i;
  }
}

void SetGoodRuns_fwd(Int_t per)
{
  for(Int_t i=0;i<nRuns;i++) GoodRuns[i] = 0;
  if (per == 0) {
    for(Int_t i=0;i<nGoodRunsLHC16r_fwd;i++)
      GoodRuns[GoodRunsLHC16r_fwd[i]-FirstRun] = GoodRunsLHC16r_fwd[i];
  } else if (per==1) {
    for(Int_t i=0;i<nGoodRunsLHC16s_fwd;i++)
      GoodRuns[GoodRunsLHC16s_fwd[i]-FirstRun] = GoodRunsLHC16r_fwd[i];
  }
}

void SetGoodRuns_mid()
{
  for(Int_t i=0;i<nRuns;i++) GoodRuns[i] = 0;
  for(Int_t i=0;i<nGoodRunsLHC16r_mid;i++)
    GoodRuns[GoodRunsLHC16r_mid[i]-FirstRun] = GoodRunsLHC16r_mid[i];
  for(Int_t i=0;i<nGoodRunsLHC16s_mid;i++)
    GoodRuns[GoodRunsLHC16s_mid[i]-FirstRun] = GoodRunsLHC16r_mid[i];
}
