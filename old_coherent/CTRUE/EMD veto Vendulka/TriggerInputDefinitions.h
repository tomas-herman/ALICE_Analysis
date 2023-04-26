

UChar_t inputId_0VBA; // 0VBA: >=1 V0A cell fired in BB timing gate
UChar_t inputId_0VBC; // 0VBC: >=1 V0C cell fired in BB timing gate
UChar_t inputId_0UBA; // 0UBA: >=1 ADA cell fired in BB timing gate
UChar_t inputId_0UBC; // 0UBC: >=1 ADC cell fired in BB timing gate
UChar_t inputId_0SH1; // 0SH1: FO fired in SPD
UChar_t inputId_0STG; 
UChar_t inputId_1ZED;
UChar_t inputId_0MUL;
UChar_t inputId_0OM2;
UChar_t inputId_0VOM;

// --------------------------------------------------

void Set2018PbPb()
{
  
  // --------------------------------
  // defines variables with the index of the
  // different L0/L1 trigger elements according to
  // Run 296786 (LHC18r), checked with 295666 (LHC18q)
  // --------------------------------
  
  inputId_0VBA = 1; // 0VBA: >=1 V0A cell fired in BB timing gate
  inputId_0VBC = 2; // 0VBC: >=1 V0C cell fired in BB timing gate
  inputId_0UBA = 4; // 0UBA: >=1 ADA cell fired in BB timing gate
  inputId_0UBC = 5; // 0UBC: >=1 ADC cell fired in BB timing gate
  inputId_0VOM = 7;
  inputId_0STG = 16;  
  inputId_0MUL = 19;  
  inputId_1ZED = 15;}


// --------------------------------------------------

void Set2015PbPb()
{
  
  // --------------------------------
  // defines variables with the index of the
  // different L0/L1 trigger elements according to
  // Run 245064
  // --------------------------------
  
  inputId_0VBA = 1; // 0VBA: >=1 V0A cell fired in BB timing gate
  inputId_0VBC = 2; // 0VBC: >=1 V0C cell fired in BB timing gate
  inputId_0UBA = 7; // 0UBA: >=1 ADA cell fired in BB timing gate
  inputId_0UBC = 8; // 0UBC: >=1 ADC cell fired in BB timing gate
  inputId_0SH1 = 13; // 0SH1: >=2 outer FO hits
  inputId_1ZED = 15;

}

void SetXeXe()
{
  
  // --------------------------------
  // defines variables with the index of the
  // different L0/L1 trigger elements according to
  // Run 280235
  // --------------------------------
  
  inputId_0VBA = 1; // 0VBA: >=1 V0A cell fired in BB timing gate
  inputId_0VBC = 2; // 0VBC: >=1 V0C cell fired in BB timing gate
  inputId_0UBA = 6; // 0UBA: >=1 ADA cell fired in BB timing gate
  inputId_0UBC = 7; // 0UBC: >=1 ADC cell fired in BB timing gate
  inputId_0SH1 = 8; // 0SH1: >=2 outer and >=2 inner FO hits 
  inputId_1ZED = 15;
  
  inputId_0OM2 = 9;
}
