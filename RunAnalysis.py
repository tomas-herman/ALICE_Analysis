# // Running macro to analyse J/Psi UPC
# //-----------------------------------------------------------------------------------

# ////////////////////////////////////////
# // Including headers and functions 
# ////////////////////////////////////////

import os

# import time
# time.sleep(1500)

# ////////////////////////////////////////
# // Rapidity range options
# ////////////////////////////////////////

# ---0------1-----2------3------4------5------6------7------8------9----
# {-4.0, -4.00, -3.75, -3.50, -3.25, -3.00, -2.75, -4.00, -3.50, -3.00};
# {-2.5, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50, -3.50, -3.00, -2.50};

# ////////////////////////////////////////
# // Computing Lumi
# ////////////////////////////////////////

# ####### 18q #######
# os.system(('root -l -q LuminosityCalculationStandAlone.C++g\(\\"18q\\",\\"CMUP6\\"\)' ) )
# os.system(('root -l -q LuminosityCalculationStandAlone.C++g\(\\"18q\\",\\"CMUP11\\"\)' ) )

# # ####### 18r #######
# os.system(('root -l -q LuminosityCalculationStandAlone.C++g\(\\"18r\\",\\"CMUP6\\"\)' ) )
# os.system(('root -l -q LuminosityCalculationStandAlone.C++g\(\\"18r\\",\\"CMUP11\\"\)' ) )

# ####### 15o #######
# os.system(('root -l -q LuminosityCalculationStandAlone.C++g\(\\"15o\\",\\"CMUP10\\"\)' ) )
# os.system(('root -l -q LuminosityCalculationStandAlone.C++g\(\\"15o\\",\\"CMUP11\\"\)' ) )

# ////////////////////////////////////////
# // Computing Efficiency
# ////////////////////////////////////////
# for i in [0,1,2,3,4,5,6]:
# for i in [0,7,8,9]:
# for i in [0]:

  # ####### CMUP11 #######
  # os.system(('root -l -q Efficiency.C++g\(%d,\\"ADvetoOn\\",\\"CMUP11\\"\)' % (i)))

  # ####### CMUP6 #######
  # // Without AD Veto
  # os.system(('root -l -q Efficiency.C++g\(%d,\\"ADvetoOff\\",\\"CMUP6\\"\)' % (i)))
  # // With AD Veto
  # os.system(('root -l -q Efficiency.C++g\(%d,\\"ADvetoOn\\",\\"CMUP6\\"\)' % (i)))

# //////////////////////////////////////////////////////////////////////////////////////
# // Getting mass fit parameters from MC and creating PDFs for Pt fit
# //////////////////////////////////////////////////////////////////////////////////////
# for i in [0,1,2,3,4,5,6]:
# for i in [0,7,8,9]:
# for i in [0]:

  # ####### CMUP11 #######
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP11\\",\\"LoadRatio\\,\\"Compute_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Create_PDFs\\"\)' % (i) ) )

  # ####### CMUP6 #######
  # // Without AD Veto
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Compute_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Create_PDFs\\"\)' % (i) ) )

  # // With AD Veto
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Compute_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Create_PDFs\\"\)' % (i) ) )

# ////////////////////////////////////////
# // Computing Psi' to J/Psi CS ratio
# ////////////////////////////////////////
# for i in [0]:

  # ####### CMUP11 #######
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP11\\",\\"RatioOnly\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )

  # ####### CMUP6 #######
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"RatioOnly\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"RatioOnly\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )

# ////////////////////////////////////////
# // Getting disoc parameters
# ////////////////////////////////////////
# for i in [0,7,8,9]:
# for i in [0]:

  # ####### CMUP11 #######
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP11\\",\\"LoadRatio\\,\\"Load_MC_fit_param\\",\\"Compute_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )

  # ####### CMUP6 #######
  # // Without AD Veto
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Compute_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"0N0N\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Compute_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"0NXN\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Compute_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"XN0N\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Compute_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"XNXN\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Compute_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )

  # // With AD Veto
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Compute_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )

# ////////////////////////////////////////////////////////////////////////////////
# // Computing yields and feed down and incoherent contributions
# ////////////////////////////////////////////////////////////////////////////////
# for i in [0,1,2,3,4,5,6]:
# for i in [0,7,8,9]:
for i in [0]:

# ####### CMUP11 #######
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP11\\",\\"LoadRatio\\,\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )

# ####### CMUP6 #######
# // Without AD Veto
  os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"0N0N\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"0NXN\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"XN0N\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"XNXN\\",\\"ADvetoOff\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )

# // With AD Veto
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"all\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"0N0N\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"0NXN\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"XN0N\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )
  # os.system(('root -l -q AnalysisJPsi.C++g\(%d,\\"XNXN\\",\\"ADvetoOn\\",\\"CMUP6\\",\\"LoadRatio\\",\\"Load_MC_fit_param\\",\\"Load_disoc_fit_param\\",\\"Load_PDFs\\"\)' % (i) ) )