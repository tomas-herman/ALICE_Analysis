

import ROOT
from array import array


binsize = 10

NRGBs = 5

stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
yellow   = [ 0.00, 0.00, 0.87, 1.00, 0.51 ]
green = [ 0.00, 0.81, 1.00, 0.20, 0.00 ]
blue  = [ 0.51, 1.00, 0.12, 0.00, 0.00 ]

stopsArray = array('d', stops)
yellowArray = array('d', yellow)
greenArray = array('d', green)
blueArray = array('d', blue)
FI = ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, yellowArray, greenArray, blueArray, binsize)
ROOT.gStyle.SetNumberContours(binsize);

# Canvas for Testing
cTest = ROOT.TCanvas('ColorTest','Test',1920,1080)
cTest.cd()

form1 = ROOT.TFormula( 'form1', 'abs(sin(x)/x)' )
sqroot = ROOT.TF1( 'sqroot', 'x*gaus(0) + [3]*form1', 0, 10 )

htest=[]
for i in range(0, binsize):
    sqroot.SetParameters( 10, 4, i, 20 )
    htest.append(ROOT.TH1D( f'test{i}', 'Test histogram with random numbers', 200, 0, 10 ))
    htest[i].FillRandom( 'sqroot', 10000 )
    htest[i].SetLineColor(FI+i) # This is the point! FI(ColorTable) + some int will generate color automatically!
    htest[i].SetLineWidth(4)
    htest[i].Draw("SAME")

# func = ROOT.TF2( 'func', 'x+y', 0, 1, 0, 1)
# hist = ROOT.TH2D("hist", "hist", 5, 0, 5, 5, 0, 5)
# for i in range (0, 5):
# 	for j in range (0,5):
# 		hist.FillRandom('func', 100)

# ROOT.gStyle.SetOptStat(0)
# hist.Draw("colz")

cTest.SaveAs("test.png")

