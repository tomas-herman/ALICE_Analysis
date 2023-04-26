# // Computing cross section
# //-----------------------------------------------------------------------------------

# ////////////////////////////////////////
# // Including headers and functions 
# ////////////////////////////////////////

import os
import math
import ROOT
import numpy as np
from sympy import *
from uncertainties import ufloat

# ////////////////////////////////////////
# // Define global variabels
# ////////////////////////////////////////

# -------------0-----1-----2-----3-----4-----5-----6-----7-----8-----9---
y_bins_min = [4.00, 4.00, 3.75, 3.50, 3.25, 3.00, 2.75, 4.00, 3.50, 3.00]
y_bins_max = [2.50, 3.75, 3.50, 3.25, 3.00, 2.75, 2.50, 3.50, 3.00, 2.50]

rapidity_range = [7,8,9]
# rapidity_range = [1,2,3,4,5,6]

# From 2018 J/Psi Analysis Luminosity Note: https://alice-notes.web.cern.ch/node/869
V0A = ufloat(0.0385, 0.0015)
ADA = ufloat(0.0009, 0.0002)
ADC = ufloat(0.0007, 0.0002)

# Updated 2018 values computed with CTRUE
pA = 0.0237
pC = 0.0237
# From AN by Igor: https://alice-notes.web.cern.ch/system/files/notes/analysis/1062/2020-04-17-ALICE_AN_ZDC_efficiency.pdf (https://alice-notes.web.cern.ch/node/1062)
eA = 0.94
eC = 0.94

# EMD correction 
F_XN0N = 0.88
F_XNXN = 0.84
# F_all = 1.014

lumi_18q_CMUP6  = 215.996
lumi_18q_CMUP11 = 216.502

lumi_18r_CMUP6  = 316.284
lumi_18r_CMUP11 = 316.286

lumi_15o_CMUP10 = 23.367
lumi_15o_CMUP11 = 192.437
# lumi_15o_CMUP10 = 0
# lumi_15o_CMUP11 = 0

lumi_CMUP6 = ufloat(lumi_18q_CMUP6+lumi_18r_CMUP6,(lumi_18q_CMUP6+lumi_18r_CMUP6)*0.0)
lumi_CMUP11 = ufloat(lumi_18q_CMUP11+lumi_18r_CMUP11+lumi_15o_CMUP11+lumi_15o_CMUP10,(lumi_18q_CMUP11+lumi_18r_CMUP11+lumi_15o_CMUP11+lumi_15o_CMUP10)*0.0)
veto = 1-V0A
veto_ADveto = (1-V0A)*(1-ADA)*(1-ADC)
veto_CMUP11 = 0.95
branching = 0.05961

# ////////////////////////////////////////
# // Function to read input values
# ////////////////////////////////////////

def Read(name, ADveto, y_min, y_max, zdc_class, trigger):
	# Open file
	file = open("Cross_section/"+trigger+"/"+ADveto+"_"+y_min+"_"+y_max+"_"+zdc_class+".txt", "r")

	# Read each line and get the desired values
	lines = file.readlines()
	for line in lines:
		if (name+" = " in line):
			return(float(line[len(name+" = "):]))

# ////////////////////////////////////////
# // Function to compute CS
# ////////////////////////////////////////

def CrossSection(ADveto, y_min, y_max, zdc_class, trigger):
	# Define variabels
	N_all = []
	N_zdc = []

	CS_all = []
	CS_not_corrected = []
	CS_corrected = []

	# Correct yields for ZDC class migration due to pilup and efficiency
	for zdc in ["all","0N0N","0NXN","XN0N","XNXN"]:
		# Read variabels
		N_JPsi = ufloat(Read("N_JPsi", ADveto, y_min, y_max, zdc, trigger), Read("N_JPsi_err", ADveto, y_min, y_max, zdc, trigger))
		# print (y_min+"-"+y_max, zdc+": {:.6f}".format(N_JPsi / ( (1+f_I+f_D) * eff * branching * veto_ADveto * lumi_CMUP6 * (float(y_min)-float(y_max)) ) /1000))
		# Compute CS
		if (zdc == "all"):
			N_all.append(N_JPsi)
		else:
			N_zdc.append(N_JPsi)

	M = np.array([ [(1-pA)*(1-pC), 	(1-eC)*(1-pA)*(1-pC), 							(1-eA)*(1-pA)*(1-pC), 							(1-eA)*(1-eC)*(1-pA)*(1-pC)], 
	               [pC*(1-pA), 		1-((1-eC)*(1-pA)*(1-pC)+eC*pA+(1-eC)*pA*pC), 	(1-eA)*pC*eC, 									(1-eA)*(1-pA)*(eC+(1-eC)*pC)],
	               [pA*(1-pC), 		(1-eC)*pA*eA, 									1-((1-eA)*(1-pA)*(1-pC)+eA*pC+(1-eA)*pA*pC), 	(1-eC)*(1-pC)*(eA+(1-eA)*pA)], 
	               [pC*pA, 			eC*pA+(1-eC)*pA*pC, 							eA*pC+(1-eA)*pA*pC, 							1-((1-eA)*(1-pA)*(eC+(1-eC)*pC)+(1-eC)*(1-pC)*(eA+(1-eA)*pA)+(1-eA)*(1-eA)*(1-eC)*(1-pA)*(1-pC))] ]) 

	Mi = np.linalg.inv(M)

	# print (M)
	# print (Mi)

	# N_corrected = (np.dot(Mi,N_zdc))

	# Compute the cross section for all zdc classes
	for zdc in ["all","0N0N","0NXN","XN0N","XNXN"]:
		# Read variabels
		f_D = ufloat(Read("f_D", ADveto, y_min, y_max, zdc, trigger), Read("f_D_err", ADveto, y_min, y_max, zdc, trigger))
		f_I = ufloat(Read("f_I", ADveto, y_min, y_max, zdc, trigger), Read("f_I_err", ADveto, y_min, y_max, zdc, trigger))
		eff = Read("eff", ADveto, y_min, y_max, zdc, trigger)
		# print (y_min+"-"+y_max, zdc+": {:.6f}".format(f_D))
		# Compute CS
		if (zdc == "all"):
			N_JPsi = N_all[0]
			# if ("ADvetoOn" == ADveto):
			# 	if ("CMUP6" == trigger):
			# 		CS_all.append(N_JPsi / ( (1+f_I+f_D) * eff * branching * veto_ADveto * lumi_CMUP6 * (float(y_min)-float(y_max)) ) /1000)
			# 	if ("CMUP11" == trigger):
			# 		CS_all.append(N_JPsi / ( (1+f_I+f_D) * eff * branching * veto_CMUP11 * lumi_CMUP11 * (float(y_min)-float(y_max)) ) /1000)
			if ("ADvetoOff" == ADveto):
				if ("CMUP6" == trigger):
					CS_all.append(N_JPsi / ( (1+f_I+f_D) * eff * veto * branching * lumi_CMUP6 * (float(y_min)-float(y_max)) ) /1000)
		else:
			if (zdc == "0N0N"): N_JPsi = N_zdc[0]
			if (zdc == "0NXN"): N_JPsi = N_zdc[1]
			if (zdc == "XN0N"): N_JPsi = N_zdc[2]/F_XN0N
			if (zdc == "XNXN"): N_JPsi = N_zdc[3]/F_XNXN

			# if ("ADvetoOn" == ADveto):
				# if ("CMUP6" == trigger):
					# need to use proper pilep and EMD vetos # CS_corrected.append(N_JPsi / ( (1+f_I+f_D) * eff * branching * veto_ADveto * lumi_CMUP6 * (float(y_min)-float(y_max)) ) /1000)
				# if ("CMUP11" == trigger):
					# need to use proper pilep and EMD vetos # CS_corrected.append(N_JPsi / ( (1+f_I+f_D) * eff * branching * veto_CMUP11 * lumi_CMUP11 * (float(y_min)-float(y_max)) ) /1000)
			if ("ADvetoOff" == ADveto):
				if ("CMUP6" == trigger):
					CS_not_corrected.append(N_JPsi / ( (1+f_I+f_D) * eff * veto * branching * lumi_CMUP6 * (float(y_min)-float(y_max)) ) /1000)

	CS_corrected = (np.dot(Mi,CS_not_corrected))

	if ("all" == zdc_class):
		# print (y_min+"-"+y_max, zdc_class+": {:.6f}".format(CS_all[0]))
		return(CS_all[0]*F_all)
	if ("0N0N" == zdc_class):
		# print (y_min+"-"+y_max, zdc_class+": {:.6f}".format(CS_corrected[0]))
		return(CS_corrected[0])
	if ("0NXN" == zdc_class):
		# print (y_min+"-"+y_max, zdc_class+": {:.6f}".format(CS_corrected[1]))
		return(CS_corrected[1])
	if ("XN0N" == zdc_class):
		# print (y_min+"-"+y_max, zdc_class+": {:.6f}".format(CS_corrected[2]))
		return(CS_corrected[2])
	if ("XNXN" == zdc_class):
		# print (y_min+"-"+y_max, zdc_class+": {:.6f}".format(CS_corrected[3]))
		return(CS_corrected[3])			

# ////////////////////////////////////////
# // Functions to plot
# ////////////////////////////////////////
# Define plotting function 
def PlotCS_axis(zoom):
	# Define variabels
	CS_val = []
	CS_err = []
	x_val = []
	x_err = []

	# Compute the cross section
	for i in rapidity_range: 
		CS_val.append(0)
		CS_err.append(0)
		x_val.append(0) 
		x_err.append(0)

	# Create cross section graph
	CS_graph = ROOT.TGraphErrors(3, np.array(x_val), np.array(CS_val), np.array(x_err), np.array(CS_err))

	if (zoom==0):
		# Zoom on ALICE forward region
		CS_graph.GetXaxis().SetLimits(-4.2,-2.2)
		CS_graph.GetYaxis().SetRangeUser(0,5)
	if (zoom==1):
		# Zoom on 0NXN, XN0N and XNXN
		CS_graph.GetXaxis().SetLimits(-4.2,-2.2)
		CS_graph.GetYaxis().SetRangeUser(0,0.45)
	if (zoom==2):
		# Zoom on 0N0N
		CS_graph.GetXaxis().SetLimits(-4.2,-2.2)
		CS_graph.GetYaxis().SetRangeUser(1.4,3.8)	
	if (zoom==3):
		# Zoom on ALICE forward+central region + LHCb forward region
		CS_graph.GetXaxis().SetLimits(-4.7,0.9)
		CS_graph.GetYaxis().SetRangeUser(1,6)	

	CS_graph.GetXaxis().SetTitle("#it{y}")
	CS_graph.GetYaxis().SetTitle("d#it{#sigma}/d#it{y} (mb)")
	CS_graph.GetYaxis().SetTitleOffset(1.1);
	CS_graph.SetLineColor(1)

	return CS_graph

def PlotCS(zdc_class, trigger, ADveto, opt):
	# Define variabels
	CS_val = []
	CS_err = []
	x_val = []
	x_err = []

	# Compute the cross section
	for i in rapidity_range: 
		CS = CrossSection(ADveto, "%.2f" %y_bins_min[i], "%.2f" %y_bins_max[i], zdc_class, trigger)
		CS_val.append(CS.nominal_value)
		# print (CS.nominal_value)
		CS_err.append(CS.std_dev)
		x_val.append( -(y_bins_max[i]+ (1.5/2/len(rapidity_range)) ) ) 
		x_err.append(1.5/2/len(rapidity_range) )

	# Create cross section graph
	CS_graph = ROOT.TGraphErrors(len(rapidity_range), np.array(x_val), np.array(CS_val), np.array(x_err), np.array(CS_err))

	if (zdc_class == "all"):
		if ("ADvetoOn" == ADveto):
			if ("CMUP6" == trigger):
				CS_graph.SetLineColor(ROOT.kRed+1)
				CS_graph.SetMarkerColor(ROOT.kRed+1)
			if ("CMUP11" == trigger):	
				CS_graph.SetLineColor(ROOT.kGreen+2)
				CS_graph.SetMarkerColor(ROOT.kGreen+2)
		if ("ADvetoOff" == ADveto):
			CS_graph.SetLineColor(ROOT.kGreen+2)
			CS_graph.SetMarkerColor(ROOT.kGreen+2)
	elif (zdc_class =="0N0N"):
		if (opt == 4 or opt == 5 or opt == 6 or opt == 7):
			if ("ADvetoOn" == ADveto):
				CS_graph.SetLineColor(ROOT.kRed+1)
				CS_graph.SetMarkerColor(ROOT.kRed+1)
			if ("ADvetoOff" == ADveto):
				CS_graph.SetLineColor(ROOT.kGreen+2)
				CS_graph.SetMarkerColor(ROOT.kGreen+2)
		else:
			CS_graph.SetLineColor(ROOT.kRed+1)
			CS_graph.SetMarkerColor(ROOT.kRed+1)
	elif (zdc_class =="0NXN"):
		if (opt == 4 or opt == 5 or opt == 6 or opt == 7):
			if ("ADvetoOn" == ADveto):
				CS_graph.SetLineColor(ROOT.kRed+1)
				CS_graph.SetMarkerColor(ROOT.kRed+1)
			if ("ADvetoOff" == ADveto):
				CS_graph.SetLineColor(ROOT.kGreen+2)
				CS_graph.SetMarkerColor(ROOT.kGreen+2)
		else:
			CS_graph.SetLineColor(ROOT.kGreen+2)
			CS_graph.SetMarkerColor(ROOT.kGreen+2)
	elif (zdc_class =="XN0N"):
		if (opt == 4 or opt == 5 or opt == 6 or opt == 7):
			if ("ADvetoOn" == ADveto):
				CS_graph.SetLineColor(ROOT.kRed+1)
				CS_graph.SetMarkerColor(ROOT.kRed+1)
			if ("ADvetoOff" == ADveto):
				CS_graph.SetLineColor(ROOT.kGreen+2)
				CS_graph.SetMarkerColor(ROOT.kGreen+2)
		else:
			CS_graph.SetLineColor(ROOT.kBlue+1)
			CS_graph.SetMarkerColor(ROOT.kBlue+1)
	elif (zdc_class =="XNXN"):
		if (opt == 4 or opt == 5 or opt == 6 or opt == 7):
			if ("ADvetoOn" == ADveto):
				CS_graph.SetLineColor(ROOT.kRed+1)
				CS_graph.SetMarkerColor(ROOT.kRed+1)
			if ("ADvetoOff" == ADveto):
				CS_graph.SetLineColor(ROOT.kGreen+2)
				CS_graph.SetMarkerColor(ROOT.kGreen+2)
		else:
			CS_graph.SetLineColor(ROOT.kMagenta+1)
			CS_graph.SetMarkerColor(ROOT.kMagenta+1)
	return CS_graph

# Define function to plot published values
def PlotPublished():
	# Define values for CS and rapidity
	n 	=  6
	x 	= (-3.875, -3.625, -3.375, -3.125, -2.875, -2.625)
	y 	= (1.621, 1.936, 2.376, 2.830, 3.014, 3.585)
	# Define rapidity range
	xh_stat = (0.000, 0.000, 0.000, 0.000, 0.000, 0.000)
	xl_stat = (0.000, 0.000, 0.000, 0.000, 0.000, 0.000)
	xh_syst = (0.125, 0.125, 0.125, 0.125, 0.125, 0.125)
	xl_syst = (0.125, 0.125, 0.125, 0.125, 0.125, 0.125)
	# Define CS errors
	yh_stat = (0.061, 0.042, 0.040, 0.047, 0.061, 0.141)
	yl_stat = (0.061, 0.042, 0.040, 0.047, 0.061, 0.141)
	yh_syst = (0.135, 0.166, 0.212, 0.253, 0.259, 0.298)
	yl_syst = (0.148, 0.190, 0.229, 0.280, 0.294, 0.368)

	# Graph with stat error
	CS_graph_published_stat = ROOT.TGraphAsymmErrors(n, np.array(x), np.array(y), np.array(xl_stat), np.array(xh_stat), np.array(yl_stat), np.array(yh_stat))
	CS_graph_published_stat.SetMarkerColor(ROOT.kBlue+1)
	CS_graph_published_stat.SetLineColor(ROOT.kBlue+1)
	# Graph with syst errors
	CS_graph_published_syst = ROOT.TGraphAsymmErrors(n, np.array(x), np.array(y), np.array(xl_syst), np.array(xh_syst), np.array(yl_syst), np.array(yh_syst))
	CS_graph_published_syst.SetLineColor(ROOT.kBlue+1)
	CS_graph_published_syst.SetFillStyle(1)

	return CS_graph_published_stat, CS_graph_published_syst

# Define function to plot published values from LHCB
def PlotPublishedLHCb():
	# Define values for CS and rapidity
	n 	=  5
	x 	= (-2.25, -2.75, -3.25, -3.75, -4.25)
	y 	= (3.0, 2.60, 2.28, 1.73, 1.10)
	# Define rapidity range
	xh_stat = (0.00, 0.00, 0.00, 0.00, 0.00)
	xl_stat = (0.00, 0.00, 0.00, 0.00, 0.00)
	xh_syst = (0.25, 0.25, 0.25, 0.25, 0.25)
	xl_syst = (0.25, 0.25, 0.25, 0.25, 0.25)
	# Define CS errors
	yh_stat = (0.4, 0.19, 0.15, 0.15, 0.22)
	yl_stat = (0.4, 0.19, 0.15, 0.15, 0.22)
	yh_syst = (0.3, 0.25, 0.21, 0.17, 0.13)
	yl_syst = (0.3, 0.25, 0.21, 0.17, 0.13)

	# Graph with stat error
	CS_graph_published_stat = ROOT.TGraphAsymmErrors(n, np.array(x), np.array(y), np.array(xl_stat), np.array(xh_stat), np.array(yl_stat), np.array(yh_stat))
	CS_graph_published_stat.SetMarkerColor(ROOT.kGreen+2)
	CS_graph_published_stat.SetLineColor(ROOT.kGreen+2)
	# Graph with syst errors
	CS_graph_published_syst = ROOT.TGraphAsymmErrors(n, np.array(x), np.array(y), np.array(xl_syst), np.array(xh_syst), np.array(yl_syst), np.array(yh_syst))
	CS_graph_published_syst.SetLineColor(ROOT.kGreen+2)
	CS_graph_published_syst.SetFillStyle(1)

	return CS_graph_published_stat, CS_graph_published_syst	

# Define function to plot published values from LHCB
def PlotPublishedCentral():
	# Define values for CS and rapidity
	n 	=  5
	x 	= (-0.475, -0.25, 0.00, 0.25, 0.475)
	y 	= (3.85008, 4.17624, 4.06709, 4.17624, 3.85008)
	# Define rapidity range
	xh_stat = (0.00, 0.00, 0.00, 0.00, 0.00)
	xl_stat = (0.00, 0.00, 0.00, 0.00, 0.00)
	xh_syst = (0.225, 0.10, 0.15, 0.10, 0.225)
	xl_syst = (0.225, 0.10, 0.15, 0.10, 0.225)
	# Define CS errors
	yh_stat = (0.109666, 0.113591, 0.114712, 0.113591, 0.109666)
	yl_stat = (0.109666, 0.113591, 0.114712, 0.113591, 0.109666)
	yh_syst = (0.273222, 0.296368, 0.288622, 0.296368, 0.273222)
	yl_syst = (0.273222, 0.296368, 0.288622, 0.296368, 0.273222)

	# Graph with stat error
	CS_graph_published_stat = ROOT.TGraphAsymmErrors(n, np.array(x), np.array(y), np.array(xl_stat), np.array(xh_stat), np.array(yl_stat), np.array(yh_stat))
	CS_graph_published_stat.SetMarkerColor(ROOT.kBlue+1)
	CS_graph_published_stat.SetLineColor(ROOT.kBlue+1)
	# Graph with syst errors
	CS_graph_published_syst = ROOT.TGraphAsymmErrors(n, np.array(x), np.array(y), np.array(xl_syst), np.array(xh_syst), np.array(yl_syst), np.array(yh_syst))
	CS_graph_published_syst.SetLineColor(ROOT.kBlue+1)
	CS_graph_published_syst.SetFillStyle(1)

	return CS_graph_published_stat, CS_graph_published_syst	

# ////////////////////////////////////////
# // Compute and plot cross section
# ////////////////////////////////////////
def Run(trigger, ADveto, opt):
	# ---------------------------
	# Set general options
	ROOT.gStyle.SetOptTitle(0)
	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetPaintTextFormat("4.3f")
	ROOT.gStyle.SetFrameLineWidth(1)
	ROOT.gStyle.SetLabelSize(0.045,"xyz")
	ROOT.gStyle.SetTitleSize(0.05,"xyz")
	ROOT.gStyle.SetTextSize(0.04)
	ROOT.gStyle.SetMarkerStyle(20)
	ROOT.gStyle.SetMarkerSize(0.8)

	# ---Legend - definition
	leg_CS = ROOT.TLegend(0.13,0.6,0.5,0.82)
	leg_CS.SetFillStyle(0)
	leg_CS.SetBorderSize(0)
	leg_CS.SetTextSize(0.04)

	# ---------------------------
	# Define Canvas
	cCS = ROOT.TCanvas("cCS","cCS",1000,800)
	cCS.SetLeftMargin(0.11);
	cCS.SetRightMargin(0.05);
	cCS.Draw() 
	cCS.cd()

	# ---Fill cross section canvas
	# ------Plot axis
	if (opt == 3 or opt == 5 or opt == 6 or opt == 7):
		CS_graph0 = PlotCS_axis(1)
		CS_graph0.Draw("same ap")
	if (opt == 4):
		CS_graph0 = PlotCS_axis(2)
		CS_graph0.Draw("same ap")	
	if (opt == 1 or opt == 2):
		CS_graph0 = PlotCS_axis(0)
		CS_graph0.Draw("same ap")
	if (opt == 0):
		CS_graph0 = PlotCS_axis(3)
		CS_graph0.Draw("same ap")		

# opt0 call once with ad veto
	if (opt==0):

		CS_graph_published_stat, CS_graph_published_syst = PlotPublished()
		CS_graph_published_stat.Draw("same p")
		CS_graph_published_syst.Draw("same e2")
		leg_CS.AddEntry(CS_graph_published_stat,"ALICE Coherent J/#psi Published","pf")

		Central_CS_graph_published_stat, Central_CS_graph_published_syst = PlotPublishedCentral()
		Central_CS_graph_published_stat.Draw("same p")
		Central_CS_graph_published_syst.Draw("same e2")
		# leg_CS.AddEntry(Central_CS_graph_published_stat,"ALICE Coherent J/#psi Published","pf")

		LHCb_CS_graph_published_stat, LHCb_CS_graph_published_syst = PlotPublishedLHCb()
		LHCb_CS_graph_published_stat.Draw("same p")
		LHCb_CS_graph_published_syst.Draw("same e2")
		leg_CS.AddEntry(LHCb_CS_graph_published_stat,"LHCb Coherent J/#psi Preliminary","pf")

		CS_graph1a = PlotCS("all", trigger, ADveto, opt)
		CS_graph1a.Draw("same p")
		leg_CS.AddEntry(CS_graph1a,"ALICE J/#psi, "+trigger+", "+ADveto+", ZDC: all", "lep")

# opt1 call once with ad veto
	if (opt==1):

		# CS_graph1b = PlotCS("all", trigger, ADvetoOff, opt)
		# CS_graph1b.Draw("same p")
		# leg_CS.AddEntry(CS_graph1b,"ALICE J/#psi, ADvetoOff, ZDC: all", "lep")

		CS_graph1b = PlotCS("all", "CMUP11", "ADvetoOn", opt)
		CS_graph1b.Draw("same p")
		leg_CS.AddEntry(CS_graph1b,"ALICE J/#psi, CMUP11, ADvetoOn, ZDC: all", "lep")

		CS_graph_published_stat, CS_graph_published_syst = PlotPublished()
		CS_graph_published_stat.Draw("same p")
		CS_graph_published_syst.Draw("same e2")
		leg_CS.AddEntry(CS_graph_published_stat,"ALICE Coherent J/#psi Published","pf")

		CS_graph1a = PlotCS("all", "CMUP6", "ADvetoOn", opt)
		CS_graph1a.Draw("same p")
		leg_CS.AddEntry(CS_graph1a,"ALICE J/#psi, CMUP6, ADvetoOn, ZDC: all", "lep")

# opt2 call once with once wihtou AD veto
	if (opt==2):
		CS_graph2 = PlotCS("0N0N", trigger, ADveto, opt)
		CS_graph2.Draw("same p")
		leg_CS.AddEntry(CS_graph2,"ALICE J/#psi, "+ADveto+", ZDC: 0N0N", "lep")

		CS_graph3 = PlotCS("0NXN", trigger, ADveto, opt)
		CS_graph3.Draw("same p")
		leg_CS.AddEntry(CS_graph3,"ALICE J/#psi, "+ADveto+", ZDC: 0NXN", "lep")

		CS_graph4 = PlotCS("XN0N", trigger, ADveto, opt)
		CS_graph4.Draw("same p")
		leg_CS.AddEntry(CS_graph4,"ALICE J/#psi, "+ADveto+", ZDC: XN0N", "lep")

		CS_graph5 = PlotCS("XNXN", trigger, ADveto, opt)
		CS_graph5.Draw("same p")
		leg_CS.AddEntry(CS_graph5,"ALICE J/#psi, "+ADveto+", ZDC: XNXN", "lep")


# opt3 call once with once wihtou AD veto
	if (opt==3):
		CS_graph3 = PlotCS("0NXN", trigger, ADveto, opt)
		CS_graph3.Draw("same p")
		leg_CS.AddEntry(CS_graph3,"ALICE J/#psi, "+ADveto+", ZDC: 0NXN", "lep")

		CS_graph4 = PlotCS("XN0N", trigger, ADveto, opt)
		CS_graph4.Draw("same p")
		leg_CS.AddEntry(CS_graph4,"ALICE J/#psi, "+ADveto+", ZDC: XN0N", "lep")

		CS_graph5 = PlotCS("XNXN", trigger, ADveto, opt)
		CS_graph5.Draw("same p")
		leg_CS.AddEntry(CS_graph5,"ALICE J/#psi, "+ADveto+", ZDC: XNXN", "lep")

# opt4 call once with AD veto
	if (opt == 4):
		CS_graph2a = PlotCS("0N0N", trigger, "ADvetoOn", opt)
		CS_graph2a.Draw("same p")
		leg_CS.AddEntry(CS_graph2a,"ALICE J/#psi, ADvetoOn, ZDC: 0N0N", "lep")

		CS_graph2b = PlotCS("0N0N", trigger, "ADvetoOff", opt)
		CS_graph2b.Draw("same p")
		leg_CS.AddEntry(CS_graph2b,"ALICE J/#psi, ADvetoOff, ZDC: 0N0N", "lep")

# opt5 call once with AD veto
	if (opt == 5):
		CS_graph3a = PlotCS("0NXN", trigger, "ADvetoOn", opt)
		CS_graph3a.Draw("same p")
		leg_CS.AddEntry(CS_graph3a,"ALICE J/#psi, ADvetoOn, ZDC: 0NXN", "lep")

		CS_graph3b = PlotCS("0NXN", trigger, "ADvetoOff", opt)
		CS_graph3b.Draw("same p")
		leg_CS.AddEntry(CS_graph3b,"ALICE J/#psi, ADvetoOff, ZDC: 0NXN", "lep")

# opt6 call once with AD veto
	if (opt == 6):
		CS_graph4a = PlotCS("XN0N", trigger, "ADvetoOn", opt)
		CS_graph4a.Draw("same p")
		leg_CS.AddEntry(CS_graph4a,"ALICE J/#psi, ADvetoOn, ZDC: XN0N", "lep")

		CS_graph4b = PlotCS("XN0N", trigger, "ADvetoOff", opt)
		CS_graph4b.Draw("same p")
		leg_CS.AddEntry(CS_graph4b,"ALICE J/#psi, ADvetoOff, ZDC: XN0N", "lep")

# opt7 call once with AD veto
	if (opt == 7):
		CS_graph5a = PlotCS("XNXN", trigger, "ADvetoOn", opt)
		CS_graph5a.Draw("same p")
		leg_CS.AddEntry(CS_graph5a,"ALICE J/#psi, ADvetoOn, ZDC: XNXN", "lep")

		CS_graph5b = PlotCS("XNXN", trigger, "ADvetoOff", opt)
		CS_graph5b.Draw("same p")
		leg_CS.AddEntry(CS_graph5b,"ALICE J/#psi, ADvetoOff, ZDC: XNXN", "lep")

	# ---Description
	text1 = ROOT.TLatex (0.13,0.85,"Preliminary #bf{ALICE, Pb+Pb #rightarrow Pb+Pb+J/#psi, #sqrt{#it{s}_{NN}} = 5.02 TeV}")
	text1.SetNDC()
	text1.Draw()

	# ---Legend - draw
	leg_CS.Draw()

	# ---Save as pdf
	if (opt==0):
		option_name = "all_my_vs_published"
	if (opt==1):
		option_name = "all_CMUP6_vs_CMUP11"	
	if (opt==2):
		if (ADveto):
			option_name = "ZDC_AD"
		else:	
			option_name = "ZDC_!AD"
	if (opt==3):
		if (ADveto):
			option_name = "ZDC_zoom_AD"
		else:	
			option_name = "ZDC_zoom_!AD"
	if (opt==4):
		option_name = "0N0N_AD_vs_!AD"
	if (opt==5):
		option_name = "0NXN_AD_vs_!AD"
	if (opt==6):
		option_name = "XN0N_AD_vs_!AD"	
	if (opt==7):
		option_name = "XNXN_AD_vs_!AD"

	CS_name = "Cross_section/Plots/%s_CS_%.2f_%.2f.png" %(str(option_name), y_bins_min[0], y_bins_max[0])
	cCS.SaveAs(CS_name)

# ////////////////////////////////////////
# // Run
# ////////////////////////////////////////

# Run("CMUP6", "ADvetoOn", 0)
# Run("CMUP6", "ADvetoOn", 1)
# Run("CMUP6", "ADvetoOn", 2)

Run("CMUP6", "ADvetoOff", 2)
Run("CMUP6", "ADvetoOff", 3)

# Run("CMUP6", "ADvetoOn", 3)
# Run("CMUP6", "ADvetoOn", 4)
# Run("CMUP6", "ADvetoOn", 5)
# Run("CMUP6", "ADvetoOn", 6)
# Run("CMUP6", "ADvetoOn", 7)