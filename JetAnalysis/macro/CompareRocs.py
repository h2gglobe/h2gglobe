import ROOT
import sys

from ROOT import TFile, TCanvas, TStyle, TApplication, TLegend, TGraph, TGraphErrors

# ROOT Setup
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(0)
ROOT.gStyle.SetHistLineWidth(4)


fin1=TFile.Open("OutputTMVAtrainingM124/TMVA_vbf_6var.root")
fin2=TFile.Open("OutputTMVAtrainingM124/TMVA_vbf_6var_diphopt.root")
fin3=TFile.Open("OutputTMVAtrainingM124/TMVA_vbf_6var_phopt.root")
fin4=TFile.Open("OutputTMVAtrainingM124/TMVA_vbf_6var_diphopt_phopt.root")
 
roc1=fin1.Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS")
roc2=fin2.Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS")
roc3=fin3.Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS")
roc4=fin4.Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS")

nre = 1

roc1.Rebin(nre)
roc2.Rebin(nre)
roc3.Rebin(nre)
roc4.Rebin(nre)

roc1.Scale(1./nre)
roc2.Scale(1./nre)
roc3.Scale(1./nre)
roc4.Scale(1./nre)

roc1.SetLineColor(1)
roc2.SetLineColor(2)
roc3.SetLineColor(3)
roc4.SetLineColor(4)


roc1.SetLineWidth(2)
roc2.SetLineWidth(2)
roc3.SetLineWidth(2)
roc4.SetLineWidth(3)



#nominal eff for standard analysis 
g = TGraphErrors()
g.SetPoint(0,0.553,0.969)
#g.SetPoint(0,0.606,0.947)
g.SetMarkerStyle(29)
g.SetLineColor(ROOT.kOrange+7)
g.SetMarkerColor(ROOT.kOrange+7)
g.SetMarkerSize(2)

g1 = TGraphErrors()
g1.SetPoint(0,0.343,0.994)
#g1.SetPoint(0,0.376,0.989)
g1.SetMarkerStyle(34)
g1.SetLineColor(6)
g1.SetMarkerColor(6)
g1.SetMarkerSize(1.5)


leg=TLegend(0.12,0.12,0.55,0.4)
leg.SetFillColor(0)
leg.SetBorderSize(1)
leg.AddEntry(roc1,"6 var","L")
leg.AddEntry(roc2,"6 var + p^{T}_{#gamma#gamma}/m_{#gamma#gamma}","L")
leg.AddEntry(roc3,"6 var + p^{T}_{#gamma 1}/m_{#gamma#gamma} + p^{T}_{#gamma 2}/m_{#gamma#gamma}","L")
leg.AddEntry(roc4,"6 var + p^{T}_{#gamma#gamma}/m_{#gamma#gamma} + p^{T}_{#gamma 1}/m_{#gamma#gamma} + p^{T}_{#gamma 2}/m_{#gamma#gamma}","L")
leg.AddEntry(g,"cut based VBF categories","P")
leg.AddEntry(g1,"cut based VBF cat1","P")


c = TCanvas("c","c")
c.SetGridx()
c.SetGridy()
ROOT.gStyle.SetOptStat(0)
roc1.GetYaxis().SetRangeUser(0.95,1.01)
roc1.GetXaxis().SetRangeUser(0.0 ,0.80)
roc1.GetXaxis().SetTitle("signal eff")
roc1.GetYaxis().SetTitle("bkg rejection")
roc1.Draw("l")
roc2.Draw("lsame")
roc3.Draw("lsame")
roc4.Draw("l")

g.Draw("psame")
g1.Draw("psame")
leg.Draw("same");
c.Update()

s = raw_input("Print canvas? ")

if s == "yes":
    print "Saving plot ..."
    c.SaveAs("compareRocs.png")
else:
    print "ciao ciao"
