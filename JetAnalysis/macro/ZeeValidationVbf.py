import ROOT
import sys

from ROOT import TFile, TCanvas, TStyle, TApplication, TLegend, TGraph, TGraphErrors

# ROOT Setup
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

fin=TFile.Open("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/zee_validation/zeevalidation-27-07-2012/histograms_CMS-HGG_test.root")

savePlots = "true"

vars = ["_Mgg0",
        "LeadJPt",
        "SubJPt",        
        "_dEta",
        "_Zep",
        "_Mjj",
        "_dPhi",
        "_LeadPhoPtOverM",
        "_SubPhoPtOverM",
        "_DiPhoPtOverM",
        "vbf_mva"]

hdata     = {}
hmc       = {}
hdata_seq = {}
hmc_seq   = {}
canvas    = {}

for i,v in enumerate(vars):
    if i < (len(vars)-1):
        hdata[i] = fin.Get("cut_VBF"+v+"_nminus1_cat0_Data")
        hmc[i]   = fin.Get("cut_VBF"+v+"_nminus1_cat0_tot")
        hdata_seq[i] = fin.Get("cut_VBF"+v+"_sequential_cat0_Data")
        hmc_seq[i]   = fin.Get("cut_VBF"+v+"_sequential_cat0_tot")
    else:
        hdata[i] = fin.Get(v+"_cat0_Data")
        hmc[i]   = fin.Get(v+"_cat0_tot")
        hdata_seq[i] = fin.Get(v+"_cat0_Data")
        hmc_seq[i]   = fin.Get(v+"_cat0_tot")
        hdata[i].Rebin(5)
        hmc[i].Rebin(5)
        hmc[i].GetXaxis().SetRangeUser(-1,1);
                
    hdata[i].SetMarkerStyle(20)
    hdata[i].SetMarkerColor(1)
    hdata[i].SetLineColor(1)
    hdata_seq[i].SetMarkerStyle(20)
    hdata_seq[i].SetMarkerColor(1)
    hmc[i].SetFillStyle(3000)
    hmc[i].SetFillColor(ROOT.kPink+1)
    hmc_seq[i].SetFillStyle(3000)
    hmc_seq[i].SetFillColor(ROOT.kPink+1)
    

for n in range(len(hmc)):
    cname = "c%d" % n
    canvas[n] = TCanvas(cname,cname,800,500)
    canvas[n].Divide(2,1)
    canvas[n].cd(1)
    hmc[n].GetYaxis().SetTitleOffset(1.4)
    hmc[n].DrawNormalized("histo")
    hdata[n].DrawNormalized("esame")
    canvas[n].cd(2)
    if n < (len(vars)-1):
        hmc_seq[n].GetYaxis().SetTitleOffset(1.4)
        hmc_seq[n].DrawNormalized("histo")
        hdata_seq[n].DrawNormalized("esame")
    if savePlots == "true":
        canvas[n].SaveAs(vars[n].replace("_","")+".png")
        
        
s = raw_input("OK?")

if s == "yes":
    print "plots saved " 
else:
    print "ciao ciao"
