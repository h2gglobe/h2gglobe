#!/bin/env python

from ROOT import TFile, TF1, TCanvas, kBlack, kBlue, kRed, kGreen, kCyan, kMagenta, TGraph, TF2, TF12, kFullCircle, gROOT
from sys import argv

def eff(fin,name,var="pt",rebin=4):

    num=fin.Get("%s_rv_cat0_tot" % var )
    den=fin.Get("%s_cat0_tot"    % var )

    if rebin > 0:
        num.Rebin(rebin)
        den.Rebin(rebin)
    
    print num.Integral() / den.Integral()

    eff = num.Clone(name)
    eff.Divide( num, den, 1., 1., "B" )
    return eff

def prob(fin,name,var="pt",profile=True,rebin=4):

    prb = fin.Get("vtxprob_%s_cat0_tot" % var)

    if rebin > 0:
        prb.RebinX(rebin)
    
    if profile:
        return prb.ProfileX(name)

    return prb.Clone(name)
    
gROOT.LoadMacro("../Macros/rootglobestyle.C")
from ROOT import setTDRStyle, gStyle
setTDRStyle()
gStyle.SetOptTitle(0)

## file="histograms_vtxeff.root"
file1=argv[1]
file2=argv[2]

fin1=TFile.Open(file1)
fin2=TFile.Open(file2)
eff1 = eff(fin1,"eff1")
eff2 = eff(fin2,"eff2")

eff2.SetLineColor(kRed)
eff2.SetMarkerColor(kRed)

c = TCanvas("c","c")
eff1.Draw()
eff2.Draw("hist same")

prob1 = prob(fin1,"prob1")
prob2 = prob(fin2,"prob2")

prob2.SetLineColor(kRed)
prob2.SetMarkerColor(kRed)

c1 = TCanvas("c1","c1")
eff1.Draw()
prob1.Draw("hist same")

c2 = TCanvas("c2","c2")
eff2.Draw()
prob2.Draw("hist same")
