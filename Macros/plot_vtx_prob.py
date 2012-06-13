#!/bin/env python

from ROOT import TFile, TF1, TCanvas, kBlack, kBlue, kRed, kGreen, kCyan, kMagenta, kWhite, TGraph, TF2, TF12, kFullCircle, gROOT, TLegend
from sys import argv

def eff(fin,name,var="pt",rebin=5):

    num=fin.Get("%s_rv_cat0_tot" % var )
    den=fin.Get("%s_cat0_tot"    % var )

    if rebin > 0:
        num.Rebin(rebin)
        den.Rebin(rebin)
    
    print num.Integral() / den.Integral()

    eff = num.Clone(name)
    eff.Divide( num, den, 1., 1., "B" )
    return eff

def prob(fin,name,var="pt",profile=True,rebin=5):

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


file=argv[1]

fin=TFile.Open(file)
eff1 = eff(fin,"eff1")
eff2 = eff(fin,"eff2","nvtx",0)

prob1 = prob(fin,"prob1")
prob2 = prob(fin,"prob2","nvtx",True,0)

prob1.SetLineColor(kRed)
prob1.SetMarkerColor(kRed)

prob2.SetLineColor(kRed)
prob2.SetMarkerColor(kRed)

c1 = TCanvas("eff_vs_pt","eff_vs_pt")
eff1.SetTitle("True vertex eff.;p_{T}(#gamma #gamma) (GeV);#varepsilon #Delta z < 1 cm")
prob1.SetTitle("Aveage vertex prob.;p_{T}(#gamma #gamma) (GeV);#varepsilon #Delta z < 1 cm")
eff1.GetYaxis().SetRangeUser(0.5, 1.1)
eff1.Draw("e0")
prob1.Draw("hist same")

leg1 = TLegend(0.5,0.2,0.8,0.4)
leg1.SetShadowColor(kWhite), leg1.SetLineColor(kWhite), leg1.SetFillColor(kWhite)
leg1.AddEntry(eff1,"","le")
leg1.AddEntry(prob1,"","l")
leg1.Draw("same")

c2 = TCanvas("eff_vs_nvtx","eff_vs_nvtx")
eff2.SetTitle("True vertex eff.;N_{vtx};#varepsilon #Delta z < 1 cm")
prob2.SetTitle("Aveage vertex prob.;N_{vtx};#varepsilon #Delta z < 1 cm")
eff2.GetYaxis().SetRangeUser(0.5, 1.1)
eff2.GetXaxis().SetRangeUser(0, 35)
eff2.Draw("e0")
prob2.Draw("hist same")

leg2 = TLegend(0.52,0.62,0.82,0.82)
leg2.SetShadowColor(kWhite), leg2.SetLineColor(kWhite), leg2.SetFillColor(kWhite)
leg2.AddEntry(eff2,"","le")
leg2.AddEntry(prob2,"","l")
leg2.Draw("same")

for c in c1, c2:
    for fmt in "C", "png", "pdf":
        c.SaveAs( "%s.%s" % (c.GetName(), fmt) )
