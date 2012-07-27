#!/bin/env python

import sys
import ROOT

def readPlot(fin, name, cat=0, which=["_sequential","_nminus1"], samples=["diphojet_8TeV","ggh_m125_8TeV","vbf_m125_8TeV",]):

    ret = []
    for p in which:
        hists = []
        for s in samples:
            nam = "%s%s_cat%d_%s" % ( name, p, cat, s )
            h = fin.Get(nam)
            ## print nam, h
            hists.append(h)
        ret.append(hists)
            
    return ret


def selectionControlPlots(fname):
    vars=[("vbf_mva",               [""]),
          ("cut_VBFLeadJPt",        ["_sequential","_nminus1"]),
          ("cut_VBFSubJPt",         ["_sequential","_nminus1"]),
          ("cut_VBF_dEta",          ["_sequential","_nminus1"]),
          ("cut_VBF_Zep",           ["_sequential","_nminus1"]),
          ("cut_VBF_Mjj",           ["_sequential","_nminus1"]),
          ("cut_VBF_dPhi",          ["_sequential","_nminus1"]),
          ("cut_VBF_DiPhoPtOverM",  ["_sequential","_nminus1"]),
          ("cut_VBF_LeadPhoPtOverM",["_sequential","_nminus1"]),
          ("cut_VBF_SubPhoPtOverM", ["_sequential","_nminus1"]),
          ]
    
    
    fin = ROOT.TFile.Open(fname)
    
    objs=[]
    for var,which in vars:
        plots = readPlot(fin,var,which=which)
        
        c = ROOT.TCanvas(var,var,800,400)
        objs.append(c)
        np = len(plots)
        if np > 1:
            c.Divide(2,np/2)
        for ip in range(np):
            ps=plots[ip]
            
            ps[0].SetLineColor(ROOT.kRed)
            ps[1].SetLineColor(ROOT.kMagenta)

            if np>1:
                c.cd(ip+1)

            p0=ps[0].Clone("%s_%d" % (var,ip) )
            objs.append(p0)
            p0.Reset("ICE")
            p0.SetEntries(0)
            p0.Draw()
            if var == "vbf_mva":
                p0.GetXaxis().SetRangeUser(-1,1)
            for p in ps:
                if var == "vbf_mva":
                    p.SetBinContent(1,0.)
                    p.Draw("hist same")
                else:
                    p.DrawNormalized("hist same")
            
            maximum=0.
            for p in ROOT.gPad.GetListOfPrimitives():
                if p.ClassName().startswith("TH1"):
                    pmax = p.GetMaximum()
                    if pmax>maximum:
                        maximum=pmax
            p0.GetYaxis().SetRangeUser(0,maximum*1.5)
            print var,p0,maximum, p0.GetYaxis().GetXmax()
            ROOT.gPad.Modified()
            ROOT.gPad.Update()
            
            for fmt in "C","png":
                c.SaveAs("%s.%s" % (c.GetName(),fmt))
                
    return objs

def eventYield(filenames,categories=[5,6],procs=["ggh_m125_8TeV","vbf_m125_8TeV","diphojet_8TeV"]):
    files = [ ROOT.TFile.Open(f) for f in filenames ]

    
    for fin in files:
        print fin.GetName()
        for cat in categories:
            plots=readPlot(fin, "mass", cat=cat, which=[""], samples=procs)[0]
            print "cat%d " % cat-1,
            for iproc in range(len(procs)):
                ## print "%1.3g" % plots[iproc].Integral(18,28), 115-135
                print "%1.3g" % plots[iproc].Integral(18,28),
            print


if __name__ == "__main__":
    
    
    objs=selectionControlPlots(sys.argv[1])

    ## eventYield(sys.argv[1:])
