import ROOT
import sys

from ROOT import TFile, TCanvas, TStyle, TApplication, TLegend, TGraph, TGraphErrors

# ROOT Setup
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)

vars = ["leadPt",
        "subleadPt",
        "diphoPtOverM",
        "leadJPt",
        "subleadJPt",
        "MJJ",
        "deltaEtaJJ",
        "Zep",        
        "deltaPhiJJ",
        "deltaPhiJJGamGam",
        "deltaPhiGamGam",
        "thetaJ1",
        "thetaJ2",
        "etaJJ"]

samples = ["ggH_m125_8TeV","ggHj_m125_8TeV","ggHjj_m125_8TeV"]

fin = ["root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/vbf_spin_studies/CMS-HGG_powheg_ggH.root",
       "root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/vbf_spin_studies/CMS-HGG_powheg_ggHj.root",
       "root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/vbf_spin_studies/CMS-HGG_powheg_ggHjj.root"]

bins=[50, 50, 50, 50, 50, 50,50,50,25,50,50,50,50,50]
low =[0, 0, 0 , 0 , 0 , 0 , 0, 0, 0, 0, 0, 0, 0, 0 ,0,0]
high=[300, 300, 2, 300, 300, 1000, 10, 10, 3.5, 3.5, 3.5, 3.5, 3.5,5]

file=[]
canvas=[]

for isample,f in enumerate(fin):
    file = TFile.Open(f)
    tree = file.Get(samples[isample])
    for jvar,v in enumerate(vars):
        hname = ('h'+v+'_'+samples[isample]).replace('_m125_8TeV','')
        #print hname
        h = ROOT.TH1F(hname,hname,bins[jvar],low[jvar],high[jvar])
        h.SetLineWidth(2)
        todraw = ' abs(' + v + ' )' 
        #print todraw
        tree.Project(hname,todraw,'category==4 || category==5 ')
        #tree.Project(hname,todraw,'category==6')
            
        if (isample < 1):
            cname = hname.replace('h','c')
            #print cname
            canvas.append(TCanvas(cname,cname,500,500))
            max=h.GetMaximum()
            h.GetYaxis().SetRangeUser(0,max*1.5)
            h.DrawNormalized()
            if (jvar < 1):
                legend = ROOT.TLegend(0.6,0.7,0.85,0.85)
                legend.SetFillColor(0)
                legend.AddEntry(h,samples[isample].replace('_m125_8TeV',''),"L")
        else:
            canvas[jvar].cd()
            h.SetLineColor(isample*2)
            h.DrawNormalized("same")
            if (jvar < 1):
                legend.AddEntry(h,samples[isample].replace('_m125_8TeV',''),"L")
            legend.Draw("same")
            canvas[jvar].Update()
            if (isample == 2 ):
                out = (h.GetName().replace('_ggHjj',''))+'.png'
                canvas[jvar].SaveAs(out)
                       
s = raw_input("OK?")
