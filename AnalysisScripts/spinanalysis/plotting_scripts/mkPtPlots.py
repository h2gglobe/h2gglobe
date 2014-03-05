"""
KEY: TTree	Data;1	Data
  KEY: TTree	diphojet_sherpa_8TeV;1	diphojet_sherpa_8TeV
	 KEY: TTree	gjet_20_8TeV_pf;1	gjet_20_8TeV_pf
	  KEY: TTree	gjet_40_8TeV_pf;1	gjet_40_8TeV_pf
		 KEY: TTree	qcd_30_8TeV_ff;1	qcd_30_8TeV_ff
		  KEY: TTree	qcd_30_8TeV_pf;1	qcd_30_8TeV_pf
			 KEY: TTree	ggh_m125_8TeV;1	ggh_m125_8TeV
			  KEY: TTree	vbf_m125_8TeV;1	vbf_m125_8TeV
				 KEY: TTree	wzh_m125_8TeV;1	wzh_m125_8TeV
				  KEY: TTree	tth_m125_8TeV;1	tth_m125_8TeV
					 KEY: TTree	ggh_m126_8TeV;1	ggh_m126_8TeV
					  KEY: TTree	vbf_m126_8TeV;1	vbf_m126_8TeV
						 KEY: TTree	wzh_m126_8TeV;1	wzh_m126_8TeV
						  KEY: TTree	tth_m126_8TeV;1	tth_m126_8TeV
							 KEY: TTree	gg_grav_m125_8TeV;1	gg_grav_m125_8TeV
							  KEY: TTree	qq_grav_m125_8TeV;1	qq_grav_m125_8TeV
								 KEY: TTree	gg_grav_m126_8TeV;1	gg_grav_m126_8TeV
								  KEY: TTree	qq_grav_m126_8TeV;1	qq_grav_m126_8TeV
									"""

#!/usr/bin/env python

import ROOT as r
import sys

if len(sys.argv)!=4:
	sys.exit('usage mkPtPlots.py <infile> <outfile> <sqrts>')

infile = r.TFile(sys.argv[1])
outfile = r.TFile(sys.argv[2],"RECREATE")
sqrts = int(sys.argv[3])

diphojet = infile.Get("spin_trees/diphojet_sherpa_%dTeV"%sqrts)
gjet20pf = infile.Get("spin_trees/gjet_20_%dTeV_pf"%sqrts)
gjet40pf = infile.Get("spin_trees/gjet_40_%dTeV_pf"%sqrts)
qcd30ff = infile.Get("spin_trees/qcd_30_%dTeV_ff"%sqrts)
qcd30pf = infile.Get("spin_trees/qcd_30_%dTeV_pf"%sqrts)

ggh = infile.Get("spin_trees/ggh_m125_%dTeV"%sqrts)
vbf = infile.Get("spin_trees/vbf_m125_%dTeV"%sqrts)
gggrav = infile.Get("spin_trees/gg_grav_m125_%dTeV"%sqrts)
qqgrav = infile.Get("spin_trees/qq_grav_m125_%dTeV"%sqrts)

data = infile.Get("spin_trees/Data")

dataLeadJPt = r.TH1F("dataLeadJPt","",100,0,500)
bkgLeadJPt = r.TH1F("bkgLeadJPt","",100,0,500)
gghLeadJPt = r.TH1F("gghLeadJPt","",100,0,500)
vbfLeadJPt = r.TH1F("vbfLeadJPt","",100,0,500)
gggravLeadJPt = r.TH1F("gggravLeadJPt","",100,0,500)
qqgravLeadJPt = r.TH1F("qqgravLeadJPt","",100,0,500)

dataHiggsPt = r.TH1F("dataHiggsPt","",500,0,500)
bkgHiggsPt = r.TH1F("bkgHiggsPt","",500,0,500)
gghHiggsPt = r.TH1F("gghHiggsPt","",500,0,500)
vbfHiggsPt = r.TH1F("vbfHiggsPt","",500,0,500)
gggravHiggsPt = r.TH1F("gggravHiggsPt","",500,0,500)
qqgravHiggsPt = r.TH1F("qqgravHiggsPt","",500,0,500)

bkgTrees=[]
bkgTrees.append(diphojet)
bkgTrees.append(gjet20pf)
bkgTrees.append(gjet40pf)
bkgTrees.append(qcd30ff)
bkgTrees.append(qcd30pf)

dataTrees=[]
dataTrees.append(data)

gghTrees=[]
gghTrees.append(ggh)

vbfTrees=[]
vbfTrees.append(vbf)

gggravTrees=[]
gggravTrees.append(gggrav)

qqgravTrees=[]
qqgravTrees.append(qqgrav)

def fillHistHiggsPt(hist,trees,isData=False):
  for tree in trees:
    for i in range(tree.GetEntries()):
      if i%10000==0: print i,'/',tree.GetEntries()
      tree.GetEntry(i)
      h = r.TLorentzVector(tree.higgs_px,tree.higgs_py,tree.higgs_pz,tree.higgs_E)
      if isData and h.M()>120. and h.M()<130.: continue
      hist.Fill(h.Pt(),tree.evweight)

def fillHistLeadJPt(hist,trees,isData=False):
  for tree in trees:
    for i in range(tree.GetEntries()):
      if i%10000==0: print i,'/',tree.GetEntries()
      tree.GetEntry(i)
      h = r.TLorentzVector(tree.higgs_px,tree.higgs_py,tree.higgs_pz,tree.higgs_E)
      if isData and h.M()>120. and h.M()<130.: continue
      hist.Fill(tree.myVBFLeadJPt,tree.evweight)

fillHistHiggsPt(dataHiggsPt,dataTrees)
#fillHistHiggsPt(bkgHiggsPt,bkgTrees)
fillHistHiggsPt(gghHiggsPt,gghTrees)
fillHistHiggsPt(vbfHiggsPt,vbfTrees)
fillHistHiggsPt(gggravHiggsPt,gggravTrees)
fillHistHiggsPt(qqgravHiggsPt,qqgravTrees)

fillHistLeadJPt(dataLeadJPt,dataTrees)
#fillHistLeadJPt(bkgLeadJPt,bkgTrees)
fillHistLeadJPt(gghLeadJPt,gghTrees)
fillHistLeadJPt(vbfLeadJPt,vbfTrees)
fillHistLeadJPt(gggravLeadJPt,gggravTrees)
fillHistLeadJPt(qqgravLeadJPt,qqgravTrees)

dataHiggsPt.SetLineColor(r.kBlack)
bkgHiggsPt.SetFillColor(r.kGreen-7)
gghHiggsPt.SetLineColor(r.kBlue)
vbfHiggsPt.SetLineColor(r.kRed)
gggravHiggsPt.SetLineColor(r.kMagenta)
qqgravHiggsPt.SetLineColor(r.kGreen+2)

dataLeadJPt.SetLineColor(r.kBlack)
bkgLeadJPt.SetFillColor(r.kGreen-7)
gghLeadJPt.SetLineColor(r.kBlue)
vbfLeadJPt.SetLineColor(r.kRed)
gggravLeadJPt.SetLineColor(r.kMagenta)
qqgravLeadJPt.SetLineColor(r.kGreen+2)

outfile.cd()
bkgHiggsPt.Write()
dataHiggsPt.Write()
gghHiggsPt.Write()
vbfHiggsPt.Write()
gggravHiggsPt.Write()
qqgravHiggsPt.Write();
bkgLeadJPt.Write()
dataLeadJPt.Write()
gghLeadJPt.Write()
vbfLeadJPt.Write()
gggravLeadJPt.Write()
qqgravLeadJPt.Write();

outfile.Close()
infile.Close()
raw_input("Ok?")


