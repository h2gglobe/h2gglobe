#!/usr/bin/env python

import os
import sys

xs={}
xs['7TeV'] = {'ggh':15.13,'vbf':1.222,'wzh':0.5785+0.3351,'tth':0.08632}
xs['8TeV'] = {'ggh':19.27,'vbf':1.578,'wzh':0.7046+0.4153,'tth':0.1293}
br = 8.57E-02

if len(sys.argv)!=3 and len(sys.argv)!=4 and len(sys.argv)!=5:
	print 'usage makeCosThetaPlots.py <before_file> <after_file> <sqrtS> <binedges>'
	sys.exit()

if len(sys.argv)==4: sqrtS = int(sys.argv[3])
else: sqrtS = 8

if len(sys.argv)==5: binedges = [float(x) for x in sys.argv[4].split(',')]
else: binedges = 0.2,0.375,0.55,0.75

import ROOT as r

def getWeight(type):
	if type==0 or type==4 or type==5:
		return xs['%dTeV'%sqrtS]['ggh']*br
	elif type==1:
		return xs['%dTeV'%sqrtS]['vbf']*br
	elif type==2:
		return xs['%dTeV'%sqrtS]['wzh']*br
	elif type==3:
		return xs['%dTeV'%sqrtS]['tth']*br
	else:
		print 'WARNING: type',type,'not recognised'
		return -1

def getCiCCat(lead_calo_eta,lead_r9,sublead_calo_eta,sublead_r9):
	category=-1
	if r.TMath.Max(r.TMath.Abs(lead_calo_eta),r.TMath.Abs(sublead_calo_eta))<1.444:
		if r.TMath.Min(lead_r9,sublead_r9)>0.94:
			category=0
		else:
			category=1
	else:
		if r.TMath.Min(lead_r9,sublead_r9)>0.94:
			category=2
		else:
			category=3
	assert(category!=-1)
	return category

def drawBinEdges(min,max):
	for edge in binedges:
		line2.SetLineStyle(r.kDashed)
		line2.SetLineWidth(3)
		line2.DrawLine(edge,min,edge,max)

r.gROOT.SetBatch()
r.TH1.SetDefaultSumw2()
r.gStyle.SetOptStat(0)

nbins=25

fBefore = r.TFile(sys.argv[1])
fAfter = r.TFile(sys.argv[2])
fOut = r.TFile('tempOut.root','RECREATE')

outProcs=['sm','gg_grav','qq_grav']
titles = ['0^{+} (SM)','2_{m}^{+} (gg)', '2_{m}^{+} (q#bar{q})']
colors = [r.kRed, r.kBlue, r.kGreen+1]
markers = [r.kFullCircle, r.kFullSquare, r.kFullTriangleUp]
beforeProcs = {'sm' : [0,1,2,3], 'gg_grav' : [4], 'qq_grav' : [5] }
afterProcs = {'sm' : ['ggh','vbf','wzh','tth'], 'gg_grav': ['gg_grav'], 'qq_grav': ['qq_grav'] }

histBefore={}
histAfter={}
histAfterCats={}

for proc in outProcs:
	histBefore[proc] = r.TH1F('before_%s'%proc,'',nbins,0,1)
	histAfter[proc] = r.TH1F('after_%s'%proc,'',nbins,0,1)
	histAfterCats[proc] = []
	for cat in range(4):
		histAfterCats[proc].append(r.TH1F('after_%s_cat%d'%(proc,cat),'',nbins,0,1))

treeBefore = fBefore.Get('CosThetaTree')
for ev in range(treeBefore.GetEntries()):
	treeBefore.GetEvent(ev)
	weight = getWeight(treeBefore.type)
	if treeBefore.type<4:
		histBefore['sm'].Fill(r.TMath.Abs(treeBefore.cosThetaCS),weight)
	elif treeBefore.type==4:
		histBefore['gg_grav'].Fill(r.TMath.Abs(treeBefore.cosThetaCS),weight)
	elif treeBefore.type==5:
		histBefore['qq_grav'].Fill(r.TMath.Abs(treeBefore.cosThetaCS),weight)
	else: print 'WARNING: type', type, 'not recognised'

for proc, types in afterProcs.items():
	for type in types:
		treeAfter = fAfter.Get('spin_trees/%s_m125_%dTeV'%(type,sqrtS))
		for ev in range(treeAfter.GetEntries()):
			treeAfter.GetEvent(ev)
			histAfter[proc].Fill(r.TMath.Abs(treeAfter.costheta_cs),treeAfter.evweight)
			cat = getCiCCat(treeAfter.lead_calo_eta,treeAfter.lead_r9,treeAfter.sublead_calo_eta,treeAfter.sublead_r9)
			histAfterCats[proc][cat].Fill(r.TMath.Abs(treeAfter.costheta_cs),treeAfter.evweight)

# effAcc plot
canv = r.TCanvas()
histAfterEffAcc = {}
for proc, hist in histAfterCats.items():
	histAfterEffAcc[proc] = histAfter[proc].Clone('%s_effacc'%histAfter[proc].GetName())
	histAfterEffAcc[proc].Divide(histBefore[proc])
	histAfterEffAcc[proc].Scale(1./histAfterEffAcc[proc].Integral())
	for cat in range(4):
		histAfterCats[proc][cat].Divide(histBefore[proc])
		histAfterCats[proc][cat].Scale(1./histAfterCats[proc][cat].Integral())

effAccRatio=[]
catmax=0.
for cat in range(4):
	effAccRatio.append(histAfterCats['gg_grav'][cat].Clone('gg_grav_sm_rat_cat%d'%cat))
	effAccRatio[cat].Divide(histAfterCats['sm'][cat])
	catmax = r.TMath.Max(catmax,effAccRatio[cat].GetMaximum())

text = r.TLatex()
text.SetNDC()
text.SetTextSize(0.03)
line2 = r.TLine()

catColors = [r.kRed,r.kBlue+1,r.kGreen+2,r.kMagenta]
catleg = r.TLegend(0.65,0.7,0.89,0.89)
catleg.SetFillColor(0)
for cat in range(4):
	effAccRatio[cat].GetYaxis().SetRangeUser(0.,2.)
	effAccRatio[cat].SetMarkerColor(catColors[cat])
	effAccRatio[cat].SetLineColor(catColors[cat])
	effAccRatio[cat].SetLineWidth(3)
	effAccRatio[cat].SetMarkerSize(0.1)
	effAccRatio[cat].GetXaxis().SetTitle('|cos(#theta*)|')
	effAccRatio[cat].GetYaxis().SetTitle('Eff #times Acc ratio to SM')
	catleg.AddEntry(effAccRatio[cat],"Category %d"%cat,"LEP")
	if cat==0: 
		effAccRatio[cat].Draw("LEP")
		line = r.TLine(0.,1.,1.,1.)
		line.Draw("same")
		drawBinEdges(0.,2.)
	effAccRatio[cat].Draw("LEPsame")
catleg.Draw("same")
text.DrawLatex(0.1,0.92,"CMS Preliminary Simulation #sqrt{s} = %dTeV"%sqrtS)
canv.Print('effacccats_%dTeV.pdf'%sqrtS)
canv.Print('effacccats_%dTeV.png'%sqrtS)

leg = r.TLegend(0.15,0.7,0.4,0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)
beforeMax=0.
afterMax=0.
effaccMax=0.
for p, proc in enumerate(outProcs):
	histBefore[proc].Scale(1./histBefore[proc].Integral())
	histBefore[proc].GetYaxis().SetTitle('a.u.')
	histBefore[proc].GetXaxis().SetTitle('|cos(#theta*)|')
	histBefore[proc].SetMarkerColor(colors[p])
	histBefore[proc].SetLineColor(colors[p])
	histBefore[proc].SetMarkerStyle(markers[p])
	histBefore[proc].SetMarkerSize(0.8)
	histAfter[proc].Scale(1./histAfter[proc].Integral())
	histAfter[proc].GetYaxis().SetTitle('a.u.')
	histAfter[proc].GetXaxis().SetTitle('|cos(#theta*)|')
	histAfter[proc].SetMarkerColor(colors[p])
	histAfter[proc].SetLineColor(colors[p])
	histAfter[proc].SetMarkerStyle(markers[p])
	histAfter[proc].SetMarkerSize(0.8)
	histAfterEffAcc[proc].GetYaxis().SetTitle('eff acc')
	histAfterEffAcc[proc].GetXaxis().SetTitle('|cos(#theta*)|')
	histAfterEffAcc[proc].SetMarkerColor(colors[p])
	histAfterEffAcc[proc].SetLineColor(colors[p])
	histAfterEffAcc[proc].SetMarkerStyle(markers[p])
	histAfterEffAcc[proc].SetMarkerSize(0.8)
	leg.AddEntry(histAfter[proc],titles[p],'LEP')
	if histBefore[proc].GetMaximum()>beforeMax: beforeMax=histBefore[proc].GetMaximum()
	if histAfter[proc].GetMaximum()>afterMax: afterMax=histAfter[proc].GetMaximum()
	if histAfterEffAcc[proc].GetMaximum()>effaccMax: effaccMax=histAfterEffAcc[proc].GetMaximum()

	fOut.cd()
	histBefore[proc].Write()
	histAfter[proc].Write()

# before
for p, proc in enumerate(outProcs):
	histBefore[proc].GetYaxis().SetRangeUser(0.,beforeMax*1.1)
	if p==0: 
		histBefore[proc].Draw('LEP')
		drawBinEdges(0.,beforeMax*1.1)
	histBefore[proc].Draw('LEPsame')
leg.Draw('same')
text.DrawLatex(0.1,0.92,"CMS Preliminary Simulation #sqrt{s} = %dTeV"%sqrtS)
canv.Modified()
canv.Print('before_%dTeV.pdf'%sqrtS)
canv.Print('before_%dTeV.png'%sqrtS)
# after
for p, proc in enumerate(outProcs):
	histAfter[proc].GetYaxis().SetRangeUser(0.,afterMax*1.2)
	if p==0: 
		histAfter[proc].Draw('LEP')
		drawBinEdges(0.,afterMax*1.2)
	histAfter[proc].Draw('LEPsame')
leg.SetX1NDC(0.64)
leg.SetX2NDC(0.89)
leg.SetY1NDC(0.7)
leg.SetY2NDC(0.89)
leg.Draw('same')
text.DrawLatex(0.1,0.92,"CMS Preliminary Simulation #sqrt{s} = %dTeV"%sqrtS)
canv.Modified()
canv.Print('after_%dTeV.pdf'%sqrtS)
canv.Print('after_%dTeV.png'%sqrtS)
# effacc
for p, proc in enumerate(outProcs):
	histAfter[proc].GetYaxis().SetRangeUser(0.,effaccMax*1.2)
	if p==0: 
		histAfterEffAcc[proc].Draw('LEP')
		drawBinEdges(0.,effaccMax*1.2)
	histAfterEffAcc[proc].Draw('LEPsame')
leg.SetX1NDC(0.64)
leg.SetX2NDC(0.89)
leg.SetY1NDC(0.7)
leg.SetY2NDC(0.89)
leg.Draw('same')
text.DrawLatex(0.1,0.92,"CMS Preliminary Simulation #sqrt{s} = %dTeV"%sqrtS)
canv.Modified()
canv.Print('after_effacc_%dTeV.pdf'%sqrtS)
canv.Print('after_effacc_%dTeV.png'%sqrtS)



