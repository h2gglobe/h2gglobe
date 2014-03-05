#!/usr/bin/env python

import ROOT as r
r.gROOT.SetBatch()

import sys
tf = r.TFile(sys.argv[1])

canv = r.TCanvas()

th2 = r.TH2F('h','',1,0,1.,1,0.4,1.8)
th2.SetStats(0)
th2.SetLineColor(0)
th2.GetXaxis().SetTitle("|cos(#theta*)|")
th2.GetYaxis().SetTitle("#sigma/#sigma_{SM}")
th2.GetXaxis().SetLabelSize(0.05);
th2.GetXaxis().SetTitleSize(0.05);
th2.GetYaxis().SetLabelSize(0.05);
th2.GetYaxis().SetTitleSize(0.05);
th2.GetXaxis().SetTitleOffset(1.0);
th2.GetYaxis().SetTitleOffset(1.0);

hcut_vals = [50,75,100,125,150,200,-1]
jcut_vals = [50,75,100,125,150,200,-1]
color = [r.kMagenta+2,r.kMagenta-4,r.kBlue,r.kBlue-7,r.kGreen+1,r.kOrange-3,r.kRed]

leg = r.TLegend(0.84,0.4,0.99,0.99)
leg.SetHeader('Cut (GeV)')
leg.SetFillColor(0)
leg.SetTextSize(0.04)

leg2 = r.TLegend(0.12,0.6,0.44,0.89)
leg2.SetFillColor(0)
smDum = r.TLine()
smDum.SetLineWidth(3)
smDum.SetLineStyle(r.kDashed)
gravDum = r.TLine()
gravDum.SetLineWidth(3)
leg2.AddEntry(smDum,'X#rightarrow#gamma#gamma 0^{+}','L')
leg2.AddEntry(gravDum,'X#rightarrow#gamma#gamma 2_{m}^{+} (100\%gg)','L')

line = r.TLine()
line.SetLineWidth(2)

th2.SetTitle('Higgs p_{T} cut')
th2.Draw("AXIS")
#line.DrawLine(0.0,1.0,1.0,1.0)

for p,cut in enumerate(reversed(hcut_vals)):
	SM = tf.Get('hcut%d_jcut-1ChannelCompSM'%cut)
	GravGG = tf.Get('hcut%d_jcut-1ChannelCompGravGG'%cut)

	x = r.Double(0)
	y = r.Double(0)
	for point in range(GravGG.GetN()):
		GravGG.GetPoint(point,x,y)
		GravGG.SetPoint(point,x,y*1.064)
	
	SM.SetLineStyle(r.kDashed)
	SM.SetLineColor(color[p])
	GravGG.SetLineColor(color[p])
	SM.SetLineWidth(3)
	GravGG.SetLineWidth(3)

	legname = '<'+str(cut)
	if cut<0.: legname='None'
	leg.AddEntry(GravGG,legname,'L')

	SM.Draw('Lsame')
	GravGG.Draw('Lsame')

leg.Draw('same')
leg2.Draw('same')
canv.Print('hcutDiffMu.pdf')
	
th2.SetTitle('Lead jet p_{T} cut')
th2.Draw("AXIS")
#line.DrawLine(0.0,1.0,1.0,1.0)

for p,cut in enumerate(reversed(jcut_vals)):
	SM = tf.Get('hcut-1_jcut%dChannelCompSM'%cut)
	GravGG = tf.Get('hcut-1_jcut%dChannelCompGravGG'%cut)

	x = r.Double(0)
	y = r.Double(0)
	for point in range(GravGG.GetN()):
		GravGG.GetPoint(point,x,y)
		GravGG.SetPoint(point,x,y*1.064)
	
	SM.SetLineStyle(r.kDashed)
	SM.SetLineColor(color[p])
	GravGG.SetLineColor(color[p])
	SM.SetLineWidth(3)
	GravGG.SetLineWidth(3)

	SM.Draw('Lsame')
	GravGG.Draw('Lsame')

leg.Draw('same')
leg2.Draw('same')
canv.Print('jcutDiffMu.pdf')
	
