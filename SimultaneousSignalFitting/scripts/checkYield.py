#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--infile")
parser.add_option("-p","--procs",default="ggh,vbf,wh,zh,tth")
parser.add_option("-c","--cats",default="0,1,2,3,4,5,6,7,8")
(options,args)=parser.parse_args()

import ROOT as r

r.gROOT.ProcessLine(".L lib/libSimultaneousSignalFit.so")

tf = r.TFile(options.infile)

ws = tf.Get('wsig_8TeV')
MH = ws.var('MH')
intL = ws.var('IntLumi')

MH.Print()
intL.Print()

print 'The output herein should be compared to the output of makeEffAcc.py'

tg = r.TGraph()


for p,m in enumerate(range(110,151,1)):
	MH.setConstant(False)
	MH.setVal(m)
	MH.setConstant(True)
	printLine = 'Signal m%3d.0: '%m
	sumTot = 0.
	for cat in [int(x) for x in options.cats.split(',')]:
		sumProc = 0.
		procLine = 'cat %d, mH=%3d.0: '%(cat,m)
		for proc in options.procs.split(','):
			normFunc = ws.function('hggpdfsmrel_%s_cat%d_norm'%(proc,cat))
			eaFunc = ws.function('fea_%s_cat%d'%(proc,cat))
			procLine += '   %s %.5f'%(proc,100.*eaFunc.getVal())
			sumProc += normFunc.getVal()*intL.getVal()
			sumTot += normFunc.getVal()*intL.getVal()
		print procLine
		printLine += '%3.5f '%sumProc
	printLine += 'tot=%3.5f'%sumTot
	print printLine
	effAccDenom=1.
	fbr = ws.function('fbr')
	effAccDenom *= fbr.getVal()
	xs = 0.
	for proc in options.procs.split(','):
		fxs = ws.function('fxs_%s'%proc)
		xs += fxs.getVal()
	effAccDenom *= xs
	effAccDenom *= intL.getVal()
	tg.SetPoint(p,m,100.*sumTot/effAccDenom)

r.gROOT.SetBatch()
canv = r.TCanvas()
tg.SetLineWidth(3)
tg.Draw("AL")
canv.Update()
canv.Modified()
canv.Print("effAccCrossCheck.pdf")
canv.Print("effAccCrossCheck.png")

tf.Close()
	
