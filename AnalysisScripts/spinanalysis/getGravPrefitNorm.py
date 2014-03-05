#!/usr/bin/env python
import os
import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s","--datacardSM")
parser.add_option("-g","--datacardGrav")
(options,args)=parser.parse_args()

import ROOT as r

fitSM = r.TGraphErrors()
fitGrav = r.TGraphErrors()
fitRatio = r.TGraphErrors()

#for p, mu2 in enumerate([0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]):
for p, mu2 in enumerate([1.0]):
  print '---------------------------'
  print '---------', mu2, '---------'
  print '---------------------------'
  os.system('combine %s -M GenerateOnly -m 125 -t -1 --expectSignal=%3.1f -n mu%3.1f --saveToys -s 0'%(options.datacardGrav,mu2,mu2))
  os.system('combine %s -M MaxLikelihoodFit -m 125 -t -1 --toysFile higgsCombinemu%3.1f.GenerateOnly.mH125.0.root --rMin=-3. -n GravFit%3.1f'%(options.datacardGrav,mu2,mu2))
  os.system('combine %s -M MaxLikelihoodFit -m 125 -t -1 --toysFile higgsCombinemu%3.1f.GenerateOnly.mH125.0.root --rMin=-3. -n SMFit%3.1f'%(options.datacardSM,mu2,mu2))

  smF = r.TFile('higgsCombineSMFit%3.1f.MaxLikelihoodFit.mH125.root'%mu2)
  smT = smF.Get('limit')
  smT.GetEntry(0) 
  smVal = smT.limit
  smT.GetEntry(1) 
  smLow = smVal-smT.limit
  smT.GetEntry(2) 
  smHigh = smT.limit-smVal
  smErr = (smHigh+smLow)/2.
  fitSM.SetPoint(p,mu2,smVal)
  fitSM.SetPointError(p,0.,smErr)
  smF.Close()

  gravF = r.TFile('higgsCombineGravFit%3.1f.MaxLikelihoodFit.mH125.root'%mu2)
  gravT = gravF.Get('limit')
  gravT.GetEntry(0) 
  gravVal = gravT.limit
  gravT.GetEntry(1) 
  gravLow = gravVal-gravT.limit
  gravT.GetEntry(2) 
  gravHigh = gravT.limit-gravVal
  gravErr = (gravHigh+gravLow)/2.
  fitSM.SetPoint(p,mu2,gravVal)
  fitSM.SetPointError(p,0.,gravErr)
  gravF.Close()
  
  rat = smVal/gravVal
  ratErr = rat*((gravErr/gravVal)**2+(smErr/smVal)**2)**0.5
  fitRatio.SetPoint(p,mu2,rat)
  fitRatio.SetPointError(p,0.,ratErr)

fitSM.SetLineColor(r.kRed)
fitGrav.SetLineColor(r.kBlue)

fitRatio.SetMarkerStyle(r.kFullCircle)
fitRatio.SetLineWidth(2)
fitRatio.GetXaxis().SetTitle("#mu_{2}^{true}")
fitRatio.GetYaxis().SetTitle("r_{0}")

fitR = fitRatio.Fit("pol0","S")
box = r.TLatex()
box.SetNDC()

canv = r.TCanvas()
fitRatio.Draw("ALEP")
box.DrawLatex(0.6,0.8,'Fit p0 = %3.2f +/- %3.2f'%(fitR.Value(0),fitR.ParError(0))) 
canv.Update()
raw_input("Ok?")
canv.Print("norm.pdf")

