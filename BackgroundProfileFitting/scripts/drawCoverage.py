#!/usr/bin/env python

import sys
import os

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-y","--year")
parser.add_option("-c","--cat",type="int",default=0)
parser.add_option("-v","--expVals",default=[],action="append")
parser.add_option("-f","--files",default=[],action="append")
parser.add_option("-I","--injectVal",type="float")
parser.add_option("-D","--outDir")
parser.add_option("-d","--datfile")
parser.add_option("-m","--runMasses",default=False,action="store_true",help="Default is to run mus")
(options,args)=parser.parse_args()

coverageValues=[0.1,0.5,1.,2.]

import ROOT as r
outf = r.TFile('CoverageFinalPlots.root','RECREATE')

def makePlot():
  
  if options.runMasses:
    valerr = 5 
    valTitle = "m_{H} (GeV)"
    ext = 'mass'
    xlow = 105
    xhigh = 155
  else:
    valerr = 0.25
    valTitle = "#mu_{gen}"
    ext = 'mu'
    xlow = 0.
    xhigh = 2.

  os.system('mkdir -p %s'%options.outDir)

  r.gROOT.SetBatch()

  dummyHist = r.TH1F('dummy','',len(options.expVals),xlow,xhigh)
  dummyHist.GetXaxis().SetTitle(valTitle)
  dummyHist.SetStats(0)

  muBiasBand = r.TGraphErrors()
  muBiasBand.SetLineColor(r.kGray)
  muBiasBand.SetFillColor(r.kGray)
  muBiasBand.SetFillStyle(1001)
  muBiasBand.SetMarkerColor(r.kGray)
  muPullBand = r.TGraphErrors()
  muPullBand.SetLineColor(r.kGray)
  muPullBand.SetFillColor(r.kGray)
  muPullBand.SetFillStyle(1001)
  muPullBand.SetMarkerColor(r.kGray)

  scanTypes=['Fab','Paul','Chi2','AIC']
  graphCol=[r.kBlue,r.kRed,r.kGreen+1,r.kMagenta]
  graphFil=[3144,3190,3002,3002]
  label=['Hgg Nominal Pol','Envelope (no pen)','Envelope (1/dof pen)','Envelope (2/dof pen)']
  plotTypes=['Coverage','Bias','Pull']

  valToFileDict={}
  for i, expM in enumerate(options.expVals):
    valToFileDict[float(expM)] = options.files[i]

  truth_mods=[]
  dummyFile = r.TFile.Open(options.files[0])
  for key in dummyFile.GetListOfKeys():
    name = key.GetName()
    truth = name.split('_mu')[0]
    if truth not in truth_mods: truth_mods.append(truth)
  dummyFile.Close()
  print truth_mods

  canv = r.TCanvas()

  siglines=[]
  sigboxes=[]
  for i, sig in enumerate(coverageValues):
    coverage = 1.-r.TMath.Prob(sig*sig,1)
    siglines.append(r.TLine(xlow,coverage,xhigh,coverage))
    siglines[i].SetLineWidth(2)
    siglines[i].SetLineStyle(2)
    sigboxes.append(r.TPaveText(1.02*xlow,coverage+0.1,1.07*xhigh,coverage))
    sigboxes[i].SetFillColor(0)
    sigboxes[i].SetLineColor(0)
    sigboxes[i].SetTextSize(0.04)
    sigboxes[i].AddText('%3.1f#sigma'%(sig))

  for truth in truth_mods:
    # coverage
    graphDict={}
    legCov = r.TLegend(0.3,0.75,0.59,0.89)
    legCov.SetFillColor(0)
    legCov.SetLineColor(0)
    legCov2 = r.TLegend(0.6,0.75,0.89,0.89)
    legCov2.SetFillColor(0)
    legCov2.SetLineColor(0)
    leg = r.TLegend(0.3,0.75,0.59,0.89)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg2 = r.TLegend(0.6,0.75,0.89,0.89)
    leg2.SetFillColor(0)
    leg2.SetLineColor(0)
    for ptype in plotTypes:
      graphDict[ptype] = {}
      if ptype=='Coverage':
        for c, stype in enumerate(scanTypes):
          graphDict[ptype][stype] = []
          for v, cov in enumerate(coverageValues):
            graphDict[ptype][stype].append(r.TGraphErrors())
            graphDict[ptype][stype][v].SetLineColor(graphCol[c])
            graphDict[ptype][stype][v].SetMarkerColor(graphCol[c])
            graphDict[ptype][stype][v].SetFillColor(graphCol[c])
            graphDict[ptype][stype][v].SetFillStyle(graphFil[c])
          if c<len(scanTypes)/2:
            leg.AddEntry(graphDict[ptype][stype][0],label[c],"F")
            legCov.AddEntry(graphDict[ptype][stype][0],label[c],"L")
          else:
            leg2.AddEntry(graphDict[ptype][stype][0],label[c],"F")
            legCov2.AddEntry(graphDict[ptype][stype][0],label[c],"L")
      else:
        for c, stype in enumerate(scanTypes):
          graphDict[ptype][stype] = r.TGraphErrors()
          graphDict[ptype][stype].SetLineColor(graphCol[c])
          graphDict[ptype][stype].SetMarkerColor(graphCol[c])
          graphDict[ptype][stype].SetFillColor(graphCol[c])
          graphDict[ptype][stype].SetFillStyle(graphFil[c])
    
    val_arr = valToFileDict.keys()
    val_arr.sort()
    for p, val in enumerate(val_arr):
      f = r.TFile.Open(valToFileDict[val])
      for c, stype in enumerate(scanTypes):
        hist = f.Get('%s_mu%s'%(truth,stype))
        if options.runMasses:
          graphDict['Bias'][stype].SetPoint(p,val,(hist.GetMean()-options.injectVal)/hist.GetRMS())
          graphDict['Bias'][stype].SetPointError(p,0,hist.GetMeanError()/hist.GetRMS())
        else:
          graphDict['Bias'][stype].SetPoint(p,val,(hist.GetMean()-val)/hist.GetRMS())
          graphDict['Bias'][stype].SetPointError(p,0,hist.GetMeanError()/hist.GetRMS())
        histPull = f.Get('%s_mu%sPull'%(truth,stype))
        graphDict['Pull'][stype].SetPoint(p,val,histPull.GetMean())
        graphDict['Pull'][stype].SetPointError(p,0,histPull.GetMeanError())
        muBiasBand.SetPoint(p,val-5,0)
        muBiasBand.SetPointError(p,valerr,0.2)
        muPullBand.SetPoint(p,val-5,0)
        muPullBand.SetPointError(p,valerr,0.14)
        for v, cov in enumerate(coverageValues):
          graph = f.Get('%s_mu%sCov%3.1f'%(truth,stype,cov))
          x = r.Double(0.)
          y = r.Double(0.)
          graph.GetPoint(0,x,y)
          yerr = graph.GetErrorY(0)
          graphDict['Coverage'][stype][v].SetPoint(p,val,y)
          graphDict['Coverage'][stype][v].SetPointError(p,valerr,yerr)
   
    muBiasBand.SetPoint(len(options.expVals),155,0)
    muBiasBand.SetPointError(len(options.expVals),valerr,0.2)
    muPullBand.SetPoint(len(options.expVals),155,0)
    muPullBand.SetPointError(len(options.expVals),valerr,0.14)

    if options.runMasses:
      dummyHist.SetTitle('Category=%d  Truth=%s  #mu_{gen}=%3.1f'%(options.cat,truth,options.injectVal))
    else:
      dummyHist.SetTitle('Category=%d  Truth=%s  m_{H}=%3.0f'%(options.cat,truth,options.injectVal))
    dummyHist.GetYaxis().SetRangeUser(0.,1.3)
    dummyHist.GetYaxis().SetTitle("Coverage");
    dummyHist.Draw()
    for sigline in siglines: sigline.Draw("same")
    for sigbox in sigboxes: sigbox.Draw("same")
    for c, stype in enumerate(scanTypes):
      for v, cov in enumerate(coverageValues):
        graphDict['Coverage'][stype][v].Draw("Psame")
        if options.runMasses: graphDict['Coverage'][stype][v].SetName('%s_cat%d_mass_%s_%s_cov%3.1f_coverage'%(options.year,options.cat,truth,stype,cov))
        else: graphDict['Coverage'][stype][v].SetName('%s_cat%d_mu_%s_%s_cov%3.1f_coverage'%(options.year,options.cat,truth,stype,cov))
        outf.cd()
        graphDict['Coverage'][stype][v].Write()
    leg.Draw("same")
    leg2.Draw("same")
    canv.Print('%s/%s_%s_coverage.pdf'%(options.outDir,ext,truth))
    canv.Print('%s/%s_%s_coverage.png'%(options.outDir,ext,truth))
    
    dummyHist.GetYaxis().SetRangeUser(-0.5,0.5)
    dummyHist.GetYaxis().SetTitle("#mu bias (syst/stat)")
    dummyHist.Draw()
    muBiasBand.Draw("E3same")
    dummyHist.Draw("AXISsame")
    for c, stype in enumerate(scanTypes):
      graphDict['Bias'][stype].Draw("E3same")
      if options.runMasses: graphDict['Bias'][stype].SetName('%s_cat%d_mass_%s_%s_bias'%(options.year,options.cat,truth,stype))
      else: graphDict['Bias'][stype].SetName('%s_cat%d_mu_%s_%s_bias'%(options.year,options.cat,truth,stype))
      outf.cd()
      graphDict['Bias'][stype].Write()
    leg.Draw("same")
    leg2.Draw("same")
    #canv.Print('%s/%s_%s_bias.pdf'%(options.outDir,ext,truth))
    #canv.Print('%s/%s_%s_bias.png'%(options.outDir,ext,truth))
  
    dummyHist.GetYaxis().SetTitle("#mu pull (fit-gen)/err")
    dummyHist.Draw()
    muPullBand.Draw("E3same")
    dummyHist.Draw("AXISsame")
    for c, stype in enumerate(scanTypes):
      graphDict['Pull'][stype].Draw("E3same")
      if options.runMasses: graphDict['Pull'][stype].SetName('%s_cat%d_mass_%s_%s_pull'%(options.year,options.cat,truth,stype))
      else: graphDict['Pull'][stype].SetName('%s_cat%d_mu_%s_%s_pull'%(options.year,options.cat,truth,stype))
      outf.cd()
      graphDict['Pull'][stype].Write()
    leg.Draw("same")
    leg2.Draw("same")
    canv.Print('%s/%s_%s_pull.pdf'%(options.outDir,ext,truth))
    canv.Print('%s/%s_%s_pull.png'%(options.outDir,ext,truth))
  
if not options.datfile:
  makePlot()
else:
  f = open(options.datfile)
  for line in f.readlines():
    if line.startswith('#') or line.startswith('\n'): continue
    if line.startswith('year'): 
      options.year = line.split('=')[1].strip('\n')
      continue
    opts = line.split()
    if opts[0].split('=')[1]=='mu':
      options.runMasses=False
    elif opts[0].split('=')[1]=='mass':
      options.runMasses=True
    else:
      sys.exit('Invalid options')
    options.cat = int(opts[1].split('=')[1])
    options.outDir = opts[2].split('=')[1]
    options.injectVal = float(opts[3].split('=')[1])
    options.expVals=[]
    options.files=[]
    for val in opts[4].split('=')[1].split(','):
      options.expVals.append(float(val))
    for file in opts[5].split('=')[1].split(','):
      options.files.append(file)
    makePlot()


  

