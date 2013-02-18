#!/usr/bin/env python

import os
import numpy

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f","--file",dest="files",default=[],action="append",help="Add a file")
parser.add_option("-c","--color",dest="colors",default=[],action="append",help="Set color")
parser.add_option("-s","--style",dest="styles",default=[],action="append",help="Set style")
parser.add_option("-w","--width",dest="widths",default=[],action="append",help="Set width")
parser.add_option("-n","--name",dest="names",default=[],action="append",help="Set name")
parser.add_option("-t","--text",dest="text",type="string",default="",help="Add Text")
parser.add_option("-e","--expected",dest="expected",type="string",default="",help="Expected only")
parser.add_option("","--limit",dest="limit",default=False,action="store_true",help="Do limit plot")
parser.add_option("","--pval",dest="pval",default=False,action="store_true",help="Do p-value plot")
parser.add_option("","--maxlh",dest="maxlh",default=False,action="store_true",help="Do best fit mu plot")
parser.add_option("-v","--verbose",dest="verbose",default=False,action="store_true")
(options,args)=parser.parse_args()

if not options.pval and not options.maxlh and not options.limit:
  sys.exit('No option set. Must set either --pval or --maxlh or --limit')

if options.verbose:
  print options.files
  print options.colors
  print options.styles
  print options.widths
  print options.names

import ROOT as r

r.gROOT.SetBatch()

mg = r.TMultiGraph()
if options.pval: leg = r.TLegend(0.6,0.12,0.89,0.4)
elif options.maxlh: leg = r.TLegend(0.6,0.8,0.89,0.89)
elif options.limit: leg = r.TLegend(0.6,0.7,0.89,0.89)
leg.SetFillColor(0)
canv = r.TCanvas()
canv.SetGrid(True)
if options.pval: canv.SetLogy()

# sigma lines
sigmas=[1,2,3,4,5]
lines=[]
labels=[]
for i,sig in enumerate(sigmas):
  y = r.RooStats.SignificanceToPValue(sig)
  lines.append(r.TLine(110,y,150,y))
  lines[i].SetLineWidth(2)
  lines[i].SetLineStyle(2)
  lines[i].SetLineColor(r.kRed)
  labels.append(r.TLatex(110 + 2, y * 1.1, "%d #sigma" % (i+1)))
  labels[i].SetTextAlign(11);
  
for k, f in enumerate(options.files):
  tf = r.TFile(f)
  tree = tf.Get('limit')
  values=[]
  for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    values.append([tree.mh, tree.limit])
  
  values.sort(key=lambda x: x[0])
  if options.pval: graph = r.TGraph()
  elif options.maxlh: graph = r.TGraphAsymmErrors()
  elif options.limit:
    graph = r.TGraphAsymmErrors()
    exp = r.TGraphAsymmErrors()
    oneSigma = r.TGraphAsymmErrors()
    twoSigma = r.TGraphAsymmErrors()

  point_counter=0
  for j in range(len(values)):
    if options.pval: graph.SetPoint(j,values[j][0],values[j][1])
    elif options.maxlh:
      if (j%4==0):
        mh = values[j][0]
        fit = values[j][1]
        low = values[j+1][1]
        high = values[j+2][1]
        graph.SetPoint(point_counter,mh,fit)
        graph.SetPointError(point_counter,0,0,abs(fit-low),abs(high-fit))
        point_counter+=1
        if options.verbose: print mh, fit, low, high
    elif options.limit:
      if (j%6==0):
        mh = values[j][0]
        down95 = values[j][1]
        down68 = values[j+1][1]
        median = values[j+2][1]
        up68 = values[j+3][1]
        up95 = values[j+4][1]
        obs = values[j+5][1]
        graph.SetPoint(point_counter,mh,obs)
        exp.SetPoint(point_counter,mh,median)
        oneSigma.SetPoint(point_counter,mh,median)
        oneSigma.SetPointError(point_counter,0,0,abs(median-down68),abs(up68-median))
        twoSigma.SetPoint(point_counter,mh,median)
        twoSigma.SetPointError(point_counter,0,0,abs(median-down95),abs(up95-median))
        point_counter+=1
        if options.verbose: print mh, median, down68, up68, down95, up95, obs

  if options.pval:
    if len(options.colors)>0:
      graph.SetLineColor(int(options.colors[k]))
    if len(options.styles)>0:
      graph.SetLineStyle(int(options.styles[k]))
    if len(options.widths)>0:
      graph.SetLineWidth(int(options.widths[k]))
    if len(options.names)>0:
      if options.names[k]!="-1": leg.AddEntry(graph,options.names[k],'L')
  elif options.maxlh:
    graph.SetMarkerStyle(21)
    graph.SetMarkerSize(0.5)
    graph.SetLineColor(1)
    graph.SetLineWidth(2)
    graph.SetFillColor(r.kGreen)
    leg.AddEntry(graph,'68% CL Band','F')
  elif options.limit:
    graph.SetMarkerStyle(21)
    graph.SetMarkerSize(0.5)
    graph.SetLineColor(1)
    graph.SetLineWidth(2)
    exp.SetLineColor(1)
    exp.SetLineStyle(2)
    oneSigma.SetLineStyle(2)
    twoSigma.SetLineStyle(2)
    oneSigma.SetFillColor(r.kYellow)
    twoSigma.SetFillColor(r.kGreen)
    if not options.expected: leg.AddEntry(graph,'Observed','L')
    #leg.AddEntry(exp,'Expected','L')
    leg.AddEntry(oneSigma,'Expected #pm 1#sigma','FL') 
    leg.AddEntry(twoSigma,'Expected #pm 2#sigma','FL') 
    
  if options.limit:
    mg.Add(twoSigma)
    mg.Add(oneSigma)
    mg.Add(exp)

  if not options.expected: mg.Add(graph)

# Draw dummy hist
dummyHist = r.TH1D("dummy","",1,110,150)
dummyHist.GetXaxis().SetTitle('m_{H} (GeV)')
if options.pval: dummyHist.GetYaxis().SetTitle('Local p-value')
elif options.maxlh: dummyHist.GetYaxis().SetTitle('Best fit #sigma/#sigma_{SM}')
elif options.limit: dummyHist.GetYaxis().SetTitle("\sigma(H#rightarrow #gamma #gamma)_{95%%CL} / \sigma(H#rightarrow #gamma #gamma)_{SM}")
dummyHist.SetTitleSize(.05,"X")
dummyHist.SetTitleOffset(0.75,"X")
dummyHist.SetTitleSize(.05,"Y")
dummyHist.SetTitleOffset(0.75,"Y")
mg.Draw("A")
dummyHist.SetMinimum(mg.GetYaxis().GetXmin())
dummyHist.SetMaximum(mg.GetYaxis().GetXmax())
if options.limit: dummyHist.GetYaxis().SetRangeUser(0.,4.)
dummyHist.SetLineColor(0)
dummyHist.SetStats(0)
dummyHist.Draw("AXIS")


# Draw Graph
if options.pval:
  mg.Draw("L")
elif options.maxlh:
  mg.Draw("3")
  mg.Draw("LPX")
elif options.limit:
  mg.Draw("3")
  mg.Draw("LPX")
  
dummyHist.Draw("AXIGSAME")

# Draw sigma lines
if options.pval:
  for i, line in enumerate(lines):
    line.Draw()
    labels[i].Draw()

# Draw line at 1 for maxlh and limit
if options.maxlh or options.limit:
  l = r.TLine(110,1.,150,1.)
  l.SetLineColor(r.kRed)
  l.SetLineWidth(2)
  l.Draw()

# Draw text 
lat = r.TLatex()
lat.SetNDC()
lat.DrawLatex(0.12,0.92,"CMS Preliminary")
lat.SetTextSize(0.03)
lat.DrawLatex(0.75,0.94,options.text)

# Draw legend
leg.Draw()

canv.RedrawAxis()

if options.pval:
  canv.Print('pval.pdf')
  canv.Print('pval.png')
elif options.maxlh:
  canv.Print('maxlh.pdf')
  canv.Print('maxlh.png')
elif options.limit:
  canv.Print('limit.pdf')
  canv.Print('limit.png')
