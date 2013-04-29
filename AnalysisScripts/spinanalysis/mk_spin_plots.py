#!/usr/bin/env python
import array

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--datfile",dest="datfile",help="Read from datfile")
parser.add_option("-i","--inputfile",dest="inputfiles",default=[],action="append")
parser.add_option("-M","--Method",dest="methods",default=[],action="append")
parser.add_option("-q","--qqbarPoints",dest="qqbarPoints",default="0.0,0.25,0.5,0.75,1.0")
parser.add_option("-c","--cosThetaBoundaries",dest="cTbounds",default="0.0,0.2,0.375,0.55,0.75,1.")
parser.add_option("-u","--unblind",action="store_true",default=False)
parser.add_option("-b","--isBatch",action="store_true",default=False)
(options,args)=parser.parse_args()

import ROOT as r
if options.isBatch: r.gROOT.SetBatch()

outf = r.TFile('HggSpinStats.root','RECREATE')
canv = r.TCanvas()
pt = r.TPaveText(0.1,0.91,0.45,0.99,"NDC");
pt.SetTextAlign(12);
pt.SetTextSize(0.04);
pt.SetFillColor(0);
pt.AddText("CMS Expected");
pt.SetBorderSize(0);
pt2 = r.TPaveText(0.55,0.91,0.9,0.99,"NDC");
pt2.SetTextAlign(32);
pt2.SetTextSize(0.04);
pt2.SetFillColor(0);
#pt2.AddText(" #sqrt{s} = 7 TeV, L = 5.051 fb^{-1}; #sqrt{s} = 8 TeV, L = 30.0 fb^{-1}");
pt2.AddText(" #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
pt2.SetBorderSize(0);

def doChannelCompatiblity():
  
  print 'Plotting Channel Compatibility ....'

  cosThetaBoundaries = array.array('f',[])
  for b in options.cTbounds.split(','):
    cosThetaBoundaries.append(float(b))

  dummyHist = r.TH1F('dummy','dummy',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  dummyHist.GetXaxis().SetTitle('|cos#theta*|')
  dummyHist.GetYaxis().SetTitle('#sigma/#sigma_{SM}')

  graphSM = r.TGraph()
  graphGrav = r.TGraph()
  graphSME = r.TGraphAsymmErrors()
  graphGravE = r.TGraphAsymmErrors()
  graphData = r.TGraphAsymmErrors()

  for f in options.inputfiles:
    if 'Grav' in f:
      gravFile = r.TFile(f)
    elif 'SM' in f:
      smFile = r.TFile(f)
    elif 'Data' in f:
      dataFile = r.TFile(f)

  sm_fit_nominal = smFile.Get('fit_nominal')
  sm_fit_alternate = smFile.Get('fit_alternate')
  grav_fit_nominal = gravFile.Get('fit_nominal')
  grav_fit_alternate = gravFile.Get('fit_alternate')
  data_fit_nominal = dataFile.Get('fit_nominal')
  data_fit_alternate = dataFile.Get('fit_alternate')

  print 'Spin0:'
  p=0
  for i in range(sm_fit_alternate.floatParsFinal().getSize()):
    chFit = sm_fit_alternate.floatParsFinal().at(i)
    if 'ChannelCompatibilityCheck' in chFit.GetName():
      graphSM.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphSME.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphSME.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
      print '\t Cat%d  %1.3f'%(p,chFit.getVal())
      p+=1
  print 'Spin2:'
  p=0
  for i in range(grav_fit_alternate.floatParsFinal().getSize()):
    chFit = grav_fit_alternate.floatParsFinal().at(i)
    if 'ChannelCompatibilityCheck' in chFit.GetName():
      graphGrav.SetPoint(p,dummyHist.GetBinCenter(p+1),1./chFit.getVal())
      graphGravE.SetPoint(p,dummyHist.GetBinCenter(p+1),1./chFit.getVal())
      graphGravE.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
      print '\t Cat%d  %1.3f'%(p,chFit.getVal())
      p+=1
  if options.unblind:
    print 'Data:'
    p=0
    for i in range(data_fit_alternate.floatParsFinal().getSize()):
      chFit = data_fit_alternate.floatParsFinal().at(i)
      if 'ChannelCompatibilityCheck' in chFit.GetName():
        graphData.SetPoint(p,dummyHist.GetBinCenter(p+1),1./chFit.getVal())
        graphData.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
        print '\t Cat%d  %1.3f'%(p,chFit.getVal())
        p+=1

  graphData.SetMarkerStyle(r.kFullCircle)
  graphData.SetMarkerColor(r.kBlack)
  graphData.SetLineColor(r.kBlack)
  graphData.SetLineWidth(2)

  graphSM.SetMarkerStyle(r.kFullSquare)
  graphSM.SetMarkerColor(r.kRed)
  graphSM.SetLineColor(r.kRed)
  graphSM.SetLineWidth(2)
  graphSM.SetFillColor(r.kRed)
  graphSM.SetFillStyle(3004)
  graphSME.SetLineColor(r.kRed)
  graphSME.SetLineWidth(2)
  graphSME.SetFillColor(r.kRed)
  graphSME.SetFillStyle(3004)
  graphGrav.SetMarkerColor(r.kBlue)
  graphGrav.SetMarkerStyle(r.kFullTriangleUp)
  graphGrav.SetLineColor(r.kBlue)
  graphGrav.SetLineWidth(2)

  leg = r.TLegend(0.11,0.75,0.4,0.89)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  leg.AddEntry(graphSM,"X#rightarrow#gamma#gamma 0^{+}","LPF");
  leg.AddEntry(graphGrav,"X#rightarrow#gamma#gamma 2^{+}_{m}","LP");
  if options.unblind: leg.AddEntry(graphData,"Observed","EP")
  th2 = r.TH2F('h','',1,0,1.,1,ymin,ymax)
  th2.SetStats(0)
  th2.SetLineColor(0)
  th2.GetXaxis().SetTitle("|cos(#theta*)|")
  th2.GetXaxis().SetTitleOffset(1.2)
  th2.GetYaxis().SetTitle("#sigma(X#rightarrow#gamma#gamma 2^{+}_{m})/#sigma(X#rightarrow#gamma#gamma 0^{+})")
  th2.Draw("AXIS")
  graphSME.Draw("E3same")
  graphSM.Draw("LPsame")
  graphGrav.Draw("LPsame")
  if options.unblind: graphData.Draw("EPsame")
  leg.Draw("same")
  f = r.TF1('f','0.',0.,1.)
  f.SetLineColor(r.kBlack)
  f.SetLineWidth(2)
  f.SetLineStyle(r.kDashed)
  f.Draw("same")
  pt.Draw('same')
  pt2.Draw('same')
  th2.Draw("AXISGsame")
  canv.Update()
  if not options.isBatch: raw_input("Looks ok?")
  canv.Update()
  canv.Print('modIndep.pdf')
  canv.Print('modIndep.png')
  
  graphSM.SetName('ChannelCompSM')
  graphSME.SetName('ChannelCompSMerr')
  graphGrav.SetName('ChannelCompGrav')
  graphData.SetName('ChannelCompData')
  outf.cd()
  graphSM.Write()
  graphSME.Write()
  graphGrav.Write()
  if options.unblind: graphData.Write()
  canv.SetName('ChannelComp')
  canv.Write()

def doSeparation():
  
  print 'Plotting Separation ...'

  testStatSM = r.TH1F('testStatSMsep','testStatSM',600,-15,15)
  testStatGRAV = r.TH1F('testStatGRAVsep','testStatGRAV',600,-15,15)
  testStatData = r.TH1F('testStatDatasep','testStatData',600,-15,15)
  testStatSM.SetStats(0)
  testStatGRAV.SetStats(0)
  testStatData.SetStats(0)

  for f in options.inputfiles:
    tf = r.TFile(f)
    tree = tf.Get('q')
    for e in range(tree.GetEntries()):
      tree.GetEntry(e)
      if tree.type>0: testStatSM.Fill(2*tree.q)
      if tree.type<0: testStatGRAV.Fill(2*tree.q)
      if tree.type==0: testStatData.Fill(2*tree.q) 
  SMprobHist = testStatSM.Integral(0,testStatSM.FindBin(testStatGRAV.GetMean()))/testStatSM.Integral();
  GRAVprobHist = testStatGRAV.Integral(testStatGRAV.FindBin(testStatSM.GetMean()),testStatGRAV.GetNbinsX())/testStatGRAV.Integral();
  SMsigmaHist = r.RooStats.PValueToSignificance(SMprobHist) 
  GRAVsigmaHist = r.RooStats.PValueToSignificance(GRAVprobHist)
  
  print '\tProb( q > median(2) | 0 ) = %1.4f = %1.2f'%(SMprobHist,SMsigmaHist)
  print '\tProb( q < median(0) | 2 ) = %1.4f = %1.2f'%(GRAVprobHist,GRAVsigmaHist)
  
  testStatSM.Rebin(10)
  testStatGRAV.Rebin(10)
 
  data = r.TArrow(testStatData.GetMean(),testStatSM.GetEntries()/20.,testStatData.GetMean(),0.,0.02,"|>")
  data.SetLineWidth(4)
  testStatGRAV.SetLineColor(r.kBlue)
  testStatGRAV.SetFillColor(r.kBlue-7)
  testStatGRAV.SetFillStyle(3001)
  testStatSM.SetLineColor(r.kRed)
  testStatSM.SetFillColor(r.kRed-7)
  testStatSM.SetFillStyle(3004)
  testStatGRAV.SetTitle("")
  testStatSM.SetTitle("")
  testStatSM.GetXaxis().SetTitle("-2 #times ln(L_{0^{+}}/L_{2^{+}_{m}})")
  testStatGRAV.GetXaxis().SetTitle("-2 #times ln(L_{0^{+}}/L_{2^{+}_{m}})")
  testStatSM.GetYaxis().SetTitle("Number of toys")
  testStatGRAV.GetYaxis().SetTitle("Number of toys")
  testStatGRAV.GetYaxis().SetTitleOffset(1.2)
  testStatSM.GetYaxis().SetTitleOffset(1.2)
  testStatGRAV.GetYaxis().SetRangeUser(0.,1.1*max(testStatSM.GetMaximum(),testStatGRAV.GetMaximum()))
  testStatSM.GetYaxis().SetRangeUser(0.,1.1*max(testStatSM.GetMaximum(),testStatGRAV.GetMaximum()))
  testStatGRAV.GetXaxis().SetTitleOffset(1.1)
  testStatSM.GetXaxis().SetTitleOffset(1.1)
  leg = r.TLegend(0.11,0.75,0.4,0.89);
  leg.SetLineColor(0);
  leg.SetFillColor(0);
  leg.AddEntry(testStatSM,"X#rightarrow#gamma#gamma 0^{+}","f")
  leg.AddEntry(testStatGRAV,"X#rightarrow#gamma#gamma 2^{+}_{m}","f")
  if options.unblind: leg.AddEntry(data,"Observed","l")

  testStatGRAV.Draw("HIST")
  testStatSM.Draw("HISTSAME")
  if options.unblind: data.Draw()
  leg.Draw("SAME")
  pt.Draw("same")
  pt2.Draw("same")
  pt3 = r.TPaveText(0.7,0.85,0.89,0.89,"NDC")
  pt3.AddText("N = %d toys"%testStatSM.GetEntries())
  pt3.SetBorderSize(0);
  pt3.SetFillColor(0);
  pt3.Draw("same");
  pt4 = r.TPaveText(0.6,0.7,0.89,0.85,"NDC");
  pt4.AddText("p (q>med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma"%GRAVsigmaHist)
  pt4.AddText("p (q<med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma"%SMsigmaHist)
  pt4.SetBorderSize(0);
  pt4.SetFillColor(0);
  pt4.Draw("same");
  canv.Update()
  if not options.isBatch: raw_input("Looks ok?")
  canv.Update()
  canv.Print('separation.pdf')
  canv.Print('separation.png')

  outf.cd()
  testStatSM.Write()
  testStatGRAV.Write()
  testStatData.SetName('testStatDatasep')
  if options.unblind: testStatData.Write()
  canv.SetName('Separation')
  canv.Write()

def getInterval(h,clevel):
  i=0
  while h.Integral(0,i)/h.Integral()<((1-clevel)/2.):
    low = h.GetBinCenter(i)
    i+=1
  i= h.GetNbinsX()
  while h.Integral(i,h.GetNbinsX())/h.Integral()<((1-clevel)/2.):
    high = h.GetBinCenter(i)
    i-=1
  return low,high

def doqqbar():
  
  print 'Plotting fqqbar ...'
  
  qqbarPoints = array.array('f',[])
  for b in options.qqbarPoints.split(','):
    qqbarPoints.append(float(b))

  grSM = r.TGraph()
  grSM68 = r.TGraphAsymmErrors()
  grSM95 = r.TGraphAsymmErrors()
  grGRAV = r.TGraph()
  grData = r.TGraph()

  for p, filename in enumerate(options.inputfiles):
    fqq = qqbarPoints[p]
    testStatSM = r.TH1F('testStatSM_fqq%1.2f'%fqq,'testStatSM',600,-15,15)
    testStatGRAV = r.TH1F('testStatGRAV_fqq%1.2f'%fqq,'testStatGRAV',600,-15,15)
    testStatData = r.TH1F('testStatData_fqq%1.2f'%fqq,'testStatGRAV',600,-15,15)
    testStatSM.SetStats(0)
    testStatGRAV.SetStats(0)
    testStatData.SetStats(0)
    tf = r.TFile.Open(filename)
    tree = tf.Get('q')
    for e in range(tree.GetEntries()):
      tree.GetEntry(e)
      if tree.type>0: testStatSM.Fill(2*tree.q)
      if tree.type<0: testStatGRAV.Fill(2*tree.q)
      if tree.type==0: testStatData.Fill(2*tree.q)
    sm68 = getInterval(testStatSM,0.68)
    sm95 = getInterval(testStatSM,0.95)
    grSM.SetPoint(p,fqq*100.,testStatSM.GetMean())
    grSM68.SetPoint(p,fqq*100.,testStatSM.GetMean())
    grSM95.SetPoint(p,fqq*100.,testStatSM.GetMean())
    grSM68.SetPointError(p,0.,0.,abs(testStatSM.GetMean()-sm68[0]),abs(sm68[1]-testStatSM.GetMean()))
    grSM95.SetPointError(p,0.,0.,abs(testStatSM.GetMean()-sm95[0]),abs(sm95[1]-testStatSM.GetMean()))
    grGRAV.SetPoint(p,fqq*100.,testStatGRAV.GetMean())
    grData.SetPoint(p,fqq*100.,testStatData.GetMean())
    outf.cd()
    testStatSM.Write()
    testStatGRAV.Write()
    if options.unblind: testStatData.Write()
    print '\t %3.1f  %1.3f  %1.3f'%(fqq*100,testStatSM.GetMean(),testStatGRAV.GetMean())

  grData.SetMarkerStyle(r.kFullCircle)
  grData.SetLineWidth(2)
  
  grSM.SetMarkerStyle(r.kFullSquare)
  grSM.SetMarkerColor(r.kRed)
  #grSM.SetLineStyle(r.kDashed)
  grSM.SetLineWidth(2)
  grSM.SetLineColor(r.kRed)

  grSM68.SetLineColor(r.kYellow)
  grSM68.SetFillColor(r.kYellow)
  grSM95.SetLineColor(r.kGreen)
  grSM95.SetFillColor(r.kGreen)

  grGRAV.SetMarkerStyle(r.kFullTriangleUp)
  grGRAV.SetMarkerColor(r.kBlue)
  #grGRAV.SetLineStyle(r.kDashed)
  grGRAV.SetLineWidth(2)
  grGRAV.SetLineColor(r.kBlue)
 
  dummyHist = r.TH1F("d",";f_{q#bar{q}} (%);-2#times ln (L_{0^{+}}/L_{2^{+}_{m}}) ",100,0,100)
  dummyHist.SetMinimum(ymin)
  dummyHist.SetMaximum(ymax)
  dummyHist.SetStats(0)
  dummyHist.Draw("AXIS")

  leg = r.TLegend(0.6,0.65,0.89,0.89);
  leg.SetLineColor(0);
  leg.SetFillColor(0);
  leg.AddEntry(grSM,"X#rightarrow#gamma#gamma 0^{+}","lp");
  leg.AddEntry(grSM68,"#pm 1#sigma expected","f");
  leg.AddEntry(grSM95,"#pm 2#sigma expected","f");
  leg.AddEntry(grGRAV,"X#rightarrow#gamma#gamma 2^{+}_{m}","lp");
  if options.unblind: leg.AddEntry(grData,"Observed","lp")

  grSM95.Draw("E3same")
  grSM68.Draw("E3same")
  grGRAV.Draw("LPsame")
  grSM.Draw("LPsame")
  if options.unblind: grData.Draw("LPsame")
  leg.Draw("SAME")
  f = r.TF1('f','0.',0.,100.)
  f.SetLineColor(r.kBlack)
  f.SetLineWidth(2)
  f.SetLineStyle(r.kDashed)
  f.Draw("same")
  pt.Draw("same")
  pt2.Draw("same")
  dummyHist.Draw("AXISGsame")
  canv.Update()
  if not options.isBatch: raw_input("Looks ok?")
  canv.Update()
  canv.Print('fqqbar.pdf')
  canv.Print('fqqbar.png')

  grSM.SetName('fqqSM')
  grGRAV.SetName('fqqGRAV')
  grSM95.SetName('fqqSM95')
  grSM68.SetName('fqqSM68')
  grData.SetName('fqqData')
  outf.cd()
  grSM.Write()
  grGRAV.Write()
  grSM95.Write()
  grSM68.Write()
  if options.unblind: grData.Write()
  canv.SetName('fqq')
  canv.Write()
    
if options.datfile:
  df = open(options.datfile)
  for line in df.readlines():
    opts = line.split(':')
    print opts
    meth = opts[0]
    options.inputfiles = opts[1].split(',')
    options.qqbarPoints = opts[2]
    options.cosThetaBoundaries=opts[3]
    if len(opts)>5:
      ymin=float(opts[4])
      ymax=float(opts[5])
    if meth=='ChannelCompatibility': doChannelCompatiblity()
    if meth=='Separation': doSeparation()
    if meth=='qqbar': doqqbar()

else:
  if len(options.methods)==0: options.methods=['ChannelCompatibility','Separation','qqbar']
  if 'ChannelCompatibility' in options.methods: doChannelCompatiblity()
  if 'Separation' in options.methods: doSeparation()
  if 'qqbar' in options.methods: doqqbar()

outf.Close()
