#!/usr/bin/env python
import array
import os
import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--datfile",dest="datfile",help="Read from datfile")
parser.add_option("-i","--inputfile",dest="inputfiles",default=[],action="append")
parser.add_option("-M","--Method",dest="methods",default=[],action="append")
parser.add_option("-q","--qqbarPoints",dest="qqbarPoints",default="0.0,0.25,0.5,0.75,1.0")
parser.add_option("-c","--cosThetaBoundaries",dest="cTbounds",default="0.0,0.2,0.375,0.55,0.75,1.")
parser.add_option("-D","--outDir",default="./")
parser.add_option("-u","--unblind",action="store_true",default=False)
parser.add_option("-b","--isBatch",action="store_true",default=False)
(options,args)=parser.parse_args()

import ROOT as r
if options.isBatch: r.gROOT.SetBatch()
r.TH1.SetDefaultSumw2()

os.system('mkdir -p %s'%options.outDir)

if not os.path.isfile('extractSignificanceStats.C'):
  sys.exit('Can\'t fine file - extractSignificanceStats.C')
r.gROOT.ProcessLine('.L extractSignificanceStats.C+')
from ROOT import extractSignificanceStats
from ROOT import setTDRStyle

outf = r.TFile('HggSpinStats.root','RECREATE')
canv = r.TCanvas()
pt = r.TPaveText(0.1,0.91,0.45,0.99,"NDC");
pt.SetTextAlign(12);
pt.SetTextSize(0.04);
pt.SetFillColor(0);
pt.AddText("CMS Preliminary");
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
  dummyHist.GetYaxis().SetTitle('#mu')

  graphSM = r.TGraph()
  graphGravGG = r.TGraph()
  graphGravQQ = r.TGraph()
  graphGrav50GG50QQ = r.TGraph()
  graphData = r.TGraph()
  graphSME = r.TGraphAsymmErrors()
  graphGravGGE = r.TGraphAsymmErrors()
  graphGravQQE = r.TGraphAsymmErrors()
  graphGrav50GG50QQE = r.TGraphAsymmErrors()
  graphDataE = r.TGraphAsymmErrors()

  histSM = r.TH1F('histSM','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histGravGG = r.TH1F('histGravGG','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histGravQQ = r.TH1F('histGravQQ','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histGrav50GG50QQ = r.TH1F('histGrav50GG50QQ','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histData = r.TH1F('histData','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histSME = r.TH1F('histSME','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histGravGGE = r.TH1F('histGravGGE','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histGravQQE = r.TH1F('histGravQQE','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histGrav50GG50QQE = r.TH1F('histGrav50GG50QQE','',len(cosThetaBoundaries)-1,cosThetaBoundaries)
  histDataE = r.TH1F('histDataE','',len(cosThetaBoundaries)-1,cosThetaBoundaries)


  for f in options.inputfiles:
    if 'GravGG.' in f or 'GRAVGG.' in f:
      gravGGFile = r.TFile(f)
    elif 'GravQQ.' in f or 'GRAVQQ.' in f:
      gravQQFile = r.TFile(f)
    elif 'Grav50GG50QQ.' in f or 'GRAV50GG50QQ.' in f:
      grav50GG50QQFile = r.TFile(f)
    elif 'SM.' in f:
      smFile = r.TFile(f)
    elif 'Data.' in f:
      dataFile = r.TFile(f)

  sm_fit_nominal = smFile.Get('fit_nominal')
  sm_fit_alternate = smFile.Get('fit_alternate')
  grav_gg_fit_nominal = gravGGFile.Get('fit_nominal')
  grav_gg_fit_alternate = gravGGFile.Get('fit_alternate')
  grav_qq_fit_nominal = gravQQFile.Get('fit_nominal')
  grav_qq_fit_alternate = gravQQFile.Get('fit_alternate')
  grav_ggqq_fit_nominal = grav50GG50QQFile.Get('fit_nominal')
  grav_ggqq_fit_alternate = grav50GG50QQFile.Get('fit_alternate')
  data_fit_nominal = dataFile.Get('fit_nominal')
  data_fit_alternate = dataFile.Get('fit_alternate')

  pointsData=[]
  pointsSM=[]
  pointsGRAVGG=[]
  pointsGRAVQQ=[]
  pointsGRAV50GG50QQ=[]

  print 'Spin0:'
  p=0
  for i in range(sm_fit_alternate.floatParsFinal().getSize()):
    chFit = sm_fit_alternate.floatParsFinal().at(i)
    if 'ChannelCompatibilityCheck' in chFit.GetName():
      graphSM.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphSME.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphSME.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
      histSM.SetBinContent(p+1,chFit.getVal())
      histSME.SetBinContent(p+1,chFit.getVal())
      histSME.SetBinError(p+1,(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.)
      pointsSM.append([chFit.getVal(),(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.])
      print '\t Cat%d  %1.3f  %1.3f'%(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      p+=1
  print 'Spin2 (gg):'
  p=0
  for i in range(grav_gg_fit_alternate.floatParsFinal().getSize()):
    chFit = grav_gg_fit_alternate.floatParsFinal().at(i)
    if 'ChannelCompatibilityCheck' in chFit.GetName():
      graphGravGG.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphGravGGE.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphGravGGE.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
      histGravGG.SetBinContent(p+1,chFit.getVal())
      histGravGGE.SetBinContent(p+1,chFit.getVal())
      histGravGGE.SetBinError(p+1,(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.)
      pointsGRAVGG.append([chFit.getVal(),(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.])
      print '\t Cat%d  %1.3f  %1.3f'%(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      p+=1
  print 'Spin2 (qq):'
  p=0
  for i in range(grav_qq_fit_alternate.floatParsFinal().getSize()):
    chFit = grav_qq_fit_alternate.floatParsFinal().at(i)
    if 'ChannelCompatibilityCheck' in chFit.GetName():
      graphGravQQ.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphGravQQE.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphGravQQE.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
      histGravQQ.SetBinContent(p+1,chFit.getVal())
      histGravQQE.SetBinContent(p+1,chFit.getVal())
      histGravQQE.SetBinError(p+1,(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.)
      pointsGRAVQQ.append([chFit.getVal(),(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.])
      print '\t Cat%d  %1.3f  %1.3f'%(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      p+=1
  print 'Spin2 (0.5gg+0.5qq):'
  p=0
  for i in range(grav_ggqq_fit_alternate.floatParsFinal().getSize()):
    chFit = grav_ggqq_fit_alternate.floatParsFinal().at(i)
    if 'ChannelCompatibilityCheck' in chFit.GetName():
      graphGrav50GG50QQ.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphGrav50GG50QQE.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      graphGrav50GG50QQE.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
      histGrav50GG50QQ.SetBinContent(p+1,chFit.getVal())
      histGrav50GG50QQE.SetBinContent(p+1,chFit.getVal())
      histGrav50GG50QQE.SetBinError(p+1,(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.)
      pointsGRAV50GG50QQ.append([chFit.getVal(),(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.])
      print '\t Cat%d  %1.3f  %1.3f'%(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
      p+=1
  if options.unblind:
    print 'Data:'
    p=0
    for i in range(data_fit_alternate.floatParsFinal().getSize()):
      chFit = data_fit_alternate.floatParsFinal().at(i)
      if 'ChannelCompatibilityCheck' in chFit.GetName():
        graphData.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
        graphDataE.SetPoint(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
        graphDataE.SetPointError(p,0.,0.,abs(chFit.getErrorLo()),abs(chFit.getErrorHi()))
        histData.SetBinContent(p+1,chFit.getVal())
        histDataE.SetBinContent(p+1,chFit.getVal())
        histDataE.SetBinError(p+1,(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.)
        pointsData.append([chFit.getVal(),(abs(chFit.getErrorLo())+abs(chFit.getErrorHi()))/2.])
        print '\t Cat%d  %1.3f  %1.3f'%(p,dummyHist.GetBinCenter(p+1),chFit.getVal())
        p+=1

  chi2SM=0.
  chi2GRAVGG=0.
  chi2GRAVQQ=0.
  chi2GRAV50GG50QQ=0.
  for p,val in enumerate(pointsData):
    chi2SM += ((val[0]-pointsSM[p][0])**2)/val[1]**2
    chi2GRAVGG += ((val[0]-pointsGRAVGG[p][0])**2)/val[1]**2
    chi2GRAVQQ += ((val[0]-pointsGRAVQQ[p][0])**2)/val[1]**2
    chi2GRAV50GG50QQ += ((val[0]-pointsGRAV50GG50QQ[p][0])**2)/val[1]**2

  if options.unblind:
    print 'Chi2 Comp DATA->SM:            ', chi2SM, ' pval = ', r.TMath.Prob(chi2SM,len(pointsData))
    print 'Chi2 Comp DATA->GRAVGG:        ', chi2GRAVGG, ' pval = ', r.TMath.Prob(chi2GRAVGG,len(pointsData))
    print 'Chi2 Comp DATA->GRAVQQ:        ', chi2GRAVQQ, ' pval = ', r.TMath.Prob(chi2GRAVQQ,len(pointsData))
    print 'Chi2 Comp DATA->GRAV50GG50QQ:  ', chi2GRAV50GG50QQ, ' pval = ', r.TMath.Prob(chi2GRAV50GG50QQ,len(pointsData))

  
  graphData.SetMarkerStyle(r.kFullCircle)
  graphData.SetMarkerColor(r.kBlack)
  graphData.SetLineColor(r.kBlack)
  graphData.SetLineWidth(2)

  graphDataE.SetMarkerStyle(r.kFullCircle)
  graphDataE.SetMarkerColor(r.kBlack)
  graphDataE.SetLineColor(r.kBlack)
  graphDataE.SetLineWidth(2)

  histData.SetMarkerStyle(r.kFullCircle)
  histData.SetMarkerColor(r.kBlack)
  histData.SetLineColor(r.kBlack)
  histData.SetLineWidth(2)

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
  graphGravGG.SetMarkerColor(r.kBlue)
  graphGravGG.SetMarkerStyle(r.kFullTriangleUp)
  graphGravGG.SetLineColor(r.kBlue)
  graphGravGG.SetLineWidth(2)
  graphGravQQ.SetMarkerColor(r.kGreen+1)
  graphGravQQ.SetMarkerStyle(r.kFullTriangleDown)
  graphGravQQ.SetLineColor(r.kGreen+1)
  graphGravQQ.SetLineWidth(2)
  graphGrav50GG50QQ.SetMarkerColor(r.kMagenta+1)
  graphGrav50GG50QQ.SetMarkerStyle(r.kFullDiamond)
  graphGrav50GG50QQ.SetMarkerSize(1.2)
  graphGrav50GG50QQ.SetLineColor(r.kMagenta+1)
  graphGrav50GG50QQ.SetLineWidth(2)

  histSM.SetMarkerStyle(r.kFullSquare)
  histSM.SetMarkerColor(r.kRed)
  histSM.SetLineColor(r.kRed)
  histSM.SetLineWidth(2)
  histSME.SetLineColor(r.kRed)
  histSME.SetLineWidth(2)
  histSME.SetFillColor(r.kRed)
  histSME.SetFillStyle(3004)
  histGravGG.SetMarkerColor(r.kBlue)
  histGravGG.SetMarkerStyle(r.kFullTriangleUp)
  histGravGG.SetLineColor(r.kBlue)
  histGravGG.SetLineWidth(2)
  histGravQQ.SetMarkerColor(r.kGreen+1)
  histGravQQ.SetMarkerStyle(r.kFullTriangleDown)
  histGravQQ.SetLineColor(r.kGreen+1)
  histGravQQ.SetLineWidth(2)
  histGrav50GG50QQ.SetMarkerColor(r.kMagenta+1)
  histGrav50GG50QQ.SetMarkerStyle(r.kFullDiamond)
  histGrav50GG50QQ.SetMarkerSize(1.2)
  histGrav50GG50QQ.SetLineColor(r.kMagenta+1)
  histGrav50GG50QQ.SetLineWidth(2)

  leg = r.TLegend(0.11,0.65,0.4,0.89)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  leg.AddEntry(graphSM,"X#rightarrow#gamma#gamma 0^{+}","LPF");
  leg.AddEntry(graphGravGG,"X#rightarrow#gamma#gamma 2^{+}_{m}(100%gg)","LP");
  leg.AddEntry(graphGravQQ,"X#rightarrow#gamma#gamma 2^{+}_{m}(100%qq)","LP");
  leg.AddEntry(graphGrav50GG50QQ,"X#rightarrow#gamma#gamma 2^{+}_{m}(50%gg,50%qq)","LP");
  if options.unblind: leg.AddEntry(graphData,"Observed","EP")
  th2 = r.TH2F('h','',1,0,1.,1,ymin,ymax)
  th2.SetStats(0)
  th2.SetLineColor(0)
  th2.GetXaxis().SetTitle("|cos(#theta*)|")
  #th2.GetYaxis().SetTitle("#sigma(X#rightarrow#gamma#gamma 2^{+}_{m})/#sigma(X#rightarrow#gamma#gamma 0^{+})")
  th2.GetYaxis().SetTitle("#sigma/#sigma_{SM}")
  th2.GetXaxis().SetLabelSize(0.05);
  th2.GetXaxis().SetTitleSize(0.05);
  th2.GetYaxis().SetLabelSize(0.05);
  th2.GetYaxis().SetTitleSize(0.05);
  th2.GetXaxis().SetTitleOffset(1.0);
  th2.GetYaxis().SetTitleOffset(1.0);
  th2.Draw("AXIS")
  #graphSME.Draw("E3same")
  graphSM.Draw("LPsame")
  graphGravGG.Draw("LPsame")
  graphGravQQ.Draw("LPsame")
  graphGrav50GG50QQ.Draw("LPsame")
  if options.unblind: graphDataE.Draw("PEsame")
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
  canv.Print(options.outDir+'/modIndep.pdf')
  canv.Print(options.outDir+'/modIndep.png')
  
  canv.Clear()
  th2.Draw("AXIS")
  histSME.Draw("E3same")
  histSM.Draw("LPsame")
  #histSM.Draw("Psame")
  histGravGG.Draw("LPsame")
  histGravQQ.Draw("LPsame")
  histGrav50GG50QQ.Draw("LPsame")
  #histGrav.Draw("Psame")
  if options.unblind: 
    histData.Draw("HISTsame")
    histData.Draw("Psame")
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
  canv.Print(options.outDir+'/modIndepHist.pdf')
  canv.Print(options.outDir+'/modIndepHist.png')
  
  graphSM.SetName('ChannelCompSM')
  graphSME.SetName('ChannelCompSMerr')
  graphGravGG.SetName('ChannelCompGravGG')
  graphGravQQ.SetName('ChannelCompGravQQ')
  graphGrav50GG50QQ.SetName('ChannelCompGrav50GG50QQ')
  graphData.SetName('ChannelCompData')
  outf.cd()
  graphSM.Write()
  graphSME.Write()
  graphGravGG.Write()
  graphGravQQ.Write()
  graphGrav50GG50QQ.Write()
  histSM.Write()
  histSME.Write()
  histGravGG.Write()
  histGravQQ.Write()
  histGrav50GG50QQ.Write()
  if options.unblind: 
    graphData.Write()
    histData.Write()
  canv.SetName('ChannelComp')
  canv.Write()

def doSeparation():
  
  print '---------------------------'
  print 'Plotting Separation ...'
  print '---------------------------'

  testStatSM = r.TH1F('testStatSMsep','testStatSM',600,-10,10)
  testStatGRAV = r.TH1F('testStatGRAVsep','testStatGRAV',600,-10,10)
  testStatData = r.TH1F('testStatDatasep','testStatData',600,-10,10)
  testStatSM.SetStats(0)
  testStatGRAV.SetStats(0)
  testStatData.SetStats(0)

  for f in options.inputfiles:
    tf = r.TFile(f)
    tree = tf.Get('q')
    for e in range(tree.GetEntries()):
      tree.GetEntry(e)
      if tree.type<0: testStatSM.Fill(-2.0*tree.q)
      if tree.type>0: testStatGRAV.Fill(-2.0*tree.q)
      if tree.type==0: testStatData.Fill(-2.0*tree.q) 
  SMprobHist = testStatSM.Integral(0,testStatSM.FindBin(testStatGRAV.GetMean()))/testStatSM.Integral();
  GRAVprobHist = testStatGRAV.Integral(testStatGRAV.FindBin(testStatSM.GetMean()),testStatGRAV.GetNbinsX())/testStatGRAV.Integral();
  SMsigmaHist = r.RooStats.PValueToSignificance(SMprobHist) 
  GRAVsigmaHist = r.RooStats.PValueToSignificance(GRAVprobHist)
  
  print '\tProb( q < median(2) | 0 ) = %1.4f = %1.2f'%(SMprobHist,SMsigmaHist)
  print '\tProb( q > median(0) | 2 ) = %1.4f = %1.2f'%(GRAVprobHist,GRAVsigmaHist)
  print '\tExpected CLs to exclude spin-2 = %1.4f'%(GRAVprobHist/0.5)
  
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
  #testStatSM.GetXaxis().SetTitle("-2 #times ln(L_{0^{+}}/L_{2^{+}_{m}})")
  #testStatGRAV.GetXaxis().SetTitle("-2 #times ln(L_{0^{+}}/L_{2^{+}_{m}})")
  testStatSM.GetXaxis().SetTitle("-2 #times ln(L_{2^{+}_{m}(gg)}/L_{0^{+}})")
  testStatGRAV.GetXaxis().SetTitle("-2 #times ln(L_{2^{+}_{m}(gg)}/L_{0^{+}})")
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
  
  canv.Clear()
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
  pt4.AddText("p (q<med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma"%GRAVsigmaHist)
  pt4.AddText("p (q>med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma"%SMsigmaHist)
  pt4.SetBorderSize(0);
  pt4.SetFillColor(0);
  pt4.Draw("same");
  canv.Update()
  if not options.isBatch: raw_input("Looks ok?")
  canv.Update()
  canv.Print(options.outDir+'/separation.pdf')
  canv.Print(options.outDir+'/separation.png')

  outf.cd()
  testStatSM.Write()
  testStatGRAV.Write()
  testStatData.SetName('testStatDatasep')
  canv.SetName('Separation')
  canv.Write()
  canv.Clear()
  if options.unblind: testStatData.Write()

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
  
  print '---------------------------'
  print 'Plotting fqqbar ...'
  print '---------------------------'
  
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
      if tree.type<0: testStatSM.Fill(-2.0*tree.q)
      if tree.type>0: testStatGRAV.Fill(-2.0*tree.q)
      if tree.type==0: testStatData.Fill(-2.0*tree.q)
    sm68 = getInterval(testStatSM,(1.-r.TMath.Prob(1,1)))
    sm95 = getInterval(testStatSM,(1.-r.TMath.Prob(4,1)))
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

  grSM68.SetLineColor(r.kGreen)
  grSM68.SetFillColor(r.kGreen)
  grSM95.SetLineColor(r.kYellow)
  grSM95.SetFillColor(r.kYellow)

  grGRAV.SetMarkerStyle(r.kFullTriangleUp)
  grGRAV.SetMarkerColor(r.kBlue)
  #grGRAV.SetLineStyle(r.kDashed)
  grGRAV.SetLineWidth(2)
  grGRAV.SetLineColor(r.kBlue)

  canv.Clear()
  dummyHist = r.TH1F("d",";f_{q#bar{q}} (%);-2#times ln (L_{2^{+}_{m}}/L_{0^{+}}) ",100,0,100)
  dummyHist.SetMinimum(ymin)
  dummyHist.SetMaximum(ymax)
  dummyHist.SetStats(0)
  dummyHist.GetXaxis().SetLabelSize(0.05);
  dummyHist.GetXaxis().SetTitleSize(0.05);
  dummyHist.GetYaxis().SetLabelSize(0.05);
  dummyHist.GetYaxis().SetTitleSize(0.05);
  dummyHist.GetXaxis().SetTitleOffset(0.9);
  dummyHist.GetYaxis().SetTitleOffset(0.9);
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
  canv.Print(options.outDir+'/fqqbar.pdf')
  canv.Print(options.outDir+'/fqqbar.png')

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

def makeStandardChannelCompatibilityPlot(file,xlow,xhigh,outextension):
  
  if not os.path.isfile('plotSingleMassSignalStrength.py'):
    sys.exit('Can\'t find fine plotSingleMassSignalStrength.py')
  os.system('./plotSingleMassSignalStrength.py -i %s -l %.1f -u %.1f -o %s -O %s'%(file,xlow,xhigh,outextension,options.outDir))

def doChannelCompatiblityStandard():

  assert(len(options.inputfiles)==2)
  makeStandardChannelCompatibilityPlot(options.inputfiles[0],ymin,ymax,options.inputfiles[1])

def makeCombineStyleSeparationPlot(ifile,leg,ofile,rebin,xlow,xhigh):

  print 'Plotting combine style separartion'
  setTDRStyle()
  extractSignificanceStats(options.unblind,leg,options.outDir+'/'+ofile,ifile,rebin,xlow,xhigh)
  if not options.isBatch: raw_input('Looks ok?')

def doSeparationComb():
  makeCombineStyleSeparationPlot(options.inputfiles[0],"2_{m}^{+} (gg)","sep_2mpgg",640,-10.,10)

def doqqbarComb():
  qqbarPoints = array.array('f',[])
  for b in options.qqbarPoints.split(','):
    qqbarPoints.append(float(b))

  for p, filename in enumerate(options.inputfiles):
    fqq = qqbarPoints[p]
    print '----------------------'
    print 'Plotting separation for fqq = ',fqq
    print '----------------------'
    fname = filename.split('{')[0]
    info = filename.split('{')[1].split('}')[0]
    rebin = int(info.split(';')[0])
    low = float(info.split(';')[1])
    high = float(info.split(';')[2])
    makeCombineStyleSeparationPlot(fname,'2_{m}^{+} (%1.0f%% qq)'%(fqq*100.),'sep_fqq%3.2f'%(fqq),rebin,low,high)
    
if options.datfile:
  df = open(options.datfile)
  for line in df.readlines():
    if line.startswith('#'): continue
    opts = line.split(':')
    print opts
    meth = opts[0]
    options.inputfiles = opts[1].split(',')
    options.qqbarPoints = opts[2]
    options.cTbounds=opts[3]
    if len(opts)>5:
      ymin=float(opts[4])
      ymax=float(opts[5])
    if meth=='ChannelCompatibility': doChannelCompatiblity()
    if meth=='Separation': doSeparation()
    if meth=='qqbar': doqqbar()
    if meth=='SeparationComb': doSeparationComb()
    if meth=='qqbarComb': doqqbarComb()
    if meth=='ChannelCompatibilityStandard': doChannelCompatiblityStandard()

else:
  if len(options.methods)==0: options.methods=['ChannelCompatibility','Separation','qqbar','SeparationComb','qqbarComb']
  if 'ChannelCompatibility' in options.methods: doChannelCompatiblity()
  if 'Separation' in options.methods: doSeparation()
  if 'qqbar' in options.methods: doqqbar()
  if 'SeparationComb' in option.methods: doSeparationComb()
  if 'qqbarComb' in options.methods: doSeparationComb()

outf.Close()
