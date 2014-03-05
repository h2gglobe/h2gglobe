#!/usr/bin/env python
print "Loading Root..."

#import pdb; pdb.set_trace()
import os
import sys
import array
import math

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-i","--infile")
parser.add_option("-t","--title",default="")
parser.add_option("-l","--xlow",type="float",default=-15.)
parser.add_option("-u","--xhigh",type="float",default=15.)
parser.add_option("-r","--rearrange",action="store_true",default=False)
parser.add_option("-D","--drawLabel",action="store_true",default=False)
parser.add_option("-L","--drawLines",action="store_true",default=False)
parser.add_option("-o","--outname",default="")
parser.add_option("-O","--outDir",default="./")
parser.add_option("-S","--sqrtS",default="comb")
(options,args)=parser.parse_args()

import ROOT as r
r.gStyle.SetOptStat(000000)
r.gStyle.SetCanvasBorderMode(0);
r.gStyle.SetCanvasColor(r.kWhite);
r.gROOT.SetBatch(1)
r.gStyle.SetPadLeftMargin(0.28);
r.gStyle.SetPadRightMargin(0.05);

print "Setting Initial Parameters."
pt = r.TPaveText(0.1,0.91,0.45,0.99,"NDC");
pt.SetTextAlign(12);
pt.SetTextSize(0.03);
pt.SetFillColor(0);
pt.AddText("CMS Preliminary");
pt.SetBorderSize(0);
pt2 = r.TPaveText(0.7,0.90,0.9,0.99,"NDC");
pt2.SetTextAlign(12);
pt2.SetTextSize(0.028);
pt2.SetFillColor(0);
if options.sqrtS=="comb":
	pt2.AddText("#splitline{#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}}{#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}}");
elif options.sqrtS=="7":
	pt2.AddText(" #sqrt{s} = 7 TeV, L = 5.1 fb^{-1}");
elif options.sqrtS=="8":
	pt2.AddText(" #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}");
else:
	sys.exit("Invalid sqrtS"+sqrtS)
pt2.SetBorderSize(0);
pt3 = r.TPaveText(0.4,0.92,0.8,0.99,"NDC")
pt3.SetTextSize(0.045)
pt3.SetTextFont(42)
pt3.AddText(options.title)
pt3.SetFillColor(0)
pt3.SetBorderSize(0)

can = r.TCanvas("Plots","Plots",800,800)
leg = r.TLegend(0.75, 0.7, 0.94, 0.89)

leg.SetFillColor(0)
leg.SetBorderSize(1)
mytext = r.TLatex()
mytext.SetTextSize(0.04)
mytext.SetTextAlign(32)
mytext.SetNDC()
rMin=options.xlow
rMax=options.xhigh
labelMap={}
Label1=r.TText()
Label1.SetNDC()
Label1.SetText(0.25,0.125,"0 #leq |cos(#theta_{CS})| < 0.2")
labelMap['SpinCat0'] = "0 #leq |cos(#theta_{CS})| < 0.2"
Label1.SetTextFont(62)
Label1.SetTextSize(.03)
Label1.SetTextAngle(90)
Label2=r.TText()
Label2.SetNDC()
Label2.SetText(0.25,0.275,"0.2 #leq |cos(#theta_{CS})| < 0.375")
labelMap['SpinCat1'] = "0.2 #leq |cos(#theta_{CS})| < 0.375"
Label2.SetTextFont(62)
Label2.SetTextSize(.03)
Label2.SetTextAngle(90)
Label3=r.TText()
Label3.SetNDC()
Label3.SetText(0.25,0.43,"0.375 #leq |cos(#theta_{CS})| < 0.55")
labelMap['SpinCat2'] = "0.375 #leq |cos(#theta_{CS})| < 0.55"
Label3.SetTextFont(62)
Label3.SetTextSize(.03)
Label3.SetTextAngle(90)
Label4=r.TText()
Label4.SetNDC()
Label4.SetText(0.25,0.6,"0.55 #leq |cos(#theta_{CS})| < 0.75")
labelMap['SpinCat3'] = "0.55 #leq |cos(#theta_{CS})| < 0.75"
Label4.SetTextFont(62)
Label4.SetTextSize(.03)
Label4.SetTextAngle(90)
Label5=r.TText()
Label5.SetNDC()
Label5.SetText(0.25,0.75,"0.75 #leq |cos(#theta_{CS})| #leq 1.0")
labelMap['SpinCat4'] = "0.75 #leq |cos(#theta_{CS})| #leq 1.0"
Label5.SetTextFont(62)
Label5.SetTextSize(.03)
Label5.SetTextAngle(90)
#splitline=5
splitline=9
myline = r.TLine(rMin,splitline,rMax,splitline)
myline.SetLineColor(r.kBlack)
myline.SetLineWidth(3)
myline.SetLineStyle(7)
cosline0 = r.TLine(rMin,4,rMax,4)
cosline0.SetLineColor(r.kBlack)
cosline0.SetLineWidth(3)
cosline0.SetLineStyle(7)
cosline1 = r.TLine(rMin,8,rMax,8)
cosline1.SetLineColor(r.kBlack)
cosline1.SetLineWidth(3)
cosline1.SetLineStyle(7)
cosline2 = r.TLine(rMin,12,rMax,12)
cosline2.SetLineColor(r.kBlack)
cosline2.SetLineWidth(3)
cosline2.SetLineStyle(7)
cosline3 = r.TLine(rMin,16,rMax,16)
cosline3.SetLineColor(r.kBlack)
cosline3.SetLineWidth(3)
cosline3.SetLineStyle(7)

dumbline = r.TLine()
dumbline.SetLineColor(r.kBlack)
dumbline.SetLineWidth(1)
#BinLabel = ["Untagged 0","Untagged 1","Untagged 2","Untagged 3","Di-jet tight","Di-jet loose","Muon","Electron","MET","Untagged 0","Untagged 1","Untagged 2","Untagged 3","Di-jet"]
intlumi=[5.1,19.6]
Energy=[7,8]
#intlumi=[12.3]
#Energy=[8]
title=options.title

file = r.TFile.Open(options.infile)
normal=file.Get("fit_nominal")
rFit=normal.floatParsFinal().find("r")
CombinedBestFitObserved=rFit.getVal()
CombinedBestFitErrorUp=rFit.getAsymErrorHi()
CombinedBestFitErrorDown=rFit.getAsymErrorLo()
if math.fabs(CombinedBestFitErrorDown)<0.000001: CombinedBestFitErrorDown=-rFit.getAsymErrorHi()
altername=file.Get("fit_alternate")
Channels=[]
BestFitObserved=[]
BestFitErrorUp=[]
BestFitErrorDown=[]
fitResults=[]
for i in range(altername.floatParsFinal().getSize()):
  #print i,altername.floatParsFinal().at(i).GetName(),altername.floatParsFinal().at(i).getVal(),altername.floatParsFinal().at(i).getAsymErrorHi(),altername.floatParsFinal().at(i).getAsymErrorLo()
  chrFit=altername.floatParsFinal().at(i)
  if chrFit.GetName().find("ChannelCompatibilityCheck")!=-1:
    Channels.append(chrFit.GetName())
    BestFitObserved.append(chrFit.getVal())
    BestFitErrorUp.append(chrFit.getAsymErrorHi())
    if math.fabs(chrFit.getAsymErrorLo())<0.000001: BestFitErrorDown.append(math.fabs(chrFit.getAsymErrorHi()))
    else: BestFitErrorDown.append(math.fabs(chrFit.getAsymErrorLo()))
    name = chrFit.GetName().split('r_')[1]
    if name.find('spinCat')==-1: 
      spinCat=-1
    else:
      spinCat = int(name[name.find('spinCat')+7])
    if name.find('kinCat')==-1:
      kinCat=-1
    else:
      kinCat = int(name[name.find('kinCat')+6])
    fitResults.append([spinCat,kinCat,chrFit.getVal(),math.fabs(chrFit.getAsymErrorLo()),chrFit.getAsymErrorHi()])

fitResults.sort(key=lambda x: x[1])
fitResults.sort(key=lambda x: x[0])
dummyHist = r.TH2F("dummy",";Best Fit #mu;",1,rMin,rMax,len(Channels),0,len(Channels))
ChannelNum=array.array("f",[x + 0.5 for x in range(len(Channels))])
graph = r.TGraphAsymmErrors()
p=0
for values in fitResults:
  print values[0], values[1], values[2], values[3], values[4]
  dummyHist.GetYaxis().SetBinLabel(p+1,'KinCat%d'%(values[1]))
  if values[1]==-1:
    dummyHist.GetYaxis().SetBinLabel(p+1,labelMap['SpinCat%d'%(values[0])])
  graph.SetPoint(p,values[2],p+0.5)
  graph.SetPointError(p,values[3],values[4],0.,0.)
  p+=1
#graph=TGraphAsymmErrors(len(Channels),array.array("f",BestFitObserved),ChannelNum,array.array("f",BestFitErrorDown),array.array("f",BestFitErrorUp),array.array("f",[0]*len(Channels)),array.array("f",[0]*len(Channels)))
graph.SetLineColor(r.kRed)
graph.SetLineWidth(1)
graph.SetMarkerStyle(21)
graph.SetMarkerSize(2)
graph.SetMarkerColor(r.kBlack)
graph.SetFillColor(r.kWhite)
#leg.AddEntry(graph,"BestFit Category")
#for i in range(len(Channels)):
  #print Channels[i],BinLabel[i]
  #dummyHist.GetYaxis().SetBinLabel(i+1,'Cat%d'%i);
dummyHist.SetTitleSize(0.045,"X")
dummyHist.SetTitleOffset(0.95,"X")
dummyHist.SetLabelSize(0.045,"X")
dummyHist.SetLabelFont(62,"Y")
dummyHist.SetLabelSize(0.035,"Y")
dummyHist.SetTitleOffset(1.1,"Y")
dummyHist.SetLineColor(r.kBlack)
dummyHist.SetFillColor(r.kGreen)
#leg.AddEntry(dummyHist,"68% Combined")
limitTree = file.Get('limit')
limitTree.GetEntry(0)
Mass = limitTree.mh
Chi2 = limitTree.limit
Chi2Prob = r.TMath.Prob(Chi2,len(fitResults))
print 'Mass: %5.1f  Chi2: %4.2f  Chi2Prob: %4.2f'%(Mass,Chi2,Chi2Prob)
leg.AddEntry(graph,"Event Class")
leg.AddEntry(dummyHist,"Combined","f")
leg.AddEntry(r.NULL,"m_{H} = %.1f GeV" %Mass,"")
leg.AddEntry(r.NULL,"#mu = %.2f#pm%.2f" %(CombinedBestFitObserved,(CombinedBestFitErrorUp+r.TMath.Abs(CombinedBestFitErrorDown))/2.),"")
leg.AddEntry(r.NULL,"#chi^{2} = %3.2f"%Chi2,"")
leg.AddEntry(r.NULL,"#chi^{2} p-val = %3.2f"%Chi2Prob,"")
dummyHist.Draw()
BestFitBand=r.TBox(CombinedBestFitObserved+CombinedBestFitErrorDown,0,CombinedBestFitObserved+CombinedBestFitErrorUp,len(Channels))
#BestFitBand.SetFillStyle(3013)
BestFitBand.SetFillColor(r.kGreen)
BestFitBand.SetLineStyle(0)
BestFitBand.Draw()
BestFitLine=r.TLine(CombinedBestFitObserved,0,CombinedBestFitObserved,len(Channels))
BestFitLine.SetLineWidth(1)
BestFitLine.SetLineColor(r.kBlack)
BestFitLine.Draw()
graph.Draw("P SAME")
#if len(intlumi)==2: mytext.DrawLatex(0.70,0.85,"#splitline{CMS preliminary}{#splitline{#sqrt{s} = %i TeV, L = %.1f fb^{-1}}{#sqrt{s} = %i TeV, L = %.1f fb^{-1}}}" %(int(Energy[0]),float(intlumi[0]),int(Energy[1]),float(intlumi[1])))
#else: mytext.DrawLatex(0.70,0.85,"#splitline{CMS preliminary}{#sqrt{s} = %i TeV, L = %.1f fb^{-1}}" %(int(Energy[0]),float(intlumi[0])))
can.RedrawAxis()
#Label1.Draw()
#Label2.Draw()
#myline.Draw()
if options.drawLabel:
  Label1.Draw()
  Label2.Draw()
  Label3.Draw()
  Label4.Draw()
  Label5.Draw()
if options.drawLines:
  cosline0.Draw()
  cosline1.Draw()
  cosline2.Draw()
  cosline3.Draw()
pt3.Draw()
leg.Draw()
pt.Draw()
pt2.Draw()
#dumbline.DrawLineNDC(0.7+0.0305,0.78-0.06,0.7+0.0305,0.78-0.037)
can.SaveAs("%s/BestFit%5.1f%s.pdf"%(options.outDir,Mass,options.outname))
can.SaveAs("%s/BestFit%5.1f%s.png"%(options.outDir,Mass,options.outname))
graph.Clear()
dummyHist.Clear()
BestFitBand.Clear()
BestFitLine.Clear()
leg.Clear()
r.gDirectory.Delete("*")


print "Done!"
