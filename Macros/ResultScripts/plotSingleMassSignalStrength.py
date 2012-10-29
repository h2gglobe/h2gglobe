print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import os
import sys
import array
import math

gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);
gROOT.SetBatch(1)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,600)
leg = TLegend(0.65, 0.75, 0.89, 0.89)
leg.SetFillColor(0)
leg.SetBorderSize(1)
mytext = TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
rMin=-15
rMax=15
#intlumi=[5.1,3.8]
#Energy=[7,8]
intlumi=[12.2]
Energy=[8]

#Masses=[x * 0.1 for x in range(1100,1501,5)]
#Masses=[123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5]
Masses=[125.0]
directory = sys.argv[1]

for Mass in Masses:
  file = TFile.Open(directory+"/higgsCombineSignalStrength.ChannelCompatibilityCheck.mH"+str(Mass)+".root")
  normal=file.Get("fit_nominal")
  rFit=normal.floatParsFinal().find("r")
  CombinedBestFitObserved=rFit.getVal()
  CombinedBestFitErrorUp=rFit.getAsymErrorHi()
  CombinedBestFitErrorDown=rFit.getAsymErrorLo()
  if math.fabs(CombinedBestFitErrorDown)<0.000001: CombinedBestFitErrorDown=rFit.getAsymErrorHi()
  altername=file.Get("fit_alternate")
  Channels=[]
  BestFitObserved=[]
  BestFitErrorUp=[]
  BestFitErrorDown=[]
  for i in range(altername.floatParsFinal().getSize()):
    print i,altername.floatParsFinal().at(i).GetName(),altername.floatParsFinal().at(i).getVal(),altername.floatParsFinal().at(i).getAsymErrorHi(),altername.floatParsFinal().at(i).getAsymErrorLo()
    chrFit=altername.floatParsFinal().at(i)
    if chrFit.GetName().find("ChannelCompatibilityCheck")!=-1:
      Channels.append(chrFit.GetName())
      BestFitObserved.append(chrFit.getVal())
      BestFitErrorUp.append(chrFit.getAsymErrorHi())
      if math.fabs(chrFit.getAsymErrorLo())<0.000001: BestFitErrorDown.append(math.fabs(chrFit.getAsymErrorHi()))
      else: BestFitErrorDown.append(math.fabs(chrFit.getAsymErrorLo()))
  dummyHist = TH2F("dummy",";Best Fit #sigma/#sigma_{SM};",1,rMin,rMax,len(Channels),0,len(Channels))
  ChannelNum=array.array("f",[x + 0.5 for x in range(len(Channels))])
  graph=TGraphAsymmErrors(len(Channels),array.array("f",BestFitObserved),ChannelNum,array.array("f",BestFitErrorDown),array.array("f",BestFitErrorUp),array.array("f",[0]*len(Channels)),array.array("f",[0]*len(Channels)))
  graph.SetLineColor(kRed)
  graph.SetLineWidth(1)
  graph.SetMarkerStyle(21)
  graph.SetMarkerSize(2)
  graph.SetMarkerColor(kBlack)
  graph.SetFillColor(kWhite)
  leg.AddEntry(graph,"BestFit Category")
  for i in range(len(Channels)):
    dummyHist.GetYaxis().SetBinLabel(i+1, Channels[i][-4:]);
  dummyHist.SetTitleSize(0.06,"X")
  dummyHist.SetTitleOffset(0.7,"X")
  dummyHist.SetLabelSize(0.07,"Y")
  dummyHist.SetLineColor(kBlack)
  dummyHist.SetFillColor(kGreen)
  leg.AddEntry(dummyHist,"68% Combined")
  dummyHist.Draw()
  BestFitBand=TBox(CombinedBestFitObserved+CombinedBestFitErrorDown,0,CombinedBestFitObserved+CombinedBestFitErrorUp,len(Channels))
  #BestFitBand.SetFillStyle(3013)
  BestFitBand.SetFillColor(kGreen)
  BestFitBand.SetLineStyle(0)
  BestFitBand.Draw()
  BestFitLine=TLine(CombinedBestFitObserved,0,CombinedBestFitObserved,len(Channels))
  BestFitLine.SetLineWidth(1)
  BestFitLine.SetLineColor(kBlack)
  BestFitLine.Draw()
  graph.Draw("P SAME")
  leg.Draw()
  if len(intlumi)==2: mytext.DrawLatex(0.13,0.8,"#splitline{CMS preliminary}{#splitline{#splitline{M_{H} = %.1f (GeV/c^{2})}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}}}" %(Mass,int(Energy[0]),float(intlumi[0]),int(Energy[1]),float(intlumi[1])))
  else: mytext.DrawLatex(0.13,0.8,"#splitline{CMS preliminary}{#splitline{M_{H} = %.1f (GeV/c^{2})}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}}" %(Mass,int(Energy[0]),float(intlumi[0])))
  can.RedrawAxis()
  can.SaveAs("BestFit"+str(Mass)+"GeV.pdf")
  graph.Clear()
  dummyHist.Clear()
  BestFitBand.Clear()
  BestFitLine.Clear()
  leg.Clear()
  gDirectory.Delete("*")


print "Done!"
