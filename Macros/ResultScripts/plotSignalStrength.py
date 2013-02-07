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
can.SetGrid(True)
leg = TLegend(0.65, 0.80, 0.89, 0.89)
leg.SetFillColor(0)
leg.SetBorderSize(1)
mytext = TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
Minimum=-2
Maximum=2.5
#intlumi=[5.1,12.2]
#Energy=[7,8]
intlumi=[5.2]
Energy=[8]

Masses=array.array("f",[x * 0.1 for x in range(1100,1501,5)])
directory = sys.argv[1]

BestFitObserved=[]
BestFitErrorUp=[]
BestFitErrorDown=[]

multigraph = TMultiGraph()
#multigraph.SetTitle(";M_{H} (GeV/c^{2});Best Fit #sigma/#sigma_{SM}")
#multigraph.SetMinimum(Minimum)
#multigraph.SetMaximum(Maximum)
dummyHist = TH1D("dummy","",1,110,150)
dummyHist.SetTitle(";M_{H} (GeV/c^{2});Best Fit #sigma/#sigma_{SM}")
dummyHist.SetTitleSize(.06,"X")
dummyHist.SetTitleOffset(0.65,"X")
dummyHist.SetTitleSize(.06,"Y")
dummyHist.SetTitleOffset(0.6,"Y")
dummyHist.SetMinimum(Minimum)
dummyHist.SetMaximum(Maximum)

for Mass in Masses:
  file = TFile.Open(directory+"/higgsCombineSignalStrength.ChannelCompatibilityCheck.mH"+str(Mass).replace(".0","")+".root")
  print file.GetName()
  normal=file.Get("fit_nominal")
  rFit=normal.floatParsFinal().find("r")
  #print "Value: %f ErrorDown: %f ErrorUp: %f" %(rFit.getVal(),rFit.getAsymErrorLo(),rFit.getAsymErrorHi())
  BestFitObserved.append(rFit.getVal())
  BestFitErrorUp.append(rFit.getAsymErrorHi())
  if math.fabs(rFit.getAsymErrorLo())<0.000001: BestFitErrorDown.append(math.fabs(rFit.getAsymErrorHi()))
  else: BestFitErrorDown.append(math.fabs(rFit.getAsymErrorLo()))

graph=TGraph(len(BestFitObserved),Masses,array.array("f",BestFitObserved))
graph.SetMarkerStyle(21)
graph.SetMarkerSize(0.5)
graph.SetLineWidth(1)
graph.SetFillColor(kGreen)
graph.SetLineColor(kBlack)
graphsigma=TGraphAsymmErrors(len(BestFitObserved),Masses,array.array("f",BestFitObserved),array.array("f",[0]*len(Masses)),array.array("f",[0]*len(Masses)),array.array("f",BestFitErrorDown),array.array("f",BestFitErrorUp))
graphsigma.SetLineColor(kGreen)
graphsigma.SetFillColor(kGreen)
multigraph.Add(graphsigma)
multigraph.Add(graph)
leg.AddEntry(graph,"68% CL Band","F")

dummyHist.Draw("AXIS")
multigraph.Draw("3LP")
dummyHist.Draw("AXIGSAME")
leg.Draw("")

line = TLine(110, 1, 150, 1)
line.SetLineWidth(2)
line.SetLineColor(kRed)
line.Draw()
can.RedrawAxis();
can.SetGrid(True)

if len(intlumi)==2: mytext.DrawLatex(0.13,0.805,"#splitline{CMS preliminary}{#splitline{#sqrt{s} = %i TeV L = %.1f fb^{-1}}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}}" %(int(Energy[0]),float(intlumi[0]),int(Energy[1]),float(intlumi[1])))
else: mytext.DrawLatex(0.13,0.83,"#splitline{CMS preliminary}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}" %(int(Energy[0]),float(intlumi[0])))

can.SaveAs(sys.argv[2])
print "Done!"
