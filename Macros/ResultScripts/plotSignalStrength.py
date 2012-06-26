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
Minimum=-3
Maximum=3
#intlumi=[5.1,3.8]
#Energy=[7,8]
intlumi=[3.8]
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
dummyHist.SetTitleOffset(0.5,"Y")
dummyHist.SetMinimum(Minimum)
dummyHist.SetMaximum(Maximum)

for Mass in Masses:
  file = TFile.Open(directory+"/higgsCombineSignalStrength.ChannelCompatibilityCheck.mH"+str(Mass)+".root")
  normal=file.Get("fit_nominal")
  #print "Value: %f ErrorDown: %f ErrorUp: %f" %(normal.floatParsFinal().find("r").getVal(),normal.floatParsFinal().find("r").getAsymErrorLo(),normal.floatParsFinal().find("r").getAsymErrorHi())
  BestFitObserved.append(normal.floatParsFinal().find("r").getVal())
  BestFitErrorUp.append(normal.floatParsFinal().find("r").getAsymErrorHi())
  if math.fabs(normal.floatParsFinal().find("r").getAsymErrorLo())<0.000001: BestFitErrorDown.append(normal.floatParsFinal().find("r").getAsymErrorHi())
  else: BestFitErrorDown.append(normal.floatParsFinal().find("r").getAsymErrorLo())

graph=TGraph(len(BestFitObserved),Masses,array.array("f",BestFitObserved))
graph.SetLineColor(kBlack)
graph.SetMarkerStyle(21)
graph.SetMarkerSize(0.5)
graph.SetLineWidth(1)
graph.SetFillColor(kGreen)
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

#if len(intlumi)==2: mytext.DrawLatex(0.18,0.82,"#splitline{CMS preliminary}{#splitline{#sqrt{s} = %i TeV L = %.1f fb^{-1}}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}}" %(int(Energy[0]),float(intlumi[0]),int(Energy[1]),float(intlumi[1])))
#else: mytext.DrawLatex(0.18,0.8,"#splitline{CMS preliminary}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}" %(int(Energy[0]),float(intlumi[0])))

can.SaveAs(sys.argv[2])
print "Done!"
