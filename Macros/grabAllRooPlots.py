import ROOT
import sys

grepmode=False
if len(sys.argv) < 2: 
	sys.exit("grabAllRooPlots.py usage:  python grabAllRooPlots.py FileName.root (optional -grep=expr , only plots containing expr in name are used)")

fileName = sys.argv[1]

if len(sys.argv) == 3:
  if "-grep=" in sys.argv[2]:
	grepmode = True
	grep	 = ((sys.argv[2]).split("="))[-1]

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetStyle("Plain")
F = ROOT.TFile(fileName)
keys = F.GetListOfKeys()
plots = []

if grepmode: print "grepping ", grep
for K in keys:
	obj = K.ReadObj()
	if grepmode:
	  if "plot" in obj.GetName() and grep in obj.GetName(): plots.append(obj.Clone())
	else:
	  if "plot" in obj.GetName(): plots.append(obj.Clone())

nPlots = len(plots)

if (nPlots >1):
  
  if int(nPlots**0.5)*int(nPlots**0.5) == nPlots: 
  	can = ROOT.TCanvas("c","Fit Plots",int(nPlots**0.5)*1200,int(nPlots**0.5)*900) # each plot needs to be pretty big
	can.Divide(int(nPlots**0.5),int(nPlots**0.5))

  else:
    sqrtnPlots = int(nPlots**0.5)+1
    can = ROOT.TCanvas("c","Fit Plots",sqrtnPlots*1200,sqrtnPlots*900) # each plot needs to be pretty big
    can.Divide(sqrtnPlots,sqrtnPlots)
else: 
  can = ROOT.TCanvas("c","Fit Plots",1200,900)
for i,P in enumerate(plots):
  can.cd(i+1)
  P.DrawClonePad()

can.SaveAs("allThePlots_%s.eps"%fileName)
