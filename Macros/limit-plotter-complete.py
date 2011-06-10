# Original Author - Nicholas Wardle
# Run this with limit-plotter-complete.py METHOD model
# 	model  = sm (standard model) or ff (fermiophobic)
#	METHOD = ProfileLikelihood, Bayesian, Frequentist


import ROOT
import array,sys,numpy
ROOT.gROOT.ProcessLine(".L medianCalc.C++")
from ROOT import medianCalc
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)

if len(sys.argv) != 3: sys.exit("Need More Arguments - run with limit-plotter-complete.py METHOD model")
#-------------------------------------------------------------------------
intlumi = str(204)
Method = sys.argv[1]
#Method = "Bayesian"

EXPName = Method+"/expected"+Method

if Method == "ProfileLikelihood":
  OBSName = Method+"/higgsCombineTest."+Method
if Method == "Bayesian":
  OBSName = Method+"/higgsCombineOBSERVED.MarkovChainMC"

EXPmasses = [110,115,120,125,130,135,140]
OBSmasses = numpy.arange(110,140.5,0.5)
#-------------------------------------------------------------------------

if sys.argv[2] == "sm":
	from theory_sm import *
elif sys.argv[2] == "ff":
	from theory_ff import *
else: sys.exit("choose either sm or ff model")

###### Nick's python for limits plus additional stuff to plot SM lines
ROOT.gROOT.ProcessLine( \
   "struct Entry{	\
    double r;		\
   };"
)
from ROOT import Entry
def getOBSERVED(file):
  try:
   tree = file.Get("limit")
  except:
   return -1
  br = tree.GetBranch("limit")
  c = Entry()
  br.SetAddress(ROOT.AddressOf(c,'r'))
  tree.GetEntry(0)
  return c.r

EXPfiles = [ROOT.TFile(EXPName+".mH%d.root"%m) for m in EXPmasses]

# Get the observed limits - Currently only does up to 1 decimal mass points
OBSfiles = []
for m in OBSmasses:
  if int(m)==m:
    OBSfiles.append(ROOT.TFile(OBSName+".mH%d.root"%m))
  else:
    OBSfiles.append(ROOT.TFile(OBSName+".mH%.1f.root"%m))

obs = [getOBSERVED(O) for O in OBSfiles]
#-------------------------------------------------------------------------

# Set-up the GRAPHS
leg=ROOT.TLegend(0.15,0.6,0.4,0.85)
leg.SetFillColor(0)
leg.SetBorderSize(0)

graph68  = ROOT.TGraphAsymmErrors()
graph95  = ROOT.TGraphAsymmErrors()
graphMed = ROOT.TGraphAsymmErrors()
graphObs = ROOT.TGraphAsymmErrors()
graphOne = ROOT.TGraphAsymmErrors()

leg.AddEntry(graphObs,"Observed Limit","L")
leg.AddEntry(graphMed,"Expected Limit","L")
leg.AddEntry(graph68,"#pm 1#sigma Expected","F")
leg.AddEntry(graph95,"#pm 2#sigma Expected","F")

MG = ROOT.TMultiGraph()

#EXPECTED
sm = 1
for i,mass,f in zip(range(len(EXPfiles)),EXPmasses,EXPfiles):
  
  median = array.array('d',[0])
  up68   = array.array('d',[0])
  dn68   = array.array('d',[0])
  up95   = array.array('d',[0])
  dn95   = array.array('d',[0])

  for j,mm in enumerate(allMasses): 
	if mm==mass: sm = xSec[j]*br[j]
  
  tree = f.Get("limit")
  medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95)
  graph68.SetPoint(i,float(mass),median[0]*sm)
  graph95.SetPoint(i,float(mass),median[0]*sm)
  graphMed.SetPoint(i,float(mass),median[0]*sm)
  graphOne.SetPoint(i,float(mass),1.*sm)

  diff95_up = abs(median[0] - up95[0])*sm
  diff95_dn = abs(median[0] - dn95[0])*sm
  diff68_up = abs(median[0] - up68[0])*sm 
  diff68_dn = abs(median[0] - dn68[0])*sm 

  
  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)
  graphMed.SetPointError(i,0,0,0,0)
  graphOne.SetPointError(i,0,0,0,0)

#OBSERVED
for i,mass in zip(range(len(OBSfiles)),OBSmasses):

  if obs[i] ==-1: continue
  for j,mm in enumerate(allMasses): 
	if mm==mass: sm = xSec[j]*br[j]
  graphObs.SetPoint(i,float(mass),obs[i]*sm)
  graphObs.SetPointError(i,0,0,0,0)

#-------------------------------------------------------------------------

xSecErrPlus = array.array('d',[0.0]*len(allMasses))
xSecErrMinus = array.array('d',[0.0]*len(allMasses))

xSec10 = array.array('d',[0.0]*len(allMasses))
xSec20 = array.array('d',[0.0]*len(allMasses))
xSec30 = array.array('d',[0.0]*len(allMasses))
xSec40 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus10 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus20 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus30 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus40 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus10 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus20 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus30 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus40 = array.array('d',[0.0]*len(allMasses))
dPlus = array.array('d',[0.0]*len(allMasses))
dPlus10 = array.array('d',[0.0]*len(allMasses))
dPlus20 = array.array('d',[0.0]*len(allMasses))
dPlus30 = array.array('d',[0.0]*len(allMasses))
dPlus40 = array.array('d',[0.0]*len(allMasses))
dMinus = array.array('d',[0.0]*len(allMasses))
dMinus10 = array.array('d',[0.0]*len(allMasses))
dMinus20 = array.array('d',[0.0]*len(allMasses))
dMinus30 = array.array('d',[0.0]*len(allMasses))
dMinus40 = array.array('d',[0.0]*len(allMasses))

for i in range(len(allMasses)):
  xSecErrPlus[i] =  xSec[i] +  xSec[i]*xSecErrPlusPercent[i]/100.;
  xSecErrMinus[i] = xSec[i] -  xSec[i]*xSecErrMinusPercent[i]/100.;  
  xSec[i]=xSec[i]*br[i];
  xSecErrPlus[i] = xSecErrPlus[i]*br[i];
  xSecErrMinus[i]= xSecErrMinus[i]*br[i];
  xSec10[i]=10.*xSec[i];
  xSec20[i]=20.*xSec[i];
  xSec30[i]=30.*xSec[i];
  xSec40[i]=40.*xSec[i];
  
  xSecErrPlus10[i]=10.*xSecErrPlus[i]
  xSecErrPlus20[i]=20.*xSecErrPlus[i]
  xSecErrPlus30[i]=30.*xSecErrPlus[i]
  xSecErrPlus40[i]=40.*xSecErrPlus[i]
  xSecErrMinus10[i]=10.*xSecErrMinus[i]
  xSecErrMinus20[i]=20.*xSecErrMinus[i]
  xSecErrMinus30[i]=30.*xSecErrMinus[i]
  xSecErrMinus40[i]=40.*xSecErrMinus[i]
  dPlus[i]  =  abs(xSecErrPlus[i]-xSec[i])
  dMinus[i] =  abs(xSec[i] - xSecErrMinus[i])
  dPlus10[i]  =  abs(xSecErrPlus10[i]-xSec10[i])
  dMinus10[i] =  abs(xSec10[i] - xSecErrMinus10[i])
  dPlus20[i]  =  xSecErrPlus20[i]-xSec20[i]
  dMinus20[i] =  xSec20[i] - xSecErrMinus20[i]
  dPlus30[i]  =  xSecErrPlus30[i]-xSec30[i]
  dMinus30[i] =  xSec30[i] - xSecErrMinus30[i]
  dPlus40[i]  =  xSecErrPlus40[i]-xSec40[i]
  dMinus40[i] =  xSec40[i] - xSecErrMinus40[i]


  
  myGraphXSecSM   = ROOT.TGraphAsymmErrors()
  myGraphXSec10SM = ROOT.TGraphAsymmErrors()
  myGraphXSec20SM = ROOT.TGraphAsymmErrors()
  myGraphXSec30SM = ROOT.TGraphAsymmErrors()
  myGraphXSec40SM = ROOT.TGraphAsymmErrors()
  
  

# DONT put in the theory bands, error will be in limit itself  
  for i in range(len(allMasses)):
    myGraphXSecSM.SetPoint(i,allMasses[i],xSec[i])
    #myGraphXSecSM.SetPointError(i,0,0,dMinus[i],dPlus[i])
    myGraphXSec10SM.SetPoint(i,allMasses[i],xSec10[i])
    #myGraphXSec10SM.SetPointError(i,0,0,  dMinus10[i],dPlus10[i])
    myGraphXSec20SM.SetPoint(i,allMasses[i],xSec20[i])
    #myGraphXSec20SM.SetPointError(i,0,0,  dMinus20[i],dPlus20[i])
    myGraphXSec30SM.SetPoint(i,allMasses[i],xSec30[i])
    #myGraphXSec30SM.SetPointError(i,0,0,  dMinus30[i],dPlus30[i])
    myGraphXSec40SM.SetPoint(i,allMasses[i],xSec40[i])
    #myGraphXSec40SM.SetPointError(i,0,0,  dMinus40[i],dPlus40[i])

#---------------------------------------------------------------------------

FILLSTYLE=1001

FILLCOLOR_95=ROOT.kYellow-4
FILLCOLOR_68=ROOT.kGreen+2
FILLCOLOR_T=ROOT.kAzure+7

graph95.SetFillColor(FILLCOLOR_95)
graph95.SetFillStyle(FILLSTYLE)
graph68.SetFillColor(FILLCOLOR_68)
graph68.SetFillStyle(FILLSTYLE)
graphMed.SetLineStyle(2)
graphMed.SetLineColor(2)
graphMed.SetLineWidth(3)
graphObs.SetLineWidth(3)
graphOne.SetLineWidth(3)

myGraphXSecSM.SetLineStyle(2)
myGraphXSecSM.SetLineColor(FILLCOLOR_T)
myGraphXSecSM.SetLineWidth(4)
myGraphXSecSM.SetFillColor(FILLCOLOR_T)
myGraphXSecSM.SetFillStyle(FILLSTYLE)
myGraphXSec10SM.SetLineStyle(2)
myGraphXSec10SM.SetLineColor(FILLCOLOR_T)
myGraphXSec10SM.SetLineWidth(4)
myGraphXSec10SM.SetFillColor(FILLCOLOR_T)
myGraphXSec10SM.SetFillStyle(FILLSTYLE)
myGraphXSec20SM.SetLineStyle(2)
myGraphXSec20SM.SetLineColor(FILLCOLOR_T)
myGraphXSec20SM.SetLineWidth(4)
myGraphXSec20SM.SetFillColor(FILLCOLOR_T)
myGraphXSec20SM.SetFillStyle(FILLSTYLE)
myGraphXSec30SM.SetLineStyle(2)
myGraphXSec30SM.SetLineColor(FILLCOLOR_T)
myGraphXSec30SM.SetLineWidth(4)
myGraphXSec30SM.SetFillColor(FILLCOLOR_T)
myGraphXSec30SM.SetFillStyle(FILLSTYLE)
myGraphXSec40SM.SetLineStyle(2)
myGraphXSec40SM.SetLineColor(FILLCOLOR_T)
myGraphXSec40SM.SetLineWidth(4)
myGraphXSec40SM.SetFillColor(FILLCOLOR_T)
myGraphXSec40SM.SetFillStyle(FILLSTYLE)


# use 95 as the default guy
MG.Add(graph95)
MG.Add(graph68)
MG.Add(graphMed)
MG.Add(myGraphXSecSM)
MG.Add(myGraphXSec10SM)
MG.Add(myGraphXSec20SM)
MG.Add(myGraphXSec30SM)
MG.Add(myGraphXSec40SM)
MG.Add(graphObs)

# -------------------------------------
C = ROOT.TCanvas("#int L = %s"%intlumi,"#int L = %s"%intlumi,1600,1100)
#C.SetLogy()
C.SetGrid(True)
MG.Draw("AL3")
MG.GetXaxis().SetTitle("m_{H}(GeV/c^{2})")
MG.GetXaxis().SetRangeUser(min(OBSmasses),max(OBSmasses))
MG.GetYaxis().SetRangeUser(0.0,1.6)
MG.GetYaxis().SetTitle("\sigma_{H}xBR(H#rightarrow #gamma #gamma) - 95% CL")
MG.SetTitle("#int L = %s"%intlumi)
mytext = ROOT.TLatex()
mytext.SetTextSize(0.04)
mytext.DrawLatex(135,0.35,"10xSM")
mytext.DrawLatex(135,0.65,"20xSM")
mytext.DrawLatex(135,1.,"30xSM")
mytext.DrawLatex(135,1.3,"40xSM")
mytext.DrawLatex(132,1.45,"#int L = %d pb^{-1}"%int(intlumi) )
leg.Draw()
C.SaveAs("limit_%s_%s.pdf"%(sys.argv[2],Method))
C.SaveAs("limit_%s_%s.gif"%(sys.argv[2],Method))
C.SaveAs("limit_%s_%s.eps"%(sys.argv[2],Method))
C.SaveAs("limit_%s_%s.ps"%(sys.argv[2],Method))
