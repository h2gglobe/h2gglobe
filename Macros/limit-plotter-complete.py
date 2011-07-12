# Original Authors - Nicholas Wardle, Nancy Marinelli, Doug Berry

# Run this with limit-plotter-complete.py METHOD model
# 	model  = sm (standard model) or ff (fermiophobic)
#	METHOD = ProfileLikelihood, Bayesian, Frequentist

# Standard Imports and calculators
import ROOT
import array,sys,numpy
ROOT.gROOT.ProcessLine(".L medianCalc.C++")
from ROOT import medianCalc
from ROOT import FrequentistLimits
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)

#-------------------------------------------------------------------------
# Configuration for the Plotter
intlumi = str(1.09)
EXPmasses = [110,115,120,125,130,135,140]
OBSmasses = numpy.arange(110,140.5,0.5)
theorySMScales = [5,10,15,20]


FILLSTYLE=1001
FILLCOLOR_95=ROOT.kYellow-4
FILLCOLOR_68=ROOT.kGreen+2
FILLCOLOR_T=ROOT.kAzure+7			# Theory lines color
#-------------------------------------------------------------------------
# UserInput
doRatio = False
if len(sys.argv) < 3: sys.exit("Need More Arguments - run with limit-plotter-complete.py METHOD model")
elif len(sys.argv) ==4 and sys.argv[3] == "-doRatio": doRatio = True
 
Method = sys.argv[1]
EXPName = Method+"/expected"+Method

if Method == "ProfileLikelihood":
  OBSName = Method+"/higgsCombineTest."+Method
if Method == "Bayesian":
  OBSName = Method+"/higgsCombineOBSERVED.MarkovChainMC"
if Method == "Frequentist":
  OBSName = Method+"/higgsCombineOBSERVED.Frequentist"

if Method == "Frequentist": EXPmasses = OBSmasses[:]

if sys.argv[2] == "sm":
	from theory_sm import *
	extraString = "SM"
elif sys.argv[2] == "ff":
	from theory_ff import *
	extraString = "FP"
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

if Method=="Frequentest":
  EXPfiles = [ROOT.TFile(EXPName+".mH%d.quant0.500.root"%m) for m in EXPmasses]
else:
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

  if not doRatio:
    for j,mm in enumerate(allMasses): 
	if mm==mass: sm = xSec[j]*br[j]
  
  tree = f.Get("limit")
  medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95)
  graph68.SetPoint(i,float(mass),median[0]*sm)
  graph95.SetPoint(i,float(mass),median[0]*sm)
  graphMed.SetPoint(i,float(mass),median[0]*sm)
  graphOne.SetPoint(i,float(mass),1.*sm)
  
  if Method == "Frequentist":
    diff95_up = abs(median[0] - FrequentistLimits(f.GetName().replace("0.500.root","0.975.root")))*sm
    diff95_dn = abs(median[0] - FrequentistLimits(f.GetName().replace("0.500.root","0.027.root")))*sm
    diff68_up = abs(median[0] - FrequentistLimits(f.GetName().replace("0.500.root","0.840.root")))*sm
    diff68_dn = abs(median[0] - FrequentistLimits(f.GetName().replace("0.500.root","0.160.root")))*sm
  else:
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
  if not doRatio:
    for j,mm in enumerate(allMasses): 
	if mm==mass: sm = xSec[j]*br[j]
  graphObs.SetPoint(i,float(mass),obs[i]*sm)
  graphObs.SetPointError(i,0,0,0,0)


#-------------------------------------------------------------------------
# Construct the theory bands

theoryArrays = []
theoryPlusArrays = []
theoryMinusArrays = []

dPlusArrays = []
dMinusArrays = []

xSecErrPlus = array.array('d',[0.0]*len(allMasses))
xSecErrMinus = array.array('d',[0.0]*len(allMasses))
dPlus = array.array('d',[0.0]*len(allMasses))

for th in theorySMScales:
  theoryArrays.append(array.array('d',[0.0]*len(allMasses)))
  theoryPlusArrays.append(array.array('d',[0.0]*len(allMasses)))
  theoryMinusArrays.append(array.array('d',[0.0]*len(allMasses)))
  dPlusArrays.append(array.array('d',[0.0]*len(allMasses)))
  dMinusArrays.append(array.array('d',[0.0]*len(allMasses)))


for i in range(len(allMasses)):
    xSecErrPlus[i] =  xSec[i] +  xSec[i]*xSecErrPlusPercent[i]/100.;
    xSecErrMinus[i] = xSec[i] -  xSec[i]*xSecErrMinusPercent[i]/100.;  

    if doRatio: xSec[i]=1.;
    else: xSec[i]=xSec[i]*br[i];

    xSecErrPlus[i] = xSecErrPlus[i]*br[i];
    xSecErrMinus[i]= xSecErrMinus[i]*br[i];

    for j,th in enumerate(theorySMScales):
      theoryArrays[j][i]=th*xSec[i];
      theoryPlusArrays[j][i]=th*xSecErrPlus[i];
      theoryMinusArrays[j][i]=th*xSecErrMinus[i];
      dPlusArrays[j][i]=abs(th*xSecErrPlus[i]-th*xSec[i])
      dMinusArrays[j][i]=abs(th*xSecErrMinus[i]-th*xSec[i])
    
  
myGraphXSecSM   = ROOT.TGraphAsymmErrors()
myGraphXSecSMScales = []
for th in theorySMScales:
  myGraphXSecSMScales.append(ROOT.TGraphAsymmErrors()) 
  
for i in range(len(allMasses)):
    myGraphXSecSM.SetPoint(i,allMasses[i],xSec[i])
    for j in range(len(theorySMScales)):
      myGraphXSecSMScales[j].SetPoint(i,allMasses[i],theoryArrays[j][i])

#---------------------------------------------------------------------------
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

MG.Add(graph95)
MG.Add(graph68)
MG.Add(graphMed)
MG.Add(myGraphXSecSM)

if not doRatio:
 for j in range(len(theorySMScales)):
  myGraphXSecSMScales[j].SetLineStyle(2)
  myGraphXSecSMScales[j].SetLineColor(FILLCOLOR_T)
  myGraphXSecSMScales[j].SetLineWidth(4)
  myGraphXSecSMScales[j].SetFillColor(FILLCOLOR_T)
  myGraphXSecSMScales[j].SetFillStyle(FILLSTYLE)
  MG.Add(myGraphXSecSMScales[j])

# use 95 as the default guy
MG.Add(graphObs)

# -------------------------------------
C = ROOT.TCanvas("#int L = %s"%intlumi,"#int L = %s"%intlumi,1600,1100)
#C.SetLogy()
C.SetGrid(True)
MG.Draw("AL3")
MG.GetXaxis().SetTitle("m_{H}(GeV/c^{2})")
MG.GetXaxis().SetRangeUser(min(OBSmasses),max(OBSmasses))
if doRatio:
 MG.GetYaxis().SetRangeUser(0.0,10)
 MG.GetYaxis().SetTitle("\sigma_{H}xBR(H#rightarrow #gamma #gamma) / %s - 95%% CL"%extraString)
else: 
 MG.GetYaxis().SetRangeUser(0.0,0.6)
 MG.GetYaxis().SetTitle("\sigma_{H}xBR(H#rightarrow #gamma #gamma) - 95%% CL")

MG.SetTitle("#int L = %s"%intlumi)
mytext = ROOT.TLatex()
mytext.SetTextSize(0.04)
SMEnd = myGraphXSecSM.Eval(max(OBSmasses))
print SMEnd
if not doRatio:
 mytext.DrawLatex(max(OBSmasses)+.3,SMEnd,"SM")
 for th in theorySMScales:
  mytext.DrawLatex(max(OBSmasses)+0.3,th*SMEnd,"%dx%s"%(th,extraString))

mytext.SetNDC()
mytext.DrawLatex(0.1,0.93,"CMS Preliminary 2011, #int L = %.2f fb^{-1}"%float(intlumi) )
leg.Draw()

#Make a bunch of extensions to the plots
C.SaveAs("limit_%s_%s.pdf"%(sys.argv[2],Method))
C.SaveAs("limit_%s_%s.gif"%(sys.argv[2],Method))
C.SaveAs("limit_%s_%s.eps"%(sys.argv[2],Method))
C.SaveAs("limit_%s_%s.ps"%(sys.argv[2],Method))

OUTTgraphs = ROOT.TFile("outputTGraphAsymmErrors_%s.root"%sys.argv[2],"RECREATE")
graphMed.SetName("median")
graphMed.Write()
graphObs.SetName("observed")
graphObs.Write()
OUTTgraphs.Write()

