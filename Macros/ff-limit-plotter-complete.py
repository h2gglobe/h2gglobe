import ROOT
import array
ROOT.gROOT.ProcessLine(".L medianCalc.C++")
from ROOT import medianCalc
ROOT.gROOT.SetStyle("Plain")
# Get The ROOT files which contain the limit plotting.
# {"Mass":File}
intlumi = str(188.5)


###### Nick's python for limits plus additional stuff to plot SM lines


files = [
#	  ROOT.TFile("all_higgsCombineTest.mH100.root")
	 ROOT.TFile("all_higgsCombineTest.mH110.root")
	, ROOT.TFile("all_higgsCombineTest.mH120.root")
	, ROOT.TFile("all_higgsCombineTest.mH130.root")
	, ROOT.TFile("all_higgsCombineTest.mH140.root")
	]

masses = [#100
       110
       ,120
       ,130
       ,140]

obs = [
#3.95013
2.65206
,2.65206
,2.607
,4.07718
]

leg=ROOT.TLegend(0.55,0.6,0.85,0.85)
#leg=ROOT.TLegend(0.65,0.7,0.95,0.955)
#leg=ROOT.TLegend(0.15,0.6,0.4,0.85)
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

for i,mass,f in zip(range(len(files)),masses,files):
  
  median = array.array('d',[0])
  up68   = array.array('d',[0])
  dn68   = array.array('d',[0])
  up95   = array.array('d',[0])
  dn95   = array.array('d',[0])
  ff = array.array('d',[2.7455*0.6027E-01,2.2849*0.2237E-01, 1.91777*0.1073E-01, 1.66784*0.5623E-02])

  
  tree = f.Get("limit")
  medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95)
  graph68.SetPoint(i,float(mass),median[0]*ff[i])
  graph95.SetPoint(i,float(mass),median[0]*ff[i])
  graphMed.SetPoint(i,float(mass),median[0]*ff[i])
  #graphObs.SetPoint(i,float(mass),obs[i])
  graphObs.SetPoint(i,float(mass),median[0]*ff[i])
  graphOne.SetPoint(i,float(mass),1.*ff[i])

  diff95_up = abs(median[0] - up95[0])*ff[i]
  diff95_dn = abs(median[0] - dn95[0])*ff[i]
  diff68_up = abs(median[0] - up68[0])*ff[i] 
  diff68_dn = abs(median[0] - dn68[0])*ff[i] 

  
  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)
  graphObs.SetPointError(i,0,0,0,0)
  graphMed.SetPointError(i,0,0,0,0)
  graphOne.SetPointError(i,0,0,0,0)

#-------------------------------------------------------------------------

allMasses = array.array('d',[110,110.5,111,111.5,112,112.5,113,113.5,114,114.5,115,115.5,116,116.5,117,117.5,118,118.5,119,119.5,120,120.5,121,121.5,122,122.5,123,123.5,124,124.5,125,125.5,126,126.5,127,127.5,128,128.5,129,129.5,130,130.5,131,131.5,132,132.5,133,133.5,134,134.5,135,135.5,136,136.5,137,137.5,138,138.5,139,139.5])


 # VBF+WH+ZH
xSec = array.array('d',[2.7455,2.72068,2.69586,2.67104,2.64622,2.6214,2.59658,2.57176,2.54694,2.52212,2.4973,2.47606,2.45482,2.43358,2.41234,2.3911,2.36986,2.34862,2.32738,2.30614,2.2849,2.26638,2.24786,2.22934,2.21082,2.1923,2.17378,2.15526,2.13674,2.11822,2.0997,2.08299,2.06628,2.04957,2.03286,2.01615,1.99944,1.98273,1.96602,1.94931,1.9326,1.91777,1.90294,1.88811,1.87328,1.85845,1.84362,1.82879,1.81396,1.79913,1.7843,1.77136,1.75842,1.74548,1.73254,1.7196,1.70666,1.69372,1.68078,1.66784])

xSecErrPlusPercent = array.array('d',[4.1,4.12,4.14,4.16,4.18,4.2,4.22,4.24,4.26,4.28,4.3,4.25,4.2,4.15,4.1,4.05,4,3.95,3.9,3.85,3.8,3.79,3.78,3.77,3.76,3.75,3.74,3.73,3.72,3.71,3.7,3.71,3.72,3.73,3.74,3.75,3.76,3.77,3.78,3.79,3.8,3.83,3.86,3.89,3.92,3.95,3.98,4.01,4.04,4.07,4.1,4.09,4.08,4.07,4.06,4.05,4.04,4.03,4.02,4.01])
  

xSecErrMinusPercent = array.array('d',[4.5,4.52,4.54,4.56,4.58,4.6,4.62,4.64,4.66,4.68,4.7,4.64,4.58,4.52,4.46,4.4,4.34,4.28,4.22,4.16,4.1,4.12,4.14,4.16,4.18,4.2,4.22,4.24,4.26,4.28,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.25,4.2,4.15,4.1,4.05,4,3.95,3.9,3.85,3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98])

br = array.array('d',[0.6027E-01,0.5720E-01,0.5432E-01,0.5160E-01,0.4905E-01,0.4665E-01,0.4439E-01,0.4226E-01,0.4026E-01,0.3837E-01,0.3658E-01,0.3490E-01,0.3331E-01,0.3181E-01,0.3039E-01,0.2905E-01,0.2778E-01,0.2658E-01,0.2544E-01,0.2436E-01,0.2334E-01,0.2237E-01,0.2145E-01,0.2058E-01,0.1975E-01,0.1896E-01,0.1821E-01,0.1750E-01,0.1682E-01,0.1617E-01,0.1556E-01,0.1497E-01,0.1441E-01,0.1387E-01,0.1336E-01,0.1287E-01,0.1240E-01,0.1196E-01,0.1153E-01,0.1112E-01,0.1073E-01,0.1035E-01,0.9994E-02,0.9649E-02,0.9318E-02,0.9000E-02,0.8694E-02,0.8401E-02,0.8119E-02,0.7847E-02,0.7586E-02,0.7335E-02,0.7093E-02,0.6860E-02,0.6635E-02,0.6418E-02,0.6209E-02,0.6007E-02,0.5811E-02,0.5623E-02])


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
#  print "debug %f %f %f" %(,xSec10[i],dPlus10,dMinus10)  
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
  
  
  
  for i in range(len(allMasses)):
    myGraphXSecSM.SetPoint(i,allMasses[i],xSec[i])
    myGraphXSecSM.SetPointError(i,0,0,dMinus[i],dPlus[i])
    myGraphXSec10SM.SetPoint(i,allMasses[i],xSec10[i])
    myGraphXSec10SM.SetPointError(i,0,0,  dMinus10[i],dPlus10[i])
    myGraphXSec20SM.SetPoint(i,allMasses[i],xSec20[i])
    myGraphXSec20SM.SetPointError(i,0,0,  dMinus20[i],dPlus20[i])
    myGraphXSec30SM.SetPoint(i,allMasses[i],xSec30[i])
    myGraphXSec30SM.SetPointError(i,0,0,  dMinus30[i],dPlus30[i])
    myGraphXSec40SM.SetPoint(i,allMasses[i],xSec40[i])
    myGraphXSec40SM.SetPointError(i,0,0,  dMinus40[i],dPlus40[i])






#---------------------------------------------------------------------------
graph95.SetFillColor(ROOT.kYellow-4)
graph95.SetFillStyle(3001)
graph68.SetFillColor(ROOT.kGreen+2)
graph68.SetFillStyle(3001)
graphMed.SetLineStyle(2)
graphMed.SetLineColor(2)
graphMed.SetLineWidth(3)
graphObs.SetLineWidth(3)
graphOne.SetLineWidth(3)
graphOne.SetLineColor(ROOT.kAzure+7)

myGraphXSecSM.SetLineStyle(2)
myGraphXSecSM.SetLineColor(ROOT.kAzure+7)
myGraphXSecSM.SetLineWidth(4)
myGraphXSecSM.SetFillColor(ROOT.kAzure+7)
myGraphXSecSM.SetFillStyle(3003)
myGraphXSec10SM.SetLineStyle(2)
myGraphXSec10SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec10SM.SetLineWidth(4)
myGraphXSec10SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec10SM.SetFillStyle(3003)
myGraphXSec20SM.SetLineStyle(2)
myGraphXSec20SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec20SM.SetLineWidth(4)
myGraphXSec20SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec20SM.SetFillStyle(3003)
myGraphXSec30SM.SetLineStyle(2)
myGraphXSec30SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec30SM.SetLineWidth(4)
myGraphXSec30SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec30SM.SetFillStyle(3003)
myGraphXSec40SM.SetLineStyle(2)
myGraphXSec40SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec40SM.SetLineWidth(4)
myGraphXSec40SM.SetFillColor(ROOT.kAzure+7)
myGraphXSec40SM.SetFillStyle(3003)



# use 95 as the default guy

MG.Add(graph95)
MG.Add(graph68)
#MG.Add(graphOne)
MG.Add(graphMed)
MG.Add(graphObs)
MG.Add(myGraphXSecSM)
MG.Add(myGraphXSec10SM)
MG.Add(myGraphXSec20SM)
MG.Add(myGraphXSec30SM)
MG.Add(myGraphXSec40SM)

# -------------------------------------
C = ROOT.TCanvas("#int L = %s"%intlumi,"#int L = %s"%intlumi,1600,1100)
#C.SetLogy()
C.SetGrid(True)
MG.Draw("AL3")
MG.GetXaxis().SetTitle("m_{H}(GeV/c^{2})")
MG.GetXaxis().SetRangeUser(108,160)
MG.GetYaxis().SetRangeUser(0.0,3.)
MG.GetYaxis().SetTitle("\sigma_{H}xBR(H#rightarrow #gamma #gamma) - 95% CL")
MG.SetTitle("#int L = %s"%intlumi)
mytext = ROOT.TLatex()
mytext.SetTextSize(0.03)
mytext.DrawLatex(125,0.38,"10xFF")
mytext.DrawLatex(125,0.7,"20xFF")
mytext.DrawLatex(125,1.05,"30xFF")
mytext.DrawLatex(125,1.4,"40xFF")
leg.Draw()
C.SaveAs("ff-limit.pdf")
C.SaveAs("ff-limit.gif")
  
	

