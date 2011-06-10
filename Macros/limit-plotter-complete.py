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
	 ROOT.TFile("higgsCombineTEST.ProfileLikelihood.mH105.1439293056.root")
	, ROOT.TFile("higgsCombineTEST.ProfileLikelihood.mH110.350559743.root")
	, ROOT.TFile("higgsCombineTEST.ProfileLikelihood.mH120.1452375510.root")
	, ROOT.TFile("higgsCombineTEST.ProfileLikelihood.mH130.-1837255987.root")
        , ROOT.TFile("higgsCombineTEST.ProfileLikelihood.mH140.244931859.root")
	]

masses = [105
       ,110
       ,120
       ,130
       ,140]

obs = [
6.08,
8.223,
6.051,
3.637,
11.5
]

#leg=ROOT.TLegend(0.55,0.6,0.85,0.85)
#leg=ROOT.TLegend(0.65,0.7,0.95,0.955)
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

for i,mass,f in zip(range(len(files)),masses,files):
  
  median = array.array('d',[0])
  up68   = array.array('d',[0])
  dn68   = array.array('d',[0])
  up95   = array.array('d',[0])
  dn95   = array.array('d',[0])
  sm = array.array('d',[24.9582*0.00178,22.7112*0.00198,19.0125*0.00225,16.1292*0.00226,13.8456*0.00194])

  
  tree = f.Get("limit")
  medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95)
  graph68.SetPoint(i,float(mass),median[0]*sm[i])
  graph95.SetPoint(i,float(mass),median[0]*sm[i])
  graphMed.SetPoint(i,float(mass),median[0]*sm[i])
  #graphObs.SetPoint(i,float(mass),obs[i])
  graphObs.SetPoint(i,float(mass),median[0]*sm[i])
  graphOne.SetPoint(i,float(mass),1.*sm[i])

  diff95_up = abs(median[0] - up95[0])*sm[i]
  diff95_dn = abs(median[0] - dn95[0])*sm[i]
  diff68_up = abs(median[0] - up68[0])*sm[i] 
  diff68_dn = abs(median[0] - dn68[0])*sm[i] 

  
  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)
  graphObs.SetPointError(i,0,0,0,0)
  graphMed.SetPointError(i,0,0,0,0)
  graphOne.SetPointError(i,0,0,0,0)

#-------------------------------------------------------------------------

allMasses = array.array('d',[105,105.5,106,106.5,107,107.5,108,108.5,109,109.5,110,110.5,111,111.5,112,112.5,113,113.5,114,114.5,115,115.5,116,116.5,117,117.5,118,118.5,119,119.5,120,120.5,121,121.5,122,122.5,123,123.5,124,124.5,125,125.5,126,126.5,127,127.5,128,128.5,129,129.5,130,130.5,131,131.5,132,132.5,133,133.5,134,134.5,135,135.5,136,136.5,137,137.5,138,138.5,139,139.5,140])
  
xSec = array.array('d',[24.9582,24.7335,24.5088,24.2841,24.0594,23.8347,23.61,23.3853,23.1606,22.9359,22.7112,22.5139,22.3165,22.1192,21.9219,21.7246,21.5272,21.3299,21.1326,20.9352,20.7379,20.5654,20.3928,20.2203,20.0477,19.8752,19.7027,19.5301,19.3576,19.185,19.0125,18.8608,18.7092,18.5576,18.4059,18.2542,18.1026,17.951,17.7993,17.6476,17.496,17.3593,17.2226,17.086,16.9493,16.8126,16.6759,16.5392,16.4026,16.2659,16.1292,16.0095,15.8898,15.7702,15.6505,15.5308,15.4111,15.2914,15.1718,15.0521,14.9324,14.8237,14.715,14.6064,14.4977,14.389,14.2803,14.1716,14.063,13.9543,13.8456])  


xSecErrPlusPercent = array.array('d',[20.8,20.76,20.72,20.68,20.64,20.6,20.56,20.52,20.48,20.44,20.4,20.36,20.32,20.28,20.24,20.2,20.16,20.12,20.08,20.04,20,19.97,19.94,19.91,19.88,19.85,19.82,19.79,19.76,19.73,19.7,19.68,19.66,19.64,19.62,19.6,19.58,19.56,19.54,19.52,19.5,19.47,19.44,19.41,19.38,19.35,19.32,19.29,19.26,19.23,19.2,19.17,19.14,19.11,19.08,19.05,19.02,18.99,18.96,18.93,18.9,18.89,18.88,18.87,18.86,18.85,18.84,18.83,18.82,18.81,18.8])

xSecErrMinusPercent = array.array('d',[15.5,15.48,15.46,15.44,15.42,15.4,15.38,15.36,15.34,15.32,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.3,15.28,15.26,15.24,15.22,15.2,15.18,15.16,15.14,15.12,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.1,15.09,15.08,15.07,15.06,15.05,15.04,15.03,15.02,15.01,15,14.99,14.98,14.97,14.96,14.95,14.94,14.93,14.92,14.91,14.9])


br = array.array('d',[0.00178,0.0018003,0.00182029,0.00183998,0.00185938,0.00187849,0.00189733,0.00191589,0.00193419,0.00195222,0.00197,0.00198737,0.00200441,0.00202113,0.00203754,0.00205365,0.00206948,0.00208501,0.00210028,0.00211527,0.00213,0.00214331,0.0021563,0.00216899,0.00218138,0.00219348,0.0022053,0.00221685,0.00222815,0.0022392,0.00225,0.00225571,0.00226125,0.00226662,0.00227182,0.00227687,0.00228177,0.00228652,0.00229114,0.00229563,0.0023,0.00229526,0.00229072,0.00228635,0.00228215,0.00227811,0.00227422,0.00227047,0.00226686,0.00226337,0.00226,0.00224526,0.00223124,0.00221791,0.0022052,0.00219308,0.00218151,0.00217044,0.00215986,0.00214972,0.00214,0.00211438,0.00209031,0.00206765,0.00204629,0.00202612,0.00200703,0.00198895,0.0019718,0.0019555,0.00194])                      

xSec10 = array.array('d',[0.0]*len(allMasses))
xSec20 = array.array('d',[0.0]*len(allMasses))
xSec30 = array.array('d',[0.0]*len(allMasses))
xSec40 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus = array.array('d',[0.0]*len(allMasses))
xSecErrPlus10 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus20 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus30 = array.array('d',[0.0]*len(allMasses))
xSecErrPlus40 = array.array('d',[0.0]*len(allMasses))
xSecErrMinus = array.array('d',[0.0]*len(allMasses))
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
#    myGraphXSecSM.SetPointError(i,0,0,dMinus[i],dPlus[i])
    myGraphXSec10SM.SetPoint(i,allMasses[i],xSec10[i])
#    myGraphXSec10SM.SetPointError(i,0,0,  dMinus10[i],dPlus10[i])
    myGraphXSec20SM.SetPoint(i,allMasses[i],xSec20[i])
#    myGraphXSec20SM.SetPointError(i,0,0,  dMinus20[i],dPlus20[i])
    myGraphXSec30SM.SetPoint(i,allMasses[i],xSec30[i])
#    myGraphXSec30SM.SetPointError(i,0,0,  dMinus30[i],dPlus30[i])
    myGraphXSec40SM.SetPoint(i,allMasses[i],xSec40[i])
#    myGraphXSec40SM.SetPointError(i,0,0,  dMinus40[i],dPlus40[i])






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

myGraphXSecSM.SetLineStyle(2)
myGraphXSecSM.SetLineColor(ROOT.kAzure+7)
myGraphXSecSM.SetLineWidth(4)
#myGraphXSecSM.SetFillColor(ROOT.kAzure+7)
#myGraphXSecSM.SetFillStyle(3003)
myGraphXSec10SM.SetLineStyle(2)
myGraphXSec10SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec10SM.SetLineWidth(4)
#myGraphXSec10SM.SetFillColor(ROOT.kAzure+7)
#myGraphXSec10SM.SetFillStyle(3003)
myGraphXSec20SM.SetLineStyle(2)
myGraphXSec20SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec20SM.SetLineWidth(4)
#myGraphXSec20SM.SetFillColor(ROOT.kAzure+7)
#myGraphXSec20SM.SetFillStyle(3003)
myGraphXSec30SM.SetLineStyle(2)
myGraphXSec30SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec30SM.SetLineWidth(4)
#myGraphXSec30SM.SetFillColor(ROOT.kAzure+7)
#myGraphXSec30SM.SetFillStyle(3003)
myGraphXSec40SM.SetLineStyle(2)
myGraphXSec40SM.SetLineColor(ROOT.kAzure+7)
myGraphXSec40SM.SetLineWidth(4)
#myGraphXSec40SM.SetFillColor(ROOT.kAzure+7)
#myGraphXSec40SM.SetFillStyle(3003)



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
MG.GetXaxis().SetRangeUser(100,140.5)
MG.GetYaxis().SetRangeUser(0.0,1.6)
MG.GetYaxis().SetTitle("\sigma_{H}xBR(H#rightarrow #gamma #gamma) - 95% CL")
MG.SetTitle("#int L = %s"%intlumi)
mytext = ROOT.TLatex()
mytext.SetTextSize(0.04)
mytext.DrawLatex(135,0.35,"10xSM")
mytext.DrawLatex(135,0.65,"20xSM")
mytext.DrawLatex(135,1.,"30xSM")
mytext.DrawLatex(135,1.3,"40xSM")
leg.Draw()
C.SaveAs("limit.pdf")
C.SaveAs("limit.gif")
  
	

