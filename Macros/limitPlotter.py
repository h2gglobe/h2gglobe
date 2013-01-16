#!/usr/bin/env python
# Original Authors - Nicholas Wardle, Nancy Marinelli, Doug Berry

# Major cleanup from limit-plotter-complete.py
#-------------------------------------------------------------------------
# UserInput
from optparse import OptionParser
parser=OptionParser()
parser.add_option("-M","--Method",dest="Method")
parser.add_option("-r","--doRatio",action="store_true")
parser.add_option("-s","--doSmooth",action="store_true")
parser.add_option("-b","--bayes",dest="bayes")
parser.add_option("-o","--outputLimits",dest="outputLimits")
parser.add_option("-e","--expectedOnly",action="store_true")
parser.add_option("-p","--path",dest="path",default="",type="str")
parser.add_option("-v","--verbose",dest="verbose",action="store_true")
parser.add_option("-a","--append",dest="append",default="",help="Append string to filename")
parser.add_option("","--sideband",dest="sideband",default=False,action="store_true")
parser.add_option("","--addline",action="append",type="str",help="add lines to the plot file.root:color:linestyle:legend entry", default = [])
parser.add_option("","--show",action="store_true")
parser.add_option("","--pval",action="store_true")
parser.add_option("","--addtxt",action="append",type="str", help="Add lines of text under CMS Preliminary",default=[])
parser.add_option("","--square",dest="square",help="Make square plots",action="store_true")
parser.add_option("","--nogrid",dest="nogrid",help="Remove grid from plots",action="store_true")
parser.add_option("","--nobox",dest="nobox",action="store_true",default=False,help="Don't draw box around text")
(options,args)=parser.parse_args()

# Standard Imports and calculators
import ROOT
import array,sys,numpy
ROOT.gROOT.ProcessLine(".x tdrstyle.cc")

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

#-------------------------------------------------------------------------
# Configuration for the Plotter
OBSmasses = []
EXPmasses = []

#OBSmassesT = [110,115,120,125,130,135,140,145,150] 
OBSmassesT = numpy.arange(110,150.1,1)
EXPmassesT = numpy.arange(110,150.1,1)
epsilon = 0.001  # make this smaller than your smallest step size

for m in OBSmassesT:
        #if "%.1f"%m=="%d.0"%(m+epsilon):continue   # sigh!
    OBSmasses.append(m)
    EXPmasses.append(m)

# Plotting Styles --------------------------------------------------------
OFFSETLOW=0
OFFSETHIGH=0
FONTSIZE=0.035
FILLSTYLE=1001
SMFILLSTYLE=3244
FILLCOLOR_95=ROOT.kYellow
FILLCOLOR_68=ROOT.kGreen
RANGEYABS=[0.0,0.6]
RANGEYRAT=[0.0,4]
RANGEMU = [-4,3.0]
MINPV = 1.0*10E-8
MAXPV = 1.0
Lines = [1.,2.,3.,4.,5.]
MINMH=int(min(EXPmasses))
MAXMH=int(max(EXPmasses))

if options.show : ROOT.gROOT.SetBatch(False)
if options.addline and not options.pval : sys.exit("Cannot addlines unless running in pvalue")

# ------------------------------------------------------------------------
# SM Signal Normalizer
if not options.doRatio:
    ROOT.gROOT.ProcessLine(".L ../libLoopAll.so")
    signalNormalizer = ROOT.Normalization_8TeV()
extraString = "SM"
# ------------------------------------------------------------------------
if options.pval:
     EXPmasses=[]
     options.doRatio=True

if options.Method=="MaxLikelihoodFit": 
    options.doRatio = True

if not options.doRatio and options.Method != "Frequentist": 
    ROOT.gROOT.ProcessLine(".L medianCalc.C+g")
    from ROOT import medianCalc
    from ROOT import FrequentistLimits
    

if options.bayes:
  BayesianFile =ROOT.TFile(options.bayes) 
  bayesObs = BayesianFile.Get("observed")
 
Method = options.Method
if not options.path: options.path=Method


EXPName = options.path+"/higgsCombineEXPECTED."+Method
if Method =="MaxLikelihoodFit": EXPName = options.path+"/higgsCombineTest."+Method
if Method == "Asymptotic" or Method == "AsymptoticNew":  EXPName = options.path+"/higgsCombineTest."+Method  # everyhting contained here
if Method == "ProfileLikelihood" or Method=="Asymptotic" or Method=="AsymptoticNew" or Method=="MaxLikelihoodFit":
  OBSName = options.path+"/higgsCombineTest."+Method
if Method == "Bayesian":
  OBSName = options.path+"/higgsCombineOBSERVED.MarkovChainMC"
if Method == "HybridNew":
  OBSName = options.path+"/higgsCombineOBSERVED.HybridNew"


if Method == "HybridNew" or Method == "Asymptotic" or Method == "AsymptoticNew": EXPmasses = OBSmasses[:]
# Done with setip
# --------------------------------------------------------------------------
ROOT.gROOT.ProcessLine( \
   "struct Entry{   \
    double r;       \
   };"
)
from ROOT import Entry
def getOBSERVED(file,entry=0):
  try:
   tree = file.Get("limit")
  except:
   return -1
  br = tree.GetBranch("limit")
  c = Entry()
  br.SetAddress(ROOT.AddressOf(c,'r'))
  tree.GetEntry(entry)
  return c.r

if Method=="HybridNew":
  EXPfiles=[]
  EXPmasses = OBSmasses[:]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon):    # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.quant0.500.root"%m))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.quant0.500.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()
  
elif Method=="Asymptotic" or Method=="AsymptoticNew" or Method=="MaxLikelihoodFit":
  EXPfiles=[]
  EXPmasses = OBSmasses[:]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.root"%(m+epsilon)))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()

else:
  EXPfiles=[]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.root"%(m+epsilon)))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()

# Get the observed limits - Currently only does up to 1 decimal mass points
OBSfiles = []
if not options.expectedOnly:
  for m in OBSmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      OBSfiles.append(ROOT.TFile(OBSName+".mH%d.root"%(m+epsilon)))
    else:
      OBSfiles.append(ROOT.TFile(OBSName+".mH%.1f.root"%m))
    if options.verbose: print "observed MH - ", m, "File - ", OBSfiles[-1].GetName()

  if Method == "Asymptotic" or Method =="AsymptoticNew" :  obs = [getOBSERVED(O,5) for O in OBSfiles] # observed is last entry in these files
  else: obs = [getOBSERVED(O) for O in OBSfiles]

else:
  obs = [0 for O in OBSmasses]
  OBSfiles = obs[:]

# -------------------------------------------------------------------------------------------------------------------------------------------
# Set-up the GRAPHS

graph68  = ROOT.TGraphAsymmErrors()
graph95  = ROOT.TGraphAsymmErrors()
graphMed = ROOT.TGraphAsymmErrors()
graphObs = ROOT.TGraphAsymmErrors()
graphOne = ROOT.TGraphAsymmErrors()
dummyGraph= ROOT.TGraphAsymmErrors()

graph68up = ROOT.TGraphErrors()
graph68dn = ROOT.TGraphErrors()
graph95up = ROOT.TGraphErrors()
graph95dn = ROOT.TGraphErrors()
graphmede = ROOT.TGraphErrors()

graph68.SetLineColor(1)
graph95.SetLineColor(1)
graph68.SetLineStyle(2)
graph95.SetLineStyle(2)
graph68.SetLineWidth(2)
graph95.SetLineWidth(2)


MG = ROOT.TMultiGraph()

def MakeMlfPlot(MG):
    legend=ROOT.TLegend(0.55,0.19,0.89,0.3)
    legend.SetFillColor(10)
    legend.SetTextFont(42)
    #legend.SetTextSize(FONTSIZE)
    graph68.SetLineStyle(1)
    legend.AddEntry(graph68,"#pm 1#sigma Uncertainty","F")

    if options.square : c = ROOT.TCanvas("c","c",600,600)
    else :c = ROOT.TCanvas("c","c",800,600)

    dhist = ROOT.TH1F("dh","dh",100,MINMH,MAXMH)
    dhist.GetYaxis().SetTitleOffset(1.2)
    dhist.GetXaxis().SetTitleOffset(1.2)
    dhist.GetYaxis().SetTitleSize(0.04)
    dhist.GetXaxis().SetTitleSize(0.04)
    dhist.GetYaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetRangeUser(MINMH,MAXMH)
    dhist.GetYaxis().SetRangeUser(RANGEMU[0],RANGEMU[1])
    dhist.GetXaxis().SetTitle("m_{H} (GeV)")
    dhist.GetYaxis().SetTitle("Best fit #sigma/#sigma_{SM}")
    dhist.Draw("AXIS")

    MG.Draw("L3")

    # ------------------------------------------------------------------------
    # Additional Lines stored in --addline -----------------------------------
    for lineF in options.addline:

        # Parse the string, should be file.root:color:linestyle:legend entry    
        vals = lineF.split(":")
        ftmp = ROOT.TFile(vals[0])
        grext = ftmp.Get("observed")
        grext.SetLineColor(int(vals[1]))
        grext.SetLineStyle(int(vals[2]))
        grext.SetLineWidth(2)
        legend.AddEntry(grext,vals[3],"L")
        grext.Draw("same")
    # ------------------------------------------------------------------------
        

    c.Update()
    text = ROOT.TLatex()
    text.SetTextColor(ROOT.kRed)
    text.SetTextSize(FONTSIZE)
    text.SetTextFont(42)

    graphOne.Draw("L")
    c.SetGrid(not options.nogrid)
    if not options.nogrid: dhist.Draw("AXIGSAME")

    mytext= ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)
    mytext.SetTextFont(42)
    mytext.SetNDC()
    
    
    box = ROOT.TPave(0.19,0.17,0.4,0.3,2,"NDC")
    box.SetLineColor(1)
    box.SetFillColor(0)
    box.SetShadowColor(0)
    if not options.nobox: box.Draw()
    mytext.DrawLatex(0.2,0.26,"CMS Preliminary")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.2,0.25-(t+1)*0.04,"%s"%(lineT))
    legend.Draw()
    ROOT.gPad.RedrawAxis();
    
    if options.show:raw_input("Looks Ok?")
    c.SaveAs("maxlhplot.pdf")
    c.SaveAs("maxlhplot.png")

#-------------------------------------------------------------------------
def MakePvalPlot(MG):

    legend=ROOT.TLegend(0.55,0.17,0.89,0.35)
    legend.SetFillColor(10)
    legend.SetTextFont(42)
    #legend.SetTextSize(FONTSIZE)
    if not options.expectedOnly: legend.AddEntry(graphObs,"Observed","L")

    if options.square : c = ROOT.TCanvas("c","c",600,600)
    else :c = ROOT.TCanvas("c","c",800,600)

    dhist = ROOT.TH1F("dh","dh",100,MINMH,MAXMH)
    dhist.GetYaxis().SetTitleOffset(1.5)
    dhist.GetXaxis().SetTitleOffset(1.2)
    dhist.GetYaxis().SetTitleSize(0.04)
    dhist.GetXaxis().SetTitleSize(0.04)
    dhist.GetYaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetRangeUser(MINMH,MAXMH)
    dhist.GetYaxis().SetRangeUser(MINPV,MAXPV)
    dhist.GetXaxis().SetTitle("m_{H} (GeV)")
    dhist.GetYaxis().SetTitle("Local p-value")
    dhist.Draw("AXIS")

    MG.Draw("L3")

    # ------------------------------------------------------------------------
    # Additional Lines stored in --addline -----------------------------------
    for lineF in options.addline:

        # Parse the string, should be file.root:color:linestyle:legend entry    
        vals = lineF.split(":")
        ftmp = ROOT.TFile(vals[0])
        grext = ftmp.Get("observed")
        grext.SetLineColor(int(vals[1]))
        grext.SetLineStyle(int(vals[2]))
        grext.SetLineWidth(2)
        legend.AddEntry(grext,vals[3],"L")
        grext.Draw("same")
    # ------------------------------------------------------------------------
        

    c.Update()
    text = ROOT.TLatex()
    text.SetTextColor(ROOT.kRed)
    text.SetTextSize(FONTSIZE)
    text.SetTextFont(42)
    
        
    Vals=[ROOT.RooStats.SignificanceToPValue(L) for L in Lines]
    TLines=[ROOT.TLine(MINMH,V,MAXMH,V) for V in Vals]

    for j,TL in enumerate(TLines):
        TL.SetLineStyle(1)
        TL.SetLineColor(2)
        TL.SetLineWidth(1)
        TL.Draw("same")
        text.DrawLatex(MAXMH+0.2,Vals[j]*0.88,"%d #sigma"%Lines[j])

    c.SetGrid(not options.nogrid)
    c.SetLogy();
    if not options.nogrid: dhist.Draw("AXIGSAME")

    mytext= ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)
    mytext.SetTextFont(42)
    mytext.SetNDC()

    box = ROOT.TPave(0.19,0.17,0.4,0.3,2,"NDC")
    box.SetLineColor(1)
    box.SetFillColor(0)
    box.SetShadowColor(0)
    if not options.nobox: box.Draw()
    mytext.DrawLatex(0.2,0.26,"CMS Preliminary")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.2,0.25-(t+1)*0.04,"%s"%(lineT))
    legend.Draw()
    ROOT.gPad.RedrawAxis();
    
    if options.show:raw_input("Looks Ok?")
    c.SaveAs("pvaluesplot.pdf")
    c.SaveAs("pvaluesplot.png")
#-------------------------------------------------------------------------

def MakeLimitPlot(MG):

    leg=ROOT.TLegend(0.53,0.66,0.89,0.90)
    leg.SetFillColor(10)

    # Different entries for the different methods
    LegendEntry = ""
    if Method == "ProfileLikelihood": LegendEntry = "PL"
    if Method == "Bayesian": LegendEntry = "Bayesian"
    if Method == "HybridNew": LegendEntry = "CLs"
    if Method == "Asymptotic": LegendEntry = "CLs (Asymptotic)"
    if Method == "AsymptoticNew": LegendEntry = "CLs (Asymptotic)"

    if not options.expectedOnly: leg.AddEntry(graphObs,"Observed","L")
    if options.bayes and not options.expectedOnly: leg.AddEntry(bayesObs,"Observed Bayesian Limit","L")
    leg.AddEntry(graph68,"Expected #pm 1#sigma","FL")
    leg.AddEntry(graph95,"Expected #pm 2#sigma","FL")

    leg.SetTextFont(42)
    #leg.SetTextSize(FONTSIZE)

    if options.square : C = ROOT.TCanvas("c","c",600,600)
    else: C = ROOT.TCanvas("c","c",700,600)

    C.SetGrid(not options.nogrid)
    dummyHist = ROOT.TH1D("dummy","",1,min(OBSmasses)-OFFSETLOW,max(OBSmasses)+OFFSETHIGH)
    dummyHist.SetTitleSize(0.04,"XY")
    dummyHist.Draw("AXIS")
    MG.Draw("L3")
    if not options.nogrid: dummyHist.Draw("AXIGSAME")

    dummyHist.GetXaxis().SetTitle("m_{H} (GeV)")
    dummyHist.GetXaxis().SetRangeUser(min(OBSmasses)-OFFSETLOW,max(OBSmasses)+OFFSETHIGH)
    if options.doRatio:
     dummyHist.GetYaxis().SetRangeUser(RANGEYRAT[0],RANGEYRAT[1])
     dummyHist.GetYaxis().SetNdivisions(5,int("%d"%(RANGEYRAT[1]-RANGEYRAT[0])),0)
     dummyHist.GetYaxis().SetNdivisions(510)
     dummyHist.GetYaxis().SetTitle("\sigma(H#rightarrow #gamma #gamma)_{95%%CL} / \sigma(H#rightarrow #gamma #gamma)_{%s}"%extraString)
    else: 
     dummyHist.GetYaxis().SetRangeUser(RANGEYABS[0],RANGEYABS[1])
     dummyHist.GetYaxis().SetNdivisions(5,int("%d"%(RANGEYABS[1]-RANGEYABS[0])),0)
     dummyHist.GetYaxis().SetTitle("\sigma #times BR(H#rightarrow #gamma #gamma)_{95%CL} (pb)")

    dummyHist.GetYaxis().SetTitleOffset(1.3)
    dummyHist.GetXaxis().SetTitleOffset(1.25)

    MG.SetTitle("")
    mytext = ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)

    mytext.SetNDC()
    mytext.SetTextFont(42)
    mytext.SetTextSize(FONTSIZE)
    
    box = ROOT.TPave(0.19,0.76,0.42,0.89,2,"NDC")
    box.SetLineColor(1)
    box.SetFillColor(0)
    box.SetShadowColor(0)
    if not options.nobox: box.Draw()
    mytext.DrawLatex(0.2,0.85,"CMS Preliminary")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.2,0.84-(t+1)*(0.04),"%s"%lineT)
  
    leg.Draw()
    ROOT.gPad.RedrawAxis();
    if options.show:raw_input("Looks Ok?")

    #Make a bunch of extensions to the plots
    outputname="limit"
    if options.doSmooth: outputname+="_smooth"
    outputname+="_"+Method
    if options.doRatio: outputname+="_ratio"
    if options.append!="": outputname+="_"+options.append
    types=[".pdf",".png",".gif",".eps",".ps"]
    for type in types: C.SaveAs(outputname+type)

#-------------------------------------------------------------------------


#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#EXPECTED + Bands
for i,mass,f in zip(range(len(EXPfiles)),EXPmasses,EXPfiles):
  if options.pval: continue
  sm = 1.
  median = array.array('d',[0])
  up68   = array.array('d',[0])
  dn68   = array.array('d',[0])
  up95   = array.array('d',[0])
  dn95   = array.array('d',[0])

  if not options.doRatio: sm = signalNormalizer.GetBR(mass)*signalNormalizer.GetXsection(mass) 
  if Method == "Asymptotic" or Method=="AsymptoticNew":   
      median[0] = getOBSERVED(f,2)
      up95[0]   = getOBSERVED(f,4)
      dn95[0]   = getOBSERVED(f,0)
      up68[0]   = getOBSERVED(f,3)
      dn68[0]   = getOBSERVED(f,1)

  elif Method=="MaxLikelihoodFit":
      median[0] = getOBSERVED(f,0)
      up95[0]   = median[0]
      dn95[0]   = median[0]
      up68[0]   = getOBSERVED(f,2)
      dn68[0]   = getOBSERVED(f,1)

  else:
    tree = f.Get("limit")
    medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95)

  graph68.SetPoint(i,float(mass),median[0]*sm)
  graph95.SetPoint(i,float(mass),median[0]*sm)
  graphMed.SetPoint(i,float(mass),median[0]*sm)
  graphOne.SetPoint(i,float(mass),1.*sm)
  
  if Method == "HybridNew":

      up95[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.975.root"))
      dn95[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.025.root"))
      up68[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.840.root"))
      dn68[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.160.root"))

  
  diff95_up = abs(median[0] - up95[0])*sm
  diff95_dn = abs(median[0] - dn95[0])*sm
  diff68_up = abs(median[0] - up68[0])*sm
  diff68_dn = abs(median[0] - dn68[0])*sm
  
  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)
  graphMed.SetPointError(i,0,0,0,0)
  graphOne.SetPointError(i,0,0,0,0)

  if options.doSmooth:  # Always fit the absolute not the ratio
    sm=1.   

    graphmede.SetPoint(i,float(mass),median[0]*sm)
    graph68up.SetPoint(i,float(mass),up68[0]*sm)
    graph68dn.SetPoint(i,float(mass),dn68[0]*sm)
    graph95up.SetPoint(i,float(mass),up95[0]*sm)
    graph95dn.SetPoint(i,float(mass),dn95[0]*sm)

# Smooth the Bands set -doSmooth
# Since i always fitted to the Absolute, need to see if i want the Ratio instead
if options.doSmooth:
 fitstring = "[0] + [1]*x*x + [2]*x*x*x +[3]*x*x*x*x + [4]*x"
 medfunc = ROOT.TF1("medfunc",fitstring,109.,150.);
 up68func = ROOT.TF1("up68func",fitstring,109.,150.);
 dn68func = ROOT.TF1("dn68func",fitstring,109.,150.);
 up95func = ROOT.TF1("up95func",fitstring,109.,150.);
 dn95func = ROOT.TF1("dn95func",fitstring,109.,150.);

 graphmede.Fit(medfunc,"R,M,EX0","Q")
 graph68up.Fit(up68func,"R,M,EX0","Q")
 graph68dn.Fit(dn68func,"R,M,EX0","Q")
 graph95up.Fit(up95func,"R,M,EX0","Q")
 graph95dn.Fit(dn95func,"R,M,EX0","Q")
 
 newCanvas = ROOT.TCanvas()
 graphmede.SetMarkerSize(0.8)
 graphmede.SetMarkerStyle(21)
 graphmede.Draw("A")
 newCanvas.SaveAs("smoothTest.pdf")

 for i,mass in zip(range(len(EXPmasses)),EXPmasses):

  sm=1.0
  if not options.doRatio:
   sm = signalNormalizer.GetBR(mass)*signalNormalizer.GetXsection(mass)

  mediansmooth = medfunc.Eval(mass)

  graphMed.SetPoint(i,mass,mediansmooth*sm)
  graph68.SetPoint(i,mass,mediansmooth*sm)
  graph95.SetPoint(i,mass,mediansmooth*sm)

  diff95_up = abs(mediansmooth - up95func.Eval(mass))*sm
  diff95_dn = abs(mediansmooth - dn95func.Eval(mass))*sm
  diff68_up = abs(mediansmooth - up68func.Eval(mass))*sm
  diff68_dn = abs(mediansmooth - dn68func.Eval(mass))*sm

  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)

# OBSERVED -------- easy as that !
for i,mass in zip(range(len(OBSfiles)),OBSmasses):

    sm = 1.;
    if obs[i] ==-1: continue
    if not options.doRatio: sm = signalNormalizer.GetBR(M)*signalNormalizer.GetXsection(M)
    graphObs.SetPoint(i,float(mass),obs[i]*sm)
    graphObs.SetPointError(i,0,0,0,0)

# Finally setup the graphs and plot
graph95.SetFillColor(FILLCOLOR_95)
graph95.SetFillStyle(FILLSTYLE)
graph68.SetFillColor(FILLCOLOR_68)
graph68.SetFillStyle(FILLSTYLE)
graphMed.SetLineStyle(2)
graphMed.SetLineColor(2)
graphMed.SetMarkerColor(2)
graphMed.SetLineWidth(3)
graphObs.SetLineWidth(3)

if options.bayes:
 bayesObs.SetLineWidth(3)
 bayesObs.SetLineColor(4)
 bayesObs.SetMarkerColor(4)
 bayesObs.SetLineStyle(7)

graphOne.SetLineWidth(3)
graphOne.SetLineColor(ROOT.kRed)
graphOne.SetMarkerColor(ROOT.kRed)
graphObs.SetMarkerStyle(20)
graphObs.SetMarkerSize(2.0)
graphObs.SetLineColor(1)

graphMed.SetLineStyle(2)
graphMed.SetLineColor(ROOT.kBlack)
if not options.pval:MG.Add(graph95)
if not options.pval:MG.Add(graph68)
if not options.pval:MG.Add(graphMed)

if not options.expectedOnly:
  MG.Add(graphObs)
  if options.bayes:
   MG.Add(bayesObs)

if not options.pval: MG.Add(graphOne)

# Plot -------------------------------------
if options.pval: MakePvalPlot(MG)
elif Method=="MaxLikelihoodFit":  MakeMlfPlot(MG)
else:MakeLimitPlot(MG)
# ------------------------------------------
if options.outputLimits:
  print "Writing Limits To ROOT file --> ",options.outputLimits
  OUTTgraphs = ROOT.TFile(options.outputLimits,"RECREATE")
  graphObs.SetName("observed")
  graphObs.Write()

  if not options.pval:
   graphMed.SetName("median")
   graphMed.Write()
   graph68.SetName("sig1")
   graph68.Write()
   graph95.SetName("sig2")
   graph95.Write()

  OUTTgraphs.Write()

