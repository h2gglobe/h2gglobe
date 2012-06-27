# Simple script to make Effiency X Acceptance plot from Binned Baseline/Massfac analysis
# run with python makeEffAcc.py CMS-HGG.root
import ROOT
import sys
import numpy

#### ROOT.gROOT.ProcessLine(".L Normalization_7TeV.C++")
### ROOT.gROOT.ProcessLine(".L Normalization_8TeV.C++")
### from ROOT import GetBR
### from ROOT import GetXsection
### GetProcXsection = GetXsection

ROOT.gSystem.Load("../libLoopAll.so")
from ROOT import Normalization_8TeV
norm = Normalization_8TeV()
GetBR = lambda x : norm.GetBR(float(x))
GetXsection = lambda x : norm.GetXsection(float(x))
GetProcXsection = lambda x,y : norm.GetXsection(x,y)

# ROOT Setup
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(0)

# Global Setup, Modify with each Reload
##### NCAT = 5
##### lumi = 5089
#lumi=3770
#systematics = ["vtxEff","idEff","E_scale","E_res","triggerEff","regSig","phoIdMva"] # These are the main contributions to eff*Acc
systematics = ["vtxEff","idEff","E_scale","E_res","triggerEff"] # These are the main contributions to eff*Acc
## Masses = range(110,152,2) 
Masses = range(110,151,1) 
# -------------------------------------------------------------

f = ROOT.TFile(sys.argv[1])
NCAT = 6
# Get The lumi from the workspace!
ws = f.Get("cms_hgg_workspace")
lRRV = ws.var("IntLumi")
lumi = lRRV.getVal()

# Some helpful output
print "File - ", sys.argv[1]
printLine = "Data:      "
Sum = 0
for i in range(NCAT):
  h = f.Get("th1f_data_mass_cat%d"%i)
  print "%d   %4.0f    %4.0f" % (i, h.Integral(1,160), h.Integral(21,100) )
  Sum+=h.Integral()
  printLine+="%3.0f"%(h.Integral())+" "
printLine+="tot=%d"%Sum

print printLine

efficiency=ROOT.TGraphAsymmErrors()
central=ROOT.TGraphAsymmErrors()
efficiencyup=ROOT.TGraphAsymmErrors()
efficiencydn=ROOT.TGraphAsymmErrors()
centralsmooth=ROOT.TGraphAsymmErrors()

fitstring = "[0] + [1]*x + [2]*x*x"
cenfunc = ROOT.TF1("cenfunc",fitstring,109.75,140.25)
upfunc = ROOT.TF1("upfunc",fitstring,109.75,140.25)
dnfunc = ROOT.TF1("dnfunc",fitstring,109.75,140.25)


for point,M in enumerate(Masses):
  printLine = "Sigl M%d: "%M
  Sum = 0
  for i in range(NCAT):
    if int(M)==M:
     h =   f.Get("th1f_sig_ggh_mass_m%d_cat%d"%(int(M),i))
     hvb = f.Get("th1f_sig_vbf_mass_m%d_cat%d"%(int(M),i))
     hvh = f.Get("th1f_sig_wzh_mass_m%d_cat%d"%(int(M),i))
     htt = f.Get("th1f_sig_tth_mass_m%d_cat%d"%(int(M),i))

     ggh = h.Integral()
     vbf = hvb.Integral()
     tth = htt.Integral()
     wzh = hvh.Integral()
     
     h.Add(hvb)
     h.Add(hvh)
     h.Add(htt)

     print "%d %3.1f   %.5f    %.5f    %.5f     %.5f" % ( i, M, 100*ggh/(GetBR(M)*( GetProcXsection(M,"ggh")*(0.975))*lumi)
							,100*vbf/(GetBR(M)*(GetProcXsection(M,"vbf"))*lumi)
							,100*wzh/(GetBR(M)*(GetProcXsection(M,"wzh"))*lumi)
							,100*tth/(GetBR(M)*(GetProcXsection(M,"tth"))*lumi)
							)
     
    else:
     h = f.Get("th1f_sig_mass_m%.1f_cat%d"%(M,i))
    Sum+=h.Integral()
    printLine+="%3.5f"%h.Integral()+" "
  printLine+="tot=%3.5f"%Sum
  
  sm =GetBR(M)*( GetProcXsection(M,"ggh")*0.975 + GetProcXsection(M,"vbf") + GetProcXsection(M,"wzh") + GetProcXsection(M,"tth") )
  effAcc = 100*Sum/(sm*lumi) # calculate Efficiency at mH
  centralsmooth.SetPoint(point,M,effAcc)
  central.SetPoint(point,M,effAcc)
  efficiency.SetPoint(point,M,effAcc)
  delUp = 0
  delDown = 0
  for s in systematics:
   syssumup=0
   syssumdn=0
   for i in range(NCAT):
    if int(M)==M:
     hup =   f.Get("th1f_sig_ggh_mass_m%d_cat%d_%sUp01_sigma"%(int(M),i,s))
     hupvb = f.Get("th1f_sig_vbf_mass_m%d_cat%d_%sUp01_sigma"%(int(M),i,s))
     hupvh = f.Get("th1f_sig_wzh_mass_m%d_cat%d_%sUp01_sigma"%(int(M),i,s))
     huptt = f.Get("th1f_sig_tth_mass_m%d_cat%d_%sUp01_sigma"%(int(M),i,s))
     hup.Add(hupvb)
     hup.Add(hupvh)
     hup.Add(huptt)

     hdn =   f.Get("th1f_sig_ggh_mass_m%d_cat%d_%sDown01_sigma"%(int(M),i,s))
     hdnvb = f.Get("th1f_sig_vbf_mass_m%d_cat%d_%sDown01_sigma"%(int(M),i,s))
     hdnvh = f.Get("th1f_sig_wzh_mass_m%d_cat%d_%sDown01_sigma"%(int(M),i,s))
     hdntt = f.Get("th1f_sig_tth_mass_m%d_cat%d_%sDown01_sigma"%(int(M),i,s))
     hdn.Add(hdnvb)
     hdn.Add(hdnvh)
     hdn.Add(hdntt)

    else:
     hup = f.Get("th1f_sig_mass_m%.1f_cat%d_%sUp01_sigma"%(M,i,s))
     hdn = f.Get("th1f_sig_mass_m%.1f_cat%d_%sDown01_sigma"%(M,i,s))
    syssumup+=hup.Integral()
    syssumdn+=hdn.Integral()

   # We make 3-sigma templates so need to scale back by 1/3
   delUp+=abs(syssumup-Sum)/3
   delDown+=abs(syssumdn-Sum)/3
   
  delUp=100*(delUp**0.5)/(sm*lumi)
  delDown=100*(delDown**0.5)/(sm*lumi)
  efficiencyup.SetPoint(point,M,delUp)
  efficiencydn.SetPoint(point,M,delDown)
  centralsmooth.SetPointError(point,0,0,0,0)
  efficiency.SetPointError(point,0,0,delDown,delUp)

  print printLine

centralsmooth.Fit(cenfunc,"R,0,EX0","")
efficiencyup.Fit(upfunc,"R,0,EX0","")
efficiencydn.Fit(dnfunc,"R,0,EX0","")

for point,M in enumerate(Masses):
 central.SetPoint(point,M,cenfunc.Eval(M))
 efficiency.SetPoint(point,M,cenfunc.Eval(M))

leg=ROOT.TLegend(0.46,0.16,0.79,0.39)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.AddEntry(central,"Higgs Signal #varepsilon #times Acc","L")
leg.AddEntry(efficiency,"#pm 1 #sigma syst. error","F")

mytext = ROOT.TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()

listy = []

MG=ROOT.TMultiGraph()
can = ROOT.TCanvas()
efficiency.SetFillColor(9)
efficiency.SetLineWidth(2)
central.SetLineWidth(2)
MG.Add(efficiency)
MG.Add(central)
MG.Draw("AL3")
MG.GetXaxis().SetTitle("m_{H} GeV")
MG.GetYaxis().SetTitle("Efficiency #times Acceptance - %")
mytext.DrawLatex(0.15,0.8,"CMS Simulation")
can.Update()
leg.Draw("same")
print "Int Lumi from workspace ", lumi
raw_input("Looks OK?")

can.Update()
print "Saving plot as effAcc_vs_mass.pdf"
can.SaveAs("effAcc_vs_mass.pdf")
can.SaveAs("effAcc_vs_mass.png")

  

	
