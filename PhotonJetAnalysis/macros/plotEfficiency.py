import sys
import os
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)
gStyle.SetOptTitle(0)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);
gStyle.SetErrorX(0.5)

#input files
mc   = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v6_mc_nvtxRw/histograms_phojetanalysis.root")
data = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v6_data/histograms_phojetanalysis.root")

# dir to save plots
outdir = sys.argv[1]
try:
    os.mkdir(outdir)
except:
    pass

leg = TLegend(0.65, 0.15, 0.89, 0.35)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
#leg.SetLineColor(kWhite)

nre = 5
cats=["","conv_","tkcone_"]
for cat in cats:
    cname = cat+"efficiencyVsPt"
    canvas = TCanvas(cname,cname,750,600)
    canvas.SetGridx()
    canvas.SetGridy()    
    canvas.cd()
    hdatad = data.Get("pt_"+cat+"cat0_tot")
    hdatan = data.Get("pt_rv_"+cat+"cat0_tot")
    hdatan.Rebin(nre)
    hdatad.Rebin(nre)    
    hdatan.Divide(hdatan,hdatad,1,1,"B")
    hdatan.SetMarkerStyle(20)
    hmcd = mc.Get("pt_"+cat+"cat0_tot")
    hmcn = mc.Get("pt_rv_"+cat+"cat0_tot")
    hmcn.Rebin(nre)
    hmcd.Rebin(nre)    
    hmcn.Divide(hmcn,hmcd,1,1,"B")
    hmcn.GetXaxis().SetRangeUser(0,250.)
    hmcn.GetYaxis().SetRangeUser(0.5,1.1)
    hmcn.SetYTitle("efficiency")
    hmcn.SetXTitle("#gamma+jet p^{T} (GeV)")
    hmcn.SetFillColor(kCyan-10)
    hmcn.SetLineColor(kBlue+2)
    hmcn.Draw("e1 p")
    hdatan.Draw("e1 same")
    if (cat==""):
        leg.AddEntry(hmcn,"MC","LP")
        leg.AddEntry(hdatan,"data","LP")
    leg.Draw("same")
    
    canvas.SaveAs(outdir+"/"+cname+".png",".png")
    canvas.SaveAs(outdir+"/"+cname+".pdf",".pdf")
    canvas.SaveAs(outdir+"/"+cname+".C",".C")
        
raw_input("OK?")
