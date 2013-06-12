import sys
import os
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);

#input files
#mc   = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v6_mc_nvtxRw/histograms_phojetanalysis.root")
#data = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v6_data/histograms_phojetanalysis.root")
#mc   = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v7_mc_pt0to20/histograms_phojetanalysis.root")
#data = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v7_data_pt0to20/histograms_phojetanalysis.root")
mc   = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v7_mc_minpt20/histograms_phojetanalysis.root")
data = TFile("root://eoscms/eos/cms/store/cmst3/user/malberti/HIGGS/PhotonPlusJet/test_v7_data_minpt20/histograms_phojetanalysis.root")


# dir to save plots
outdir = sys.argv[1]
try:
    os.mkdir(outdir)
except:
    pass

#comparison of event and kinematic variables
evtvars=["nvtx","photon_pt","photon_eta","jet_pt","jet_eta","pt"]
for i in range(len(evtvars)):
    cname = evtvars[i]
    canvas = TCanvas(cname,cname,750,600)
    canvas.cd()
    hname = evtvars[i]+"_cat0_tot"
    hdata = data.Get(hname)
    hdata.SetMarkerStyle(20)
    hmc = mc.Get(hname)
    hmc.SetTitle(evtvars[i])
    hmc.DrawNormalized("histo",hdata.GetSumOfWeights())
    hdata.DrawNormalized("esame", hdata.GetSumOfWeights())
    canvas.SaveAs(outdir+"/"+cname+".png",".png")
    canvas.SaveAs(outdir+"/"+cname+".pdf",".pdf")
    canvas.SaveAs(outdir+"/"+cname+".C",".C")


raw_input("continue?")

#compare fraction of converted photons
vars=["photon_pt","photon_eta","pt"]
labels=["photon pT","photon eta","(gamma+jet) pT"]
nre =[10,1,10] 
for i in range(len(vars)):
    hdatad = data.Get(vars[i]+"_cat0_tot")
    hmcd = mc.Get(vars[i]+"_cat0_tot")
    hdatad.Rebin(nre[i])
    hmcd.Rebin(nre[i])
    
    cname = vars[i]+"_frac_goodconv"
    canvas = TCanvas(cname,cname,750,600)
    canvas.cd()
    hdatan = data.Get(vars[i]+"_conv_cat0_tot")    
    hdatan.Rebin(nre[i])
    hdatan.Divide(hdatan,hdatad)
    hdatan.SetMarkerStyle(20)
    hmcn = mc.Get(vars[i]+"_conv_cat0_tot")
    hmcn.Rebin(nre[i])
    hmcn.Divide(hmcn,hmcd)
    hmcn.SetXTitle(labels[i])
    hmcn.GetYaxis().SetRangeUser(0,1.0)
    hmcn.SetYTitle("fraction of good conv")
    hmcn.SetTitle("fraction of photons with good conv")
    hmcn.Draw("histo")
    hdatan.Draw("esame")
    print outdir+"/"+cname+".png"
    canvas.SaveAs(outdir+"/"+cname+".png",".png")
    canvas.SaveAs(outdir+"/"+cname+".pdf",".pdf")
    canvas.SaveAs(outdir+"/"+cname+".C",".C")
    
    cname2 = vars[i]+"_frac_tkcone"
    canvas2 = TCanvas(cname2,cname2,750,600)
    canvas2.cd()
    hdatan = data.Get(vars[i]+"_tkcone_cat0_tot")
    hdatan.Rebin(nre[i])
    hdatan.Divide(hdatan,hdatad)
    hdatan.SetMarkerStyle(20)
    hmcn = mc.Get(vars[i]+"_tkcone_cat0_tot")
    hmcn.Rebin(nre[i])
    hmcn.Divide(hmcn,hmcd)
    hmcn.GetYaxis().SetRangeUser(0,1.0)
    hmcn.SetXTitle(labels[i])
    hmcn.SetYTitle("fraction of photons with ntk(0.05) > 0")
    hmcn.SetTitle("fraction of photons with ntk(0.05) > 0 and no good conversions")
    hmcn.Draw("histo")
    hdatan.Draw("esame")
    print outdir+"/"+cname2+".png"
    canvas2.SaveAs(outdir+"/"+cname2+".png",".png")
    canvas2.SaveAs(outdir+"/"+cname2+".pdf",".pdf")
    canvas2.SaveAs(outdir+"/"+cname2+".C",".C")
    
raw_input("continue?")

#comparison of vertex related variables
colors = [kGreen+2, kRed+2]
fillstyle = [3002,3005]
leg = TLegend(0.65, 0.70, 0.89, 0.89)
leg.SetFillColor(0)
leg.SetBorderSize(1)
leg.SetLineColor(kWhite)
mytext = TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
legends = ["Jet Tagged Vertex","Non-Jet Tagged Vertex"]
variables = ["sumpt2","ptbal","ptasym","pulltoconv","nconv","sumpt2in","sumpt2out","nchin","vtxmva","evtmva"]
#variables = ["sumpt2","sumpt2in","sumpt2out","nchin"]
histnames       = []
histnamesConv   = []
histnamesTkCone = []
histnamesMc       = []
histnamesConvMc   = []
histnamesTkConeMc = []
for var in variables:
    histnames.append([var+"_rv_cat0_tot",var+"_wv_cat0_tot"])
    histnamesConv.append([var+"_rv_conv_cat0_tot",var+"_wv_conv_cat0_tot"])
    histnamesTkCone.append([var+"_rv_tkcone_cat0_tot",var+"_wv_tkcone_cat0_tot"])

    #histnames.append([var+"_rv_cat0_QCD_Pt80to170_8TeV",var+"_wv_cat0_QCD_Pt80to170_8TeV"])
    #histnamesConv.append([var+"_rv_conv_cat0_QCD_Pt80to170_8TeV",var+"_wv_conv_cat0_QCD_Pt80to170_8TeV"])
    #histnamesTkCone.append([var+"_rv_tkcone_cat0_QCD_Pt80to170_8TeV",var+"_wv_tkcone_cat0_QCD_Pt80to170_8TeV"])

    #histnamesMc.append([var+"_rv_cat0_QCD_Pt250to350_8TeV",var+"_wv_cat0_QCD_Pt250to350_8TeV"])
    #histnamesConvMc.append([var+"_rv_conv_cat0_QCD_Pt250to350_8TeV",var+"_wv_conv_cat0_QCD_Pt250to350_8TeV"])
    #histnamesTkConeMc.append([var+"_rv_tkcone_cat0_QCD_Pt250to350_8TeV",var+"_wv_tkcone_cat0_QCD_Pt250to350_8TeV"])

    #histnames.append([var+"_rv_cat0_G_Pt80to120_8TeV",var+"_wv_cat0_G_Pt80to120_8TeV"])
    #histnamesConv.append([var+"_rv_conv_cat0_G_Pt80to120_8TeV",var+"_wv_conv_cat0_G_Pt80to120_8TeV"])
    #histnamesTkCone.append([var+"_rv_tkcone_cat0_G_Pt80to120_8TeV",var+"_wv_tkcone_cat0_G_Pt80to120_8TeV"])
    #histnamesMc.append([var+"_rv_cat0_G_Pt470to800_8TeV",var+"_wv_cat0_G_Pt470to800_8TeV"])
    #histnamesConvMc.append([var+"_rv_conv_cat0_G_Pt470to800_8TeV",var+"_wv_conv_cat0_G_Pt470to800_8TeV"])
    #histnamesTkConeMc.append([var+"_rv_tkcone_cat0_G_Pt470to800_8TeV",var+"_wv_tkcone_cat0_G_Pt470to800_8TeV"])


#print histnames
    
histos = [histnames,histnamesConv,histnamesTkCone]
#histosMc = [histnamesMc,histnamesConvMc,histnamesTkConeMc]
histosMc = histos
suffix = ["","_conv","_tksincone"]

for k,h in enumerate(histos):
    for i in range(len(variables)):
        cname = variables[i]+suffix[k]
        canvas = TCanvas(cname,cname,750,600)
        if (variables[i]=="pulltoconv" or variables[i]=="vtxmva" or variables[i]=="evtmva"):
            canvas.SetLogy()
        canvas.cd()
        for j in range(len(h[i])):
        
            hdata = data.Get(h[i][j])
            if (hdata.Integral()>0):
                hdata.Scale(1/hdata.Integral())
            hdata.SetLineColor(colors[j])
            hdata.SetMarkerStyle(20)
            hdata.SetMarkerColor(colors[j])
            #hmc = mc.Get(h[i][j])
            hmc = mc.Get(histosMc[k][i][j])
            if (hmc.Integral()>0):
                hmc.Scale(1/hmc.Integral())
            hmc.SetLineColor(colors[j])
            hmc.SetFillColor(colors[j])
            hmc.SetFillStyle(fillstyle[j])
            
            if (k==0 and i==0):
                leg.AddEntry(hmc,legends[j],"F")
                
            if (j==0):
                hmc.SetTitle(variables[i])
                hmc.SetMinimum(0.0001)
                hmctemp = mc.Get(h[i][1])
                maxy = max(hmc.GetMaximum(), hmctemp.GetMaximum()/hmctemp.Integral()) 
                hmc.SetMaximum(maxy*1.5)
                hmc.SetYTitle("a.u.")
                if (variables[i] == "nchin"):
                    hmc.GetXaxis().SetRangeUser(0,6)
                if (variables[i] == "pulltoconv"):
                    hmc.GetXaxis().SetRangeUser(0,20)
                hmc.Draw("histo")
                hdata.Draw("esame")
            else:
                hmc.Draw("histosame")
                hdata.Draw("esame")
                         
        leg.Draw("same")
        canvas.SaveAs(outdir+"/"+cname+".png",".png")
        canvas.SaveAs(outdir+"/"+cname+".pdf",".pdf")
        canvas.SaveAs(outdir+"/"+cname+".C",".C")







            
raw_input("OK?")
