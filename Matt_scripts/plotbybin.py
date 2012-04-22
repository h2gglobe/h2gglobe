import sys,os
import ROOT
import numpy
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i","--input",dest="tfileName")
parser.add_option("-m","--binmass",dest="mass")
(options,args)=parser.parse_args()

lumistring = "5.09 fb^{-1}"
sigscale = 5.
plotOutDir="plots"

def checkBins(h1,h2):
  if h1.GetNbinsX()!=h2.GetNbinsX():
    sys.exit(h1.GetName()+" and " +h2.GetName()+" have different number of bins")

def checkMass(m):
  if m!=110 or m!=115 or m!=120 or m!=125 or m!=130 or m!=135 or m!=140 or m!=150:
    sys.exit("%d is not a valid mass"%m)

def getMin(hists):
  smallest=10000000
  for h in hists:
    low = h.GetMinimum()
    if low<smallest:
      smallest=low
  if smallest-2*h.GetBinError(h.GetMinimumBin())<0:
    return 0
  else:
    return smallest-2*h.GetBinError(h.GetMinimumBin())

def getMax(hists):
  biggest=0.
  for h in hists:
    high = h.GetMaximum()
    if high>biggest:
      biggest=high
  return biggest+4*h.GetBinError(h.GetMaximumBin())

tFile = ROOT.TFile(options.tfileName,"UPDATE")

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

if not os.path.isdir('plots'):
  os.makedirs('plots')
binM=int(options.mass)

temp = tFile.Get("data_%d.0"%binM)
nbins = temp.GetNbinsX()

if binM==110:
  mlow  = 0.0
  mhigh = 2.5
elif binM==140:
  mlow  = -2.5
  mhigh = 5.0
elif binM==150:
  mlow  = -5.0
  mhigh = 0.1
else:
  mlow  = -2.5
  mhigh = 2.5

nMH = len(numpy.arange(mlow,mhigh,0.5))

h_data = [ROOT.TH1F("data_binM%d_bin%d"%(binM,b+1),"",nMH,0,nMH) for b in range(nbins) ]
h_bkg = [ROOT.TH1F("bkg_binM%d_bin%d"%(binM,b+1),"",nMH,0,nMH) for b in range(nbins) ]
h_bkg1sig = [ROOT.TH1F("bkg1sig_binM%d_bin%d"%(binM,b+1),"",nMH,0,nMH) for b in range(nbins) ]
h_bkg2sig = [ROOT.TH1F("bkg2sig_binM%d_bin%d"%(binM,b+1),"",nMH,0,nMH) for b in range(nbins) ]
h_bkgML = [ROOT.TH1F("bkgML_binM%d_bin%d"%(binM,b+1),"",nMH,0,nMH) for b in range(nbins) ]
h_sbML = [ROOT.TH1F("sbML_binM%d_bin%d"%(binM,b+1),"",nMH,0,nMH) for b in range(nbins) ]
h_sig = [ROOT.TH1F("sig_binM%d_bin%d"%(binM,b+1),"",nMH,0,nMH) for b in range(nbins) ]

point=0
for mdiff in numpy.arange(mlow,mhigh,0.5):
  mH = binM+mdiff
  if mH<110 or mH>150: continue
  
  dataH = tFile.Get("data_%3.1f"%mH)
  bkgH = tFile.Get("bkg_%3.1f"%mH)
  bkgH1sig = tFile.Get("bkg1_%3.1f"%mH)
  bkgH2sig = tFile.Get("bkg2_%3.1f"%mH)
  bkgmlH = tFile.Get("bkgml_%3.1f"%mH)
  sbmlH = tFile.Get("sbml_%3.1f"%mH)
  sigH = tFile.Get("sig_%3.1f"%mH)

  checkBins(dataH,bkgH)
  checkBins(dataH,bkgH1sig)
  checkBins(dataH,bkgH2sig)
  checkBins(dataH,bkgmlH)
  checkBins(dataH,sbmlH)
  checkBins(dataH,sigH)

  for b in range(dataH.GetNbinsX()):
    
    h_data[b].SetBinContent(point+1,dataH.GetBinContent(b+1))
    h_data[b].SetBinError(point+1,dataH.GetBinError(b+1))
    h_data[b].GetXaxis().SetBinLabel(point+1,"%3.1f"%mH)
    h_bkg[b].SetBinContent(point+1,bkgH.GetBinContent(b+1))
    h_bkg[b].SetBinError(point+1,bkgH.GetBinError(b+1))
    h_bkg[b].GetXaxis().SetBinLabel(point+1,"%3.1f"%mH)
    h_bkg1sig[b].SetBinContent(point+1,bkgH1sig.GetBinContent(b+1))
    h_bkg1sig[b].SetBinError(point+1,bkgH1sig.GetBinError(b+1))
    h_bkg1sig[b].GetXaxis().SetBinLabel(point+1,"%3.1f"%mH)
    h_bkg2sig[b].SetBinContent(point+1,bkgH2sig.GetBinContent(b+1))
    h_bkg2sig[b].SetBinError(point+1,bkgH2sig.GetBinError(b+1))
    h_bkg2sig[b].GetXaxis().SetBinLabel(point+1,"%3.1f"%mH)
    h_bkgML[b].SetBinContent(point+1,bkgmlH.GetBinContent(b+1))
    h_bkgML[b].SetBinError(point+1,bkgmlH.GetBinError(b+1))
    h_bkgML[b].GetXaxis().SetBinLabel(point+1,"%3.1f"%mH)
    h_sbML[b].SetBinContent(point+1,sbmlH.GetBinContent(b+1))
    h_sbML[b].SetBinError(point+1,sbmlH.GetBinError(b+1))
    h_sbML[b].GetXaxis().SetBinLabel(point+1,"%3.1f"%mH)
    h_sig[b].SetBinContent(point+1,sigH.GetBinContent(b+1))
    h_sig[b].SetBinError(point+1,sigH.GetBinError(b+1))
    h_sig[b].GetXaxis().SetBinLabel(point+1,"%3.1f"%mH)

  point=point+1

for b in range(nbins):

  h_data[b].SetMarkerStyle(20)
  h_bkg[b].SetLineColor(1); h_bkg[b].SetLineWidth(2)
  h_bkg1sig[b].SetLineColor(1); h_bkg1sig[b].SetLineWidth(2); h_bkg1sig[b].SetFillStyle(1001); h_bkg1sig[b].SetFillColor(3)
  h_bkg2sig[b].SetLineColor(1); h_bkg2sig[b].SetLineWidth(2); h_bkg2sig[b].SetFillStyle(1001); h_bkg2sig[b].SetFillColor(5)
  h_bkgML[b].SetLineColor(4); h_bkgML[b].SetLineWidth(3)
  h_sbML[b].SetLineColor(2); h_sbML[b].SetLineWidth(3)
  h_sig[b].SetLineColor(2); h_sig[b].SetLineStyle(2); h_sig[b].SetLineWidth(2)

  sb = ROOT.TCanvas("sb","sb",1600,1600)
  up = ROOT.TPad("u","u",0.01,0.3,0.99,0.99)
  #up.SetLogy()
  #dp.SetLogy()
  dp = ROOT.TPad("d","d",0.01,0.01,0.99,0.3)
  up.SetNumber(1)
  dp.SetNumber(2)
  up.Draw()
  dp.Draw()
  sb.cd(1)
  min=getMin([h_data[b],h_bkg[b],h_bkg1sig[b],h_bkg2sig[b],h_bkgML[b],h_sbML[b]])
  max=getMax([h_data[b],h_bkg[b],h_bkg1sig[b],h_bkg2sig[b],h_bkgML[b],h_sbML[b]])

  leg = ROOT.TLegend(0.6,0.5,0.89,0.89)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetLineColor(0)
  leg.SetLineStyle(0)
  leg.AddEntry(h_data[b],"Data","lep")
  leg.AddEntry(h_sig[b],"Higgs Signal (x5)","l")
  leg.AddEntry(h_bkg[b],"Background (datacard)","l")
  leg.AddEntry(h_bkg1sig[b],"Bkg #pm 1#sigma (datacard)","f")
  leg.AddEntry(h_bkg2sig[b],"Bkg #pm 2#sigma (datacard)","f")
  leg.AddEntry(h_bkgML[b],"Background (ML fit)","l")
  leg.AddEntry(h_sbML[b],"S+B (ML fit)","l")

  if b<5:
    min = getMin([h_data[b],h_bkg2sig[b],h_bkg1sig[b],h_bkg[b],h_bkgML[b],h_sbML[b]])
    max = getMax([h_data[b],h_bkg2sig[b],h_bkg1sig[b],h_bkg[b],h_bkgML[b],h_sbML[b]])
  else: 
    min = getMin([h_data[b],h_bkg2sig[b],h_bkg1sig[b],h_bkg[b],h_bkgML[b],h_sbML[b],h_sig[b]])
    max = getMax([h_data[b],h_bkg2sig[b],h_bkg1sig[b],h_bkg[b],h_bkgML[b],h_sbML[b],h_sig[b]])
  
  #h_data[b].SetTitle("Bin %d (binning mass %d) "%(b+1,binM))
  h_data[b].GetXaxis().SetTitle("m_{H} (GeV)")
  h_data[b].GetYaxis().SetRangeUser(min,max)
  h_data[b].Draw("9le1p")
  h_bkg2sig[b].Draw("9sameE2")
  h_bkg1sig[b].Draw("9sameE2")
  h_bkg[b].Draw("9samehist")
  h_bkgML[b].Draw("9samehist")
  h_sbML[b].Draw("9samehist")
  h_sig[b].Draw("9samehist")
  h_data[b].Draw("9samele1p") 
  leg.Draw("same")
  mytext = ROOT.TLatex();mytext.SetTextSize(0.03);mytext.SetNDC();mytext.DrawLatex(0.1,0.94,"CMS preliminary,  #sqrt{s} = 7 TeV ");mytext.SetTextSize(0.04)
  mytext.DrawLatex(0.5,0.94,"#int L = %s"%(lumistring))
  if b+1==nbins:
    mytext.DrawLatex(0.15,0.8,"Bin VBF")
    up.SaveAs("plots/plotbybin_m%d_bVBF.pdf"%(binM))
    up.SaveAs("plots/plotbybin_m%d_bVBF.png"%(binM))
  else:
    mytext.DrawLatex(0.15,0.8,"Bin %d"%(b+1))
    up.SaveAs("plots/plotbybin_m%d_b%d.pdf"%(binM,b+1))
    up.SaveAs("plots/plotbybin_m%d_b%d.png"%(binM,b+1))

  sb.cd(2)

  h_sig[b].Sumw2()
  sob=h_sig[b].Clone(); sobml=h_sig[b].Clone(); sobsbml=h_sig[b].Clone()
  sob.Divide(h_bkg[b])
  sobml.Divide(h_bkgML[b])
  sobsbml.Divide(h_sbML[b])
  sob.GetXaxis().SetLabelSize(0.1)
  sob.GetYaxis().SetLabelSize(0.08)
  sob.SetLineStyle(1)
  sob.SetLineColor(1)
  sobml.SetLineStyle(1)
  sobml.SetLineColor(2)
  sobsbml.SetLineStyle(1)
  sobsbml.SetLineColor(4)
  leg3 = ROOT.TLegend(0.11,0.5,0.4,0.88);leg3.SetFillColor(0);leg3.SetBorderSize(0)
  leg3.AddEntry(sob,"S/B (x5 SM)","lep")
  leg3.AddEntry(sobml,"S/B_{ML} (x5 SM)","lep")
  leg3.AddEntry(sobsbml,"S/S+B_{ML} (x5 SM)","lep")
  min=getMin([sob,sobml,sobsbml])
  max=getMax([sob,sobml,sobsbml])
  sob.GetYaxis().SetRangeUser(min,max)
  sob.Draw("9ep")
  sobml.Draw("9epsame")
  sobsbml.Draw("9epsame")
  leg3.Draw("same")
  if b+1==nbins:
    sb.SaveAs("plots/frac_bin_m%d_bVBF.pdf"%(binM))
    sb.SaveAs("plots/frac_bin_m%d_bVBF.png"%(binM))
  else:
    sb.SaveAs("plots/frac_bin_m%d_b%d.pdf"%(binM,b+1))
    sb.SaveAs("plots/frac_bin_m%d_b%d.png"%(binM,b+1))

  h_data[b].Write()
  h_bkg[b].Write()
  h_bkg1sig[b].Write()
  h_bkg2sig[b].Write()
  h_bkgML[b].Write()
  h_sbML[b].Write()
  h_sig[b].Write()
  
tFile.Close()
    

