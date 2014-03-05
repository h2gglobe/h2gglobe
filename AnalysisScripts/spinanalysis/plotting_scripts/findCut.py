import ROOT as r
import sys

r.gStyle.SetOptStat(0)

infile = r.TFile(sys.argv[1])

dataHiggsPt = infile.Get("dataHiggsPt")
bkgHiggsPt = infile.Get("bkgHiggsPt")
gghHiggsPt = infile.Get("gghHiggsPt")
vbfHiggsPt = infile.Get("vbfHiggsPt")
gggravHiggsPt = infile.Get("gggravHiggsPt")
qqgravHiggsPt = infile.Get("qqgravHiggsPt")

def calcPurity(hist):
  purHist = hist.Clone('%s_pur'%hist.GetName())
  for bin in range(1,hist.GetNbinsX()+1):
    purHist.SetBinContent(bin,hist.Integral(0,bin)/hist.Integral())
  return purHist

def calcRoc(shist,bhist):
  assert(shist.GetNbinsX()==bhist.GetNbinsX())
  assert(shist.GetBinLowEdge(1)==bhist.GetBinLowEdge(1))
  assert(shist.GetBinLowEdge(shist.GetNbinsX()+1)==bhist.GetBinLowEdge(bhist.GetNbinsX()+1))
  
  gr = r.TGraph()
  p=0
  for bin in range(1,shist.GetNbinsX()+1):
    gr.SetPoint(p,bhist.GetBinContent(bin),shist.GetBinContent(bin))
    p+=1
  gr.SetLineColor(shist.GetLineColor())
  return gr 

dataHiggsPtPur = calcPurity(dataHiggsPt)
bkgHiggsPtPur = calcPurity(bkgHiggsPt)
gghHiggsPtPur = calcPurity(gghHiggsPt)
vbfHiggsPtPur = calcPurity(vbfHiggsPt)
gggravHiggsPtPur = calcPurity(gggravHiggsPt)
qqgravHiggsPtPur = calcPurity(qqgravHiggsPt)

dataVbfHiggsPtRoc = calcRoc(dataHiggsPtPur,vbfHiggsPtPur)
bkgVbfHiggsPtRoc = calcRoc(bkgHiggsPtPur,vbfHiggsPtPur)
gghVbfHiggsPtRoc = calcRoc(gghHiggsPtPur,vbfHiggsPtPur)
gggravVbfHiggsPtRoc = calcRoc(gggravHiggsPtPur,vbfHiggsPtPur)
qqgravVbfHiggsPtRoc = calcRoc(qqgravHiggsPtPur,vbfHiggsPtPur)

bkgHiggsPtPur.SetLineColor(r.kOrange+1)
bkgHiggsPtPur.SetLineStyle(7)
bkgHiggsPtPur.SetFillStyle(0)

dataHiggsPtPur.SetLineWidth(2)
bkgHiggsPtPur.SetLineWidth(2)
gghHiggsPtPur.SetLineWidth(2)
vbfHiggsPtPur.SetLineWidth(2)
gggravHiggsPtPur.SetLineWidth(2)
qqgravHiggsPtPur.SetLineWidth(2)

dataVbfHiggsPtRoc.SetLineWidth(2)
bkgVbfHiggsPtRoc.SetLineWidth(2)
gghVbfHiggsPtRoc.SetLineWidth(2)
gggravVbfHiggsPtRoc.SetLineWidth(2)
qqgravVbfHiggsPtRoc.SetLineWidth(2)
canv = r.TCanvas()

leg = r.TLegend(0.7,0.11,0.89,0.5)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(dataHiggsPtPur,"Data","L")
leg.AddEntry(bkgHiggsPtPur,"Bkg","L")
leg.AddEntry(gghHiggsPtPur,"ggh","L")
leg.AddEntry(vbfHiggsPtPur,"vbf","L")
leg.AddEntry(gggravHiggsPtPur,"gg_grav","L")
leg.AddEntry(qqgravHiggsPtPur,"qq_grav","L")
dataHiggsPtPur.Draw()
bkgHiggsPtPur.Draw("Lsame")
gghHiggsPtPur.Draw("same")
vbfHiggsPtPur.Draw("same")
gggravHiggsPtPur.Draw("same")
qqgravHiggsPtPur.Draw("same")
leg.Draw()
raw_input('Ok?')

leg = r.TLegend(0.7,0.11,0.89,0.5)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(dataVbfHiggsPtRoc,"Data","L")
leg.AddEntry(bkgVbfHiggsPtRoc,"Bkg","L")
leg.AddEntry(gghVbfHiggsPtRoc,"ggh","L")
leg.AddEntry(gggravVbfHiggsPtRoc,"gg_grav","L")
leg.AddEntry(qqgravVbfHiggsPtRoc,"qq_grav","L")
dataVbfHiggsPtRoc.Draw("AL")
bkgVbfHiggsPtRoc.Draw("Lsame")
gghVbfHiggsPtRoc.Draw("Lsame")
gggravVbfHiggsPtRoc.Draw("Lsame")
qqgravVbfHiggsPtRoc.Draw("Lsame")
leg.Draw()
canv.Modified()
canv.Update()
raw_input('Ok?')

dataHiggsPt
bkgHiggsPt 
gghHiggsPt.Scale(1./gghHiggsPt.Integral())
vbfHiggsPt.Scale(1./vbfHiggsPt.Integral()) 
gggravHiggsPt.Scale(1./gggravHiggsPt.Integral())
qqgravHiggsPt.Scale(1./qqgravHiggsPt.Integral())

gghHiggsPt.Draw()
vbfHiggsPt.Draw("same")
gggravHiggsPt.Draw("same")
qqgravHiggsPt.Draw("same")
canv.Modified()
canv.Update()
raw_input('Ok?')

