import sys
import ROOT
from ROOT import TColor
ROOT.gROOT.SetBatch(1)
intlumi="3.06"

def Setup(gr,col,style,width,leg,name,MG):
	gr.SetLineColor(col)
	gr.SetLineWidth(width)
	#gr.SetLineColor(0)
	gr.SetLineStyle(style)
	gr.SetMarkerColor(col)
	gr.SetMarkerSize(0.8)
	gr.SetMarkerStyle(20)
	fitstring = "[0] + [1]*x*x + [2]*x*x*x +[3]*x*x*x*x + [4]*x"
	grSmooth = ROOT.TF1("grSmooth",fitstring,109.,150.)
	#gr.Fit(grSmooth,"R,M,EXO","Q")
	#grSmooth.SetLineWidth(2)
	#grSmooth.SetLineColor(col)
	MG.Add(gr)
	leg.AddEntry(gr,name,"L")
	return grSmooth

def fitSetup(gr,col,style,MG):
	gr.SetLineWidth(2)
	gr.SetLineColor(col)
	gr.SetLineStyle(style)
	gr.SetMarkerColor(col)
	gr.SetMarkerSize(0.8)
	gr.SetMarkerStyle(20)
	gr.SetFillColor(9)
	MG.Add(gr)

MG=ROOT.TMultiGraph()
fitMG=ROOT.TMultiGraph()

files = [
"Jun15UnblindOldWeights3fb/MassFac3fbOldWeights/exppval_comb.root",
"Jun15UnblindOldWeights3fb/MassFac3fbOldWeights/exppval_7TeV_5fb.root",
"Jun15UnblindOldWeights3fb/MassFac3fbOldWeights/exppval_8TeV_3fb.root"
]

fitfiles = [
#"/vols/cms02/nw709/hgg/src_cvs/dec14/CMSSW_4_2_8/src/HiggsAnalysis/CombinedLimit/appmva/bdtcount/limit/bf.root",
]

medfiles = [
]

onepoint = ROOT.TGraphErrors()
onepoint.SetPoint(0,124.0,0.0012012)
onepoint.SetPointError(0,0,0.0003465)
onepoint.SetMarkerColor(38)
onepoint.SetMarkerStyle(20)
onepoint.SetMarkerSize(1.5)

FILES = []
MEDFILES = []
GR=[]
fitFILES = []
fitGR=[]

colors=[1,ROOT.kBlue+2,ROOT.kRed]
styles=[9,9,9]
widths=[3,2,2]
colors_med=[38]
styles_med=[7]
widths_med=[2,2,2]

fitcolors=[1,2,3,4,6,7,9]
fitstyles=[1,1,1,1,1,1,1,1]

names=[
#"Local P-value (Bin "+str(int(bin)-1)+")"
#"1xSM Higgs \sqrt{s}=7TeV L=5.1fb^{-1}",
#"Basline",
#"MassFactorized MVA",
#"Sideband MVA"
"#splitline{Combined}{(1xSM Higgs Expected)}",
"\sqrt{s}=7TeV L=5.1fb^{-1}",
"\sqrt{s}=8TeV L=3.1fb^{-1}"
]

names_med=["1xSM Higgs Median Expected"]

# significance lines
#Lines = [1.,2.,3.,4.]
Lines = [1.,2.,3.]#,4.,5.,6.]
Vals=[ROOT.RooStats.SignificanceToPValue(L) for L in Lines]
TLines=[ROOT.TLine(110,V,150,V) for V in Vals]

fitLines = [1.]
fitTLines=[ROOT.TLine(110,fl,150,fl) for fl in fitLines]

legend=ROOT.TLegend(0.11,0.135,0.35,0.4)
legend.SetFillColor(0)
legend.SetBorderSize(0)

for i in range(len(medfiles)):
	MEDFILES.append(ROOT.TFile(medfiles[i]))
	GR.append(MEDFILES[i].Get("median"))
	Setup(GR[-1],colors_med[i],styles_med[i],widths_med[i],legend,names_med[i],MG)

GRsmooth=[]
for i in range(len(files)):
	FILES.append(ROOT.TFile(files[i]))
	ggr =FILES[i].Get("observed")
	print ggr.GetName() 
	print ROOT.RooStats.PValueToSignificance(ggr.Eval(124.0))
	print ggr.Eval(124.0)
	GR.append(FILES[i].Get("observed"))
	
	GRsmooth.append(Setup(GR[-1],colors[i],styles[i],widths[i],legend,names[i],MG))

#MG.Add(onepoint)
#legend.AddEntry(onepoint,"Observed p-value (Ensemble)","PEl")

for i in range(len(fitfiles)):
	fitFILES.append(ROOT.TFile(fitfiles[i]))
	fitGR.append(fitFILES[i].Get("bestfit"))
	fitSetup(fitGR[i],fitcolors[i],fitstyles[i],fitMG)

mytext= ROOT.TLatex()
mytext.SetTextSize(0.045)

ROOT.gROOT.SetStyle("Plain")
c = ROOT.TCanvas()

"""
ROOT.gROOT.SetStyle("Plain")
c = ROOT.TCanvas("c","c",1600,1800)
# Make 2 pads for the signal strength fit.
up = ROOT.TPad("u","u",0.01,0.7,0.99,0.99);
dp = ROOT.TPad("d","d",0.01,0.01,0.99,0.7);
up.SetNumber(1);
dp.SetNumber(2);
up.Draw();
dp.Draw();


c.cd(1)
fitMG.Draw("AL")#E4")
#fitMGLINE = ROOT.TGraph()
#for m in range(110,151): fitMGLINE.SetPoint(m-110,float(m),fitMG.Eval(float(m)))
fitMG.Draw("LP")
fitMG.GetXaxis().SetRangeUser(110,150)
fitMG.GetYaxis().SetRangeUser(0,3)
fitMG.GetXaxis().SetTitle("")
fitMG.GetYaxis().SetTitle("#mu")
fitMG.GetYaxis().SetTitleSize(0.15)
fitMG.GetYaxis().SetTitleOffset(0.12)
#fitMG.GetYaxis().SetLabelSize(0.1)
#fitMG.GetXaxis().SetLabelSize(0.1)
for j,TL in enumerate(fitTLines):
	TL.SetLineStyle(1)
	TL.SetLineColor(1)
	TL.SetLineWidth(2)
	TL.Draw("same")
#	text.DrawLatex(151,fitLines[j],"%dxSM"%Lines[j])
#

fitMG.GetYaxis().SetNdivisions(3)
fitMG.GetYaxis().SetLabelSize(0.11);
fitMG.GetXaxis().SetLabelSize(0.11);
up.SetGrid(True)
#dp.SetLogy()
#
c.cd(2)
"""
MG.Draw("AL")#LP")
MG.GetXaxis().SetRangeUser(110,150)
MG.GetYaxis().SetRangeUser(1.e-4,1.)
for i, grSmo in enumerate(GRsmooth):
	grSmo.SetLineColor(colors[i])
	#grSmo.Draw("same")

MG.GetXaxis().SetTitle("m_{H} (GeV)")
MG.GetYaxis().SetTitle("p-value")
MG.GetYaxis().SetTitleOffset(0.6)
MG.GetYaxis().SetTitleSize(0.05)
MG.GetXaxis().SetTitleSize(0.05)
MG.GetYaxis().SetLabelSize(0.05)
MG.GetXaxis().SetLabelSize(0.05)
#onepoint.Draw("PE2L")
text = ROOT.TLatex()

#text.SetNDC()
for j,TL in enumerate(TLines):
	TL.SetLineStyle(2)
	TL.SetLineColor(1)
	TL.SetLineWidth(1)
	TL.Draw("same")
	text.DrawLatex(150.5,Vals[j]*0.88,"%d #sigma"%Lines[j])
#mytext.SetFillColor(8)	
mytext.SetNDC()
#mytext.SetFillColor(0)
#mytext.SetLineColor(0)
legend.Draw()

Box = ROOT.TBox(137,0.0015,147,0.0055)
Box.SetFillColor(10)
Box.SetFillStyle(1001)
#Box.Draw()
#mytext
mytext.DrawLatex(0.55,0.30,"CMS preliminary")
mytext.DrawLatex(0.55,0.22,"#splitline{#sqrt{s} = 7 TeV L = 5.1 fb^{-1}}{#sqrt{s} = 8 TeV L = 3.1 fb^{-1}}")
#mytext.DrawLatex(0.55,0.22,"#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 3.1 fb^{-1}}")
c.SetGridx()
c.SetLogy()


#raw_input() 
#dp.SaveAs("pvals_bin"+bin+".pdf")
c.SaveAs(sys.argv[1]+"/pvalues.pdf")
# dp.SaveAs(sys.argv[1]+"/pvalues_nobf.pdf")
c.SaveAs(sys.argv[1]+"/pvalues.png")
c.SaveAs(sys.argv[1]+"/pvalues.C")
# dp.SaveAs(sys.argv[1]+"/pvalues_nobf.png")

