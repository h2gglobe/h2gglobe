import sys
import ROOT
ROOT.gROOT.SetBatch(1)
intlumi="5.09"

def Setup(gr,col,style,leg,name,MG):
	gr.SetLineWidth(2)
	gr.SetLineColor(col)
	gr.SetLineStyle(style)
	gr.SetMarkerColor(col)
	gr.SetMarkerSize(0.8)
	gr.SetMarkerStyle(20)
	leg.AddEntry(gr,name,"L")
	MG.Add(gr)

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
"pval.root"
]
medfiles = [
#"pvals_median.root"
]

onepoint = ROOT.TGraphErrors()
onepoint.SetMarkerColor(38)
onepoint.SetMarkerStyle(20)
onepoint.SetMarkerSize(1.5)

fitfiles = [
"maxlh.root"
]

FILES = []
MEDFILES = []
GR=[]
fitFILES = []
fitGR=[]

colors=[1,2,4]
styles=[1,1,1]
colors_med=[38]
styles_med=[7]

fitcolors=[1]
fitstyles=[1]

names=[
"#splitline{Local P-value}{}",
#"#splitline{Local P-value}{(no di-jet)}",
#"#splitline{Local P-value}{(di-jet only)}"
]
names_med=["1xSM Higgs Median Expected"]

# significance lines
Lines = [1.,2.,3.]
Vals=[ROOT.RooStats.SignificanceToPValue(L) for L in Lines]
TLines=[ROOT.TLine(120,V,130,V) for V in Vals]

BinningBoundaries=[110,112.5,117.5,122.5,127.5,132.5,142.5,150]
BBTLines=[ROOT.TLine(BinningBoundaries[v],0.0001,BinningBoundaries[v],2.5) for v in range(3,5)]

fitLines = [1.]
fitTLines=[ROOT.TLine(120,fl,130,fl) for fl in fitLines]

legend=ROOT.TLegend(0.11,0.1,0.45,0.28)
legend.SetFillColor(0)
legend.SetBorderSize(0)

for i in range(len(medfiles)):
	MEDFILES.append(ROOT.TFile(medfiles[i]))
	GR.append(MEDFILES[i].Get("median"))
	Setup(GR[-1],colors_med[i],styles_med[i],legend,names_med[i],MG)

for i in range(len(files)):
	FILES.append(ROOT.TFile(files[i]))
	ggr =FILES[i].Get("observed")
	print ggr.GetName() 
	print ROOT.RooStats.PValueToSignificance(ggr.Eval(124.0))
	print ggr.Eval(124.0)
	GR.append(FILES[i].Get("observed"))
	onepoint.SetPoint(0,124.0,0.000265217)
	onepoint.SetPointError(0,0,0.0000240085)
	
	Setup(GR[-1],colors[i],styles[i],legend,names[i],MG)

#MG.Add(onepoint)
legend.AddEntry(onepoint,"Observed p-value (Ensemble)","PEl")

for i in range(len(fitfiles)):
	fitFILES.append(ROOT.TFile(fitfiles[i]))
	fitGR.append(fitFILES[i].Get("bestfit"))
	fitSetup(fitGR[i],fitcolors[i],fitstyles[i],fitMG)

ROOT.gROOT.SetStyle("Plain")
c = ROOT.TCanvas("c","c",1600,1800)
# Make 2 pads for the signal strength fit.
up = ROOT.TPad("u","u",0.01,0.7,0.99,0.99);
dp = ROOT.TPad("d","d",0.01,0.01,0.99,0.7);
up.SetNumber(1);
dp.SetNumber(2);
up.Draw();
dp.Draw();

mytext= ROOT.TLatex()
mytext.SetTextSize(0.045)

c.cd(1)
fitMG.Draw("AE4C")
#fitMGLINE = ROOT.TGraph()
#for m in range(110,151): fitMGLINE.SetPoint(m-110,float(m),fitMG.Eval(float(m)))
fitMG.Draw("CP")
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
MG.Draw("ALP")
MG.GetXaxis().SetRangeUser(110,150)
MG.GetYaxis().SetRangeUser(0.0001,2.5)

MG.GetXaxis().SetTitle("m_{H} (GeV)")
MG.GetYaxis().SetTitle("p-value")
MG.GetYaxis().SetTitleOffset(0.6)
MG.GetYaxis().SetTitleSize(0.05)
MG.GetXaxis().SetTitleSize(0.05)
MG.GetYaxis().SetLabelSize(0.05)
MG.GetXaxis().SetLabelSize(0.05)
onepoint.Draw("PE2L")
text = ROOT.TLatex()

#text.SetNDC()
for j,TL in enumerate(TLines):
	TL.SetLineStyle(7)
	TL.SetLineColor(2)
	TL.SetLineWidth(1)
	TL.Draw("same")
	text.DrawLatex(130.5,Vals[j]*0.88,"%d #sigma"%Lines[j])

for j,TL in enumerate(BBTLines):
  TL.SetLineStyle(7)
  TL.SetLineColor(4)
  TL.SetLineWidth(2)
  TL.Draw("same")
#mytext.SetFillColor(8)	
mytext.SetNDC()
legend.SetFillColor(10)
legend.SetFillStyle(0)
legend.Draw()

Box = ROOT.TBox(137,0.0015,147,0.0055)
Box.SetFillColor(10)
Box.SetFillStyle(1001)
Box.Draw()
#mytext
mytext.DrawLatex(0.55,0.22,"#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = %.2f fb^{-1}}"%float(intlumi))
dp.SetGrid(True)
dp.SetLogy()


#raw_input() 
c.SaveAs("pvalues.pdf")
#dp.SaveAs("pvals_bin"+bin+".pdf")
#c.SaveAs("pvaluesFit_bin"+bin+".pdf")


