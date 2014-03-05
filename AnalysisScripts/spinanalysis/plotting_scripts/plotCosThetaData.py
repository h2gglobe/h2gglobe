import ROOT as r
import sys

if len(sys.argv)!=3:
	sys.exit('usage python plotCosThetaData.py <file> <sqrtS>')

tf = r.TFile(sys.argv[1])
sqrtS = int(sys.argv[2])
r.TH1.SetDefaultSumw2()
r.gStyle.SetOptStat(0)

def fillHist(h,trees,mlow,mhigh):
	for tree in trees:
		for e in range(tree.GetEntries()):
			tree.GetEntry(e)
			if tree.higgs_mass>=mlow and tree.higgs_mass<=mhigh:
				h.Fill(r.TMath.Abs(tree.costheta_cs),tree.evweight)
	h.Scale(1./h.Integral())

dataTrees = [tf.Get('spin_trees/Data')]
smTrees = [tf.Get('spin_trees/ggh_m125_%dTeV'%sqrtS),
					 tf.Get('spin_trees/vbf_m125_%dTeV'%sqrtS),
					 tf.Get('spin_trees/wzh_m125_%dTeV'%sqrtS),
					 tf.Get('spin_trees/tth_m125_%dTeV'%sqrtS)]
ggTrees = [tf.Get('spin_trees/gg_grav_m125_%dTeV'%sqrtS)]
qqTrees = [tf.Get('spin_trees/qq_grav_m125_%dTeV'%sqrtS)]

dataHist_m120 = r.TH1F('data_m120','',25,0.,1.)
dataHist_m125 = r.TH1F('data_m125','',25,0.,1.)
dataHist_m130 = r.TH1F('data_m130','',25,0.,1.)
smHist = r.TH1F('sm','',25,0.,1.)
ggHist = r.TH1F('gg','',25,0.,1.)
qqHist = r.TH1F('qq','',25,0.,1.)

fillHist(dataHist_m120,dataTrees,120-2.5,120+2.5)
fillHist(dataHist_m125,dataTrees,125-2.5,125+2.5)
fillHist(dataHist_m130,dataTrees,130-2.5,130+2.5)
fillHist(smHist,smTrees,100.,180.)
fillHist(ggHist,ggTrees,100.,180.)
fillHist(qqHist,qqTrees,100.,180.)

smHist.SetLineColor(r.kRed)
ggHist.SetLineColor(r.kBlue)
qqHist.SetLineColor(r.kGreen+1)
dataHist_m120.SetLineColor(r.kGray+1)
dataHist_m125.SetLineColor(r.kMagenta)
dataHist_m130.SetLineColor(r.kMagenta+2)
dataHist_m120.SetMarkerColor(r.kGray+1)
dataHist_m125.SetMarkerColor(r.kMagenta)
dataHist_m130.SetMarkerColor(r.kMagenta+2)

smHist.SetLineWidth(2)
ggHist.SetLineWidth(2)
qqHist.SetLineWidth(2)
smHist.SetLineStyle(2)
ggHist.SetLineStyle(2)
qqHist.SetLineStyle(2)

dataHist_m120.SetLineWidth(2)
dataHist_m125.SetLineWidth(2)
dataHist_m130.SetLineWidth(2)
dataHist_m120.SetMarkerStyle(r.kFullCircle)
dataHist_m125.SetMarkerStyle(r.kFullSquare)
dataHist_m130.SetMarkerStyle(r.kFullTriangleUp)

smHist.GetYaxis().SetTitle('a.u.')
smHist.GetXaxis().SetTitle('|cos(#theta_{CS})|')
smHist.GetYaxis().SetRangeUser(0.,0.1)

leg = r.TLegend(0.53,0.7,0.74,0.89)
leg.SetHeader('Simulation')
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(smHist,'0^{+} (SM) m_{H}=125','L')
leg.AddEntry(ggHist,'2_{m}^{+} (gg) m_{H}=125','L')
leg.AddEntry(qqHist,'2_{m}^{+} (qq) m_{H}=125','L')
leg2 = r.TLegend(0.74,0.7,0.89,0.89)
leg2.SetHeader('Data')
leg2.SetFillColor(0)
leg2.SetLineColor(0)
leg2.AddEntry(dataHist_m120,'m_{#gamma#gamma}=120#pm2.5','LEP')
leg2.AddEntry(dataHist_m125,'m_{#gamma#gamma}=125#pm2.5','LEP')
leg2.AddEntry(dataHist_m130,'m_{#gamma#gamma}=130#pm2.5','LEP')

pt = r.TPaveText(0.1,0.91,0.45,0.99,"NDC");
pt.SetTextAlign(12);
pt.SetTextSize(0.04);
pt.SetFillColor(0);
pt.AddText("CMS Preliminary");
pt.SetBorderSize(0);
pt2 = r.TPaveText(0.75,0.90,0.9,0.99,"NDC");
pt2.SetTextAlign(32);
pt2.SetTextSize(0.04);
pt2.SetFillColor(0);
if sqrtS==7:
	pt2.AddText('#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}')
elif sqrtS==8:
	pt2.AddText('#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}')


canv = r.TCanvas()
smHist.Draw('HIST')
ggHist.Draw('HISTsame')
qqHist.Draw('HISTsame')
dataHist_m120.Draw('LEPsame')
dataHist_m125.Draw('LEPsame')
dataHist_m130.Draw('LEPsame')
leg.Draw('same')
leg2.Draw('same')
pt.Draw('same')
pt2.Draw('same')
canv.Print('cosTheta_data_%dTeV.pdf'%sqrtS)
canv.Print('cosTheta_data_%dTeV.png'%sqrtS)
raw_input('Looks ok?\n')
