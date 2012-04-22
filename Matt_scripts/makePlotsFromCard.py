import ROOT
import os,numpy,sys,math,array
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--input",dest="tfileName")
parser.add_option("-o","--output",dest="outfileName")
parser.add_option("-a","--additional",dest="addlines")
parser.add_option("-m","--mass",dest="mass",type="float")
parser.add_option("","--noVbfTag",dest="includeVBF",default=True,action="store_false")

(options,args)=parser.parse_args()
options.splitSignal=False

ROOT.gROOT.SetBatch(1)
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)

# make defaults
# PLOT OPS ----------------
lumistring = "5.09 fb^{-1}"
sigscale   = 5.
plotOutDir="./plots/"
mass=float(options.mass)

def getMLBkg(mlfile):
  mlF = open(mlfile)
  mlBkg=[]
  for l in mlF.readlines():
    bkg = l.split()[len(l.split())-1]
    if l.split()[1]=="bkg":
      mlBkg.append(float(bkg))
  return mlBkg

def getMLsb(mlfile):
  mlF = open(mlfile)
  mlSB=[]
  for l in mlF.readlines():
    sb = l.split()[len(l.split())-2]
    if l.split()[1]=="bkg":
      thisBkg=float(sb)
    else:
      thisBkg+=float(sb)
      if l.split()[1]=="vbf":
        mlSB.append(thisBkg)
  return mlSB 

def plainBin(hist):
	nb = hist.GetNbinsX()
	h2 = ROOT.TH1F(hist.GetName()+"new","",nb,0,nb)
	for i in range (1,nb+1):
		h2.SetBinContent(i,hist.GetBinContent(i))
		h2.SetBinError(i,hist.GetBinError(i))
		if (options.includeVBF):
			if hist.GetBinLowEdge(i+1) < 1.:
			  h2.GetXaxis().SetBinLabel(i,"BDT Bin %d "%(i))
			else: 
			  h2.GetXaxis().SetBinLabel(i," Di-jet ")
	h2.GetXaxis().SetNdivisions(nb)
	return h2

def plotDistributions(mass,data,signals,bkg,bkgml,sbml,errors):

	sb = ROOT.TCanvas("sb","sb",1600,1600)
	up = ROOT.TPad("u","u",0.01,0.3,0.99,0.99)
	up.SetLogy()
	#dp.SetLogy()
	dp = ROOT.TPad("d","d",0.01,0.01,0.99,0.3)
	up.SetNumber(1)
	dp.SetNumber(2)
	up.Draw()
	dp.Draw()
	sb.cd(1)
  
	if options.splitSignal: # last signal is separated off
	  for i in range(1,len(signals)-1):
		signals[0].Add(signals[i])
	else:
	  for i in range(1,len(signals)):
		signals[0].Add(signals[i])

	nbins = data.GetNbinsX()

	flatdata    = plainBin(data)
	flatsignal  = plainBin(signals[0])
	flatsignal1 = plainBin(signals[-1])

	flatbkg  = plainBin(bkg);flatbkg.SetLineColor(1);flatbkg.SetLineWidth(1)

	flatbkgml = plainBin(bkgml); flatbkgml.SetLineColor(2); flatbkgml.SetLineWidth(2)
	flatsbml  = plainBin(sbml); flatsbml.SetLineColor(4); flatsbml.SetLineWidth(2)

	fNew  = flatbkg.Clone()
	fNew2 = flatbkg.Clone()

	flatdata.SetMarkerStyle(20);flatdata.SetMarkerSize(1.0)
		
	fNew.SetLineColor(1);fNew.SetLineWidth(2);fNew.SetFillStyle(1001);fNew.SetFillColor(3)
	fNew2.SetFillColor(5);fNew2.SetLineColor(1);fNew2.SetLineWidth(2);fNew2.SetFillStyle(1001)

	fNewT = fNew.Clone();fNew2T = fNew2.Clone();fNewT.SetFillStyle(1001);fNew2T.SetFillStyle(1001)

	flatsignal.SetLineWidth(2);flatsignal.SetLineStyle(ROOT.kDashed);flatsignal.SetLineColor(ROOT.kRed);flatsignal.Scale(sigscale)
	flatsignal1.SetLineWidth(2);flatsignal1.SetLineStyle(ROOT.kDashed); flatsignal1.SetLineColor(ROOT.kGreen+4);flatsignal1.Scale(sigscale)
		
	leg = ROOT.TLegend(0.6,0.5,0.88,0.88);leg.SetFillColor(0);leg.SetBorderSize(0)

	for b in range(1,nbins+1):
		additional = errors[b-1]
  		fNew.SetBinError(b,((fNew.GetBinError(b)**2)+((fNew.GetBinContent(b)*additional)**2))**0.5)
  		fNew2.SetBinError(b,2*(((fNew2.GetBinError(b)**2)+((fNew2.GetBinContent(b)*additional)**2))**0.5))
	if (not options.includeVBF): flatdata.GetXaxis().SetTitle("Category")
	flatdata.Draw("9");fNew2.Draw("9sameE2");fNew.Draw("9sameE2");flatbkg.Draw("9samehist")
	if options.splitSignal: 
	  sigst = ROOT.THStack();sigst.Add(flatsignal);sigst.Add(flatsignal1);sigst.Draw("9samehist")
	else:  flatsignal.Draw("9samehist")
	flatdata.Draw("9sameP");flatdata.SetMinimum(1.0);flatdata.SetMaximum(20*flatdata.Integral())
	flatbkgml.Draw("9samehist"); flatsbml.Draw("9samehist")

  
	leg.AddEntry(flatdata,"Data","PLE")
	if options.splitSignal:
	  leg.AddEntry(flatsignal,"Higgs (GG,WZ,TT), m_{H}=%3.1f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	  leg.AddEntry(flatsignal1,"Higgs, m_{H}=%3.1f GeV (x%d)"%(mass,int(sigscale)) ,"L")

	else: leg.AddEntry(flatsignal,"Higgs, m_{H}=%3.1f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	leg.AddEntry(flatbkg,"Background","L");leg.AddEntry(fNewT,"\pm 1\sigma","F");leg.AddEntry(fNew2T,"\pm 2\sigma","F")
	leg.AddEntry(flatbkgml,"Background ML","L")
	leg.AddEntry(flatsbml,"S+B ML","L")
	leg.Draw()
	
	chi2bkgFlat = flatdata.Chi2Test(flatbkg)
	chi2bkgML = flatdata.Chi2Test(flatbkgml)
	chi2sbML = flatdata.Chi2Test(flatsbml)
	lat = ROOT.TPaveText(0.3,0.6,0.5,0.89,"NDC")
	lat.SetLineColor(0)
	lat.SetFillColor(0)
	lat.SetShadowColor(0)
	lat.SetTextAlign(11)
	lat.AddText("#chi^{2} (B) = %1.3f"%chi2bkgFlat)
	lat.AddText("#chi^{2} (B_{ML}) = %1.3f"%chi2bkgML)
	lat.AddText("#chi^{2} (SB_{ML}) = %1.3f"%chi2sbML)
  
	mytext = ROOT.TLatex();mytext.SetTextSize(0.03);mytext.SetNDC();mytext.DrawLatex(0.1,0.94,"CMS preliminary,  #sqrt{s} = 7 TeV ");mytext.SetTextSize(0.04)
	mytext.DrawLatex(0.5,0.94,"#int L = %s"%(lumistring))
	lat.Draw() 
	leg.Draw() 
	up.SaveAs(plotOutDir+"/model_m%3.1f.pdf"%mass);up.SaveAs(plotOutDir+"/model_m%3.1f.png"%mass)
  
  # write to file
	outF = ROOT.TFile(options.outfileName,"UPDATE")
	flatdata.SetName("data_%3.1f"%mass)
	flatbkg.SetName("bkg_%3.1f"%mass)
	fNew.SetName("bkg1_%3.1f"%mass)
	fNew2.SetName("bkg2_%3.1f"%mass)
	flatbkgml.SetName("bkgml_%3.1f"%mass)
	flatsbml.SetName("sbml_%3.1f"%mass)
	flatsignal.SetName("sig_%3.1f"%mass)
	flatdata.Write(); flatbkg.Write(); fNew.Write(); fNew2.Write(); flatbkgml.Write(); flatsbml.Write(); flatsignal.Write()
	outF.Close()

	sb.cd(2)
	flatsignal.Sumw2()
	sob = flatsignal.Clone(); sobml = flatsignal.Clone(); sobsbml = flatsignal.Clone()
	sob.GetXaxis().SetLabelSize(0.1)
	sob.GetYaxis().SetLabelSize(0.08)
	sob.Divide(flatbkg)
	sobml.Divide(flatbkgml)
	sobsbml.Divide(flatsbml)
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
	sob.Draw("9ep")
	sobml.Draw("9epsame")
	sobsbml.Draw("9epsame")
	leg3.Draw("same")
	#mytext.DrawLatex(0.2,0.6,"#int L = %s"%(lumistring))
	sb.SaveAs(plotOutDir+"/frac_model_m%3.1f.pdf"%mass);sb.SaveAs(plotOutDir+"/frac_model_m%3.1f.png"%mass)

	d = ROOT.TCanvas()
	leg2 = ROOT.TLegend(0.56,0.56,0.88,0.88)
	leg2.SetFillColor(0);leg2.SetBorderSize(0)
	if (not options.includeVBF): flatdata.GetXaxis().SetTitle("Category")
	flatdata.GetYaxis().SetTitle("Data - Background");
	datErrs = []
	for b in range(1,nbins+1): datErrs.append((flatdata.GetBinContent(b))**0.5);
	flatdata.Add(flatbkg,-1)
	for b in range(1,nbins+1): 
		flatdata.SetBinError(b,((datErrs[b-1]*datErrs[b-1]) +(fNew.GetBinError(b)*fNew.GetBinError(b)))**0.5 )
	flatbkg.Add(flatbkg,-1)

	flatdata.Draw("9");flatbkg.Draw("9samehist")
	if options.splitSignal: 
	  sigst = ROOT.THStack();sigst.Add(flatsignal);sigst.Add(flatsignal1)
	  sigst.Draw("9samehist")
	else:
	  flatsignal.Draw("9samehist")
	flatdata.Draw("9sameP");flatdata.SetMaximum(250);flatdata.SetMinimum(-100)

	leg2.AddEntry(flatdata,"Data","PLE")
	if options.splitSignal:
	  leg2.AddEntry(flatsignal,"Higgs (GG,WZ,TT), m_{H}=%3.0f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	  leg2.AddEntry(flatsignal1,"Higgs, m_{H}=%3.0f GeV (x%d)"%(mass,int(sigscale)) ,"L")

	else: leg2.AddEntry(flatsignal,"Higgs, m_{H}=%3.0f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	leg2.AddEntry(flatbkg,"Background","L")

	mytext.SetTextSize(0.04)

	leg2.Draw()
	mytext = ROOT.TLatex();mytext.SetTextSize(0.03);mytext.SetNDC();mytext.DrawLatex(0.1,0.92,"CMS preliminary,  #sqrt{s} = 7 TeV ");mytext.SetTextSize(0.04)
	mytext.DrawLatex(0.2,0.8,"#int L = %s"%(lumistring))
	d.SaveAs(plotOutDir+"/diff_model_m%3.1f.pdf"%mass);d.SaveAs(plotOutDir+"/diff_model_m%3.1f.png"%mass)
# First read the data-card
tfile = open(options.tfileName,"r")
dcardlines = tfile.readlines()
# we need observation, rate, bkg_norm and massBias lines 
infoDict = {"data":[],"rate":[],"bkg_norm":[],"massBias":[],"ncat":int(1),"bkg":[],"ggh":[],"vbf":[],"wzh":[],"tth":[]}
for dl in dcardlines:
	info = dl.split()
	if len(info)<1: continue
	if info[0]=="observation":
		infoDict["data"]=info[1:]
		infoDict["ncat"]=len(info[1:])
	elif info[0]=="rate":
		infoDict["rate"]=info[1:]
	elif info[0]=="bkg_norm":
		infoDict["bkg_norm"]=info[2:]
	elif "massBias" in info[0]:
		infoDict["massBias"].append(info[2:])

# get info from mlfit
infoDict["bkgml"] = getMLBkg(options.addlines)
infoDict["sbml"]  = getMLsb(options.addlines)

# make some bins
binEdges = [-1.+2*float(n)/infoDict["ncat"] for n in range(infoDict["ncat"])]
binEdges.append(1.)
binEdges.append(1.04)
arrbinEdges = array.array("d",binEdges)
# Now we have everyting to constuct histograms:
ncat = infoDict["ncat"]
if ncat!=len(infoDict["bkgml"]) and ncat !=len(infoDict["sbml"]): sys.exit('Wrong number of bins')
dataHist = ROOT.TH1F("dataHist","",ncat,arrbinEdges)
bkgHist  = ROOT.TH1F("bkgHist","",ncat,arrbinEdges)
gghHist  = ROOT.TH1F("gghHist","",ncat,arrbinEdges)
vbfHist  = ROOT.TH1F("vbfHist","",ncat,arrbinEdges)
wzhHist  = ROOT.TH1F("wzhHist","",ncat,arrbinEdges)
tthHist  = ROOT.TH1F("tthHist","",ncat,arrbinEdges)
bkgmlHist = ROOT.TH1F("bkgmlHist","",ncat,arrbinEdges)
sbmlHist = ROOT.TH1F("sbmlHist","",ncat,arrbinEdges)
		
# unfold the signals + bkg rates
for b in range(infoDict["ncat"]):
	infoDict["ggh"].append(infoDict["rate"][b*5+0])
	infoDict["vbf"].append(infoDict["rate"][b*5+1])
	infoDict["wzh"].append(infoDict["rate"][b*5+2])
	infoDict["tth"].append(infoDict["rate"][b*5+3])
	infoDict["bkg"].append(infoDict["rate"][b*5+4])

for b in range(infoDict["ncat"]):
	dataHist.SetBinContent(b+1,float(infoDict["data"][b]))	
	bkgHist.SetBinContent(b+1,float(infoDict["bkg"][b]))	
	gghHist.SetBinContent(b+1,float(infoDict["ggh"][b]))	
	vbfHist.SetBinContent(b+1,float(infoDict["vbf"][b]))	
	wzhHist.SetBinContent(b+1,float(infoDict["wzh"][b]))	
	tthHist.SetBinContent(b+1,float(infoDict["tth"][b]))
	bkgmlHist.SetBinContent(b+1,float(infoDict["bkgml"][b]))
	sbmlHist.SetBinContent(b+1,float(infoDict["sbml"][b]))

# final thing is to make the errors
binErrors=[0 for b in range(infoDict["ncat"])]
for mb in infoDict["massBias"]:
  
  for b in range(infoDict["ncat"]):
	errVal = mb[b*5+4]
	if errVal!="-":
		errsc = float((errVal.split("/"))[1])-1
		err = float(infoDict["bkg"][b])*errsc
		binErrors[b]+=(err*err)

for b in range(infoDict["ncat"]): 
	bkgHist.SetBinError(b+1,binErrors[b]**0.5)

#finall norm errors
normErrors = []
for b in range(infoDict["ncat"]): 
	normErrors.append(float(infoDict["bkg_norm"][b*5+4])-1)
# Can finally make plots
plotDistributions(options.mass,dataHist.Clone(),[gghHist.Clone(),wzhHist.Clone(),tthHist.Clone(),vbfHist.Clone()],bkgHist.Clone(),bkgmlHist.Clone(),sbmlHist.Clone(),normErrors)
