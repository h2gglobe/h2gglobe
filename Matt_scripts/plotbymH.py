import ROOT
import sys,os

path=sys.argv[1]
outfile=sys.argv[2]
mass=float(sys.argv[3])

mlF = open('%s/%3.1fout.txt'%(path,mass))
dcF = open('%s/mva-datacard_grad_%3.1f.txt'%(path,mass))

tf = ROOT.TFile("../CMS-HGG_jan16reload_batch_refit.root_interpolated.root")
bkgH = tf.Get("th1f_bkg_grad_%3.1f_cat0_fitsb_biascorr"%mass)
bkgH2 = bkgH.Clone()
bkgFlat = bkgH.Clone()

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

if not os.path.isdir('plots'):
  os.makedirs('plots')

bkgFlat.SetLineColor(1); bkgFlat.SetLineWidth(2)
bkgH.SetLineColor(1); bkgH.SetLineWidth(2); bkgH.SetFillStyle(1001); bkgH.SetFillColor(3)
bkgH2.SetLineColor(1); bkgH2.SetLineWidth(2); bkgH.SetFillStyle(1001); bkgH2.SetFillColor(5)

bkgML = ROOT.TH1F("BkgML","BkgML",bkgH.GetNbinsX(),0,bkgH.GetNbinsX())
sbML  = ROOT.TH1F("SBML","SBML",bkgH.GetNbinsX(),0,bkgH.GetNbinsX())
dataH = ROOT.TH1F("data","data",bkgH.GetNbinsX(),0,bkgH.GetNbinsX())
dataH.SetMarkerStyle(20)
dataH.SetMarkerSize(1)

mlBkg=[]
mlSB=[]
dcBkg=[]
dcBkgEup=[]
dcBkgEdown=[]
paulErrs=[]

dcData=[]

for l in mlF.readlines():
  bkg = l.split()[len(l.split())-1]
  sb = l.split()[len(l.split())-2]
  if l.split()[1]=="bkg":
    mlBkg.append(float(bkg))
    thisBkg=float(sb)
  else:
    thisBkg+=float(sb)
    if l.split()[1]=="vbf":
      mlSB.append(thisBkg)

for l in dcF.readlines():
  if len(l.split())>0:
    if l.split()[0]=='observation':
      for el in l.split():
        if el=='observation': continue
        else: dcData.append(float(el))
    if l.split()[0]=='rate':
      for i,el in enumerate(l.split()):
        if el=='rate': continue
        elif i%5==0: dcBkg.append(float(el))
    if l.split()[0]=='bkg_norm':
      for i,el in enumerate(l.split()):
        if el=='bkg_norm': continue
        elif el=='lnN': continue
        elif (i-1)%5==0: paulErrs.append(float(el))

if len(paulErrs)!=bkgH.GetNbinsX():
  exit("ERROR - bins not equal")

for i, err in enumerate(paulErrs):
  bkgH.SetBinError(i+1,((bkgH.GetBinError(i+1)**2)+((bkgH.GetBinContent(i+1)*(err-1))**2))**0.5)
  bkgH2.SetBinError(i+1,2*((bkgH2.GetBinError(i+1)**2)+((bkgH2.GetBinContent(i+1)*(err-1))**2))**0.5)

if len(mlBkg)!=bkgML.GetNbinsX():
  exit("ERROR - bins not equal")

for i, num in enumerate(mlBkg):
  if (bkgML.GetBinLowEdge(i+1)<7.):
    bkgML.GetXaxis().SetBinLabel(i+1,"Bin %d"%(i+1))
  else:
    bkgML.GetXaxis().SetBinLabel(i+1,"Dijet")
  bkgML.SetBinContent(i+1,num)
  
bkgML.SetLineColor(2)
bkgML.SetLineWidth(3)

for i, num in enumerate(mlSB):
  if (sbML.GetBinLowEdge(i+1)<7.):
    sbML.GetXaxis().SetBinLabel(i+1,"Bin %d"%(i+1))
  else:
    sbML.GetXaxis().SetBinLabel(i+1,"Dijet")
  sbML.SetBinContent(i+1,num)

sbML.SetLineColor(4)
sbML.SetLineWidth(3)
  
for i, num in enumerate(dcData):
  if (dataH.GetBinLowEdge(i+1)<7):
    dataH.GetXaxis().SetBinLabel(i+1,"Bin %d"%(i+1))
  else:
    dataH.GetXaxis().SetBinLabel(i+1,"Dijet")
  dataH.SetBinContent(i+1,num)

leg = ROOT.TLegend(0.6,0.5,0.89,0.89)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetLineColor(0)
leg.AddEntry(dataH,"Data","lep")
leg.AddEntry(bkgFlat,"Background (datacard)","l")
leg.AddEntry(bkgH,"Bkg #pm 1#sigma (datacard)","f")
leg.AddEntry(bkgH2,"Bkg #pm 2#sigma (datacard)","f")
leg.AddEntry(bkgML,"Background (ML fit)","l")
leg.AddEntry(sbML,"S+B (ML fit)","l")

chi2bkgFlat = dataH.Chi2Test(bkgFlat)
chi2bkgML = dataH.Chi2Test(bkgML)
chi2sbML = dataH.Chi2Test(sbML)

lat = ROOT.TPaveText(0.3,0.6,0.5,0.89,"NDC")
lat.SetLineColor(0)
lat.SetFillColor(0)
lat.SetTextAlign(11)
lat.AddText("#chi^{2} (B) = %1.3f"%chi2bkgFlat)
lat.AddText("#chi^{2} (B_{ML}) = %1.3f"%chi2bkgML)
lat.AddText("#chi^{2} (SB_{ML}) = %1.3f"%chi2sbML)

canv = ROOT.TCanvas()
dataH.SetTitle("Mass %s.0"%mass)
dataH.GetYaxis().SetRangeUser(1,20000);
dataH.Draw("9le1p")
bkgH2.Draw("9sameE2")
bkgH.Draw("9sameE2")
bkgFlat.Draw("9samehist")
bkgML.Draw("9samehist")
sbML.Draw("9samehist")
dataH.Draw("9samele1p")
leg.Draw("same")
lat.Draw("same")
canv.SetLogy()
canv.SaveAs("plots/plotbymh_%3.1f.pdf"%(mass))
canv.SaveAs("plots/plotbymh_%3.1f.png"%(mass))

outF = ROOT.TFile(outfile,"UPDATE")
dataH.SetName("data_%3.1f"%mass)
dataH.Write()
bkgFlat.SetName("bkg_ic_%3.1f"%mass)
bkgFlat.Write()
bkgH.SetName("bkg_ic_1sig_%3.1f"%mass)
bkgH.Write()
bkgH2.SetName("bkg_ic_2sig_%3.1f"%mass)
bkgH2.Write()
bkgML.SetName("bkg_ml_%3.1f"%mass)
bkgML.Write()
sbML.SetName("sb_ml_%3.1f"%mass)
sbML.Write()
outF.Close()


  
