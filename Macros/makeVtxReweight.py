#================================================================================================
# makeVtxRegweight.py
# Constructs a ROOT file with the Efficiency to select the right/wrong vertex in the H->gg system
# Adapted from original script by P. Musella
#================================================================================================

class options:
   def __init__(self):
	self.eff = "/afs/cern.ch/user/m/malberti/public/scaleFactors/vtxIdEff_vs_bosonPt_globe_Glu121_S6_PUweights_2011_0100_73500.root "
	self.eff_ratio = "/afs/cern.ch/user/m/malberti/public/scaleFactors/BDT_vtxIdScaleFactorFromZmumu_DYJetsToLL_Fall11_S6_Run2011all.root "
	self.n_categories = 8
	self.outfile = "vertex_reweighing.root"

# prevent ROOT from parsing command line
from ROOT import *
from math import sqrt

o=options()

gStyle.SetMarkerSize(1.5)
gROOT.SetBatch(True)

##
## Read files
## 
### files = []
### files.append( { "file" : o.eff, "label" : "eff" } )
### 
### ph = PlotHelper(files)
### 
### # histograms to be read
### histos_to_read = [
###     ("higgsPt",                [("Rebin",2)]),
###     ("matchVtxHiggsPt",          [("Rebin",2)]),
###     ]
### 
### # efficiencies to compute
### efficiencies = [
###     ("effMatchVsPt","matchVtxHiggsPt", "higgsPt", "Selection efficiency (matching);p_{T}(#gamma,#gamma); Efficiency (Matching)"),
### ]
### 
### # read and divide histograms
### ph.read(histos_to_read)
### ph.divide(efficiencies,[(yrange,(0.,1.05))])
### mc_eff = ph.histos["eff"]["effMatchVsPt"]

fr = TFile.Open(o.eff_ratio)
eff_ratio = fr.Get("scaleFactor")

eff_ratio_func = TF1("eff_ratio_func", "[0]+([1] + [2]*x)*(x<=20)",0,200)
eff_ratio.Fit(eff_ratio_func)

fe = TFile.Open(o.eff)
mc_eff = fe.Get("efficiency")

pass_rewei = mc_eff.Clone("pass_rewei")
pass_rewei.SetTitle("Wrong vertex reweighing; p_{T}(#gamma #gamma); Weight")
fail_rewei = mc_eff.Clone("fail_rewei") 
fail_rewei.SetTitle("Wrong vertex reweighing; p_{T}(#gamma #gamma); Weight")


lfwe = 0.
for i in range(mc_eff.GetN()):
    if mc_eff.GetX()[i] != eff_ratio.GetX()[i]:
         print "Efficeincy and ratio have different binning %d %f %f" % ( i, mc_eff.GetX()[i] , eff_ratio.GetX()[i] )
         sys.exit(1)
        
    pw  = eff_ratio.GetY()[i]
    pwe = eff_ratio.GetErrorY(i)

    x = mc_eff.GetX()[i]
    xe = mc_eff.GetErrorX(i)

    fpw  = eff_ratio_func.Eval(x)
    ## print pw - fpw
    ## pwe = 1.e-2

    eff  = mc_eff.GetY()[i]
    effe = mc_eff.GetErrorY(i)
    ## if eff < 0.95 and eff*pw <0.95:
    if eff < 0.98:
        fw  = (1. - eff*pw)  / (1. - eff)
        fwe = sqrt( (1.-fw)*(1.-fw)*pwe*pwe + fw*fw/((1-eff)*(1.-eff))*effe*effe )
        lfwe = fwe
    else:
        fw  = 1.
        try:
            fwe = sqrt( (1.-fw)*(1.-fw)*pwe*pwe + fw*fw/((1-eff)*(1.-eff))*effe*effe )
        except:
            fwe = sqrt( (1.-fw)*(1.-fw)*pwe*pwe + fw*fw/(0.02*0.02)*effe*effe )
        ## fwe = lfwe

    pass_rewei.SetPoint(i,x,pw)
    pass_rewei.SetPointError(i,xe,xe,pwe,pwe)
    fail_rewei.SetPoint(i,x,fw)
    fail_rewei.SetPointError(i,xe,xe,fwe,fwe)

    print "x = %1.1f pass_w = %1.2f +- %1.2f fail_w = %1.2f +- %1.2f " % ( x, pw, pwe, fw, fwe ) 
    
fout = TFile.Open(o.outfile, "recreate")
fout.cd()

for c in range(o.n_categories):
    pass_rewei.Clone("ratioVertex_cat%d_pass" % c).Write()
    fail_rewei.Clone("ratioVertex_cat%d_fail" % c).Write()

fout.Close()
fr.Close()

gROOT.Reset()
