import ROOT
import sys
import numpy

files = [ ROOT.TFile.Open(f) for f in sys.argv[1:] ]


def drawModel(mh,ncat=6,color=ROOT.kBlue,opt="hist",procs = [ "ggh", "vbf", "wzh", "tth" ]):

    for p in procs:
        c = ROOT.TCanvas("canv_%s_%1.5g" % (p,mh),"canv_%s_%1.5g" % (p,mh), 900, 400*ncat/3)
        c.Divide(3, ncat / 3 )
        globals()[c.GetName().replace(".","_")] = c
        for cat in range(ncat):
            c.cd(cat+1)
            h = ROOT.gDirectory.Get("th1f_sig_%s_mass_m%1.5g_cat%d" % ( p, mh, cat ))
            
            h.SetLineColor(color)
            h.SetMarkerColor(color)

            h.Draw(opt)

            print h.GetName(), h.Integral()
            

def getPoint(h):

    from math import sqrt

    qi  = ROOT.GraphToTF1("%s_func" % h.GetName(), ROOT.TGraph(h) )
    fit = ROOT.TF1( "%s_func" % h.GetName(), qi, 110., 180.,1, "GraphToTF1")
    fit.SetParameter(0,0)
    hm  = fit.GetMaximum()*0.5
    hmx = fit.GetMaximumX()
    high = fit.GetX(hm, hmx, 180)
    low  = fit.GetX(hm, 110, hmx)

    return h.Integral(), h.Integral()/sqrt(h.GetEntries()), high-low, h.GetMean()
    

def drawYields(masses,ncat=6,color=ROOT.kBlue,opt="APL",procs = [ "ggh", "vbf", "wzh", "tth" ]):
    
    ROOT.gROOT.LoadMacro("ResultScripts/GraphToTF1.C+")

    for p in procs:
        c = ROOT.TCanvas("canv_%s" % (p),"canv_%s" % (p), 900, 400*ncat/3)
        c.Divide(3, ncat / 3 )
        globals()[c.GetName().replace(".","_")] = c

        d = ROOT.TCanvas("canv_%sFWHM" % (p),"canv_%sFWHM" % (p), 900, 400*ncat/3)
        d.Divide(3, ncat / 3 )
        globals()[d.GetName().replace(".","_")] = d

        e = ROOT.TCanvas("canv_%sMean" % (p),"canv_%sMean" % (p), 900, 400*ncat/3)
        e.Divide(3, ncat / 3 )
        globals()[d.GetName().replace(".","_")] = e

        for cat in range(ncat):
            gcat = ROOT.TGraphErrors()
            gcat.SetName("%s_cat%d" % ( p, cat) )
            gcatFWHM = ROOT.TGraphErrors()
            gcatFWHM.SetName("%s_cat%d_FWMH" % ( p, cat) )
            gcatMean = ROOT.TGraphErrors()
            gcatMean.SetName("%s_cat%d_Mean" % ( p, cat) )
            globals()[gcat.GetName().replace(".","_")] = gcat
            globals()[gcatFWHM.GetName().replace(".","_")] = gcatFWHM
            globals()[gcatMean.GetName().replace(".","_")] = gcatMean
            for mh in masses:
                name = "th1f_sig_%s_mass_m%1.5g_cat%d" % ( p, mh, cat )
                h = ROOT.gDirectory.Get(name)
                try:
                    integral, integralE, fwhm, mean = getPoint(h)
                    ip = gcat.GetN()
                    gcat.SetPoint(ip,mh,integral)
                    if( abs(round(mh / 5.)*5. - mh) < 1.e-5 ):
                        gcat.SetPointError(ip,0.,integralE)
                    
                    gcatFWHM.SetPoint(ip,mh,fwhm)
                    gcatMean.SetPoint(ip,mh,mean)
                    
                except:
                    print "%s not found" % name

            
            c.cd(cat+1)
            gcat.SetLineColor(color)
            gcat.SetMarkerColor(color)
            gcat.SetMarkerStyle(ROOT.kFullCircle)
            
            gcat.Print()
            gcat.Draw(opt)

            d.cd(cat+1)
            gcatFWHM.SetLineColor(color)
            gcatFWHM.SetMarkerColor(color)
            gcatFWHM.SetMarkerStyle(ROOT.kFullCircle)

            gcatFWHM.Print()
            gcatFWHM.Draw(opt)

            e.cd(cat+1)
            gcatMean.SetLineColor(color)
            gcatMean.SetMarkerColor(color)
            gcatMean.SetMarkerStyle(ROOT.kFullCircle)

            gcatMean.Print()
            gcatMean.Draw(opt)
