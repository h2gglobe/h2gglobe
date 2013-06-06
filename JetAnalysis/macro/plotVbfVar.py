#!/bin/env python

import sys
import os
import ROOT

from math import fabs

def tdrGrid(tdrStyle, gridOn):
#def tdrGrid(TStyle tdrStyle, bool gridOn):
    tdrStyle.SetPadGridX(gridOn)
    tdrStyle.SetPadGridY(gridOn)
    return

# fixOverlay: Redraws the axis

def fixOverlay():
  gPad.RedrawAxis()
  return

def setTDRStyle():
  tdrStyle=ROOT.TStyle("tdrStyle","Style for P-TDR")

# For the canvas:
  tdrStyle.SetCanvasBorderMode(0) 
  tdrStyle.SetCanvasColor(ROOT.kWhite) 
  tdrStyle.SetCanvasDefH(600)  #Height of canvas
  tdrStyle.SetCanvasDefW(800)  #Width of canvas
  tdrStyle.SetCanvasDefX(0)    #POsition on screen
  tdrStyle.SetCanvasDefY(0) 

# For the Pad:
  tdrStyle.SetPadBorderMode(0) 
  # tdrStyle.SetPadBorderSize(Width_t size = 1) 
  tdrStyle.SetPadColor(ROOT.kWhite) 
  tdrStyle.SetPadGridX(False) 
  tdrStyle.SetPadGridY(False) 
  tdrStyle.SetGridColor(0) 
  tdrStyle.SetGridStyle(3) 
  tdrStyle.SetGridWidth(1) 

# For the frame:
  tdrStyle.SetFrameBorderMode(0) 
  tdrStyle.SetFrameBorderSize(1) 
  tdrStyle.SetFrameFillColor(0) 
  tdrStyle.SetFrameFillStyle(0) 
  tdrStyle.SetFrameLineColor(1) 
  tdrStyle.SetFrameLineStyle(1) 
  tdrStyle.SetFrameLineWidth(1) 

# For the histo:
  # tdrStyle.SetHistFillColor(1) 
  # tdrStyle.SetHistFillStyle(0) 
  tdrStyle.SetHistLineColor(1) 
  tdrStyle.SetHistLineStyle(0) 
  tdrStyle.SetHistLineWidth(1) 
  # tdrStyle.SetLegoInnerR(Float_t rad = 0.5) 
  # tdrStyle.SetNumberContours(Int_t number = 20) 

  tdrStyle.SetEndErrorSize(2) 
  #tdrStyle.SetErrorMarker(20)   # Seems to give an error
  tdrStyle.SetErrorX(0.) 
  
  tdrStyle.SetMarkerStyle(20) 

#For the fit/function:
  tdrStyle.SetOptFit(0) 
  tdrStyle.SetFitFormat("5.4g") 
  tdrStyle.SetFuncColor(2) 
  tdrStyle.SetFuncStyle(1) 
  tdrStyle.SetFuncWidth(1) 

#For the date:
  tdrStyle.SetOptDate(0) 
  # tdrStyle.SetDateX(Float_t x = 0.01) 
  # tdrStyle.SetDateY(Float_t y = 0.01) 

# For the statistics box:
  tdrStyle.SetOptFile(0) 
  tdrStyle.SetOptStat(0)  # To display the mean and RMS:   SetOptStat("mr") 
  tdrStyle.SetStatColor(ROOT.kWhite) 
  tdrStyle.SetStatFont(42) 
  tdrStyle.SetStatFontSize(0.025) 
  tdrStyle.SetStatTextColor(1) 
  tdrStyle.SetStatFormat("6.4g") 
  tdrStyle.SetStatBorderSize(1) 
  tdrStyle.SetStatH(0.1) 
  tdrStyle.SetStatW(0.15) 
  # tdrStyle.SetStatStyle(Style_t style = 1001) 
  # tdrStyle.SetStatX(Float_t x = 0) 
  # tdrStyle.SetStatY(Float_t y = 0) 

# Margins:
  tdrStyle.SetPadTopMargin(0.05) 
  tdrStyle.SetPadBottomMargin(0.13) 
  tdrStyle.SetPadLeftMargin(0.16) 
  tdrStyle.SetPadRightMargin(0.10) 

# For the Global title:
  tdrStyle.SetOptTitle(0)     # 0=No Title
  tdrStyle.SetTitleFont(42) 
  tdrStyle.SetTitleColor(1) 
  tdrStyle.SetTitleTextColor(1) 
  tdrStyle.SetTitleFillColor(10) 
  tdrStyle.SetTitleFontSize(0.07) 
  # tdrStyle.SetTitleH(0)  # Set the height of the title box
  # tdrStyle.SetTitleW(0)  # Set the width of the title box
  # tdrStyle.SetTitleX(0)  # Set the position of the title box
  # tdrStyle.SetTitleY(0.985)  # Set the position of the title box
  # tdrStyle.SetTitleStyle(Style_t style = 1001) 
  # tdrStyle.SetTitleBorderSize(2) 

# For the axis titles:
  tdrStyle.SetTitleColor(1, "XYZ") 
  tdrStyle.SetTitleFont(42, "XYZ") 
  tdrStyle.SetTitleSize(0.08, "XYZ") 
  # tdrStyle.SetTitleXSize(Float_t size = 0.02)  # Another way to set the size?
  # tdrStyle.SetTitleYSize(Float_t size = 0.02) 
  tdrStyle.SetTitleXOffset(1.5) 
  tdrStyle.SetTitleYOffset(1.5)
  # tdrStyle.SetTitleOffset(1.1, "Y")  # Another way to set the Offset

# For the axis labels:
  tdrStyle.SetLabelColor(1, "XYZ") 
  tdrStyle.SetLabelFont(42, "XYZ") 
  tdrStyle.SetLabelOffset(0.007, "XYZ") 
  tdrStyle.SetLabelSize(0.05, "XYZ") 

# For the axis:
  tdrStyle.SetAxisColor(1, "XYZ") 
  tdrStyle.SetStripDecimals(ROOT.kTRUE) 
  tdrStyle.SetTickLength(0.03, "XYZ") 
  tdrStyle.SetNdivisions(510, "XYZ") 
  tdrStyle.SetPadTickX(0)   # 0=Text labels (and tics) only on bottom, 1=Text labels on top and bottom
  tdrStyle.SetPadTickY(1) 

# Change for log plots:
  tdrStyle.SetOptLogx(0) 
  tdrStyle.SetOptLogy(0) 
  tdrStyle.SetOptLogz(0) 

# Postscript options:
  tdrStyle.SetPaperSize(20.,20.) 
  # tdrStyle.SetLineScalePS(Float_t scale = 3) 
  # tdrStyle.SetLineStyleString(Int_t i, const char* text) 
  # tdrStyle.SetHeaderPS(const char* header) 
  # tdrStyle.SetTitlePS(const char* pstitle) 

  # tdrStyle.SetBarOffset(Float_t baroff = 0.5) 
  # tdrStyle.SetBarWidth(Float_t barwidth = 0.5) 
  # tdrStyle.SetPaintTextFormat(const char* format = "g") 
  # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0) 
  # tdrStyle.SetTimeOffset(Double_t toffset) 
  # tdrStyle.SetHistMinimumZero(kTRUE) 

  #gROOT.ForceStyle()   # Try this if stuff doesn't work right
  
  # Find RooFit include
  #gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include") 

  tdrStyle.cd() 
  return


#
# Handy histogram manipulation
# 
def applyModifs(h,modifs):
    for method in modifs:
        args = None
        ret = None
        if type(method) == tuple:
            method, args = method
        if type(method) == str:
            if hasattr(h,method):
                method = getattr(h,method)
            else:
                method = globals()[method]
        if args == None:            
            try:
                ret = method(h)
            except:
                ret = method()
        else:            
            if not ( type(args) == tuple or type(args) == list ):
                args = tuple([args])
            try:
                ret = method(h,*args)
            except:
                ret = method(*args)
        if ret and ret != h:
            print "Replacing ", h, ret, method
            h = ret
            ret = None
    return h
# 
def setcolors(h,color):
    h.SetLineColor(color)
    h.SetFillColor(color)
    h.SetMarkerColor(color)
# 
def legopt(h,opt):
    h.legopt = opt
# 
def xrange(h,a,b):
    h.GetXaxis().SetRangeUser(a,b)
# 
def xtitle(h,tit):
    h.GetXaxis().SetTitle(tit)
# 
def ytitle(h,tit):
    h.GetYaxis().SetTitle( (tit % { "binw" : h.GetBinWidth(0), "xtitle" : h.GetXaxis().GetTitle() }).replace("(GeV)","") )
            
#
# Read plots from globe histogram files
#
def readPlot(fin, name, cat=0, which=["_sequential","_nminus1"], samples=["diphojet_8TeV","ggh_m125_8TeV","vbf_m125_8TeV",]):

    ret = []
    for p in which:
        hists = []
        for s in samples:
            nam = "%s%s_cat%d_%s" % ( name, p, cat, s )
            h = fin.Get(nam)
            print nam
            h.GetXaxis().SetTitle(  h.GetXaxis().GetTitle().replace("@"," ") )
            ## print nam, h
            hists.append(h)
        ret.append(hists)
            
    return ret

#
# Quick MVA control plots selection
#
def selectionControlPlots(fname):
    vars=[("vbf_mva",               [""]),
          ("cut_VBFLeadJPt",        ["_sequential","_nminus1"]),
          ("cut_VBFSubJPt",         ["_sequential","_nminus1"]),
          ("cut_VBF_dEta",          ["_sequential","_nminus1"]),
          ("cut_VBF_Zep",           ["_sequential","_nminus1"]),
          ("cut_VBF_Mjj",           ["_sequential","_nminus1"]),
          ("cut_VBF_dPhi",          ["_sequential","_nminus1"]),
          ("cut_VBF_DiPhoPtOverM",  ["_sequential","_nminus1"]),
          ("cut_VBF_LeadPhoPtOverM",["_sequential","_nminus1"]),
          ("cut_VBF_SubPhoPtOverM", ["_sequential","_nminus1"]),
          ]
    
    
    fin = ROOT.TFile.Open(fname)
    
    objs=[]
    for var,which in vars:
        plots = readPlot(fin,var,which=which)
        
        c = ROOT.TCanvas(var,var,800,400)
        objs.append(c)
        np = len(plots)
        if np > 1:
            c.Divide(2,np/2)
        for ip in range(np):
            ps=plots[ip]
            
            ps[0].SetLineColor(ROOT.kRed)
            ps[1].SetLineColor(ROOT.kMagenta)

            if np>1:
                c.cd(ip+1)

            p0=ps[0].Clone("%s_%d" % (var,ip) )
            objs.append(p0)
            p0.Reset("ICE")
            p0.SetEntries(0)
            p0.Draw()
            if var == "vbf_mva":
                p0.GetXaxis().SetRangeUser(-1,1)
            for p in ps:
                if var == "vbf_mva":
                    p.SetBinContent(1,0.)
                    p.Draw("hist same")
                else:
                    p.DrawNormalized("hist same")
            
            maximum=0.
            for p in ROOT.gPad.GetListOfPrimitives():
                if p.ClassName().startswith("TH1"):
                    pmax = p.GetMaximum()
                    if pmax>maximum:
                        maximum=pmax
            p0.GetYaxis().SetRangeUser(0,maximum*1.5)
            print var,p0,maximum, p0.GetYaxis().GetXmax()
            ROOT.gPad.Modified()
            ROOT.gPad.Update()
                        
            for fmt in "C","png":
                c.SaveAs("%s.%s" % (c.GetName(),fmt))
                
    return objs

#
# Dump event yields for different samples
#
def eventYield(filenames,categories=[5,6],procs=["ggh_m125_8TeV","vbf_m125_8TeV","diphojet_8TeV"]):
    files = [ ROOT.TFile.Open(f) for f in filenames ]

    
    for fin in files:
        print fin.GetName()
        for cat in categories:
            plots=readPlot(fin, "mass", cat=cat, which=[""], samples=procs)[0]
            print "cat%d " % cat-1,
            for iproc in range(len(procs)):
                ## print "%1.3g" % plots[iproc].Integral(18,28), 115-135
                print "%1.3g" % plots[iproc].Integral(18,28),
            print

#
# Read histograms for a given process, applying manipulators
#
def readProc(fin,name,title,style,subproc,plot,plotmodifs,category):

    names = subproc.keys()
    print fin, plot, names, name
    histos = readPlot(fin,plot,which=[""],samples=names,cat=category)[0]
    print histos
    for iplot in range(len(histos)):
        h = histos[iplot]
        hname = names[iplot]
        print h,hname
        h = applyModifs(h,subproc[hname])

    print len(histos)
    sum = histos[0].Clone(name)
    sum.SetTitle(title)
    
    for h in histos[1:]:
        sum.Add(h)
        
    sum = applyModifs(sum,plotmodifs)
    sum = applyModifs(sum,style)
    
    return sum

#
# Prepare canvas and legend
#
def makeCanvAndLeg(name,legPos):
    canv = ROOT.TCanvas(name)
    leg  = ROOT.TLegend(*legPos)

    leg.SetFillStyle(0), leg.SetLineColor(ROOT.kWhite)## , leg.SetShadowColor(ROOT.kWhite)

    return canv, leg

def makeLegend(legPos):
    leg  = ROOT.TLegend(*legPos)
    leg.SetFillStyle(0), leg.SetLineColor(ROOT.kWhite)## , leg.SetShadowColor(ROOT.kWhite)

    return leg

#
# Make THStack out of python list
#
def makeStack(name,histos):
    stk = ROOT.THStack()
    for h in histos:
        stk.Add(h)
    return stk

def stackTitles(stk):
    stk.GetHistogram().GetXaxis().SetTitle( stk.GetHists()[0].GetXaxis().GetTitle() )
    stk.GetHistogram().GetYaxis().SetTitle( stk.GetHists()[0].GetYaxis().GetTitle() )

#
# Make THStack out of python list
#
def makeEnvelope(name,histos,stPlus=None,stMinus=None):
    nominal = histos[0]
    errPlus  = nominal.Clone( "%s_ErrPlus" % name )
    errMinus = nominal.Clone( "%s_ErrMinus" % name )
    if stPlus:
        applyModifs( errPlus, stPlus )
        applyModifs( errMinus, stMinus )
    for ibin in range(nominal.GetNbinsX()):
        hist = ROOT.TH1F("hist","hist",11,-5.,5.)
        hist.Reset("ICE")
        hist.Sumw2()
        points = []
        
        plus  = nominal.GetBinContent(ibin+1)
        minus = nominal.GetBinContent(ibin+1)
        nom   = nominal.GetBinContent(ibin+1)
        err   = nominal.GetBinError(ibin+1)
        points.append( [nom,err] )
        hist.Fill(nom)
        for h in histos[1:]:
            content =  h.GetBinContent(ibin+1)
            err     =  h.GetBinError(ibin+1)
            hist.Fill(content)
            points.append( [content,err] )
            if content < minus:
                minus = content
            if content > plus:
                plus = content

        if hist.GetRMS() == 0.:
            errPlus.SetBinContent(ibin+1,plus)
            errMinus.SetBinContent(ibin+1,minus)
            continue
            
        hist2 = ROOT.TH1F("hist2","hist2",11,hist.GetMean()-5.*hist.GetRMS(),hist.GetMean()+5.*hist.GetRMS())
        hist2.Sumw2()
        for p,e in points:
            hist2.Fill(p)
            
        func = ROOT.TF1("func","[0]*exp( -0.5*pow( (x-[1])/( (x>=0)*[2] + (x<=0)*[3] ) ,2.) )",hist2.GetMean()-5.*hist2.GetRMS(),hist2.GetMean()+5.*hist2.GetRMS())
        func.SetParameters(len(histos),hist2.GetMean(),hist2.GetRMS(), hist2.GetRMS())
        ## func.SetParLimits(2,-0.5*hist.GetMean(),0.5*hist.GetMean())
        ## func.SetParLimits(3,-0.1*hist.GetMean(),0.1*hist.GetMean())
        stat = int( hist2.Fit( func, "L" ) )
        if stat == 0:
            errPlus.SetBinContent(ibin+1, max(nom,func.GetParameter(1)+fabs(func.GetParameter(2))))
            errMinus.SetBinContent(ibin+1,min(nom,func.GetParameter(1)-fabs(func.GetParameter(3))))
        else:
            errPlus.SetBinContent(ibin+1,plus)
            errMinus.SetBinContent(ibin+1,minus)
        
    return makeStack(name,[errPlus,errMinus,nominal])


def shiftHisto(offset,nominal,shifted):
    for ibin in range(offset.GetNbinsX()):
        shifted.SetBinContent(ibin+1, shifted.GetBinContent(ibin+1) + offset.GetBinContent(ibin+1) - nominal.GetBinContent(ibin+1) )
        shifted.SetBinError(ibin+1, ROOT.TMath.Sqrt( shifted.GetBinError(ibin+1)*shifted.GetBinError(ibin+1) - 
                                                     nominal.GetBinError(ibin+1)*nominal.GetBinError(ibin+1) ) )
        ## shifted.SetBinContent(ibin+1, shifted.GetBinContent(ibin+1) - nominal.GetBinContent(ibin+1) )
    return shifted

def ratioHisto(num,den,ytit):
    num.Divide(den)
    ytitle(num,ytit)
    return num

#
# Draw a THStack
#
def drawStack(stk, method, option):
    ymax = 0.
    if "DrawNormalized" in method:
        rng = None
        if "[" in method:
            rng = [ float(f) for f in method.split("DrawNormalized")[1].split("[")[1].split("]")[0].split(",") ]
            print rng
        histos = [ stk.GetStack().At(0) ]
        if "nostack" in option:
            histos = stk.GetHists()
            option = option.replace("nostack","")
        for h in histos:
            h.SetFillStyle(0)
            h.SetLineWidth(2)
            bmin = -1
            bmax = -1
            if rng:
                bmin = h.FindBin(rng[0])
                bmax = h.FindBin(rng[1])
            h.Scale(1./h.Integral(bmin,bmax))
            h.Draw("%s SAME" % option)
            ymax = max(ymax, h.GetMaximum()) 
    else:
        getattr(stk,method.split(",")[0])("%s" % option)
        ymax = stk.GetMaximum(option)
        print stk.GetName(), stk.GetHists()[0].Integral()
        
    return ymax

#
# Perform data/MC comparison for many plots and categories
#
def dataMcComparison(data, bkg, sig, plots, categories=[0], savefmts=["C","png","pdf"]):

    objs = []
    canvs = []
    # loop over categories
    for cat in categories:
        print cat
        # loop over plots
        for plot in plots:

            plotname, plotmodifs, drawopts, legPos = plot
            dm, dataopt, bkgopt, sigopt = drawopts
            
            bkghists = []
            sighists = []
            datahists = []
            
            # read background MC
            if bkg != None:
                bkgfile, bkgprocs = bkg
                bkghists = [ readProc(bkgfile,*subprocs,plot=plotname,plotmodifs=plotmodifs,category=cat) for subprocs in bkgprocs ]
                
            # read signal MC
            if sig != None:
                sigfile, sigprocs = sig
                sighists = [ readProc(sigfile,*subprocs,plot=plotname,plotmodifs=plotmodifs,category=cat) for subprocs in sigprocs ]
                
            # read data
            if data != None:
                datafile, dataprocs = data
                datahists = [ readProc(datafile,*subprocs,plot=plotname,plotmodifs=plotmodifs,category=cat) for subprocs in dataprocs ]
                
            # collect histograms
            allhists = datahists+bkghists+sighists
            objs += allhists
            
            # make empty frame histogram for convenience
            frame = allhists[0].Clone()
            frame.Reset("ICE")
            frame.SetEntries(0)
            objs.append(frame)
            ymax = 0.
            ymin = 0.

            # allocate canvas and legend and draw frame
            canv,leg = makeCanvAndLeg("%s_cat%d" % ( plotname, cat), legPos )
            objs.append(canv)
            objs.append(leg)
            canvs.append(canv)
            frame.Draw()

            # draw background first
            if len(bkghists) > 0:
                bkgstk = makeStack("bkg_%s_cat%d" % ( plotname, cat), bkghists)
                ### getattr(bkgstk,dm)("%s SAME" % bkgopt)
                ### ymax = max(ymax,bkgstk.GetMaximum(bkgopt))
                ymax = max(ymax,drawStack(bkgstk,dm,"%s SAME"%bkgopt))
                objs.append(bkgstk)
                
            # then data
            if len(datahists) > 0:
                datastk = makeStack("data_%s_cat%d" % ( plotname, cat),datahists)
                ### getattr(datastk,dm)("%s SAME" % dataopt)
                ### ymax = max(ymax,datastk.GetMaximum())
                ymax = max(ymax,drawStack(datastk,dm,"%s SAME"%dataopt))
                objs.append(datastk)

            # and finally signal
            if len(sighists) > 0:
                sigstk = makeStack("sig_%s_cat%d" % ( plotname, cat),sighists)
                ### getattr(sigstk,dm)("%s SAME" % sigopt)
                ### ymax = max(ymax,sigstk.GetMaximum(sigopt))
                ymax = max(ymax,drawStack(sigstk,dm,"%s SAME"%sigopt))
                objs.append(sigstk)

            # make legend
            for h in allhists:
                legopt = "f"
                if hasattr(h,"legopt"):
                    legopt = h.legopt
                leg.AddEntry(h,"",legopt)

            # adjust yaxis
            frame.GetYaxis().SetRangeUser(ymin,ymax*2)
            leg.Draw("same")
            canv.RedrawAxis()

            # if needed draw inset with zoom-in
            if "DrawInset" in dm:
                inset =  [ float(f) for f in dm.split("DrawInset")[1].split("[")[1].split("]")[0].split(",") ]
                rng = inset[0:2]
                pos = inset[2:]
                
                padname = "%s_cat%d_inset" % ( plotname, cat)
                pad = ROOT.TPad(padname, padname, *pos)
                objs.append(pad)
                pad.Draw("")
                pad.SetFillStyle(0)
                
                pad.cd()
                padframe = frame.Clone()
                padframe.GetXaxis().SetRangeUser(*rng)
                padframe.GetYaxis().SetRangeUser(ymin,ymax*1.2)
                padframe.Draw()
                objs.append(padframe)
                
                if len(bkghists) > 0:
                    drawStack(bkgstk,"Draw",bkgopt+" same")
                
                if len(datahists) > 0:
                    drawStack(datastk,"Draw",dataopt+" same")
                    
                if len(sighists) > 0:
                    drawStack(sigstk,"Draw",sigopt+" same")

                pad.RedrawAxis()
                
    # save plots
    for c in canvs:
        for fmt in savefmts:
            c.SaveAs("%s.%s" % (c.GetName(),fmt))
            
    return objs

if __name__ == "__main__":

    #
    # open input files
    #
    fdata = ROOT.TFile.Open(sys.argv[1])
    fbkg  = ROOT.TFile.Open(sys.argv[2])
    fsig  = ROOT.TFile.Open(sys.argv[3])

    #
    # prepare output folder
    #
    outdir = sys.argv[4]
    try:
        os.mkdir(outdir)
    except:
        pass
    os.chdir(outdir)

    #
    # draw options and style
    # 
    defdrawopt = ("Draw","e","hist","hist nostack") ## drawing method(s), drawing option data, MC background, MC signal
    mvadrawopt = ("Draw,DrawInset[0.85,1,0.45,0.48,0.9,0.92]","e","hist","hist nostack")
    ## ## mvadrawopt = ("DrawNormalized[0.9,1]","e","hist nostack","hist nostack")
    ## mvadrawopt = ("DrawNormalized[0.9,1],DrawInset[0.85,1,0.45,0.426,0.9,0.9]","e","hist nostack","hist nostack")
    ## defdrawopt = ("DrawNormalized","e","hist nostack","hist nostack")
    deflegPos  = (0.137,0.525,0.474,0.888)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    
    #
    # datasets to be plot
    #
    mcFudge = 1.7

    # data
    data = [ fdata, [ ("data", "data",
                       [("SetMarkerStyle",(ROOT.kFullCircle)),(legopt,"pl")],
                       { "Data" : [] }
                       )
                      ]
             ]
    # background
    bkg  = [ fbkg,  [ ("pp","2 Prompt #gamma",          ## process name
                       [(setcolors,591),(legopt,"f"),("Scale",mcFudge)], ## style
                       { "diphojet_8TeV"      : [],     ## subprocesses to add, with possible manpulators
                         "dipho_Box_25_8TeV"  : [],
                         "dipho_Box_250_8TeV" : [],
                         ### "gjet_20_8TeV_pp"    : [],
                         ### "gjet_40_8TeV_pp"    : []
                         }
                       ),
                      ("pp","1 Prompt #gamma 1 Fake #gamma",
                       [(setcolors,406),(legopt,"f"),("Scale",mcFudge)],
                       { ## "qcd_30_8TeV_pf"   : [("Smooth",25)],
                         ## "qcd_40_8TeV_pf"   : [("Smooth",25)],
                        "gjet_20_8TeV_pf"  : [],
                        "gjet_40_8TeV_pf"  : []
                        }
                       ),
                      ### ("dy","Z/#gamma^{*} #rightarrow ee",
                      ###  [(setcolors,606),(legopt,"f"),("Scale",mcFudge)],
                      ###  { "DYJetsToLL"  : []
                      ###   }
                      ###  ),
                      ]
             ]
    # signal
    sig   = [ fsig,  [ ("ggH", "25 x ggH H #rightarrow #gamma #gamma",
                        [(setcolors,ROOT.kMagenta),("SetLineWidth",2),("SetFillStyle",0),("Scale",25.),(legopt,"l")],
                        { "ggh_m125_8TeV" : [] }
                        ),
                       ("qqH", "25 x qqH H #rightarrow #gamma #gamma",
                        [(setcolors,ROOT.kRed+1),("SetLineWidth",2),("SetFillStyle",0),("Scale",25.),(legopt,"l")],
                        { "vbf_m125_8TeV" : [] }
                        )
                       ]
              ]
    
    # Make data/MC comparison
    objs=dataMcComparison( data = data,
                           sig = sig,
                           bkg = bkg,
                           categories = [0],
                           ## categories = [0,5,6],
                           ## categories = [5,14],
                           ## categories = [0,5,9,14],
                           plots = [ ("vbf_mva",[("SetBinContent",(1,0.)),("Rebin",4),(xrange,(-1.,1)),
                                                 ("SetBinError",(1,0.)),(xtitle,"MVA"),## (ytitle,"A.U.")
                                                 (ytitle,"Events/%(binw)1.2g")
                                                 ],
                                      mvadrawopt,deflegPos),
                                     
                                     ("cut_VBFLeadJPt_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBFSubJPt_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBF_dEta_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBF_Zep_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBF_Mjj_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBF_dPhi_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBF_DiPhoPtOverM_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBF_LeadPhoPtOverM_sequential",[],defdrawopt,deflegPos),
                                     ("cut_VBF_SubPhoPtOverM_sequential",[],defdrawopt,deflegPos),

                                     ("pho_eta",[("Rebin",4)],defdrawopt,deflegPos),
                                     ("pho_r9",[],defdrawopt,deflegPos),
                                     ("pho1_eta",[("Rebin",4)],defdrawopt,deflegPos),
                                     ("pho2_eta",[("Rebin",4)],defdrawopt,deflegPos),
                                     ("pho1_pt",[],defdrawopt,deflegPos),
                                     ("pho2_pt",[],defdrawopt,deflegPos),
                                     
                                     ("cut_VBFLeadJPt_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBFSubJPt_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBF_dEta_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBF_Zep_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBF_Mjj_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBF_dPhi_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBF_DiPhoPtOverM_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBF_LeadPhoPtOverM_nminus1",[],defdrawopt,deflegPos),
                                     ("cut_VBF_SubPhoPtOverM_nminus1",[],defdrawopt,deflegPos),
                                     
                                     ] 
                           )

    
    ## eventYield(sys.argv[1:])
