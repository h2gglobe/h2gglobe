#!/bin/env python

import sys
import os
import ROOT

#
# Handy histogram manipulation
# 
def applyModifs(h,modifs):
    for method in modifs:
        args = None
        if type(method) == tuple:
            method, args = method
        if type(method) == str:
            if hasattr(h,method):
                method = getattr(h,method)
            else:
                method = globals()[method]
        if args == None:            
            try:
                method(h)
            except:
                method()
        else:            
            if not ( type(args) == tuple or type(args) == list ):
                args = tuple([args])
            try:
                method(h,*args)
            except:
                method(*args)
                
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
    h.GetYaxis().SetTitle(tit % { "binw" : h.GetBinWidth(0) } )
            
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
    print fin, plot, names
    histos = readPlot(fin,plot,which=[""],samples=names)[0]
    for iplot in range(len(histos)):
        h = histos[iplot]
        hname = names[iplot]
        applyModifs(h,subproc[hname])

    sum = histos[0].Clone(name)
    sum.SetTitle(title)

    for h in histos[1:]:
        sum.Add(h)

    applyModifs(sum,plotmodifs)
    applyModifs(sum,style)
    
    return sum

#
# Prepare canvas and legend
#
def makeCanvAndLeg(name,legPos):
    canv = ROOT.TCanvas(name)
    leg  = ROOT.TLegend(*legPos)

    leg.SetFillStyle(0), leg.SetLineColor(ROOT.kWhite)## , leg.SetShadowColor(ROOT.kWhite)

    return canv, leg

#
# Make THStack out of python list
#
def makeStack(name,histos):
    stk = ROOT.THStack()
    for h in histos:
        stk.Add(h)
    return stk

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
        getattr(stk,method.split(",")[0])("%s SAME" % option)
        ymax = stk.GetMaximum(option)
        
    return ymax

#
# Perform data/MC comparison for many plots and categories
#
def dataMcComparison(data, bkg, sig, plots, categories=[0], savefmts=["C","png","pdf"]):

    objs = []
    canvs = []
    # loop over categories
    for cat in categories:
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
                ymax = max(ymax,drawStack(bkgstk,dm,bkgopt))
                objs.append(bkgstk)
                
            # then data
            if len(datahists) > 0:
                datastk = makeStack("data_%s_cat%d" % ( plotname, cat),datahists)
                ### getattr(datastk,dm)("%s SAME" % dataopt)
                ### ymax = max(ymax,datastk.GetMaximum())
                ymax = max(ymax,drawStack(datastk,dm,dataopt))
                objs.append(datastk)

            # and finally signal
            if len(sighists) > 0:
                sigstk = makeStack("sig_%s_cat%d" % ( plotname, cat),sighists)
                ### getattr(sigstk,dm)("%s SAME" % sigopt)
                ### ymax = max(ymax,sigstk.GetMaximum(sigopt))
                ymax = max(ymax,drawStack(sigstk,dm,sigopt))
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
                    drawStack(bkgstk,"Draw",bkgopt)
                
                if len(datahists) > 0:
                    drawStack(datastk,"Draw",dataopt)
                    
                if len(sighists) > 0:
                    drawStack(sigstk,"Draw",sigopt)

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
    mcFudge = 1.4*2.

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
                       [(setcolors,406),(legopt,"f"),("Scale",1.4)],
                       { ## "qcd_30_8TeV_pf"   : [("Smooth",25)],
                         ## "qcd_40_8TeV_pf"   : [("Smooth",25)],
                        "gjet_20_8TeV_pf"  : [],
                        "gjet_40_8TeV_pf"  : []
                        }
                       ),
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
                           plots = [ ("vbf_mva",[("SetBinContent",(1,0.)),("Rebin",5),(xrange,(-1.,1)),
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
