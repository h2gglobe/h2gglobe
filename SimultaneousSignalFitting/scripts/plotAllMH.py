#!/usr/bin/env python

import os
import sys
import numpy

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-i","--inputfile",dest="inputfile",type="str",default="dat/config.dat",help="Input config file")
parser.add_option("-p","--proc",dest="proc",type="str",help="Process")
parser.add_option("-c","--cat",dest="cat",type="int",help="Category")
parser.add_option("--runAll",dest="runAll",action="store_true",default=False,help="Do all procs all cats")
parser.add_option("-L","--mHlow",dest="mHlow",type="int",default=110,help="Lowest mH")
parser.add_option("-H","--mHhigh",dest="mHhigh",type="int",default=150,help="Highest mH")
parser.add_option("-S","--mHstep",dest="mHstep",type="float",default=1,help="Step of mH")
(options,args)=parser.parse_args()

import ROOT

def plot(proc,cat):
  if proc=='all':
    if cat==-1:
      pdf = ws.pdf('sigpdfrelAllCats_allProcs')
      if not pdf:
        print 'sigpdfrelAllCats_allProcs not found'
        return
    else:
      pdf = ws.pdf("sigpdfrelcat%d_allProcs"%cat)
      if not pdf:
        print 'sigpdfrelcat%d_allProcs not found'%cat
        return
  else:
    pdf = ws.pdf("hggpdfrel_%s_cat%d"%(proc,cat))
    if not pdf:
      print 'hggpdfrel_%s_cat%d not found'%(proc,cat)
      return

  if not pdf: sys.exit("Not found")

  plot = mass.frame(ROOT.RooFit.Range(options.mHlow,options.mHhigh))

  for m in numpy.arange(options.mHlow,options.mHhigh+options.mHstep,options.mHstep):
    mh.setVal(m)
    pdf.plotOn(plot)
    #if m==110: pdf.plotOn(plot)
    #else: pdf.plotOn(plot,ROOT.RooFit.DrawOption("same"))

  c = ROOT.TCanvas()
  plot.Draw()
  if cat==-1:
    c.Print("custom_plots/%s_catall.pdf"%(proc))
    c.Print("custom_plots/%s_catall.png"%(proc))
  else:
    c.Print("custom_plots/%s_cat%d.pdf"%(proc,cat))
    c.Print("custom_plots/%s_cat%d.png"%(proc,cat))
   

inF = ROOT.TFile(options.inputfile)
ws = inF.Get("wsig_8TeV")

mass = ws.var("CMS_hgg_mass")
mh = ws.var("MH")

ROOT.gROOT.SetBatch()
os.system("mkdir -p custom_plots")


if options.runAll:
  ROOT.gROOT.SetBatch()
  plot('all',-1)
  for cats in range(9):
    plot('all',cats)
  for procs in ['ggh','vbf','wzh','tth']:
    for cats in range(9):
      plot(procs,cats)
else:
  plot(options.proc,options.cat)


