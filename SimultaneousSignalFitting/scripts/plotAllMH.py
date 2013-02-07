#!/usr/bin/env python

import os
import sys
import ROOT

def plot(proc,cat):
  pdf = ws.pdf("hggpdfrel_%s_cat%s"%(proc,cat))

  if not pdf: sys.exit("Not found")

  plot = mass.frame(ROOT.RooFit.Range(105,155))

  for m in range(110,150):
    mh.setVal(m)
    pdf.plotOn(plot)
    #if m==110: pdf.plotOn(plot)
    #else: pdf.plotOn(plot,ROOT.RooFit.DrawOption("same"))

  c = ROOT.TCanvas()
  plot.Draw()
  c.Print("custom_plots/%s_cat%s.pdf"%(proc,cat))
  c.Print("custom_plots/%s_cat%s.png"%(proc,cat))
   

inF = ROOT.TFile(sys.argv[1])
ws = inF.Get("wsig_8TeV")

mass = ws.var("CMS_hgg_mass")
mh = ws.var("MH")

ROOT.gROOT.SetBatch()
os.system("mkdir -p custom_plots")

proc = raw_input("Which process?\n")
cat = raw_input("Which cat?\n")

if proc=='all' or cat=='all':
  ROOT.gROOT.SetBatch()
  for procs in ['ggh','vbf','wzh','tth']:
    for cats in range(9):
      plot(procs,cats)
else:
  plot(proc,cat)


