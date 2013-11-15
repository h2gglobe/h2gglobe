#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-o","--outfile",help="Outputfile")
parser.add_option("-S","--sqrts",type='int',default=8,help="Sqrt(S) COM energy for finding strings etc")

(options,args) = parser.parse_args()

files = args[:]

import sys
import ROOT

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gROOT.LoadMacro("scripts/roofit_iterate.h+")
workspace = ROOT.RooWorkspace("multipdf","multipdf")

for fi in files:

  f   = ROOT.TFile.Open(fi)
  win = f.Get("multipdf")
  ext = "%dTeV"%options.sqrts

  cat = None
  cNum = -1;

  while cat == None:
   cNum+=1
   if cNum>1000: sys.exit ("Cannot Find Category index in file %s" %fi)
   cat = win.cat("pdfindex_%d_%s"%(cNum,ext))

  print "Merging Category ",cNum

  pdf   = win.pdf("CMS_hgg_cat%d_%s_bkgshape"%(cNum,ext))
  norm  = win.var("CMS_hgg_cat%d_%s_bkgshape_norm"%(cNum,ext))
  data = win.data("roohist_data_mass_cat%d"%(cNum))

  getattr(workspace,'import')(cat)
  getattr(workspace,'import')(pdf)
  getattr(workspace,'import')(norm)
  getattr(workspace,'import')(data)


outfile   = ROOT.TFile(options.outfile,"RECREATE")
workspace.Write()
outfile.Close()  
