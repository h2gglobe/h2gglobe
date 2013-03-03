#!/usr/bin/env python

import sys
import numpy
import os

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-i","--input",dest="input",type="str",help="Input file")
parser.add_option("--runAll",dest="runAll",action="store_true",default=False,help="Do all procs all cats")
parser.add_option("-p","--proc",dest="proc",type="str",help="Process")
parser.add_option("-c","--cat",dest="cat",type="int",help="Category")
parser.add_option("-d","--dir",dest="dir",type="str",help="Dump values to file in this directory (not print to screen)")
(options,args)=parser.parse_args()

def dumpValues(ws,proc,cat):
  mh = ws.var('MH')
  mass = ws.var('CMS_hgg_mass')
  pdf = ws.pdf('hggpdfrel_%s_cat%d'%(proc,cat))
  comps = pdf.getComponents()
  it = comps.createIterator()
  comp = it.Next()
  vals=[]
  for i in range(comps.getSize()):
    name = comp.GetName()
    if 'func' in name:
      g = int(name.split('_g')[1].split('_%s'%proc)[0])
      type = name.split('_g')[0].split('func_')[1]
      for m in numpy.arange(110,151,5):
        mh.setVal(m)
        vals.append('%s_mh%d_g%d %1.5f'%(type,m,g,comp.getVal()))
    comp = it.Next()
  
  vals.sort()
  if options.dir:
    f = open('%s/initFit_%s_cat%d.dat'%(options.dir,proc,cat),'w')
  
  for s in vals:
    if options.dir:
      f.write('%s\n'%s)
    else:
      print s

  if options.dir: f.close()

import ROOT as r
r.gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so')

if options.dir:
  os.system('mkdir -p %s'%options.dir)

tf = r.TFile(options.input)
ws = tf.Get('wsig_8TeV')
if options.runAll:
  for proc in ['ggh','vbf','wzh','tth']:
    for cat in range(0,9):
      dumpValues(ws,proc,cat)
else:
  dumpValues(ws,options.proc,options.cat)
