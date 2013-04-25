#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--infile",dest="infile",help="File to search")
parser.add_option("-g","--grep",dest="grep",default=[],action="append",help="Grep string (can take multiple arguments)")
parser.add_option("-w","--workspace",dest="workspace",help="Use workspace instead")
(options,args)=parser.parse_args()

import ROOT as r

tf = r.TFile(options.infile)
if options.workspace: ws = tf.Get(options.workspace)

for key in tf.GetListOfKeys():
  if 'th1f' in key.GetName():
    for grep in options.grep:
      if grep in key.GetName():
        th1f = tf.Get(key.GetName())
        print '%60s  %8d  %10.4f'%(th1f.GetName(),th1f.GetEntries(),th1f.Integral())
        if options.workspace:
          dset = ws.data(key.GetName().replace('th1f','roohist'))
          print '%60s  %8d  %10.4f'%(dset.GetName(),dset.numEntries(),dset.sumEntries())

