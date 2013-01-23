import os
import sys
import ROOT as r

inFile = r.TFile(sys.argv[1])
#outFile = r.TFile(sys.argv[2],"RECREATE")

inWS = inFile.Get('wsig_8TeV')

#outWS = r.RooFit.RooWorkspace('cms_hgg_workspace')

for arg in inWS.allPdfs():
  print arg.GetName()


