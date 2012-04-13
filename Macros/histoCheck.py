import ROOT
import sys

filename = sys.argv[1]
tfile = ROOT.TFile(filename)

for key in tfile.GetListOfKeys():
  if 'th1f' in key.GetName() and ('BDT' in key.GetName() or 'VBF' in key.GetName()):
    th1f = tfile.Get(key.GetName())
    print "%60s %10d %4.4f"%(th1f.GetName(),th1f.GetEntries(),th1f.Integral())
