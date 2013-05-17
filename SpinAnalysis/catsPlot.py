#!/usr/bin/env python

import os
import re
import sys
import glob
import fnmatch
import ROOT

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

baseDirectory = sys.argv[1]
subDirs = [x[1] for x in os.walk(baseDirectory)]
subDirs = [baseDirectory+"/"+x for x in subDirs[0]]

home='$CMSSW_BASE/src/h2gglobe/SpinAnalysis'
startDir = os.getenv("PWD")

print "Starting dir: "+startDir

haddLine = "hadd -f %s/%s/separation.root"%(startDir,baseDirectory)

for dir in subDirs:
  #jobs = [x[1] for x in os.walk(dir)]
  #jobs = [dir+"/"+x for x in jobs[0]]
  jobs = glob.glob(dir+"/job*")
  print "Checking "+dir+" with %d jobs."%len(jobs)
  #os.system("cd "+startDir+"/"+dir)
  os.chdir(startDir+"/"+dir)

  temp = re.findall(r'\d+',dir.split("/")[-1])
  print "  This one has "+temp[0]+" categories in cosTheta*."

  if not os.path.exists("qmu.root"):
    print "\t  Unable to find qmu.root for %s"%dir
  else:
    processCommand = '%s/bin/diyPlot %s'%(home,startDir+"/"+dir+"/qmu.root --nCats "+temp[0])
    #print processCommand
    os.system(processCommand)
    os.chdir(startDir)

    haddLine += " %s/%s/%sCats_separation.root"%(startDir, dir, temp[0])

#print haddLine
if haddLine != "hadd -f %s/%s/separation.root"%(startDir,baseDirectory):
  os.system(haddLine)

  file = ROOT.TFile.Open("%s/%s/separation.root"%(startDir,baseDirectory))
  tree = file.Get("Separation")
  values = []
  for ev in tree:
    values.append((ev.UnbinnedGravSigma, ev.UnbinnedGravSigmaErr, ev.nCats))

  minVal = min([x[2] for x in values])
  maxVal = max([x[2] for x in values])

  hist = ROOT.TH1F("Separation", "Separation;nCats;p( q < median(0) | 2 )", int(len(values)+2), minVal-1.5, maxVal+1.5)
  for x in values:
    hist.SetBinContent(hist.FindBin(x[2]), x[0])
    hist.SetBinError(hist.FindBin(x[2]), x[1])
  errors = ROOT.TH1F(hist)
  errors2 = ROOT.TH1F(hist)

  errors.SetNameTitle("Separation", "Separation;nCats;p( q < median(0) | 2 )")

  canvas = ROOT.TCanvas("Sep", "Sep", 800, 600)
  #kBlue, kRed, kGreen, kCyan, kMagenta, kYellow
  errors.SetFillColor(ROOT.kBlue-10)
  hist.SetLineColor(ROOT.kBlue)
  errors2.SetLineColor(ROOT.kBlue)
  errors.Draw("E2")
  errors2.Draw("SAME E1 X0")
  hist.Draw("SAME hist")

  canvas.Update()
  canvas.SaveAs("%s/%s/Separation.pdf"%(startDir,baseDirectory))
  canvas.SaveAs("%s/%s/Separation.png"%(startDir,baseDirectory))
  canvas.SaveAs("%s/%s/Separation.C"%(startDir,baseDirectory))
  canvas.SaveAs("%s/%s/Separation.root"%(startDir,baseDirectory))
