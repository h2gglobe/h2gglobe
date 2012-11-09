#!/usr/bin/env python

import math, ROOT

fin=open("results.csv")

sumDeltaMu = 0.
sumDeltaMuSq = 0.
n = 0.

h = ROOT.TH1F("deltaMu","deltaMu;#delta #mu;n psedo-experiments",65,0.25,0.75)

for line in fin.read().split("\n"):
    line =line.lstrip().rstrip()
    if line == "":
        continue

    vals = [ float(v) for v in line.split(",") if v != "" ]

    if len(vals) != 7:
        continue
    ipart, mass1, mu1, sig1, mass2, mu2, sig2 = vals
    deltaMu = mu1 - mu2

    sumDeltaMu += deltaMu 
    sumDeltaMuSq += deltaMu*deltaMu
    n += 1.

    h.Fill(deltaMu)
    
avgJkn = sumDeltaMu / n
varJkn = (n-1)/n * ( sumDeltaMuSq - n*avgJkn*avgJkn )

print avgJkn,
print math.sqrt(varJkn)

ROOT.gStyle.SetOptFit(1)
c=ROOT.TCanvas()
h.Fit("gaus")
h.Draw()

for fmt in "C","png","pdf":
    c.SaveAs("jackknife.%s" % fmt )

