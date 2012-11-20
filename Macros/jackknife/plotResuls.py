#!/usr/bin/env python

import math, ROOT

import sys

dir=sys.argv[1]

## fin=open("results.csv")
fin=open("%s/results.csv" % dir)
## fin=open("results3.csv")
## fin=open("results_binned.csv")

fnev=open("nevents.txt")

### 125.5 1.46047 2.76678
### 123.5 1.01719 1.95587
### 1.45488
### 1.49518
### 1.01719

diffen = math.sqrt(0.547316*0.596972 + 0.532582*0.599454)

mun1 = 1.50375
mun2 = 1.01719

ntot1 = 29930
ntot2 = 29159

sumMu1 = 0.
sumMu2 = 0.
sumDeltaMu = 0.
sumDeltaMuSq = 0.
n = 0.

### h = ROOT.TH1F("deltaMu","deltaMu;#delta #mu;n psedo-experiments",65,0.25,0.75)
## h = ROOT.TH1F("deltaMu","deltaMu;#delta #mu_{(i)} - #delta #mu_{n};n psedo-experiments",65,-0.2,0.3)
h = ROOT.TH1F("deltaMu","deltaMu;#delta #mu_{(i)} / err_{#delta #mu_{(i)}}  - #delta #mu_{n} / err_{#delta #mu_{n}};n psedo-experiments",65,-0.2,0.3)
hnev1 = ROOT.TH1F("nevents1","ICHEP;n^{i}_{events} / n_{events};n psedo-experiments",20,0.996,0.998)
hnev2 = ROOT.TH1F("nevents2","ICHEP-like;n^{i}_{events} / n_{events};n psedo-experiments",20,0.996,0.998)
hnev2.SetLineColor(ROOT.kRed)

for line in fnev.read().split("\n"):
    line =line.lstrip().rstrip().replace("/",",")
    if line == "":
        continue
    vals = [ float(v) for v in line.split(",") if v != "" ]
    if len(vals) != 2:
        continue
    hnev1.Fill( vals[0] / ntot1)
    hnev2.Fill( vals[1] / ntot2)

for line in fin.read().split("\n"):
    line =line.lstrip().rstrip().replace("/",",")
    if line == "":
        continue

    vals = [ float(v) for v in line.split(",") if v != "" ]

    ## if len(vals) != 7:
    if len(vals) != 11:
        print "Pronblem with line", line
        continue
    ## ipart, mass1, mu1, sig1, mass2, mu2, sig2 = vals
    ipart, mass1, mu1, mu1ep, mu1em, sig1, mass2, mu2, mu2ep, mu2em, sig2 = vals
    ## deltaMu = (mu1 - mu2)/mu1
    diffe = math.sqrt( -(mu1ep*mu1em  + mu2ep*mu2em) )
    ## mu1 *= 1-0.016
    ## mu2 *= 1-0.016
    mu1 -= mun1
    mu2 -= mun2
    ### mu1 /= diffe
    ### mu2 /= diffe
    deltaMu = (mu1 - mu2)

    sumDeltaMu += deltaMu 
    sumDeltaMuSq += deltaMu*deltaMu
    sumMu1 += mu1
    sumMu2 += mu2
    n += 1.

    
    h.Fill(deltaMu)

deltaMun = (mun1 - mun2)
### deltaMun /= diffen
    
mun1 = 0.
mun2 = 0.

avgDelta = sumDeltaMu / n
varJkn = (n-1)/n * ( sumDeltaMuSq - n*avgDelta*avgDelta )
avgJkn = n*(mun1-mun2) - (n-1) * avgDelta

print mun1, sumMu1 / n, n*mun1 - (n-1)*sumMu1/n, mun2, sumMu2 /n, n*mun2 - (n-1)*sumMu2/n

### mun2 /= mun1
### mun1 = 1.

print avgDelta, (mun1-mun2), avgJkn,
print math.sqrt(varJkn),  deltaMun / math.sqrt(varJkn)

ROOT.gStyle.SetOptFit(1)
c=ROOT.TCanvas()
h.Fit("gaus")
h.Draw()

ROOT.gStyle.SetOptFit(1)
d=ROOT.TCanvas()
hnev1.Draw("")
hnev2.Draw("same")

for fmt in "C","png","pdf":
    c.SaveAs("%s/jackknife.%s" % (dir,fmt) )
    d.SaveAs("%s/nevents.%s"   % (dir,fmt) )

