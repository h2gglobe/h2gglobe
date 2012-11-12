#!/usr/bin/env python

import math, ROOT

fin=open("scan.csv")

##masses=[110,115,120,125,130,135,140,145,150]
masses=[110,115,120,125,130,135,140,145,150]

h = ROOT.TH2F("deltaMuVsMass","deltaMu;m_{H};#delta #mu_{(i)} / err_{#delta #mu_{(i)}}  - #delta #mu_{n} / err_{#delta #mu_{n}};n psedo-experiments",9,107.5,152.5,65,-0.3205,0.3205)

sumDeltaMu   = {}
sumDeltaMuSq = {}

for m in masses:
    sumDeltaMu[m] = 0.
    sumDeltaMuSq[m] = 0.
n = 0
sumDeltaMu[-1] = 0.
sumDeltaMuSq[-1] = 0.

for line in fin.read().split("\n"):
    line =line.lstrip().rstrip().replace("/",",")
    if line == "":
        continue

    vals = [ float(v) for v in line.split(",") if v != "" ]

    npoints = (len(vals)-1)/4

    count = False
    for i in range(npoints):
        mass1 = vals[1+4*i]
        mu1   = vals[1+4*i+1]
        mass2 = vals[1+4*i+2]
        mu2   = vals[1+4*i+3]

        if mass1 != mass2:
            continue

        if mass1 in masses:
            count = True
            deltaMu = mu1 - mu2
            h.Fill(mass1,deltaMu)

            sumDeltaMu[mass1] += deltaMu
            sumDeltaMuSq[mass1] += deltaMu*deltaMu

            sumDeltaMu[-1] += deltaMu
            sumDeltaMuSq[-1] += deltaMu*deltaMu

    if count:
        n += 1.

avg = 0
for m in masses:
    avgDelta = sumDeltaMu[m] / n
    varJkn = (n-1)/n * ( sumDeltaMuSq[m] - n*avgDelta*avgDelta )
    avg += math.sqrt(varJkn)
    
    print m, avgDelta, math.sqrt(varJkn)

avg /= len(masses)

print avg

### ROOT.gStyle.SetOptFit(1)
### c=ROOT.TCanvas()
### h.Fit("gaus")
### h.Draw()
### 
### ROOT.gStyle.SetOptFit(1)
### d=ROOT.TCanvas()
### hnev1.Draw("")
### hnev2.Draw("same")
### 
### for fmt in "C","png","pdf":
###     c.SaveAs("jackknife.%s" % fmt )
###     d.SaveAs("nevents.%s" % fmt )

