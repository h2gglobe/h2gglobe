#! /usr/bin/python

import sys, math
import ROOT

def readMus(filename):
    par = {}
    file = open(filename)
    lines = file.readlines()
    file.close()

    for l in lines:
        items = l.split("\n")[0].split()
        par[int(items[0])] = (float(items[1]), float(items[2]), float(items[3])) 
    return par

parA = readMus(sys.argv[1])
parB = readMus(sys.argv[2])

nExp = min(len(parA), len(parB))

mu = ROOT.TH1F("sigmadmu", "#delta#mu_{Average}-#delta#mu", 50, -0.1, 0.1)
dmu = ROOT.TH1F("dmu", "#delta#mu", 50, 0.28, 0.4)
mu1 = ROOT.TH1F("mu", "#mu1", 50, 1.2, 1.6)
mu2 = ROOT.TH1F("mu", "#mu2", 50, 0.9, 1.3)
mass = ROOT.TH1F("mass", "mass", 100, -0.1, 0.1)
pvalue = ROOT.TH1F("pvalue", "pvalue", 100, -0.05, 0.05)

sumMu = 0
rms = 0
for i in xrange(nExp):
    sumMu = sumMu + (parA[i][0] - parB[i][0])
    rms = rms + math.pow((parA[i][0] - parB[i][0]), 2)
    
avgMu = sumMu/float(nExp)
    
for i in xrange(389):
    mu.Fill(avgMu - (parA[i][0] - parB[i][0]))
    dmu.Fill((parA[i][0] - parB[i][0]))
    mu1.Fill(parA[i][0])
    mu2.Fill(parB[i][0])
    
c = ROOT.TCanvas("c", "c")
c.Divide(2,2)
c.cd(1)
mu1.Draw()
mu1.Fit("gaus")
c.cd(2)
mu2.Draw()
mu2.Fit("gaus")
c.cd(3)
dmu.Draw()
dmu.Fit("gaus")
c.cd(4)
mu.Draw()
mu.Fit("gaus")
#c.cd(2)
#mass.Draw()
#c.cd(3)
#pvalue.Draw()

print "Average dMu   : %1.2f" % (avgMu)
print "Average s(dMu): %1.2f" % (float(nExp-1)/math.sqrt(float(nExp))*mu.GetRMS())

raw_input()

