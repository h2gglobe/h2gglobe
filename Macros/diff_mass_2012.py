#!/usr/bin/env python

import sys
from sys import argv
import json 
from pprint import pprint 
from math import fabs,sqrt, acos
import operator,time
from numpy import uint32

import ROOT
ROOT.gROOT.SetBatch(True)

minrun=999999
maxrun=0

lastrun=999999

doVertex=False

isMVA=False


class EventInfo():
    def __init__(self): 
        self.info={}
    
    def setup(self,line):
        self.info["line"]=str(line)
        for lineitem in line.split():
            if lineitem.find(":") != -1:
                name=lineitem.split(":")[0]
                value=lineitem.split(":")[1]
                
                try: 
                    self.info[name]=float(value)
                except:
                    self.info[name]=value
            
     
def getEventKey(line):
    run=-999
    lumi=-999
    event=-999
    for lineitem in line.split():
        if lineitem.find("run") != -1:
            run= int(lineitem.split(":")[1])
        elif lineitem.find("lumi") != -1:
            lumi= int(lineitem.split(":")[1])
        elif lineitem.find("event") != -1:
            event= uint32(lineitem.split(":")[1])
        if run!=-999:
            if lumi!=-999:
                if event!=-999:
                    break
    if run!=-999:
        global minrun, maxrun
        if int(run)<minrun:
            minrun=int(run)
        if int(run)>maxrun:
            maxrun=int(run)

    return run,lumi,event 


def checkForMVA(line):
    if line.find("diphoBDT") != -1:
        global isMVA
        isMVA=True

def getUnderOverFlows(hist):
    nbins = h.GetNbinsX()
    hist.Fill(h.GetBinCenter(1),hist.GetBinContent(0))
    hist.Fill(h.GetBinCenter(nbins),hist.GetBinContent(nbins+1))
    hist.SetBinContent(0,0)
    hist.SetBinContent(nbins+1,0)


def bookHisto(histos,name,nbins,min,max,relative=False):
    ymin = -max*0.2
    ymax = max*0.2
    nybins = 2*nbins
    den = ""
    if relative:
        ymin = -0.5
        ymax = 0.5
        nybins = 100
        den = "/ %s" % name
    if name.find("mgg")!=-1 or  name.find("run")!=-1:
        #h = ROOT.TProfile(name, "%s" % ( name ) ,nbins,min,max)
        h = ROOT.TH2F(name, "%s" % ( name ) ,nbins,min,max,nybins,-0.03,0.03)
    elif name.find("eta")!=-1:
        h = ROOT.TH2F(name, "%s;%s; (%s_{A} - %s_{B}) %s" % ( name, name, name, name, den ),nbins,min,max,nybins,-3,3.)
    else:
        h = ROOT.TH2F(name, "%s;%s; (%s_{A} - %s_{B}) %s" % ( name, name, name, name, den ),nbins,min,max,nybins,ymin,ymax)
    h.relative = relative
    h.unsynch = []
    histos[name] = h
    return h

def fillHisto(h, varA, varB, ev, relative=False, otherxval=False, xval=1.):
    y=(varA-varB)
    if relative and varA != 0.:
        y /= varA
    if not otherxval:
        h.Fill(varA, y )
    else:
        h.Fill(xval, y )

    if abs(y) > 0.03:
        h.unsynch.append(ev)

histosdiffmassdiffpho = {}    
histosdiffmass = {}
histossamemass = {}
histosallmass = {}
histosmass = {}

t0=time.time()

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(111111)


fn1 = argv.pop(1)
fn2 = argv.pop(1)
        
file1 = open(fn1)
file2 = open(fn2)

events1={}
events2={}

print "starting list1"
counter=0
for line in file1.readlines():
    run,lumi,event=getEventKey(line)
    if run!=-999:
        checkForMVA(line)
        neweventinfo=EventInfo()
        neweventinfo.setup(line)
        events1[run,lumi,event]=neweventinfo
        counter=counter+1
        if counter%20000 == 0:
            print counter,time.time()-t0

print "starting list2"
counter=0
for line in file2.readlines():
    run,lumi,event=getEventKey(line)
    if run !=-999:
        neweventinfo=EventInfo()
        neweventinfo.setup(line)
        events2[run,lumi,event]=neweventinfo
        counter=counter+1
        if counter%20000 == 0:
            print counter,time.time()-t0

events1keys=events1.keys()
events2keys=events2.keys()

print "starting common keys"
print time.time()-t0
commonkeys = set(events1keys).intersection(set(events2keys))
print time.time()-t0

print "common, list1, list2",len(commonkeys),len(events1keys),len(events2keys)

# events in ICHEP only
only1 = set(events1keys) -  set(events2keys)
# events in A+B only
only2 = set(events2keys) -  set(events1keys)

#file_out = open("onlyICHEP.txt", "w")

r9 = []
sceta = []
hoe = []
sigieie = []
ecaliso = []
hcaliso = []
trckiso = []
chpfis2 = []
chpfis3 = []
phoid = []
phoeta = []
ptom = []
e = []
eerr = []
eerrsmeared = []

for i in xrange(2):
    r9.append(ROOT.TH1F("r9"+str(i), "r9"+str(i), 100, 0, 1))
    sceta.append(ROOT.TH1F("eta"+str(i), "eta"+str(i), 100, -2.5, 2.5))
    hoe.append(ROOT.TH1F("hoe"+str(i), "hoe"+str(i), 100, 0, .1))
    sigieie.append(ROOT.TH1F("sigieie"+str(i), "sigieie"+str(i), 100, 0, 0.05))
    ecaliso.append(ROOT.TH1F("ecaliso"+str(i), "ecaliso"+str(i), 100, 0, 10))
    hcaliso.append(ROOT.TH1F("hcaliso"+str(i), "hcaliso"+str(i), 100, 0, 10))
    trckiso.append(ROOT.TH1F("trckiso"+str(i), "trckiso"+str(i), 100, 0, 10))
    chpfis2.append(ROOT.TH1F("chpfiso2"+str(i), "chpfiso2"+str(i), 100, 0, 10))
    chpfis3.append(ROOT.TH1F("chpfiso3"+str(i), "chpfiso3"+str(i), 100, 0, 10))
    phoid.append(ROOT.TH1F("phoid"+str(i), "phoid"+str(i), 100, -1, 1))
    phoeta.append(ROOT.TH1F("phoeta"+str(i), "phoeta"+str(i), 100, -2.5, 2.5))
    ptom.append(ROOT.TH1F("ptom"+str(i), "ptom"+str(i), 100, 0, 1))
    e.append(ROOT.TH1F("e"+str(i), "e"+str(i), 100, 0, 100))
    eerr.append(ROOT.TH1F("eerr"+str(i), "eerr"+str(i), 100, 0, .5))
    eerrsmeared.append(ROOT.TH1F("eerrsmeared"+str(i), "eerrsmeared"+str(i), 100, 0, .5))

sigmrv = ROOT.TH1F("sigmarv"+str(i), "sigmarv"+str(i), 100, 0, 0.06)
sigmwv = ROOT.TH1F("sigmawv"+str(i), "sigmawv"+str(i), 100, 0, 0.06)
vtxprob = ROOT.TH1F("vtxprob"+str(i), "vtxprob"+str(i), 100, 0, 1)
cosdphi = ROOT.TH1F("cosdphi"+str(i), "r9"+str(i), 100, -1, 1)
diphoBDT = ROOT.TH1F("diphoBDT"+str(i), "diphoBDT"+str(i), 100, -1, 1)
mgg = ROOT.TH1F("mgg"+str(i), "mgg"+str(i), 160, 100, 180)

counter=0
for eventkey in only1:
    setICHEPOnly = events1[eventkey].info
    if (setICHEPOnly['diphoBDT'] > -0.05):
        #file_out.write(str(eventkey)+"\n")
        r9[0].Fill(setICHEPOnly['r9_1'])
        sceta[0].Fill(setICHEPOnly['sceta_1'])
        hoe[0].Fill(setICHEPOnly['hoe_1'])
        sigieie[0].Fill(setICHEPOnly['sigieie_1'])
        ecaliso[0].Fill(setICHEPOnly['ecaliso_1'])
        hcaliso[0].Fill(setICHEPOnly['hcaliso_1'])
        trckiso[0].Fill(setICHEPOnly['trckiso_1'])
        chpfis2[0].Fill(setICHEPOnly['chpfis2_1'])
        chpfis3[0].Fill(setICHEPOnly['chpfis3_1'])
        phoid[0].Fill(setICHEPOnly['phoid_1'])
        phoeta[0].Fill(setICHEPOnly['phoeta_1'])
        ptom[0].Fill(setICHEPOnly['pt_1/m'])
        e[0].Fill(setICHEPOnly['e_1'])
        eerr[0].Fill(setICHEPOnly['eerr_1'])
        eerrsmeared[0].Fill(setICHEPOnly['eerrsmeared_1'])

        r9[1].Fill(setICHEPOnly['r9_2'])
        sceta[1].Fill(setICHEPOnly['sceta_2'])
        hoe[1].Fill(setICHEPOnly['hoe_2'])
        sigieie[1].Fill(setICHEPOnly['sigieie_2'])
        ecaliso[1].Fill(setICHEPOnly['ecaliso_2'])
        hcaliso[1].Fill(setICHEPOnly['hcaliso_2'])
        trckiso[1].Fill(setICHEPOnly['trckiso_2'])
        chpfis2[1].Fill(setICHEPOnly['chpfis2_2'])
        chpfis3[1].Fill(setICHEPOnly['chpfiso3_2'])
        phoid[1].Fill(setICHEPOnly['phoid_2'])
        phoeta[1].Fill(setICHEPOnly['phoeta_2'])
        ptom[1].Fill(setICHEPOnly['pt_2/m'])
        e[1].Fill(setICHEPOnly['e_2'])
        eerr[1].Fill(setICHEPOnly['eerr_2'])
        eerrsmeared[1].Fill(setICHEPOnly['eerrsmeared_2'])

        sigmrv.Fill(setICHEPOnly['sigmrv'])
        sigmwv.Fill(setICHEPOnly['sigmwv'])
        vtxprob.Fill(setICHEPOnly['vtxprob'])
        cosdphi.Fill(setICHEPOnly['cosdphi'])
        diphoBDT.Fill(setICHEPOnly['diphoBDT'])
        mgg.Fill(setICHEPOnly['mgg'])
        counter+=1
#file_out.close()



for i in xrange(2):
    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    r9[i].Draw()
    c.SaveAs("r9"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    sceta[i].Draw()
    c.SaveAs("eta"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    hoe[i].Draw()
    c.SaveAs("hoe"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    sigieie[i].Draw()
    c.SaveAs("sigieie"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    ecaliso[i].Draw()
    c.SaveAs("ecaliso"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    hcaliso[i].Draw()
    c.SaveAs("hcaliso"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    trckiso[i].Draw()
    c.SaveAs("trckiso"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    chpfis2[i].Draw()
    c.SaveAs("chpfiso2"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    chpfis3[i].Draw()
    c.SaveAs("chpfiso3"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    phoid[i].Draw()
    c.SaveAs("phoid"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    phoeta[i].Draw()
    c.SaveAs("phoeta"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    ptom[i].Draw()
    c.SaveAs("ptom"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    e[i].Draw()
    c.SaveAs("e"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    eerr[i].Draw()
    c.SaveAs("eerr"+str(i)+".png", "png")

    c = ROOT.TCanvas("cr9"+str(i))
    c.cd()
    eerrsmeared[i].Draw()
    c.SaveAs("eerrsmeared"+str(i)+".png", "png")


c = ROOT.TCanvas("cr9"+str(i))
c.cd()
sigmrv.Draw()
c.SaveAs("sigmrv.png", "png")
c = ROOT.TCanvas("cr9"+str(i))
c.cd()
sigmwv.Draw()
c.SaveAs("sigmwv.png", "png")
c = ROOT.TCanvas("cr9"+str(i))
c.cd()
vtxprob.Draw()
c.SaveAs("vtxprob.png", "png")
c = ROOT.TCanvas("cr9"+str(i))
c.cd()
cosdphi.Draw()
c.SaveAs("cosdphi.png", "png")
c = ROOT.TCanvas("cr9"+str(i))
c.cd()
diphoBDT.Draw()
c.SaveAs("diphoBDT.png", "png")
c = ROOT.TCanvas("cr9"+str(i))
c.cd()
mgg.Draw()
c.SaveAs("mgg.png", "png")


print "Only ICHEP %d" % counter
print "Only A+B %d" % len(only2)

nruns = (maxrun - minrun) / 100

print isMVA

if isMVA:
    etaname = "sceta_"
else:
    etaname = "scEta"

print "setup time",time.time()-t0

if isMVA:
    bookHisto(histosdiffmass,"diff0mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histossamemass,"same0mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histosdiffmass,"diff1mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histossamemass,"same1mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histosdiffmass,"diff2mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histossamemass,"same2mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histosdiffmass,"diff3mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histossamemass,"same3mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histosdiffmassdiffpho,"diffdiffpho0mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histosdiffmassdiffpho,"diffdiffpho1mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histosdiffmassdiffpho,"diffdiffpho2mass_diphoBDT"    ,51,-0.12,1.02,True)
    bookHisto(histosdiffmassdiffpho,"diffdiffpho3mass_diphoBDT"    ,51,-0.12,1.02,True)
    
bookHisto(histosmass,"allmass_mgg"   ,160,100,180,True)
bookHisto(histosallmass,"all0mass_mgg"   ,160,100,180,True)
bookHisto(histosallmass,"all1mass_mgg"   ,160,100,180,True)
bookHisto(histosallmass,"all2mass_mgg"   ,160,100,180,True)
bookHisto(histosallmass,"all3mass_mgg"   ,160,100,180,True)

bookHisto(histosmass,"allmass_run"   ,nruns,minrun,maxrun,True)
bookHisto(histosallmass,"all0mass_run"   ,nruns,minrun,maxrun,True)
bookHisto(histosallmass,"all1mass_run"   ,nruns,minrun,maxrun,True)
bookHisto(histosallmass,"all2mass_run"   ,nruns,minrun,maxrun,True)
bookHisto(histosallmass,"all3mass_run"   ,nruns,minrun,maxrun,True)

bookHisto(histosdiffmass,"same0mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmass,"same0mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmass,"same0mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmass,"same0mass_r9_1"     ,51,0.938,1.02,True)
bookHisto(histosdiffmass,"same0mass_r9_2"     ,51,0.938,1.02,True)
bookHisto(histosdiffmass,"same0mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same0mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same0mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmass,"same0mass_sigmrv"      ,120,0,5,True)
bookHisto(histosdiffmass,"same0mass_sigmwv"      ,120,0,8,True)
bookHisto(histosdiffmass,"same0mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmass,"same0mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmass,"same1mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmass,"same1mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmass,"same1mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmass,"same1mass_r9_1"     ,51,0.3,0.942,True)
bookHisto(histosdiffmass,"same1mass_r9_2"     ,51,0.3,0.942,True)
bookHisto(histosdiffmass,"same1mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same1mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same1mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmass,"same1mass_sigmrv"      ,120,0,5,True)
bookHisto(histosdiffmass,"same1mass_sigmwv"      ,120,0,8,True)
bookHisto(histosdiffmass,"same1mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmass,"same1mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmass,"same2mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmass,"same2mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmass,"same2mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmass,"same2mass_r9_1"     ,51,0.938,1.02,True)
bookHisto(histosdiffmass,"same2mass_r9_2"     ,51,0.938,1.02,True)
bookHisto(histosdiffmass,"same2mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same2mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same2mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmass,"same2mass_sigmrv"      ,120,0,10,True)
bookHisto(histosdiffmass,"same2mass_sigmwv"      ,120,0,10,True)
bookHisto(histosdiffmass,"same2mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmass,"same2mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmass,"same3mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmass,"same3mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmass,"same3mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmass,"same3mass_r9_1"     ,51,0.5,1.02,True)
bookHisto(histosdiffmass,"same3mass_r9_2"     ,51,0.5,1.02,True)
bookHisto(histosdiffmass,"same3mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same3mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmass,"same3mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmass,"same3mass_sigmrv"      ,120,0,10,True)
bookHisto(histosdiffmass,"same3mass_sigmwv"      ,120,0,10,True)
bookHisto(histosdiffmass,"same3mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmass,"same3mass_%s2" % etaname      ,50,-3,3,False)

bookHisto(histosdiffmassdiffpho,"swapped0mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_r9_1"     ,51,0.938,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_r9_2"     ,51,0.938,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_sigmrv"      ,120,0,5,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_sigmwv"      ,120,0,8,True)
bookHisto(histosdiffmassdiffpho,"swapped0mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmassdiffpho,"swapped0mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmassdiffpho,"swapped1mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_r9_1"     ,51,0.3,0.942,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_r9_2"     ,51,0.3,0.942,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_sigmrv"      ,120,0,5,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_sigmwv"      ,120,0,8,True)
bookHisto(histosdiffmassdiffpho,"swapped1mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmassdiffpho,"swapped1mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmassdiffpho,"swapped2mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_r9_1"     ,51,0.938,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_r9_2"     ,51,0.938,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_sigmrv"      ,120,0,10,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_sigmwv"      ,120,0,10,True)
bookHisto(histosdiffmassdiffpho,"swapped2mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmassdiffpho,"swapped2mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmassdiffpho,"swapped3mass_mgg"   ,160,100,180,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_e_1"   ,140,30,170,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_e_2"   ,120,20,140,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_r9_1"     ,51,0.5,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_r9_2"     ,51,0.5,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_sigmrv"      ,120,0,10,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_sigmwv"      ,120,0,10,True)
bookHisto(histosdiffmassdiffpho,"swapped3mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histosdiffmassdiffpho,"swapped3mass_%s2" % etaname      ,50,-3,3,False)

bookHisto(histossamemass,"diff0mass_mgg"   ,160,100,180,True)
bookHisto(histossamemass,"diff0mass_e_1"    ,140,30,170,True)
bookHisto(histossamemass,"diff0mass_e_2"    ,120,20,140,True)
bookHisto(histossamemass,"diff0mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff0mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff0mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histossamemass,"diff0mass_sigmrv"      ,120,0,5,True)
bookHisto(histossamemass,"diff0mass_sigmwv"      ,120,0,8,True)
bookHisto(histossamemass,"diff0mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histossamemass,"diff0mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histossamemass,"diff1mass_mgg"   ,160,100,180,True)
bookHisto(histossamemass,"diff1mass_e_1"    ,140,30,170,True)
bookHisto(histossamemass,"diff1mass_e_2"    ,120,20,140,True)
bookHisto(histossamemass,"diff1mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff1mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff1mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histossamemass,"diff1mass_sigmrv"      ,120,0,6,True)
bookHisto(histossamemass,"diff1mass_sigmwv"      ,120,0,8,True)
bookHisto(histossamemass,"diff1mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histossamemass,"diff1mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histossamemass,"diff2mass_mgg"   ,160,100,180,True)
bookHisto(histossamemass,"diff2mass_e_1"    ,140,30,170,True)
bookHisto(histossamemass,"diff2mass_e_2"    ,120,20,140,True)
bookHisto(histossamemass,"diff2mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff2mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff2mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histossamemass,"diff2mass_sigmrv"      ,120,0,6,True)
bookHisto(histossamemass,"diff2mass_sigmwv"      ,120,0,10,True)
bookHisto(histossamemass,"diff2mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histossamemass,"diff2mass_%s2" % etaname      ,50,-3,3,False)
bookHisto(histossamemass,"diff3mass_mgg"   ,160,100,180,True)
bookHisto(histossamemass,"diff3mass_e_1"    ,140,30,170,True)
bookHisto(histossamemass,"diff3mass_e_2"    ,120,20,140,True)
bookHisto(histossamemass,"diff3mass_phoid_1"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff3mass_phoid_2"     ,51,-0.32,0.52,True)
bookHisto(histossamemass,"diff3mass_vtxprob"     ,26,-0.02,1.02,True)
bookHisto(histossamemass,"diff3mass_sigmrv"      ,120,0,6,True)
bookHisto(histossamemass,"diff3mass_sigmwv"      ,120,0,10,True)
bookHisto(histossamemass,"diff3mass_%s1" % etaname      ,50,-3,3,False)
bookHisto(histossamemass,"diff3mass_%s2" % etaname      ,50,-3,3,False)

counter=0
counter2=0
for eventkey in commonkeys:
    if counter%10000:
        print "common events",counter,time.time()-t0

    setA = events1[eventkey].info
    setB = events2[eventkey].info
        

    if len(setA) == 0 or len(setB) == 0:
        continue

    #setA["%s1" % etaname] = abs( setA["%s1" % etaname] )
    #setB["%s1" % etaname] = abs( setB["%s1" % etaname] )
    #setA["%s2" % etaname] = abs( setA["%s2" % etaname] )
    #setB["%s2" % etaname] = abs( setB["%s2" % etaname] )
   
    if abs( setA["%s1" % etaname] ) < 1.5 and  abs( setA["%s2" % etaname] ) < 1.5: 
        if setA["r9_1"] > 0.94 and setA["r9_2"] > 0.94:
            setA["dipho_cat"]=0
        else:
            setA["dipho_cat"]=1
    else:
        if setA["r9_1"] > 0.94 and setA["r9_2"] > 0.94:
            setA["dipho_cat"]=2
        else:
            setA["dipho_cat"]=3
    
    if abs( setB["%s1" % etaname] ) < 1.5 and  abs( setB["%s2" % etaname] ) < 1.5: 
        if setB["r9_1"] > 0.94 and setB["r9_2"] > 0.94:
            setB["dipho_cat"]=0
        else:
            setB["dipho_cat"]=1
    else:
        if setB["r9_1"] > 0.94 and setB["r9_2"] > 0.94:
            setB["dipho_cat"]=2
        else:
            setB["dipho_cat"]=3
    
    try:
        for name,hist in histosmass.iteritems():
            try:
                varname=name.split("mass_")[1]
                if name.find("run")!=-1:
                    fillHisto( hist, setA["mgg"],  setB["mgg"], eventkey, hist.relative, True, setA[varname] )
                else:
                    fillHisto( hist, setA[varname], setB[varname], eventkey, hist.relative )
            except Exception, e:
                print "fail in mass",e,eventkey
                pass
        
        for name,hist in histosallmass.iteritems():
            varname=name.split("mass_")[1]
            try:
                if setA["dipho_cat"]== int(name.split("mass_")[0][3]):
                    if varname.find("run")!=-1:
                        fillHisto( hist, setA["mgg"],  setB["mgg"], eventkey, hist.relative, True, setA[varname] )
                    else:
                        fillHisto( hist, setA[varname], setB[varname], eventkey, hist.relative )
            except Exception, e:
                print "fail in allmass",e,eventkey
                pass
    except Exception, e:
        print "fail in try",e,eventkey
        print "line A",setA["line"]
        print "line B",setB["line"]
        continue

    #sqrtdiffB1=sqrt(fabs(setB["eerr_1"]*setB["eerr_1"]-setB["eerrsmeared_1"]*setB["eerrsmeared_1"]))/setB["e_1"]    
    #histoswithcuts["sqrt_err2_minuserrsmear2_B"].Fill( setB["sceta_1"], sqrtdiffB1 )
        
    try:
        if (abs( 1. - setA["mgg"] / setB["mgg"]) < 0.03):
            if ((fabs(setA["sceta_1"] - setB["sceta_1"]) < 0.05) and (fabs(setA["sceta_2"] - setB["sceta_2"]) < 0.05) and (fabs(setA['cosdphi'] - setB['cosdphi']) < 0.05)):
                for name,hist in histosdiffmass.iteritems():
                        varname=name.split("mass_")[1]
                        try:
                            if setA["dipho_cat"]== int(name.split("mass_")[0][4]):
                                fillHisto( hist, setA[varname], setB[varname], eventkey, hist.relative )
                        except Exception, e:
                            print "fail in diffmass",e,eventkey
                            pass
            elif ((fabs(setA["sceta_2"] - setB["sceta_1"]) < 0.05) and (fabs(setA["sceta_1"] - setB["sceta_2"]) < 0.05) and (fabs(setA['cosdphi'] - setB['cosdphi']) < 0.05)):
                for name,hist in histosdiffmassdiffpho.iteritems():
                        varname=name.split("mass_")[1]
                        try:
                            #print name, name.split("mass_")[0][-1]
                            varname2 = varname
                            if (varname[-1] == "1"):
                                varname2 = varname[:-1]+"2"
                            if (varname[-1] == "2"):
                                varname2 = varname[:-1]+"1"
                                
                            if setA["dipho_cat"]== int(name.split("mass_")[0][-1]):
                                fillHisto( hist, setA[varname], setB[varname2], eventkey, hist.relative )
                        except Exception, e:
                            print "fail in diffmassdiffpho",e,eventkey
                            pass
            else:
                counter2+=1
        
        else:
            for name,hist in histossamemass.iteritems():
                varname=name.split("mass_")[1]
                try:
                    if setA["dipho_cat"]== int(name.split("mass_")[0][4]):
                        fillHisto( hist, setA[varname], setB[varname], eventkey, hist.relative )
                except Exception, e:
                    print "fail in samemass",e,eventkey
                    pass
        
    except Exception, e:
        print "fail in try",e,eventkey
        print "line A",setA["line"]
        print "line B",setB["line"]
        continue

                   
print "histos time",time.time()-t0
 
listsToDump = ["ptbal", "ptasym","sceta_1","sceta_2","p2conv"]
allhist={}
for hist in histosdiffmassdiffpho.keys():
    allhist[hist]=histosdiffmassdiffpho[hist]
for hist in histosdiffmass.keys():
    allhist[hist]=histosdiffmass[hist]
for hist in histossamemass.keys():
    allhist[hist]=histossamemass[hist]
for hist in histosallmass.keys():
    allhist[hist]=histosallmass[hist]
for hist in histosmass.keys():
    allhist[hist]=histosmass[hist]
#for name,h in allhist.iteritems():
names = allhist.keys()
names.sort()
for name in names:
    h=allhist[name]
    c = ROOT.TCanvas ( name )
    d = None
    c.SetGridx()
    c.SetGridy()
    c.SetLogz()
    c.cd()
    if h.IsA().InheritsFrom(ROOT.TProfile.Class()):
        print "TProfile?"
        h.Draw()
            
    elif h.IsA().InheritsFrom(ROOT.TH2.Class()):
        h.Draw("colz")
        d = ROOT.TCanvas ( "%s_proj" % name )
        if not "vertex" in name:
            d.SetLogy()
        prj = h.ProjectionY()
        getUnderOverFlows( prj )
        prj.DrawNormalized()
        sync = prj.Integral( prj.FindBin(-3.e-3), prj.FindBin(3.e-3 ) )
        all = prj.GetEntries()
        if all > 0:
            print "%s synch : %d / %d = %1.4f RMS = %1.4f" % ( name, sync, all, sync/all, prj.GetRMS() )
            if name in listsToDump:
                f = open("plotsmass/%s.txt" % name, "w+")
                for ev in h.unsynch:
                    f.write("run:%d\tlumi:%d\tevent:%d\n" % ev )
                f.close()
    else:
        h.Draw()
    
    c.SaveAs( "plotsmass/%s.png" % name )
    if d:
        d.SaveAs( "plotsmass/%s_proj.png" % name )
        
print "all time",time.time()-t0

print "%s %d" % ( fn1, len(events1) )
print "%s %d" % ( fn2, len(events2) )

print "Common %d" % len(commonkeys)

print "Different photons", counter2
### pprint(events2)
