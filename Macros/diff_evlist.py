#!/usr/bin/env python

# ############################################################################
#
# The script is developed for purposes of synchronizing H->gammagamma analysis
# performed in different frameworks. It takes as input two txt files
# containing event information in format <field_name>:<field_value>, 
# and outputs the following information: 
#
# 1) event counts and printouts of common and unique events
# 2) same as 1) per event category
# 3) level of synchrnoization per <field_name>
# 4) plots of differences per <field_name>
# 
#
# Original version: P.Musella
# Improved version: V.Rekovic  Jan/2013 
# (added synchronization per category)
#
#
# #############################################################################

from sys import argv
import json 
from pprint import pprint 
from math import fabs
import operator
import array

import ROOT

minrun=99999
maxrun=0

Commoncat = [0] * 10
CommonOnly1cat = [0] * 10
CommonOnly2cat = [0] * 10
Only1cat = [0] * 10
Only2cat = [0] * 10

lastrun=9999999999999

def getlist(input):
    lst = {}
    lstDuplEvents = {}
    lstRuns = {}

    global minrun, maxrun

    for line in input.split("\n"):
        vars = {}
        try:
            for i in line.replace(": ",":").replace("="," ").replace("\t"," ").split(" "):
                if i != "":
                    j = i.split(":")
                    #if j[0] == "event":
                        #globals()[j[0]] = (int(j[1])) & 0xFFFFFFFF
												#continue
                    if j[0] == "run" or j[0] == "lumi" or j[0] == "event" or j[0] == "type":
                        globals()[j[0]] = abs(int(j[1]))
                    else:
                         vars[j[0]] = float(j[1])
                    
            #if run > lastrun :
            #if run > lastrun or run == 200961 or run == 200976 or run == 201191 :
            #if run > lastrun or run == 201191 :
            if run > lastrun :
                continue
            
            if run<minrun:
                minrun=run
            if run>maxrun:
                maxrun=run
                
            if(lst.has_key((run, lumi, event))):
                #print run, lumi, event
                #print "newVars = ", "mgg:", vars["mgg"]
                #print "oldVars = ", "mgg:", lst[(run,lumi,event)]["mgg"]
                lstDuplEvents[ (run, lumi, event) ] = vars
                continue

            lst[  (run, lumi, event) ] = vars
            lstRuns[  (run) ] = vars
            
        except Exception, e:
            #print line
            print e
            
    evts = lst.keys()
    duplEvts = lstDuplEvents.keys()
    runs = lstRuns.keys()

    print " Number of unique evts    : %d" % ( len(evts) )
    print " Number of duplicated evts: %d" % ( len(duplEvts)-1 )
    print " Number of unique runs: %d" % ( len(runs) )
    return lst

def bookHisto(name,nbins,min,max,relative=False):
    ymin = -max*0.2
    ymax = max*0.2
    nybins = 2*nbins
    den = ""
    if relative:
        ymin = -0.01
        ymax = 0.01
        nybins = 100
        den = "/ %s" % name
    h = ROOT.TH2F(name, "%s;%s; (%s_{A} - %s_{B}) %s" % ( name, name, name, name, den ),nbins,min,max,nybins,ymin,ymax)
    h.relative = relative
    return h

def fillHisto(h, varA, varB, relative=False):
    y=(varA-varB)
    if relative and varA != 0.:
        y /= varA
    h.Fill(varA, y )

    
histos = {
    "mgg"                     : bookHisto("mgg",                          1000,   100,   180,  True),
    #"run"                     : bookHisto("run",                          1000,   000,   000,  True),
    #"event"                   : bookHisto("event",                        1000,   000,   000,  True),
    "rho"                     : bookHisto("rho",                          1000,   0.0,   0.6,  True),
    "pho1_ind"                : bookHisto("pho1_ind",                      100,     0,   100,  True),
    "pho1_scInd"              : bookHisto("pho1_scInd",                    100,     0,   100,  True),
    "pho1_r9"                 : bookHisto("pho1_r9",                      1000,     0,   1.0,  True),
    "pho1_scEta"              : bookHisto("pho1_scEta",                   1000,   -3.,   3.0,  True),
    "pho1_pt"                 : bookHisto("pho1_pt",                      1000,     0,   300,  True),
    "pho1_eta"                : bookHisto("pho1_eta",                     1000,   -3.,   3.0,  True),
    "pho1_phi"                : bookHisto("pho1_phi",                     1000, -3.15,  3.15,  True),
    "pho1_e"                  : bookHisto("pho1_e",                       1000,     0,  2000,  True),
    "pho1_eErr"               : bookHisto("pho1_eErr",                    1000,     0,   300,  True),
    "pho1_isConv"             : bookHisto("pho1_isConv",                     2,     0,     2,  True),
    "pho1_HoE"                : bookHisto("pho1_HoE",                      100,     0,   0.3,  True),
    "pho1_hcalIso03"          : bookHisto("pho1_hcalIso03",               1000,     0,    40,  True),
    "pho1_trkIso03"           : bookHisto("pho1_trkIso03",                1000,     0,    40,  True),
    "pho1_pfChargedIsoGood02" : bookHisto("pho1_pfChargedIsoGood02",      1000,     0,    40,  True),
    "pho1_pfChargedIsoGood03" : bookHisto("pho1_pfChargedIsoGood03",      1000,     0,    40,  True),
    "pho1_pfChargedIsoBad03"  : bookHisto("pho1_pfChargedIsoBad03",       1000,     0,    40,  True),
    "pho1_pfPhotonIso03"      : bookHisto("pho1_pfPhotonIso03",           1000,     0,    40,  True),
    "pho1_pfNeutralIso03"     : bookHisto("pho1_pfNeutralIso03",          1000,     0,    40,  True),
    "pho1_sieie"              : bookHisto("pho1_sieie",                   1000,     0,   0.5,  True),
    "pho1_sieip"              : bookHisto("pho1_sieip",                   1000,     0,   0.5,  True),
    "pho1_etaWidth"           : bookHisto("pho1_etaWidth",                1000,     0,   0.5,  True),
    "pho1_phiWidth"           : bookHisto("pho1_phiWidth",                1000,     0,   0.5,  True),
    "pho1_lambdaRatio"        : bookHisto("pho1_lambdaRatio",             1000,    -1,     1,  True),
    "pho1_s4Ratio"            : bookHisto("pho1_s4Ratio",                 1000,    -1,     1,  True),
    "pho1_ESEffSigmaRR"       : bookHisto("pho1_ESEffSigmaRR",            1000,     0,   1.0,  True),
    "pho2_ind"                : bookHisto("pho2_ind",                      100,     0,   100,  True),
    "pho2_scInd"              : bookHisto("pho2_scInd",                    100,     0,   100,  True),
    "pho2_r9"                 : bookHisto("pho2_r9",                      1000,     0,   1.0,  True),
    "pho2_scEta"              : bookHisto("pho2_scEta",                   1000,   -3.,   3.0,  True),
    "pho2_pt"                 : bookHisto("pho2_pt",                      1000,     0,   300,  True),
    "pho2_eta"                : bookHisto("pho2_eta",                     1000,   -3.,   3.0,  True),
    "pho2_phi"                : bookHisto("pho2_phi",                     1000, -3.15,  3.15,  True),
    "pho2_e"                  : bookHisto("pho2_e",                       1000,     0,  2000,  True),
    "pho2_eErr"               : bookHisto("pho2_eErr",                    1000,     0,   300,  True),
    "pho2_isConv"             : bookHisto("pho2_isConv",                     2,     0,     2,  True),
    "pho2_HoE"                : bookHisto("pho2_HoE",                      100,     0,   0.3,  True),
    "pho2_hcalIso03"          : bookHisto("pho2_hcalIso03",               1000,     0,    40,  True),
    "pho2_trkIso03"           : bookHisto("pho2_trkIso03",                1000,     0,    40,  True),
    "pho2_pfChargedIsoGood02" : bookHisto("pho2_pfChargedIsoGood02",      1000,     0,    40,  True),
    "pho2_pfChargedIsoGood03" : bookHisto("pho2_pfChargedIsoGood03",      1000,     0,    40,  True),
    "pho2_pfChargedIsoBad03"  : bookHisto("pho2_pfChargedIsoBad03",       1000,     0,    40,  True),
    "pho2_pfPhotonIso03"      : bookHisto("pho2_pfPhotonIso03",           1000,     0,    40,  True),
    "pho2_pfNeutralIso03"     : bookHisto("pho2_pfNeutralIso03",          1000,     0,    40,  True),
    "pho2_sieie"              : bookHisto("pho2_sieie",                   1000,     0,   0.5,  True),
    "pho2_sieip"              : bookHisto("pho2_sieip",                   1000,     0,   0.5,  True),
    "pho2_etaWidth"           : bookHisto("pho2_etaWidth",                1000,     0,   0.5,  True),
    "pho2_phiWidth"           : bookHisto("pho2_phiWidth",                1000,     0,   0.5,  True),
    "pho2_lambdaRatio"        : bookHisto("pho2_lambdaRatio",             1000,    -1,     1,  True),
    "pho2_s4Ratio"            : bookHisto("pho2_s4Ratio",                 1000,    -1,     1,  True),
    "pho2_ESEffSigmaRR"       : bookHisto("pho2_ESEffSigmaRR",            1000,     0,   1.0,  True),
    "mass"                    : bookHisto("mass",                         1800,   100,   180,  True),
    "rVtxSigmaMoM"            : bookHisto("rVtxSigmaMoM",                 1000,     0,   0.5,  True),
    "wVtxSigmaMoM"            : bookHisto("wVtxSigmaMoM",                 1000,     0,   0.5,  True),
    "pho1_ptOverM"            : bookHisto("pho1_ptOverM",                 1000,     0,   1.0,  True),
    "pho2_ptOverM"            : bookHisto("pho2_ptOverM",                 1000,     0,   1.0,  True),
    "vtxIndex"                : bookHisto("vtxIndex",                      100,     0,   100,  True),
    "vtxProb"                 : bookHisto("vtxProb",                      1000,     0,   1.0,  True),
    "cosDPhi"                 : bookHisto("cosDPhi",                      1000,   -1.,   1.0,  True),
    "ptBal"                   : bookHisto("ptBal",                        1000,   -31,    31,  True),
    "ptAsym"                  : bookHisto("ptAsym",                       1000,   -31,    31,  True),
    "logSPt2"                 : bookHisto("logSPt2",                      1000,    -3,     3,  True),
    "p2Conv"                  : bookHisto("p2Conv",                       1000,   -31,    31,  True),
    "nConv"                   : bookHisto("nConv",                          50,     0,    50,  True),
    "jet1_ind"                : bookHisto("jet1_ind",                      500,     0,   500,  True),
    "jet1_eta"                : bookHisto("jet1_eta",                     1000,   -6.,     6,  True),
    "jet1_pt"                 : bookHisto("jet1_pt",                      1000,     0,  2000,  True),
    "jet2_ind"                : bookHisto("jet2_ind",                      500,     0,   500,  True),
    "jet2_eta"                : bookHisto("jet2_eta",                     1000,    -6,     6,  True),
    "jet2_pt"                 : bookHisto("jet2_pt",                      1000,     0,  2000,  True),
    "dijet_dEta"              : bookHisto("dijet_dEta",                   1000,     0,    15,  True),
    "dijet_Zep"               : bookHisto("dijet_Zep",                    1000,     0,    10,  True),
    "dijet_dPhi"              : bookHisto("dijet_dPhi",                   1000,     0,  3.15,  True),
    "dijet_Mjj"               : bookHisto("dijet_Mjj",                    1000,     0,  3000,  True)
}

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(111111)

fn1 = argv.pop(1)
fn2 = argv.pop(1)
        
file1 = open(fn1)
file2 = open(fn2)

print "\n"
print "list1: %s" % fn1
print "=========================="
list1 = getlist( file1.read() )
print "\n"
print "list2: %s" % fn2
print "=========================="
list2 = getlist( file2.read() )
print "\n"

events1 = list1.keys()
events2 = list2.keys()

common = set(events1).intersection(  set(events2) )

only1 = list(set(events1) -  set(events2))
only2 = list(set(events2) -  set(events1))

only1.sort()
only2.sort()

nruns = (maxrun - minrun) / 100

for ev in common:

    setA = list1[ev]
    setB = list2[ev]

    if len(setA) == 0 or len(setB) == 0:
        continue

    if abs( 1. - setA["pho1_scEta"] / setB["pho1_scEta"] ) > 1.e-3 or abs( 1. - setA["pho2_scEta"] / setB["pho2_scEta"] ) > 1.e-3:
            fillHisto( histos["pho1_scEta"], setA["pho1_scEta"], setB["pho1_scEta"], True )
            fillHisto( histos["pho2_scEta"], setA["pho2_scEta"], setB["pho2_scEta"], True )
    else:
        for name,hist in histos.iteritems():
            try:
                fillHisto( hist, setB[name], setA[name], hist.relative )
            except Exception, e:
                pass
        
    #if setA["evcat"] != setB["evcat"]:
		#	#print "setA = ", setA
		#	#print "setB = ", setB
		#	#print "\n"
		#	print "(event, lumi, run)", ev
		#	print "setA[evcat] =", setA["evcat"]," setB[evcat] =", setB["evcat"]
		#	print "setA[diphoBDT] =", setA["diphoBDT"]," setB[diphoBDT] =", setB["diphoBDT"]
		#	for name in setB:
		#		if name == "evcat":
		#			continue
		#		try:
		#			delta=setA[name]-setB[name]
		#			if(delta != 0):
		#				print "delta(",name,") =", delta
		#		except Exception, e:
		#			pass

    #for i in range(0,9): 
    #    
		#	if setA["evcat"] == i and setB["evcat"] == i:
		#		Commoncat[i] = Commoncat[i] + 1
		#	if setA["evcat"] == i and setB["evcat"] != i:
		#		CommonOnly1cat[i] = CommonOnly1cat[i] + 1
		#	if setA["evcat"] != i and setB["evcat"] == i:
		#		CommonOnly2cat[i] = CommonOnly2cat[i] + 1


#for ev in only1:
#
#  setA = list1[ev]
#  print "Only1 " , ev , " diphoBTD:" , setA["diphoBDT"], "  mgg:" , setA["mgg"] 
#  for i in range(0,9):
#    if setA["evcat"] == i:
#      Only1cat[i] = Only1cat[i] + 1
#
#for ev in only2:
#
#  setB = list2[ev]
#  print "Only2 " , ev , " diphoBTD:" , setB["diphoBDT"], "  mgg:" , setB["mgg"] 
#  for i in range(0,9):
#    if setB["evcat"] == i:
#      Only2cat[i] = Only2cat[i] + 1



#print "\n"
#print "%s %d" % ( fn1, len(events1) )
#print "%s %d" % ( fn2, len(events2) )
print "\n"
print "Common %d" % len(common)
#print common
print "============================"
print "Only1 %d" % len(only1)
print "============================"
#print only1

#print "============================"
print "Only2 %d" % len(only2)
print "============================"
#print only2

### pprint(events2)

print "\n"
#print "------------------------------------------------------------------"
print "CAT      ", repr("ANY").rjust(5),
for i in range(0,9):
   #print " ", i,
   print repr(i).rjust(5),
print "   "
print "============================================================================"
print "Common   ",
print repr(len(common)).rjust(5),
for i in range(0,9):
   #print " ", Commoncat[i],
   print repr(Commoncat[i]).rjust(5),
print "   "
#print "Common1  ",
print "Common1  ", 
print repr(len(common)).rjust(5),
for i in range(0,9):
   print repr(CommonOnly1cat[i]).rjust(5),
print "   "
#print "Common2  ",
print "Common2  ", 
print repr(len(common)).rjust(5),
for i in range(0,9):
   #print " ", CommonOnly2cat[i],
   print repr(CommonOnly2cat[i]).rjust(5),
print "   "
print "============================================================================"
print "Only1    ",
print repr(len(only1)).rjust(5),
for i in range(0,9):
   #print " ", Only1cat[i],
   print repr(Only1cat[i]).rjust(5),
print "   "
print "Only2    ",
print repr(len(only2)).rjust(5),
for i in range(0,9):
   #print " ", Only2cat[i],
   print repr(Only2cat[i]).rjust(5),
print "   "
print "============================================================================"

print "\n"
print "Synchronization w/in .3% of the value"
print "============================================================================"

names=histos.keys()
names.sort()
for name in names:
    h = histos[name]
    c = ROOT.TCanvas ( name )
    d = None
    c.SetGridx()
    c.SetGridy()
    c.SetLogz()
    c.cd()
    if h.IsA().InheritsFrom(ROOT.TH2.Class()):
        #h.Draw("colz")
        d = ROOT.TCanvas ( "%s_proj" % name )
        #if not "vertex" in name:
            #d.SetLogy()
        prj = h.ProjectionY()
        #prj.DrawNormalized()
        sync = prj.Integral( prj.FindBin(-3.e-3), prj.FindBin(3.e-3 ) )
        all = prj.GetEntries()
        if all > 0:
            print "%s synch : %d / %d = %1.4f RMS = %1.4f" % ( name, sync, all, sync/all, prj.GetRMS() )
            if sync/all < 0.98:
                print name,"needs to be synced!!!!!!!!!!!"
        #h.Draw()
        #c.SaveAs( "plots/%s.pdf" % name )
    #else:
        #h.Draw()
        #c.SaveAs( "plots/%s.pdf" % name )
    #if d:
        #d.Draw()
        #d.SaveAs( "plots/%s_proj.pdf" % name )
        
