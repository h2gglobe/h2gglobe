#!/usr/bin/env python

from sys import argv
import json 
from pprint import pprint 
from math import fabs
import operator

import ROOT

minrun=999999
maxrun=0

def getlist(input,minmass=0,maxmass=1000.):
    lst = {}
    global minrun, maxrun
    print minmass, maxmass
    for line in input.split("\n"):
        l = [ i for i in line.replace("="," ").split(" ")  if i != "" ]
        try:
            if len(l) < 9 and not line.startswith("Passing"):
                continue
            if line.startswith("Passing"):
                run, ls, ev, vtx, cat, mass, pt = int(l[1]), int(l[2]), int(l[3]), 0, 0, float(l[5]), float(l[7])
            elif "VBF" in line:
                ## Run = 177449  LS = 402  Event = 606455702  SelVtx = 0  CAT4 = 0  ggM = 149.3 ggPt =  189.6  jetEta1 = 1.218  jetEta2 = -2.775  jetPhi1 = 0.7197  jetPhi2 = -0.5363  jetEt1 = 138.6  jetEt2 = 33.11 Mjj 496 dEtajj 3.993 Zeppenfeld 1.785 dPhijjgg 2.973 VBF itype 0
                run, ls, ev, vtx, cat, mass, pt = int(l[1]), int(l[3]), int(l[5]), int(l[7]), 4, float(l[11]), float(l[13])
            elif len(l)==14:
                type, run, ls, ev, mass, cat, vtx, pt = int(l[1]), int(l[3]), int(l[5]), int(l[7]), float(l[9]), int(l[11]), int(l[13]), 0.
            elif len(l)>13:
                type, run, ls, ev, vtx, cat, mass, pt = int(l[1]), int(l[3]), int(l[5]), int(l[7]), int(l[9]), int(l[11]), float(l[13]), float(l[15])
                if type != 0:
                    continue
            elif len(l)==12:
                type, run, ls, ev, vtx, cat, mass, pt = int(l[1]), int(l[3]), int(l[5]), int(l[7]), int(l[9]), 0, float(l[11]), 0.
            elif len(l)==10:
                type, run, ls, ev, vtx, cat, mass, pt = int(l[1]), int(l[3]), int(l[5]), int(l[7]), 0, 0, float(l[9]), 0.
            else:
                run, ls, ev, vtx, cat, mass, pt = int(l[1]), int(l[3]), int(l[5]), int(l[7]), int(l[9]), float(l[11]), float(l[13])
        except Exception, e:
            print line
            print e
        

        if mass >= minmass and mass <= maxmass and run>1:
            if run<minrun:
                minrun=run
            if run>maxrun:
                maxrun=run
            ## cat = cat % 4
            ### if pt < 40:
            ###     cat = cat + 4
            lst[ (run, ls, ev) ] = ( vtx, cat, mass, pt )
    return lst

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(111111)

fn1 = argv.pop(1)
fn2 = argv.pop(1)
        
file1 = open(fn1)
file2 = open(fn2)

## list1 = getlist( file1.read(),100,180 )
list1 = getlist( file1.read() )
## print list1

## list2 = getlist( file2.read(),100,180 )
list2 = getlist( file2.read() )
## print list2

def cmp(x,y):
    if x[0] < y[0] : return True
    elif x[0] == y[0]:
        if x[1] < y[1] : return True
        elif x[1] == y[1]:
            if x[2] < y[2] : return True
                
    return False

events1 = list1.keys()
events2 = list2.keys()

events1.sort( cmp )
events2.sort( cmp )

common = set(events1).intersection(  set(events2) )

only1 = list(set(events1) -  set(events2))
only2 = list(set(events2) -  set(events1))

only1.sort()
only2.sort()

## only1.sort(cmp) 
## only2.sort(cmp) 

diffvtx  = []
diffmass = []
diffpt   = []
diffcat  = []

diffmassplots = []
diffmassrel = []
diffmassprofiles = []
diffmassprofiles_mass = []
diffmass_scatter = []
nruns = (maxrun - minrun) / 100
## ncat = 8
ncat = 5
print nruns, maxrun, minrun
for cat in range(0,ncat):
    diffmassplots.append( ROOT.TH1F("diffmass_cat%d" % cat ,"diffmass_cat%d; #Delta m_{#gamma #gamma} (GeV/c^{2}); Events / bin" %cat , 100, -10., 10. ) )
    diffmassrel.append( ROOT.TH1F("diffrel_cat%d" % cat ,"diffrel_cat%d; #Delta m_{#gamma #gamma} / m; Events / bin" %cat , 100, -0.2, 0.2 ) )
    diffmassprofiles.append( ROOT.TProfile("diffrel_run_cat%d" % cat ,"diffmass_run_cat%d; run number; #Delta m_{#gamma #gamma} / m; Events / bin" %cat,
                                           nruns, minrun, maxrun) ) 
    diffmassprofiles_mass.append( ROOT.TProfile("diffrel_mass_cat%d" % cat ,"diffmass_mass_cat%d; m_{#gamma #gamma} (GeV/c^{2}); #Delta m_{#gamma #gamma} / m; Events / bin" %cat,
                                           80, 100., 180.) ) 
    diffmass_scatter.append( ROOT.TH2F("scatter_mass_cat%d" % cat ,"scatter_mass_cat%d; m^{1}_{#gamma #gamma} (GeV/c^{2}); m^{2}_{#gamma #gamma} (GeV/c^{2}); Events / bin" %cat,
                                           80, 100., 180., 80, 100., 180.) ) 

diffmassplots.append( ROOT.TH1F("diffmass_other","diffmass_other; #Delta m_{#gamma #gamma} (GeV/c^{2}); Events / bin", 100, -10., 10. ) )
diffmassrel.append( ROOT.TH1F("diffrel_other","diffrel_other; #Delta m_{#gamma #gamma} / m; Events / bin", 100, -0.2, 0.2 ) )
diffmassprofiles.append( ROOT.TProfile("diffrel_run_other","diffmass_run_other; run number; #Delta m_{#gamma #gamma} / m",
                                           nruns, minrun, maxrun) )
diffmassprofiles_mass.append( ROOT.TProfile("diffrel_mass_other", "diffmass_mass_other; m_{#gamma #gamma} (GeV/c^{2}); #Delta m_{#gamma #gamma} / m; Events / bin",
                                           80, 100., 180.) )
diffmass_scatter.append( ROOT.TH2F("scatter_mass_other","scatter_mass_other; m^{1}_{#gamma #gamma} (GeV/c^{2}); m^{2}_{#gamma #gamma} (GeV/c^{2}); Events / bin",
                                           80, 100., 180., 80, 100., 180.) ) 


for ev in common:
    vtx1, cat1, mass1, pt1 = list1[ ev ]
    vtx2, cat2, mass2, pt2 = list2[ ev ]

    if vtx1 != vtx2:
        diffvtx.append( (ev, list1[ ev ], list2[ ev ])  )
    if fabs(1. - mass1/mass2) > 0.003:
        diffmass.append( (ev, list1[ ev ], list2[ ev ]) )
    ## if fabs(1. - pt1/pt2) > 0.05:
    ##     diffpt.append( (ev, list1[ ev ], list2[ ev ])  )
    if cat1 != cat2:
        diffcat.append( (ev, list1[ ev ], list2[ ev ]) )

    delta_mass = mass1 - mass2
    diffmassplots[ cat1 ].Fill( delta_mass )
    diffmassrel[ cat1 ].Fill( delta_mass/mass1 )
    diffmassprofiles[ cat1 ].Fill( ev[0], delta_mass/mass1 )
    diffmassprofiles_mass[ cat1 ].Fill( mass1, delta_mass/mass1 )
    diffmass_scatter[ cat1 ].Fill( mass1, mass2 )
    ## if cat1 == cat2:
    ##     diffmassplots[ cat1 ].Fill( delta_mass )
    ##     diffmassrel[ cat1 ].Fill( delta_mass/mass1 )
    ##     diffmassprofiles[ cat1 ].Fill( ev[0], delta_mass/mass1 )
    ##     diffmassprofiles_mass[ cat1 ].Fill( mass1, delta_mass/mass1 )
    ## else:
    ##     diffmassplots[ ncat ].Fill( delta_mass )
    ##     diffmassrel[ ncat ].Fill( delta_mass/mass1 )
    ##     diffmassprofiles[ ncat ].Fill( ev[0], delta_mass/mass1 )
    ##     diffmassprofiles_mass[ ncat ].Fill( mass1, delta_mass/mass1 )
    
print "Only in %s %d" % ( fn1, len(only1) )
## for run,lumi,event in only1:
##     print "%d %d %d" % ( run,lumi,event )
print "Only in %s %d" % ( fn2, len(only2) )
## for run,lumi,event in only2:
##     print "%d %d %d" % ( run,lumi,event )
## pprint(only2)

fonly1 = open( "only_%s" % fn1.rsplit("/",1)[1], "w+"  )
fonly2 = open( "only_%s" % fn2.rsplit("/",1)[1], "w+"  )
## print >> fonly1, json.dumps( list(only1) )
## print >> fonly2, json.dumps( list(only2) )

for ev in only1:
    print >> fonly1, "%d %d %d %1.3f" % (ev[0], ev[1], ev[2], list1[ev][2] )

for ev in only2:
    print >> fonly2, "%d %d %d %1.3f" % (ev[0], ev[1], ev[2], list2[ev][2] )

fonly1.close()
fonly2.close()

print "Common %d" % len(common)

print "Different vertex %d" % len(diffvtx)
## pprint(diffvtx)

print "Different mass %d" % len(diffmass)
## pprint(diffmass)

print "Different pt %d" % len(diffpt)
## pprint(diffpt)

print "Different category %d" % len(diffcat)
## pprint(diffpt)

for p in diffmass_scatter+diffmassrel+diffmassplots+diffmassprofiles+diffmassprofiles_mass:
## for p in diffmass_scatter:
    c = ROOT.TCanvas ( p.GetName() )
    c.SetGridx()
    c.SetGridy()
    c.cd()
    if p.IsA().InheritsFrom(ROOT.TH2.Class()):
        p.Draw("colz")
    else:
        p.Draw()
    c.SaveAs( "%s.png" % p.GetName() )
    
