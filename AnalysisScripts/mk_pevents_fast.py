#!/usr/bin/env python

from os import popen
from sys import argv
import ROOT

if (len(argv)<3):
    print "Must provide eos directory to run on and an output file!\n"
    print "\tExample:\n\tpython mk_pevents_fast.py /store/group/phys_higgs/cmshgg/reduced/V13_04_00_v1/mc/Summer12_S7_8TeV/GluGluToHToGG_M-115_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1/ baseline/tmp_datafiles_baseline_oct1review.dat.pevents &"
    print "\tpython mk_pevents_fast.py /store/group/phys_higgs/cmshgg/reduced/V13_04_00_v1/mc/Summer12_S7_8TeV/ baseline/tmp_datafiles_baseline_oct1review.dat.pevents \"grep HToGG\" &"
    exit(1)

eos="/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select"

if (len(argv)==3): filelist=popen(eos+" find -f "+argv[1]+" | grep .root | grep -v pileup").readlines()
else: filelist=popen(eos+" find -f "+argv[1]+" | grep .root | grep -v pileup | "+argv[3]).readlines()

outfile=open(argv[2],'w')

tfilelist = []
for i in range(len(filelist)):
    filelist[i]=filelist[i].strip("\n")
    tfilelist.append(ROOT.TFile.Open("root://eoscms/"+filelist[i]))
    global_variables = tfilelist[i].Get("global_variables")
    global_variables.GetEntry(0)
    numevents = global_variables.tot_events
    outfile.write("root://eoscms/"+filelist[i]+"="+str(numevents)+"\n")

outfile.close()
#print "pevent file is complete eos is now closing files please suspend program."
for file in tfilelist: file.Close()
