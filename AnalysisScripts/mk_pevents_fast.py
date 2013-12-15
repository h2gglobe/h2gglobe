#!/usr/bin/env python

from os import popen
from sys import argv

if (len(argv)==1):
    print "Must provide either a configuation dat or an eos directory and an output file!\n"
    print "\tExamples:"
    print "\tpython mk_pevents_fast.py baseline/datafiles_baseline_oct1review.dat"
    print "\tpython mk_pevents_fast.py baseline/tmp_datafiles_baseline_oct1review.dat.pevents /store/group/phys_higgs/cmshgg/reduced/V13_04_00_v1/mc/Summer12_S7_8TeV/GluGluToHToGG_M-115_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1/"
    print "\tpython mk_pevents_fast.py baseline/tmp_datafiles_baseline_oct1review.dat.pevents /store/group/phys_higgs/cmshgg/reduced/V13_04_00_v1/mc/Summer12_S7_8TeV/ \"grep HToGG\""
    exit(1)

import ROOT
eos="/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select"

if len(argv)==2: filelist=popen("grep CaDir "+argv[1]+" | grep -v typ=0 | grep -v '#' | sed 's|.*CaDir=\\(\\S*\\).*|\\1|' | xargs -i "+eos+" find -f {} | grep .root").readlines()
if len(argv)==3: filelist=popen(eos+" find -f "+argv[2]+" | grep .root ").readlines()
if len(argv)==4: filelist=popen(eos+" find -f "+argv[2]+" | grep .root | "+argv[3]).readlines()

if len(argv)==2:
    if argv[1].rfind("/")==-1:outfile=open("tmp_"+argv[1]+".pevents",'w')
    else: outfile=open(argv[1][:argv[1].rfind("/")+1]+"tmp_"+argv[1][argv[1].rfind("/")+1:]+".pevents",'w')
else: outfile=open(argv[1],'w')


tfilelist = []
for i in range(len(filelist)):
#    print filelist[i]
    filelist[i]=filelist[i].strip("\n")
    try:
        tfilelist.append(ROOT.TFile.Open("root://eoscms/"+filelist[i]))
        if tfilelist[i].Get("global_variables")!=None:
            global_variables = tfilelist[i].Get("global_variables")
            global_variables.GetEntry(0)
            numevents = global_variables.processedEvents
            outfile.write("root://eoscms/"+filelist[i]+"="+str(numevents)+"\n")
    except:
        print "failed to open",filelist[i]
outfile.close()
for file in tfilelist: file.Close()
