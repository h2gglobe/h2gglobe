#!/usr/bin/env python

import sys, os, re
import subprocess, time

#----------------------------------------------------------------------

def floatEquals(x,y):
    return abs(x-y) < 1e-4

def drange(start, stop, step):
    element = start
    while element < stop:
        yield element
        element=element+step

#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------
from optparse import OptionParser

parser = OptionParser("""

%prog [options] input_file1.root [ input_file2.root ... ]
%prog [options] input_file1.log [ input_file2.log ... ]

    takes the observed and expected limit from the
    per-mass point ROOT files (output of CMS' combine program)
    and produces a CSV file (by default written to stdout)
"""
)

parser.add_option(
    "--datacard",
    help="This is the combine datacard. ",
    default=None,
    type = str,
    )

parser.add_option(
    "--massstep",
    help="combine will be called every step from lowmass to highmass.",
    type=float,
    default=1,
    )

parser.add_option(
    "--lowmass",
    help="combine will be called every step from lowmass to highmass.",
    type=float,
    default=110,
    )

parser.add_option(
    "--highmass",
    help="combine will be called every step from lowmass to highmass.",
    type=float,
    default=150,
    )

parser.add_option(
    "--label",
    help="label in directory after the date/time",
    type=str,
    default="nolabel",
    )

(options, ARGV) = parser.parse_args()

#----------------------------------------

import ROOT

errFlag = False    

lowmass=options.lowmass
highmass=options.highmass + options.massstep/2.0
massstep=options.massstep
datacard=os.getcwd()+"/"+options.datacard
print datacard


massList = drange( lowmass, highmass, massstep)

thisdir=time.strftime("combine-%d%m%y-%X-")+str(options.label)

thisdir=thisdir.replace(":", "")
print thisdir

os.system("mkdir "+thisdir)

os.system("cp "+options.datacard+" "+thisdir)

os.chdir(thisdir)

for mass in massList:
  thisCmd='combine -d '+datacard+' -M Asymptotic --minosAlgo="stepping"  -m '+str(mass)+' --minimizerStrategy 1 --saveWorkspace'
  print thisCmd
  os.system(thisCmd)
  #comb = subprocess.Popen(thisCmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
      #stderr=subprocess.PIPE)
  
  #comb.wait()
  #print "comb.stdout()"+comb.stdout()

os.chdir("../")
