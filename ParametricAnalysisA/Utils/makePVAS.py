#!/usr/bin/env python

import sys, os, re
import shutil, subprocess, time

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

python makePVAS.py  --datacard=""

    This prog takes a functional datacard and makes
    asymptotic pvalues between 110 and 150 in steps of
    1 GeV unless these values are changed by the user.
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

thisdir=time.strftime("combine-pvalue-%d%m%y-%X-")+str(options.label)

thisdir=thisdir.replace(":", "")
print thisdir

os.mkdir(thisdir)

shutil.copy(options.datacard, thisdir)

os.chdir(thisdir)

for mass in massList:
  logFname = "mass-%.1f.log" % mass
  cmdParts=[
    'combine',
    '-d '+datacard,
    '--signif', 
    #'--signif --pvalue',
    '-M ProfileLikelihood',
    '-t 0',
    '-m '+str(mass),
    '2>&1 | cat >> %s &' % logFname,
  ] 


  thisCmd=" ".join(cmdParts)
  
  os.system(thisCmd)
  #comb = subprocess.Popen(thisCmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
      #stderr=subprocess.PIPE)
  
  #comb.wait()
  #print "comb.stdout()"+comb.stdout()

os.chdir("../")
