#!/usr/bin/env python
import sys,os,commands,glob
from optparse import OptionParser

myjobsN = []
def cback(option,opt_str,value,parser):
	value = value.split(",")
	for v in value: myjobsN.append(int(v))

parser = OptionParser()
parser.add_option("","--resubFailed",action="store_true",dest="resubFailed",default=False)
parser.add_option("-q","--queue",dest="queue",default="1nh")
parser.add_option("-d","--directory",dest="directory")
parser.add_option("-j","--jobs",dest="jobs",action="callback",callback=cback,type='string')
parser.add_option("-l","--label",dest="label",default="sub")

(options,args)=parser.parse_args()

jobs = glob.glob( "%s/%s*.sh" % (options.directory,options.label ) )

subjobs = []
for n in myjobsN:
	subjobs.append("%s/%s%d.sh" % (options.directory,options.label,n ))
if len(subjobs)>0 :jobs = subjobs[:]

for j in jobs:
	os.system("chmod 775 %s"%j)
	if options.resubFailed and os.path.isfile("%s.done"%j): continue
	else:
	   os.system("rm %s.run"%j)
	   os.system("rm %s.fail"%j)
	   os.system("rm %s.done"%j)
	   os.system("rm %s.log"%j)
	   print "Submitting Job: "bsub -q %s -o %s.log $PWD/%s "%(options.queue,j,j)"
   	   os.system("bsub -q %s -o %s.log $PWD/%s "%(options.queue,j,j))
