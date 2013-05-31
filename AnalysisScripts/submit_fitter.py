#!/usr/bin/env python
import sys,os,commands,glob
from optparse import OptionParser

myjobsN = []
def cback(option,opt_str,value,parser):
	value = value.split(",")
	for v in value:
		if v!="": myjobsN.append(int(v))
			
parser = OptionParser()
parser.add_option("","--resubFailed",action="store_true",dest="resubFailed",default=False)
parser.add_option("","--submitMissing",action="store_true",dest="submitMissing",default=False)
parser.add_option("-q","--queue",dest="queue",default="1nh")
parser.add_option("-d","--directory",dest="directory")
parser.add_option("-j","--jobs",dest="jobs",action="callback",callback=cback,type='string')
parser.add_option("-l","--label",dest="label",default="sub")
parser.add_option("","--runIC",dest="runIC",default=False,action="store_true")

(options,args)=parser.parse_args()

workingDir=os.getcwd()

jobs = glob.glob( "%s/%s*.sh" % (options.directory,options.label ) )

subjobs = []
for n in myjobsN:
	subjobs.append("%s/%s%d.sh" % (options.directory,options.label,n ))
if len(subjobs)>0 :jobs = subjobs[:]

for j in jobs:
	os.system("chmod 775 %s"%j)
	if not options.submitMissing and options.resubFailed and not os.path.isfile("%s.fail"%j): continue
	elif not options.submitMissing or not os.path.isfile("%s.done"%j):
	   os.system("rm -f %s.run"%j)
	   os.system("rm -f %s.fail"%j)
	   os.system("rm -f %s.done"%j)
	   os.system("rm -f %s.log"%j)
	   if options.runIC:
	   	print "Submitting Job: qsub -q %s -o %s/%s.log -e %s/%s.err %s "%(options.queue,workingDir,j,workingDir,j,j)
   	   	os.system("qsub -q %s -o %s/%s.log -e %s/%s.err %s "%(options.queue,workingDir,j,workingDir,j,j))
	   else:
	     print "Submitting Job: bsub -q %s -o %s.log $PWD/%s" %(options.queue,j,j)
   	     os.system("bsub -q %s -o %s.log $PWD/%s "%(options.queue,j,j))
