#!/usr/bin/env python
import sys,os,commands,glob
from optparse import OptionParser
parser = OptionParser()
parser.add_option("","--resubFailed",action="store_true",dest="resubFailed",default=False)
parser.add_option("-q","--queue",dest="queue",default="1nh")
parser.add_option("-f","--taskdir",dest="taskdir")

(options,args)=parser.parse_args()

jobs = glob.glob( "%s/sub*.sh" % options.taskdir )

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
