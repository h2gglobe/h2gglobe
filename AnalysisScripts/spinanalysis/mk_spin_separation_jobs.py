#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--card")
parser.add_option("-d","--dir")
parser.add_option("-f","--fqqpoints",default="0.0,0.25,0.5,0.75,1.0")
parser.add_option("-n","--njobs",type="int")
parser.add_option("-t","--toysperjob",type="int")
parser.add_option("-q","--queue")
parser.add_option("--isGrid",action="store_true",default=False)
parser.add_option("--skipWorkspace",dest="skipWorkspace",action="store_true",default=False)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
(options,args)=parser.parse_args()

fqqpoints = [] 
for p in options.fqqpoints.split(","):
	fqqpoints.append(float(p))

import os
import sys
os.system('mkdir -p %s'%options.dir)
combfile = options.card.replace('.txt','.root')

if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/combine')):
	sys.exit('ERROR - CombinedLimit package must be installed')
if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/text2workspace.py')):
	sys.exit('ERROR - CombinedLimit package must be installed')
if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/combineCards.py')):
	sys.exit('ERROR - CombinedLimit package must be installed')

if not options.skipWorkspace:
	print 'Using combine tool to make workspace from card', options.card
	os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating --PO=muFloating --PO higgsMassRange=122,128 -o %s/%s'%(options.card,options.dir,combfile))

def writeForBatch()

def writeForGrid()


for fqq in fqqpoints:
	for job in range(options.njobs):
		f = open('%s/%s/sub_fqq%4.2f_job%d.sh'%(os.getcwd(),options.dir,fqq,job),'w')
		f.write('#!/bin/bash\n')
		if not options.isGrid:
			f.write('rm -f %s.done\n'%f.name)
			f.write('rm -f %s.fail\n'%f.name)
			f.write('touch %s.run\n'%f.name)
			f.write('mkdir scratch\n')
			f.write('cd scratch\n')
			f.write('cd %s/%s\n'%(os.getcwd(),options.dir))
			f.write('eval `scramv1 runtime -sh`\n')
			f.write('cd -\n')
			f.write('cp %s/%s/%s .\n'%(os.getcwd(),options.dir,combfile))
		f.write('if ( combine %s -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --redefineSignalPOIs r --setPhysicsModelParameters fqq=%4.2f --freezeNuisances fqq --cminDefaultMinimizerType Minuit2 -s 0 -n fqq%4.2fjob%d'%(combfile,options.toysperjob,fqq,fqq,job))
		if options.isGrid: f.write(' >& log.txt')
		f.write(' )\n')
		if options.isGrid:
			f.write('\tthen echo OK\n')
		else:
			f.write('\tthen touch %s.done\n'%f.name)
			f.write('\telse touch %s.fail\n'%f.name)
		f.write('fi\n')
		if options.isGrid:
			outFile = 'outputFile_fqq%4.2f_job%d'%(fqq,job)
			f.write('mkdir %s\n'%outFile) 
			f.write('mv higgsCombinefqq%4.2fjob%d.HybridNew.mH120.root %s/\n'%(fqq,job,outFile))
			f.write('mv log.txt %s/\n'%(outFile))
			f.write('tar cvfz %s.tgz %s/\n'%(outFile,outFile))
		else:
			f.write('rm -f %s.run\n'%f.name)
			f.write('cp higgsCombinefqq%4.2fjob%d.HybridNew.mH120.root %s/%s/\n'%(fqq,job,os.getcwd(),options.dir))
		f.close()
		os.system('chmod +x %s'%f.name)
		if not options.dryRun and not options.isGrid:
			os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))

		if options.isGrid:
			# write crab cfg
			crabF = open('%s'%f.name.replace('.sh','.cfg'),'w')
			crabF.write('[CRAB]\n')
			crabF.write('jobtype = cmssw\n')
			crabF.write('scheduler = glidein\n')
			crabF.write('\n')
			crabF.write('[CMSSW]\n')
			crabF.write('output_file = %s.tgz\n'%outFile)
			crabF.write('datasetpath=None\n')
			crabF.write('pset=None\n')
			crabF.write('total_number_of_events=1\n')
			crabF.write('number_of_jobs=1\n')
			crabF.write('\n')
			crabF.write('[USER]\n')
			crabF.write('script_exe=%s\n'%f.name)
			crabF.write('additional_input_files=combine,%s\n'%combfile)
			crabF.write('return_data=1\n')
			crabF.close()

