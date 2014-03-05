#!/usr/bin/env python

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-i","--card")
parser.add_option("-d","--dir")
#parser.add_option("-t","--ntoys",type="int")
parser.add_option("-j","--njobs",type="int")
parser.add_option("-s","--start_job",type="int",default=0)
parser.add_option("-t","--toysperjob",type="int",default=1)
parser.add_option("-f","--fqqpoints",default="0.0,0.25,0.5,0.75,1.0")
parser.add_option("-q","--queue")
parser.add_option("-v","--verbose",type="int")
parser.add_option("--skipWorkspace",dest="skipWorkspace",action="store_true",default=False)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
parser.add_option("--dataFitsOnly",action="store_true",default=False)
parser.add_option("--toysOnly",action="store_true",default=False)
parser.add_option("--runGrid",action="store_true",default=False)
(options,args)=parser.parse_args()

fqqpoints = []
model_opts = [] # [[fqq,x]]
model_opts.append([0.0,0.]) #sm
for p in options.fqqpoints.split(","):
	model_opts.append([float(p),1.])
	fqqpoints.append(float(p))

import os
import sys
import math
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

def writeScript(subfile,rootfiles,exec_line,outloc):
	subfile.write('#!/bin/bash\n')
	subfile.write('rm -f %s.done\n'%subfile.name)
	subfile.write('rm -f %s.fail\n'%subfile.name)
	subfile.write('touch %s.run\n'%subfile.name)
	subfile.write('mkdir -p scratch\n')
	subfile.write('cd scratch\n')
	subfile.write('cd %s/%s\n'%(os.getcwd(),options.dir))
	subfile.write('eval `scramv1 runtime -sh`\n')
	subfile.write('cd -\n')
	for rootfile in rootfiles:
		subfile.write('while [ ! -f %s/%s/%s ]\n'%(os.getcwd(),options.dir,rootfile))
		subfile.write('\tdo sleep 30\n')
		subfile.write('\techo \"Waiting for data best fit\"\n')
		subfile.write('done\n')
		subfile.write('cp %s/%s/%s .\n'%(os.getcwd(),options.dir,rootfile))
	subfile.write('if ( %s )\n'%exec_line)
	subfile.write('\tthen touch %s.done\n'%subfile.name)
	subfile.write('\telse touch %s.fail\n'%subfile.name)
	subfile.write('fi\n')
	subfile.write('rm -f %s.run\n'%subfile.name)
	subfile.write('mkdir -p %s/%s/%s\n'%(os.getcwd(),options.dir,outloc))
	subfile.write('mv higgsCombineJob*.root %s/%s/%s/\n'%(os.getcwd(),options.dir,outloc))

def writeScriptGrid(subfile,exec_line,outloc):
	subfile.write('#!/bin/bash\n')
	subfile.write('if [ -e %s ]; then\n'%outloc)
	subfile.write('\trm -rf %s\n'%outloc)
	subfile.write('fi\n')
	subfile.write('mkdir %s\n'%outloc)
	subfile.write('echo \"max events from CRAB: $MaxEvents\"\n')
	subfile.write('n=\"$MaxEvents\"\n')
	subfile.write('if [ \"$n\" = "" ]; then\n')
	subfile.write('\tn=\"$1\"\n')
	subfile.write('fi\n')
	subfile.write('if [ \"$n\" = "" ]; then\n')
	subfile.write('\techo \"Error: missing number of experiments"\n')
	subfile.write('\texit 2;\n')
	subfile.write('fi\n')
	subfile.write('if ( %s >& %s/%s.log )\n'%(exec_line,outloc,os.path.basename(subfile.name)))
	subfile.write('\tthen touch %s/%s.done\n'%(outloc,os.path.basename(subfile.name)))
	subfile.write('\telse touch %s/%s.fail\n'%(outloc,os.path.basename(subfile.name)))
	subfile.write('fi\n')
	subfile.write('mv higgsCombineGen*.root %s/\n'%outloc)
	subfile.write('echo \"pack the results\"\n')
	subfile.write('tar cvfz %s.tgz %s\n'%(outloc,outloc))

# __main__
if options.runGrid:
	os.system('mkdir -p %s/%s/gridJobs'%(os.getcwd(),options.dir))
	os.system('ln -s %s/%s/%s %s/%s/gridJobs/%s'%(os.getcwd(),options.dir,combfile,os.getcwd(),options.dir,combfile))
	for root, dirs, files in os.walk('%s/%s/data_fits/'%(os.getcwd(),options.dir)):
		for f in files:
			os.system('ln -s %s/%s/data_fits/%s %s/%s/gridJobs/%s'%(os.getcwd(),options.dir,f,os.getcwd(),options.dir,f))
	os.system('ln -s $CMSSW_BASE/bin/$SCRAM_ARCH/combine %s/%s/gridJobs/combine'%(os.getcwd(),options.dir))

for fqq,x in model_opts:

	outloc = 'data_fits'
	name = 'DataFit_fqq%4.2f_x%1.0f'%(fqq,x)
	bestfitname = '%s/higgsCombine%s.MultiDimFit.mH120.root'%(outloc,name)
	if not options.toysOnly:
		# require fit to data and stash in workspace
		f = open('%s/%s/sub_data_fit_fqq%4.2f_x%1.0f.sh'%(os.getcwd(),options.dir,fqq,x),'w')
		line = 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x --redefineSignalPOIs r,MH --saveWorkspace -n %s'%(combfile,fqq,x,name)
		if options.verbose: line += ' -v %d'%options.verbose
		writeScript(f,[combfile],line,outloc)
		f.close()
		os.system('chmod +x %s'%f.name)
		if not options.dryRun:
			os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))

	if not options.dataFitsOnly:
		
		if options.runGrid:
			outloc = 'outputToy_fqq%4.2f_x%1.0f'%(fqq,x)
			sub_script = open('%s/%s/gridJobs/sub_grid_fqq%4.2f_x%1.0f.sh'%(os.getcwd(),options.dir,fqq,x),'w')
			# throw toy 
			name = 'Gen_fqq%4.2f_x%1.0f'%(fqq,x)
			line = './combine %s --snapshotName MultiDimFit -M GenerateOnly --toysFrequentist --bypassFrequentistFit -t $n --expectSignal=1 --redefineSignalPOIs r -s 0 -n %s --saveToys'%(os.path.basename(bestfitname),name)
			toyname = 'higgsCombine%s.GenerateOnly.mH120.0.root'%name
			# fit toy
			line += ' \\\n && '
			if x==0: # if sm fit back with model and all other spin models
				line += './combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq%4.2f_x%1.0f -t $n --toysFile %s'%(combfile,fqq,x,name,fqq,x,toyname)
				for fqq2 in fqqpoints:
					line += ' \\\n && '
					line += './combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq%4.2f_x1 -t $n --toysFile %s'%(combfile,fqq2,name,fqq2,toyname)
			else: # otherwise its a spin model so fit with model and the sm
				line += './combine %s -M MultiDimFit --setPhysicsModelParameters fqq=0.00,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq0.00_x0 -t $n --toysFile %s'%(combfile,name,toyname)
				line += ' \\\n && '
				line += './combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq%4.2f_x%1.0f -t $n --toysFile %s'%(combfile,fqq,x,name,fqq,x,toyname)
			writeScriptGrid(sub_script,line,outloc)
			sub_script.close()
			os.system('chmod +x %s'%sub_script.name)

			grid_cfg = open('%s/%s/gridJobs/sub_grid_fqq%4.2f_x%1.0f_crab.cfg'%(os.getcwd(),options.dir,fqq,x),'w')
			grid_cfg.write('[CRAB]\n')
			grid_cfg.write('jobtype = cmssw\n')
			grid_cfg.write('scheduler = remoteGlidein\n')
			grid_cfg.write('\n')
			grid_cfg.write('[CMSSW]\n')
			grid_cfg.write('output_file = %s.tgz\n'%outloc)
			grid_cfg.write('datasetpath=None\n')
			grid_cfg.write('pset=None\n')
			grid_cfg.write('total_number_of_events=%d\n'%(options.toysperjob*options.njobs))
			grid_cfg.write('number_of_jobs=%d\n'%options.njobs)
			grid_cfg.write('\n')
			grid_cfg.write('[USER]\n')
			grid_cfg.write('ui_working_dir = crab_fqq%4.2f_x%1.0f\n'%(fqq,x))
			grid_cfg.write('script_exe = %s\n'%os.path.basename(sub_script.name))
			grid_cfg.write('additional_input_files = combine,%s,%s\n'%(combfile,os.path.basename(bestfitname)))
			grid_cfg.write('return_data = 1\n')

		else:
			for jobn in range(options.start_job,options.start_job+options.njobs):
				os.system('mkdir -p %s/%s/job%d'%(os.getcwd(),options.dir,jobn))
				sub_script = open('%s/%s/job%d/sub_job%d_fqq%4.2f_x%1.0f.sh'%(os.getcwd(),options.dir,jobn,jobn,fqq,x),'w')
				sub_script.write('#!/bin/bash\n')
				# throw toy 
				name = 'Job%d_fqq%4.2f_x%1.0f'%(jobn,fqq,x)
				line = 'combine %s --snapshotName MultiDimFit -M GenerateOnly --toysFrequentist --bypassFrequentistFit -t %d --expectSignal=1 --redefineSignalPOIs r -s 0 -n %s --saveToys'%(os.path.basename(bestfitname),options.toysperjob,name)
				toyname = 'higgsCombine%s.GenerateOnly.mH120.0.root'%name
				# fit toy
				line += ' \\\n && '
				if x==0: # if sm fit back with model and all other spin models
					line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq%4.2f_x%1.0f -t %d --toysFile %s'%(combfile,fqq,x,name,fqq,x,options.toysperjob,toyname)
					for fqq2 in fqqpoints:
						line += ' \\\n && '
						line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq%4.2f_x1 -t %d --toysFile %s'%(combfile,fqq2,name,fqq2,options.toysperjob,toyname)
				else: # otherwise its a spin model so fit with model and the sm
					line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=0.00,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq0.00_x0 -t %d --toysFile %s'%(combfile,name,options.toysperjob,toyname)
					line += ' \\\n && '
					line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %s_Fit_fqq%4.2f_x%1.0f -t %d --toysFile %s'%(combfile,fqq,x,name,fqq,x,options.toysperjob,toyname)
				
				writeScript(sub_script,[combfile,bestfitname],line,'job%d'%jobn)
				sub_script.close()
				os.system('chmod +x %s'%sub_script.name)
				if options.dryRun:
					print 'bsub -q %s -o %s.log %s'%(options.queue,sub_script.name,sub_script.name)
				else:
					os.system('bsub -q %s -o %s.log %s'%(options.queue,sub_script.name,sub_script.name))

if options.runGrid:
	print '---*** Jobs set up for running with CRAB. ***---'
	print '---*** Do the following to create and submit jobs: '
	print '\t cd %s/gridJob'%options.dir
	print '\t for f in `find sub*_crab.cfg`; do crab -create -cfg $f; done'
	print '\t for f in `find crab* -maxdepth 0 -type d`; do crab -submit -c $f; done'
	print '---*** Check job status on dashboard or with: '
	print '\t for f in `find crab* -maxdepth 0 -type d`; do crab -status -c $f; done'
	print '---*** Get output when jobs done with: '
	print '\t for f in `find crab* -maxdepth 0 -type d`; do crab -getoutput -c $f; done'


"""			
				
			#for toyn in range(options.start_toy+(jobn*options.toysperjob),options.start_toy+((jobn+1)*options.toysperjob)):
				#outloc = 'toy%d'%toyn
				#f = open('%s/%s/sub_toy%d_fqq%4.2f_x%1.0f.sh'%(os.getcwd(),options.dir,toyn,fqq,x),'w')
				for toyn in range(options.toysperjob):
					name = 'Job%d_fqq%4.2f_x%1.0f'%(toyn,fqq,x)
					line = 'combine %s --snapshotName MultiDimFit -M GenerateOnly --toysFrequentist --bypassFrequentistFit -t %d --expectSignal=1 --redefineSignalPOIs r -s 0 -n %s --saveToys'%(os.path.basename(bestfitname),toyn,name)
					if options.verbose: line += ' -v %d'%options.verbose
					line += ' \\\n && '
					toyname = 'higgsCombine%s.GenerateOnly.mH120.0.root'%name
					# if sm fit back with model and all other spin models
					if x==0:
						line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %sFit_fqq%4.2f_x%1.0f -t 1 --toysFile %s'%(combfile,fqq,x,name,fqq,x,toyname)
						if options.verbose: line += ' -v %d'%options.verbose
					for fqq2 in fqqpoints:
						line += ' \\\n && '
						line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %sFit_fqq%4.2f_x1 -t 1 --toysFile %s'%(combfile,fqq2,name,fqq2,toyname)
						if options.verbose: line += ' -v %d'%options.verbose
					writeScript(f,[combfile,bestfitname],line,outloc)
				# otherwise its a spin model so fit with model and the sm
				else:
					line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=0.00,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %sFit_fqq0.00_x0 -t 1 --toysFile %s'%(combfile,name,toyname)
					if options.verbose: line += ' -v %d'%options.verbose
					line += ' \\\n && '
					line += 'combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x --redefineSignalPOIs r,MH -n %sFit_fqq%4.2f_x%1.0f -t 1 --toysFile %s'%(combfile,fqq,x,name,fqq,x,toyname)
					if options.verbose: line += ' -v %d'%options.verbose
					writeScript(f,[combfile,bestfitname],line,outloc)
					
				f.close()
				os.system('chmod +x %s'%f.name)
				sub_script.write('bash %s >& %s.log\n'%(f.name,f.name))
			sub_script.close()
			os.system('chmod +x %s'%sub_script.name)
			if not options.dryRun:
				os.system('bsub -q %s -o %s.log %s'%(options.queue,sub_script.name,sub_script.name))
"""
