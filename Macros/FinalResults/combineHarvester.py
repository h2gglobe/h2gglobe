#!/usr/bin/env python
# vim: ai sw=2 ts=2 mouse=a number

import os
import numpy
import sys
import fnmatch

from optparse import OptionParser
from optparse import OptionGroup

parser = OptionParser()
parser.add_option("-d","--datfile",help="Pick up running options from datfile")
parser.add_option("-q","--queue",help="Which batch queue")
parser.add_option("--dryRun",default=False,action="store_true",help="Dont submit")
parser.add_option("--runLocal",default=False,action="store_true",help="Run locally")
parser.add_option("--skipWorkspace",default=False,action="store_true",help="Dont remake MultiDim workspace")
parser.add_option("--hadd",help="Trawl passed directory and hadd files. To be used when jobs are complete.")
parser.add_option("-v","--verbose",default=False,action="store_true")
#parser.add_option("--blindStd",default=False,action="store_true",help="Run standard suite of blind plots")
#parser.add_option("--unblindSimple",default=False,action="store_true",help="Run simple set of unblind plots (limit, pval, best fit mu)")
#parser.add_option("--unblindFull",default=False,action="store_true",help="Run full suite of unblind plots")
specOpts = OptionGroup(parser,"Specific options")
specOpts.add_option("--datacard")
specOpts.add_option("--files")
specOpts.add_option("--outDir")
specOpts.add_option("--method")
specOpts.add_option("--expected",type="int")
specOpts.add_option("--mhLow",type="float")
specOpts.add_option("--mhHigh",type="float")
specOpts.add_option("--mhStep",type="float")
specOpts.add_option("--muLow",type="float")
specOpts.add_option("--muHigh",type="float")
specOpts.add_option("--rvLow",type="float")
specOpts.add_option("--rvHigh",type="float")
specOpts.add_option("--rfLow",type="float")
specOpts.add_option("--rfHigh",type="float")
specOpts.add_option("--cvLow",type="float")
specOpts.add_option("--cvHigh",type="float")
specOpts.add_option("--cfLow",type="float")
specOpts.add_option("--cfHigh",type="float")
specOpts.add_option("--jobs",type="int")
specOpts.add_option("--pointsperjob",type="int",default=1)
specOpts.add_option("--expectSignal",type="float")
specOpts.add_option("--expectSignalMass",type="float")
specOpts.add_option("--splitChannels")
specOpts.add_option("--additionalOptions")
parser.add_option_group(specOpts)
(opts,args) = parser.parse_args()

if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/combine')):
	sys.exit('ERROR - CombinedLimit package must be installed')
if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/text2workspace.py')):
	sys.exit('ERROR - CombinedLimit package must be installed')
if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/combineCards.py')):
	sys.exit('ERROR - CombinedLimit package must be installed')

cwd = os.getcwd()
allowedMethods = ['Asymptotic','AsymptoticGrid','ProfileLikelihood','ChannelCompatibilityCheck','MHScan','MuScan','RVScan','RFScan','RVRFScan','MuMHScan']

def checkValidMethod():
	if opts.method not in allowedMethods: sys.exit('%s is not a valid method'%opts.method)

def configureMassFromNJobs():
	if opts.mhLow and opts.mhHigh and opts.mhStep:
		masses = numpy.arange(opts.mhLow,opts.mhHigh+opts.mhStep,opts.mhStep)
		if len(masses)<opts.jobs: sys.exit("Can't have more masses than number of jobs")
		else:
			opts.masses_per_job = [[] for x in range(opts.jobs)]
			while len(masses)!=0:
				for j in range(opts.jobs):
					if len(masses)==0: break
					opts.masses_per_job[j].append(masses[0])
					masses = numpy.delete(masses,0)
		if len(opts.masses_per_job)!=opts.jobs: sys.exit('ERROR - len job config (%d) not equal to njobs (%d)'%(len(opts.masses_per_job),opts.jobs))

def splitCard():
	if not opts.splitChannels: sys.exit('Channel splitting options not specified')
	f = open(opts.datacard)
	allCats = set()
	for line in f.readlines():
		if line.startswith('bin'):
			for el in line.split()[1:]:
				allCats.add(el)
	f.close()
	veto = ""
	for cat in allCats:
		if cat in opts.splitChannels: continue
		else: veto += "|ch1_"+cat
	veto=veto[1:]
	splitCardName = opts.datacard.replace('.txt','')
	for cat in opts.splitChannels: splitCardName += '_'+cat
	splitCardName += '.txt'
	os.system('combineCards.py --xc="%s" %s > %s'%(veto,opts.datacard,splitCardName))
	opts.datacard = splitCardName

def writePreamble(sub_file):
	sub_file.write('#!/bin/bash\n')
	sub_file.write('touch %s.run\n'%os.path.abspath(sub_file.name))
	sub_file.write('cd %s\n'%os.getcwd())
	sub_file.write('eval `scramv1 runtime -sh`\n')
	sub_file.write('cd -\n')
	sub_file.write('mkdir -p scratch\n')
	sub_file.write('cd scratch\n')
	sub_file.write('cp -p $CMSSW_BASE/bin/$SCRAM_ARCH/combine .\n')
	sub_file.write('cp -p %s .\n'%os.path.abspath(opts.datacard))
	for file in opts.files.split(','):
		sub_file.write('cp -p %s .\n'%os.path.abspath(file))

def writePostamble(sub_file, exec_line):
	sub_file.write('if ( %s ) then\n'%exec_line)
	sub_file.write('\t mv higgsCombine*.root %s\n'%os.path.abspath(opts.outDir))
	sub_file.write('\t touch %s.done\n'%os.path.abspath(sub_file.name))
	sub_file.write('else\n')
	sub_file.write('\t touch %s.fail\n'%os.path.abspath(sub_file.name))
	sub_file.write('fi\n')
	sub_file.write('rm -f %s.run\n'%os.path.abspath(sub_file.name))
	sub_file.close()
	os.system('chmod +x %s'%os.path.abspath(sub_file.name))
	if not opts.dryRun and opts.queue:
		os.system('rm -f %s.done'%os.path.abspath(sub_file.name))
		os.system('rm -f %s.fail'%os.path.abspath(sub_file.name))
		os.system('rm -f %s.log'%os.path.abspath(sub_file.name))
		os.system('bsub -q %s -o %s.log %s'%(opts.queue,os.path.abspath(sub_file.name),os.path.abspath(sub_file.name)))
	if opts.runLocal:
		os.system('bash %s'%os.path.abspath(sub_file.name))

def writeAsymptotic():
	print 'Writing Asymptotic'
	try:
		assert(opts.masses_per_job)
	except AssertionError:
		sys.exit('No masses have been defined')

	for j, mass_set in enumerate(opts.masses_per_job):
		file = open('%s/sub_job%d.sh'%(opts.outDir,j),'w')
		writePreamble(file)
		exec_line = ''
		for mass in mass_set:
			exec_line +=	'combine %s -M Asymptotic -m %6.2f --cminDefaultMinimizerType=Minuit2'%(opts.datacard,mass)
			if opts.expected: exec_line += ' --run=expected'
			if mass!=mass_set[-1]: exec_line += ' && '
		writePostamble(file,exec_line)

def writeAsymptoticGrid():
	print 'Writing AsymptoticGrid'
	
	if not os.path.exists(os.path.expandvars('$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/makeAsymptoticGrid.py')):
		sys.exit('ERROR - CombinedLimit package must be installed')
	
	try:
		assert(opts.masses_per_job)
	except AssertionError:
		sys.exit('No masses have been defined')
	
	# create specialised limit grid workspace
	if not opts.skipWorkspace:
		print 'Creating workspace for %s...'%opts.method
		ws_exec_line = 'text2workspace.py %s -o %s'%(os.path.abspath(opts.datacard),os.path.abspath(opts.datacard).replace('.txt','.root')) 
		if opts.verbose: print '\t', ws_exec_line 
		os.system(ws_exec_line)
	opts.datacard = opts.datacard.replace('.txt','.root')

	# sub jobs through combine
	for j, mass_set in enumerate(opts.masses_per_job):
		for mass in mass_set:
			os.system('python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/makeAsymptoticGrid.py -w %s -m %6.2f -n 10 -r %3.1f %3.1f --runLimit --nCPU=3 -d %s'%(opts.datacard,mass,opts.muLow,opts.muHigh,os.path.abspath(opts.outDir)))
			sub_file_name = os.path.abspath(opts.outDir+'/limitgrid_%5.2f.sh'%(mass)) 
			if not opts.dryRun and opts.queue:
				os.system('rm -f %s.log'%os.path.abspath(sub_file_name))
				os.system('bsub -q %s -n 3 -R "span[hosts=1]" -o %s.log %s'%(opts.queue,os.path.abspath(sub_file_name),os.path.abspath(sub_file_name)))
			if opts.runLocal:
				os.system('bash %s'%os.path.abspath(sub_file_name))

	# switch back
	opts.datacard = opts.datacard.replace('.root','.txt')

def writeProfileLikelhood():

	print 'Writing ProfileLikelihood'
	try:
		assert(opts.masses_per_job)
	except AssertionError:
		sys.exit('No masses have been defined')

	tempcardstore = opts.datacard
	if opts.splitChannels: splitCard()
	for j, mass_set in enumerate(opts.masses_per_job):
		file = open('%s/sub_job%d.sh'%(opts.outDir,j),'w')
		writePreamble(file)
		exec_line = ''
		for mass in mass_set:
			exec_line +=	'combine %s -M ProfileLikelihood -m %6.2f --signif --pval --cminDefaultMinimizerType=Minuit2'%(opts.datacard,mass)
			if opts.expected: exec_line += ' -t -1'
			if opts.expectSignal: exec_line += ' --expectSignal=%3.1f'%opts.expectSignal
			if opts.expectSignalMass: exec_line += ' --expectSignalMass=%6.2f'%opts.expectSignalMass
			if mass!=mass_set[-1]: exec_line += ' && '
		writePostamble(file,exec_line)
		# change back
		opts.datacard = tempcardstore

def writeChannelCompatibility():

	print 'Writing ChannelCompatibility'
	try:
		assert(opts.mh)
	except AssertionError:
		sys.exit('mh is not defined')

	file = open('%s/sub_m%6.2f.sh'%(opts.outDir,opts.mh),'w')
	writePreamble(file)
	exec_line = 'combine %s -M ChannelCompatibilityCheck -m %6.2f --rMin=-25. --saveFitResult --cminDefaultMinimizerType=Minuit2'%(opts.datacard,opts.mh)
	writePostamble(file,exec_line)

def writeMultiDimFit():

	print 'Writing MultiDim Scan'
	ws_args = { "RVRFScan" 	 : "-P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs" ,
							"RVScan"	 	 : "-P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs" ,
							"RVnpRFScan" : "-P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs" ,
							"RFScan"	 	 : "-P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs" ,
							"RFnpRVScan" : "-P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs" ,
							"MuScan"		 : "",
							"CVCFScan"	 : "-P HiggsAnalysis.CombinedLimit.HiggsCouplingsLOSM:cVcF",
							"MHScan"		 : "-P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs --PO higgsMassRange=120,130",
							"MuMHScan"	 : "-P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass"
						}
	
	combine_args = { "RVRFScan" 	 : "-P RV -P RF" ,
                   "RVScan"	 	 	 : "--floatOtherPOIs=1 -P RV" ,
                   "RVnpRFScan"	 : "--floatOtherPOIs=0 -P RV" ,
                   "RFScan"	 	 	 : "--floatOtherPOIs=1 -P RF" ,
                   "RFnpRVScan"	 : "--floatOtherPOIs=0 -P RF" ,
                   "MuScan"		 	 : "-P r",
                   "CVCFScan"	 	 : "-P CV -P CF",
                   "MHScan"		 	 : "--floatOtherPOIs=1 -P MH",
                   "MuMHScan"	 	 : "-P r -P MH"
								 }
	par_ranges = {}
	if opts.rvLow!=None and opts.rvHigh!=None and opts.rfLow!=None and opts.rfHigh!=None:
		par_ranges["RVRFScan"]		= "RV=%4.2f,%4.2f:RF=%4.2f,%4.2f"%(opts.rvLow,opts.rvHigh,opts.rfLow,opts.rfHigh)
	if opts.rvLow!=None and opts.rvHigh!=None:
		par_ranges["RVScan"]			= "RV=%4.2f,%4.2f"%(opts.rvLow,opts.rvHigh) 
	if opts.rvLow!=None and opts.rvHigh!=None:
		par_ranges["RVnpRFScan"]	= "RV=%4.2f,%4.2f"%(opts.rvLow,opts.rvHigh)
	if opts.rfLow!=None and opts.rfHigh!=None:
		par_ranges["RFScan"]			= "RF=%4.2f,%4.2f"%(opts.rfLow,opts.rfHigh)
	if opts.rfLow!=None and opts.rfHigh!=None:
		par_ranges["RFnpRVScan"]	= "RF=%4.2f,%4.2f"%(opts.rfLow,opts.rfHigh)
	if opts.muLow!=None and opts.muHigh!=None:
		par_ranges["MuScan"]			= "r=%4.2f,%4.2f"%(opts.muLow,opts.muHigh) 
	if opts.cvLow!=None and opts.cvHigh!=None and opts.cfLow!=None and opts.cfHigh!=None:
		par_ranges["CVCFScan"]		= "CV=%4.2f,%4.2f:CF=%4.2f,%4.2f"%(opts.cvLow,opts.cvHigh,opts.cfLow,opts.cfHigh)
	if opts.mhLow!=None and opts.mhHigh!=None:
		par_ranges["MHScan"]			= "MH=%6.2f,%6.2f"%(opts.mhLow,opts.mhHigh)
	if opts.muLow!=None and opts.muHigh!=None and opts.mhLow!=None and opts.mhHigh!=None:
		par_ranges["MuMHScan"]		= "r=%4.2f,%4.2f:MH=%6.2f,%6.2f"%(opts.muLow,opts.muHigh,opts.mhLow,opts.mhHigh)

	# create specialised MultiDimFit workspace
	if not opts.skipWorkspace:
		print 'Creating workspace for %s...'%opts.method
		ws_exec_line = 'text2workspace.py %s -o %s %s'%(os.path.abspath(opts.datacard),os.path.abspath(opts.datacard).replace('.txt',opts.method+'.root'),ws_args[opts.method]) 
		if opts.verbose: print '\t', ws_exec_line 
		os.system(ws_exec_line)
	opts.datacard = opts.datacard.replace('.txt',opts.method+'.root')

	# make job scripts
	for i in range(opts.jobs):
		file = open('%s/sub_m%6.2f_job%d.sh'%(opts.outDir,opts.mh,i),'w')
		writePreamble(file)
		exec_line = 'combine %s -M MultiDimFit --X-rtd ADDNLL_FASTEXIT --cminDefaultMinimizerType Minuit2 --algo=grid --saveNLL --setPhysicsModelParameterRanges %s %s --points=%d --firstPoint=%d --lastPoint=%d -n Job%d'%(opts.datacard,par_ranges[opts.method],combine_args[opts.method],opts.pointsperjob*opts.jobs,i*opts.pointsperjob,(i+1)*opts.pointsperjob-1,i)
		if opts.mh: exec_line += ' -m %6.2f'%opts.mh
		if opts.expected: exec_line += ' -t -1'
		if opts.expectSignal: exec_line += ' --expectSignal %4.2f'%opts.expectSignal
		if opts.expectSignalMass: exec_line += ' --expectSignalMass %6.2f'%opts.expectSignalMass
		if opts.additionalOptions: exec_line += ' %s'%opts.additionalOptions 
		if opts.verbose: print '\t', exec_line
		writePostamble(file,exec_line)

	opts.datacard = opts.datacard.replace(opts.method+'.root','.txt')

def run():
	os.system('mkdir -p %s'%opts.outDir)
	if opts.verbose: print 'Made directory', opts.outDir
	checkValidMethod()
	if opts.method=='Asymptotic' or opts.method=='AsymptoticGrid' or opts.method=='ProfileLikelihood':
		configureMassFromNJobs()
	if opts.method=='Asymptotic':
		writeAsymptotic()
	elif opts.method=='AsymptoticGrid':
		writeAsymptoticGrid()
	elif opts.method=='ProfileLikelihood':
		writeProfileLikelhood()
	elif opts.method=='ChannelCompatibilityCheck':
		writeChannelCompatibility()
	else:
		writeMultiDimFit()

def resetDefaultConfig():
	for opt in specOpts.option_list:
		opt_name = opt.dest.strip('--')
		if opt_name=='datacard' or opt_name=='files': continue
		else: setattr(opts,opt_name,None)

def configure(config_line):
	# could automate this but makes it easier to read and add options this way
	resetDefaultConfig()
	if opts.verbose: print config_line
	for option in config_line.split():
		if option.startswith('outDir='): opts.outDir = option.split('=')[1]
		if option.startswith('method='): opts.method = option.split('=')[1]
		if option.startswith('expected='): opts.expected = int(option.split('=')[1])
		if option.startswith('expectSignal='): opts.expectSignal = float(option.split('=')[1])
		if option.startswith('expectSignalMass='): opts.expectSignalMass = float(option.split('=')[1])
		if option.startswith('mhLow='): opts.mhLow = float(option.split('=')[1])
		if option.startswith('mhHigh='): opts.mhHigh = float(option.split('=')[1])
		if option.startswith('mhStep='): opts.mhStep = float(option.split('=')[1])
		if option.startswith('jobs='): opts.jobs = int(option.split('=')[1])
		if option.startswith('pointsperjob='): opts.pointsperjob = int(option.split('=')[1])
		if option.startswith('splitChannels='): opts.splitChannels = option.split('=')[1].split(',')
		if option.startswith('mh='): opts.mh = float(option.split('=')[1])
		if option.startswith('muLow='): opts.muLow = float(option.split('=')[1])
		if option.startswith('muHigh='): opts.muHigh = float(option.split('=')[1])
		if option.startswith('rvLow='): opts.rvLow = float(option.split('=')[1])
		if option.startswith('rvHigh='): opts.rvHigh = float(option.split('=')[1])
		if option.startswith('rfLow='): opts.rfLow = float(option.split('=')[1])
		if option.startswith('rfHigh='): opts.rfHigh = float(option.split('=')[1])
		if option.startswith('cvLow='): opts.cvLow = float(option.split('=')[1])
		if option.startswith('cvHigh='): opts.cvHigh = float(option.split('=')[1])
		if option.startswith('cfLow='): opts.cfLow = float(option.split('=')[1])
		if option.startswith('cfHigh='): opts.cfHigh = float(option.split('=')[1])
		if option.startswith('opts='): opts.additionalOptions = option.split('=')[1].replace(',',' ')
	if opts.verbose: print opts
	run()

def trawlHadd():
	list_of_dirs=set()
	for root, dirs, files in os.walk(opts.hadd):
		for x in files:
			if 'higgsCombine' in x and '.root' in x: 
				list_of_dirs.add(root)

	for dir in list_of_dirs:
		for root, dirs, files in os.walk(dir):
			list_of_files=''
			for file in fnmatch.filter(files,'higgsCombine*.root'):
				list_of_files += ' '+os.path.join(root,'%s'%file)
			print root, ' -- ', len(list_of_files.split())
			exec_line = 'hadd -f %s/%s.root%s'%(dir,os.path.basename(dir),list_of_files)
			if opts.verbose: print exec_line
			os.system(exec_line)

if opts.hadd:
	trawlHadd()
elif opts.datfile:
	datfile = open(opts.datfile)
	for line in datfile.readlines():
		line=line.strip('\n')
		if line.startswith('#') or len(line)==0: 
			continue
		if line.startswith('datacard'): 
			opts.datacard = line.split('=')[1]
			assert('.txt' in opts.datacard)
			continue
		if line.startswith('files'):
			opts.files = line.split('=')[1]
			continue
		configure(line)
else:
	# default setup here
	print 'Not yet implemented'
