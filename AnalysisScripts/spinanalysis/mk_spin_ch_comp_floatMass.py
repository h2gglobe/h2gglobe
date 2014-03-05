#!/usr/bin/env python

from optparse import OptionParser
parser=OptionParser()
parser.add_option("--cardSpin")
parser.add_option("--cardJustSM")
parser.add_option("-d","--dir")
parser.add_option("-k","--kinCats",type="int",default=4)
parser.add_option("-s","--spinCats",type="int",default=5)
parser.add_option("-g","--gridpoints",type="int")
parser.add_option("-G","--seperateGridJobs",action="store_true",default=False)
parser.add_option("--skipWorkspace",dest="skipWorkspace",action="store_true",default=False)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
parser.add_option("--dataOnly",action="store_true",default=False)
parser.add_option("--toysOnly",action="store_true",default=False)
parser.add_option("-q","--queue")
(options,args)=parser.parse_args()

import os
import sys
import time
os.system('mkdir -p %s'%options.dir)
import ROOT as r

if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/combine')):
	sys.exit('ERROR - CombinedLimit package must be installed')
if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/text2workspace.py')):
	sys.exit('ERROR - CombinedLimit package must be installed')
if not os.path.exists(os.path.expandvars('$CMSSW_BASE/bin/$SCRAM_ARCH/combineCards.py')):
	sys.exit('ERROR - CombinedLimit package must be installed')

def buildCatsMap():
	card = open(options.cardSpin)
	cats = set()
	lines = card.readlines()
	for line in lines:
		if line.startswith('bin'):
			for el in line.split()[1:]:
				cats.add(el)
			break
	line=''
	for sCat in range(options.spinCats):
		line += '--PO map='
		inc_cats=[]
		for kCat in range(options.kinCats):
			cat = kCat*options.spinCats+sCat
			for gcat in cats:
				if 'cat%d_'%cat in gcat:
					inc_cats.append(gcat)
		for cat in inc_cats:
			line += '.*'+cat+'/.*H.*,'
		line = line[:-1]
		line+=':r_spinCat%d[1.,-5.,5] '%sCat
	line = line[:-1]
	return line

def checkMulti(cardname,inc_cats):
	os.system('mv %s temp.txt'%cardname)
	newf = open(cardname,'w')
	oldf = open('temp.txt')
	for line in oldf.readlines():
		if 'discrete' in line:
			for cat in inc_cats:
				if 'pdfindex_%s'%(cat.split('cat')[1]) in line: newf.write(line)
				else: continue
		else:
			newf.write(line)
	newf.close()
	oldf.close()
	os.system('rm -f temp.txt')

def mk_spin_cat_cards():
	new_cards=[]
	card = open(options.cardSpin)
	cats = set()
	lines = card.readlines()
	for line in lines:
		if line.startswith('bin'):
			for el in line.split()[1:]:
				cats.add(el)
			break
	for sCat in range(options.spinCats):
		exc_string=''
		inc_cats=[]
		for kCat in range(options.kinCats):
			cat = kCat*options.spinCats+sCat
			for gcat in cats:
				if 'cat%d_'%cat in gcat:
					inc_cats.append(gcat)

		for cat in cats:
			if cat not in inc_cats:
				exc_string += 'ch1_%s|'%cat
		exc_string = exc_string[:-1]
		newcardname = card.name.replace('.txt','_spinCat%d.txt'%sCat)
		os.system('combineCards.py %s --xc=\"%s\" > %s'%(card.name,exc_string,newcardname))
		checkMulti(newcardname,inc_cats)
		new_cards.append(newcardname)
	return new_cards

def writeScript(subfile,rootfiles,exec_line,outloc=''):
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
		subfile.write('\techo \"Waiting for files\"\n')
		subfile.write('done\n')
		subfile.write('cp %s/%s/%s .\n'%(os.getcwd(),options.dir,rootfile))
	subfile.write('if ( %s )\n'%exec_line)
	subfile.write('\tthen touch %s.done\n'%subfile.name)
	subfile.write('\telse touch %s.fail\n'%subfile.name)
	subfile.write('fi\n')
	subfile.write('rm -f %s.run\n'%subfile.name)
	subfile.write('mkdir -p %s/%s/%s\n'%(os.getcwd(),options.dir,outloc))
	subfile.write('mv higgsCombine*.root %s/%s/%s/\n'%(os.getcwd(),options.dir,outloc))

# MAIN
# build cat map
catsMap = buildCatsMap()

# make per cat card
#per_cat_cards = mk_spin_cat_cards()
#print per_cat_cards

# make workspaces for SM fit
if not options.skipWorkspace:
	os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel %s --PO higgsMassRange=120,130 -o %s/%s'%(options.cardJustSM,catsMap,options.dir,options.cardJustSM.replace('.txt','.root'))) 
options.cardJustSM = options.cardJustSM.replace('.txt','.root')

if not options.dataOnly:
	# fit this to data (which will profile 5 signal strengths independently and the mass) to find best fit mass
	f = open('%s/%s/sub_sm_fit_to_data_multisignal.sh'%(os.getcwd(),options.dir),'w')
	line = 'combine %s -M MultiDimFit --redefineSignalPOIs '%options.cardJustSM
	for sCat in range(options.spinCats):
		line += 'r_spinCat%d,'%sCat
	line +='MH'
	line += ' --setPhysicsModelParameterRanges '
	for sCat in range(options.spinCats):
		line += 'r_spinCat%d=-5,5.:'%sCat
	line = line[:-1]
	line += ' -n SMFitToDataMultiSignal'
	writeScript(f,[options.cardJustSM],line)
	f.close()
	os.system('chmod +x %s'%f.name)
	ans = raw_input('INITIAL FIT TO DATA. Run locally? Will submit to batch otherwise (y/n). Type enter if the fit is already done and you want to go to the next step.\n')
	if ans=='y' or ans=='Y' or ans=='yes':
		print 'Fitting data..'
		os.system('bash %s >& %s.log'%(f.name,f.name))
		print 'Done'
	elif ans=='n' or ans=='N' or ans=='no':
		os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))

	# now need to get this back out
	while not os.path.exists('%s/%s/higgsCombineSMFitToDataMultiSignal.MultiDimFit.mH120.root'%(os.getcwd(),options.dir)):
		print 'Waiting for data fit...'
		print 'Fitting output should be available in %s.log'%f.name
		time.sleep(30)
	print 'Data fit file found'

	tf = r.TFile('%s/%s/higgsCombineSMFitToDataMultiSignal.MultiDimFit.mH120.root'%(os.getcwd(),options.dir))
	tree = tf.Get('limit')
	tree.GetEntry(0)
	mass = tree.MH
	tf.Close()
	print 'Best fit mass was found to be %6.2f'%mass
	print 'Will now generate spin model toys from this mass point and fit them back'

	# now generate a toy at this mass for each spin hypothesis
	# make spin workspace first
	if not options.skipWorkspace:
		os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating --PO=muFloating --PO higgsMassRange=122,128 -o %s/%s'%(options.cardSpin,options.dir,options.cardSpin.replace('.txt','.root')))
	options.cardSpin = options.cardSpin.replace('.txt','.root')

	# gen toys
	for name,fqq,x in [ ['sm',0.0,0.0] , ['gravgg',0.0,1.0] , ['grav50gg50qq',0.5,1.0] , ['gravqq',1.0,1.0] ]:
		f = open('%s/%s/sub_gen_%s.sh'%(os.getcwd(),options.dir,name),'w')
		line = 'combine %s -M GenerateOnly -m %6.2f -t -1 --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x,MH --redefineSignalPOIs r --expectSignal=1 -n _%s_asimov --saveToys -s 0'%(options.cardSpin,mass,fqq,x,name)
		writeScript(f,[options.cardSpin],line)
		f.close()
		os.system('chmod +x %s'%f.name)
		print 'Throwing Asimov for %s hypothesis - fqq=%4.2f x=%1.0f'%(name,fqq,x)
		os.system('bash %s >& %s.log'%(f.name,f.name))
		print 'Done'

	# now fit the toys back - use sm workspace for this - don't need errors so don't run grid
	for name,fqq,x in [ ['sm',0.0,0.0] , ['gravgg',0.0,1.0] , ['grav50gg50qq',0.5,1.0] , ['gravqq',1.0,1.0] ]:
		f = open('%s/%s/sub_fit_sm_to_%s.sh'%(os.getcwd(),options.dir,name),'w')
		line = 'combine %s -M MultiDimFit --redefineSignalPOIs '%options.cardJustSM
		for sCat in range(options.spinCats):
			line += 'r_spinCat%d,'%sCat
		line +='MH'
		line += ' --setPhysicsModelParameterRanges '
		for sCat in range(options.spinCats):
			line += 'r_spinCat%d=-5,5.:'%sCat
		line = line[:-1]
		line += ' -n _sm_fit_to_%s_asimov_MultiSignal'%name
		mstring = (('%6.2f'%mass).replace('.00','')).replace('.0','')
		toyfile = ('higgsCombine_%s_asimov.GenerateOnly.mH%s.0.root'%(name,mstring))
		line += (' -t -1 --toysFile %s'%toyfile) 
		writeScript(f,[options.cardJustSM,toyfile],line)
		f.close()
		os.system('chmod +x %s'%f.name)
		print 'Fitting sm to Asimov for %s hypothesis - fqq=%4.2f x=%1.0f'%(name,fqq,x)
		if not options.dryRun:
			os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))

if not options.toysOnly:
	# should also do data fit again but doing the proper profiles to get the error
	for cat in range(options.spinCats):
		if options.seperateGridJobs:
			for job in range(options.gridpoints):
				f = open('%s/%s/sub_sm_fit_to_data_spinCat%d_job%d.sh'%(os.getcwd(),options.dir,cat,job),'w')
				line = 'combine %s -M MultiDimFit -P r_spinCat%d --floatOtherPOIs=1 --algo=grid --points=%d --squareDistPoiStep --firstPoint=%d --lastPoint=%d'%(options.cardJustSM,cat,options.gridpoints,job,job+1)
				line += ' --setPhysicsModelParameterRanges '
				for sCat in range(options.spinCats):
					line += 'r_spinCat%d=-1,4.:'%sCat
				line = line[:-1]
				line += ' -n SMFitToDataSpinCat%dJob%d'%(cat,job)
				writeScript(f,[options.cardJustSM],line)
				f.close()
				os.system('chmod +x %s'%f.name)
				print 'Fitting sm to data - profiling r_spinCat%d'%cat
				if not options.dryRun:
					os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))
		else:
			f = open('%s/%s/sub_sm_fit_to_data_spinCat%d.sh'%(os.getcwd(),options.dir,cat),'w')
			line = 'combine %s -M MultiDimFit -P r_spinCat%d --floatOtherPOIs=1 --algo=grid --points=%d --squareDistPoiStep'%(options.cardJustSM,cat,options.gridpoints)
			line += ' --setPhysicsModelParameterRanges '
			for sCat in range(options.spinCats):
				line += 'r_spinCat%d=-5,5.:'%sCat
			line = line[:-1]
			line += ' -n SMFitToDataSpinCat%d'%cat
			writeScript(f,[options.cardJustSM],line)
			f.close()
			os.system('chmod +x %s'%f.name)
			print 'Fitting sm to data - profiling r_spinCat%d'%cat
			if not options.dryRun:
				os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))

"""

assert len(per_cat_cards)==options.spinCats
# now generate toy in each spin cat for different spin models
for cat in range(options.spinCats):
	card = per_cat_cards[cat]
	if not options.skipWorkspace:
		os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating --PO=muFloating --PO higgsMassRange=122,128 -o %s/%s'%(card,options.dir,card.replace('.txt','.root')))
	card = card.replace('.txt','.root')

	for label,name,fqq,x in [ ['sm','SMAsimov',0.0,0.0] , ['gravgg','GravGGAsimov',0.0,1.0] , ['grav50gg50qq','Grav50GG50QQAsimov',0.5,1.0] , ['gravqq','GravQQAsimov',1.0,1.0] ]:
		f = open('%s/%s/sub_gen_%s_cat%d.sh'%(os.getcwd(),options.dir,label,cat),'w')
		line = 'combine %s -M GenerateOnly -m %6.2f -t -1 --setPhysicsModelParameters fqq=%4.2f,x=%1.0f --freezeNuisances fqq,x,MH --redefineSignalPOIs r --expectSignal=1 -n %sCat%d --saveToys -s 0'%(card,mass,fqq,x,name,cat)
		writeScript(f,[card],line)
		f.close()
		os.system('chmod +x %s'%f.name)
		print 'Throwing Asimov'
		os.system('bash %s >& %s.log'%(f.name,f.name))
		print 'Done'

		toyfile = ('higgsCombine%sCat%d.GenerateOnly.mH%6.2f.0.root'%(name,cat,mass)).replace('.00','')
		f = open('%s/%s/sub_fit_smto%s_cat%d.sh'%(os.getcwd(),options.dir,label,cat),'w')
		line = 'combine %s -M MultiDimFit -m %6.2f --setPhysicsModelParameters fqq=0.0,x=0. --freezeNuisances fqq,x,MH --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-5.,5. -n SMFit%sCat%d -t -1 --toysFile %s'%(card,mass,name,cat,toyfile)
		writeScript(f,[card,toyfile],line)
		f.close()
		os.system('chmod +x %s'%f.name)
		if not allLocal and not allBatch: 
			ans = raw_input('FIT TO TOY. Run locally? Will submit to batch otherwise (y/n). Type Y (N) to run all local or all batch.\n')
		if ans=='Y':
			allLocal=True
			os.system('bash %s >& %s.log'%(f.name,f.name))
		elif ans=='N':
			allBatch=True
			os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))
		if ans=='y' or ans=='Y' or ans=='yes' or allLocal:
			os.system('bash %s >& %s.log'%(f.name,f.name))
		elif ans=='n' or ans=='N' or ans=='no' or allBatch:
			os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))
		else:
			print 'Nothing to do'
		
def writeGen(spinCat,fqq,x,label,name):
	card = options.card.replace('.txt','_spinCat%d.root'%spinCat)
	f = open('%s/%s/sub_gen_%s_cat%d.sh'%(os.getcwd(),options.dir,label,spinCat),'w')
	f.write('#!/bin/bash\n')
	f.write('rm -f %s.done\n'%f.name)
	f.write('rm -f %s.fail\n'%f.name)
	f.write('touch %s.run\n'%f.name)
	f.write('mkdir -p scratch\n')
	f.write('cd scratch\n')
	f.write('cd %s/%s\n'%(os.getcwd(),options.dir))
	f.write('eval `scramv1 runtime -sh`\n')
	f.write('cd -\n')
	f.write('cp %s/%s/%s .\n'%(os.getcwd(),options.dir,card))
	f.write('if ( combine %s -M GenerateOnly -m %6.2f -t -1 --setPhysicsModelParameters fqq=%3.1f,x=%3.1f --freezeNuisances fqq,x --redefineSignalPOIs r --expectSignal=1 -n %s --saveToys -s 0 )\n'%(card,options.genMass,fqq,x,name))
	f.write('\tthen touch %s.done\n'%f.name)
	f.write('\telse touch %s.fail\n'%f.name)
	f.write('fi\n')
	f.write('rm -f %s.run\n'%f.name)
	f.write('cp higgsCombine*.root %s/%s/\n'%(os.getcwd(),options.dir))
	f.close()
	os.system('chmod +x %s'%f.name)
	if not options.dryRun:
		os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))

def writeFit(spinCat,fqq,x,label,name,toylabel,toyname,isData=False):
	card = options.card.replace('.txt','_spinCat%d.root'%spinCat)
	toyfile = 'higgsCombine%s.GenerateOnly.mH%6.2f.0.root'%(toyname,options.genMass)
	f = open('%s/%s/sub_fit_%s_cat%d.sh'%(os.getcwd(),options.dir,label,spinCat),'w')
	f.write('#!/bin/bash\n')
	f.write('rm -f %s.done\n'%f.name)
	f.write('rm -f %s.fail\n'%f.name)
	f.write('touch %s.run\n'%f.name)
	f.write('mkdir -p scratch\n')
	f.write('cd scratch\n')
	f.write('cd %s/%s\n'%(os.getcwd(),options.dir))
	f.write('eval `scramv1 runtime -sh`\n')
	f.write('cd -\n')
	f.write('cp %s/%s/%s .\n'%(os.getcwd(),options.dir,card))
	if not isData:
		f.write('while [ ! -f %s/%s/sub_gen_%s_cat%d.sh.done ]\n'%(os.getcwd(),options.dir,toylabel,spinCat))
		f.write('\tdo sleep 10\n')
		f.write('\techo \"Waiting for toy\"\n')
		f.write('done\n')
		f.write('cp %s/%s/%s .\n'%(os.getcwd(),options.dir,toyfile))
	f.write('if ( combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%3.1f,x=%3.1f --freezeNuisances fqq,x,MH --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-5.,5. -m %6.2f -n %s'%(card,fqq,x,options.genMass,name))
	if not isData: f.write(' -t -1 --toysFile %s'%toyfile)
	if options.gridpoints: f.write(' --algo=grid --points %d --squareDistPoiStep'%(options.gridpoints))
	f.write(' )\n')
	f.write('\tthen touch %s.done\n'%f.name)
	f.write('\telse touch %s.fail\n'%f.name)
	f.write('fi\n')
	f.write('rm -f %s.run\n'%f.name)
	f.write('cp higgsCombine*.root %s/%s/\n'%(os.getcwd(),options.dir))
	f.close()
	os.system('chmod +x %s'%f.name)
	if not options.dryRun:
		os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))

def writeFitComb(fqq,x,label,name):
	card = options.card.replace('.txt','.root')
	f = open('%s/%s/sub_fit_%s.sh'%(os.getcwd(),options.dir,label),'w')
	f.write('#!/bin/bash\n')
	f.write('rm -f %s.done\n'%f.name)
	f.write('rm -f %s.fail\n'%f.name)
	f.write('touch %s.run\n'%f.name)
	f.write('mkdir -p scratch\n')
	f.write('cd scratch\n')
	f.write('cd %s/%s\n'%(os.getcwd(),options.dir))
	f.write('eval `scramv1 runtime -sh`\n')
	f.write('cd -\n')
	f.write('cp %s/%s/%s .\n'%(os.getcwd(),options.dir,card))
	f.write('if ( combine %s -M MultiDimFit --setPhysicsModelParameters fqq=%3.1f,x=%3.1f --freezeNuisances fqq,x,MH --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-5.,5. -m %6.2f -n %s'%(card,fqq,x,options.genMass,name))
	if options.gridpoints: f.write(' --algo=grid --points %d --squareDistPoiStep'%(options.gridpoints))
	f.write(' )\n')
	f.write('\tthen touch %s.done\n'%f.name)
	f.write('\telse touch %s.fail\n'%f.name)
	f.write('fi\n')
	f.write('rm -f %s.run\n'%f.name)
	f.write('cp higgsCombine*.root %s/%s/\n'%(os.getcwd(),options.dir))
	f.close()
	os.system('chmod +x %s'%f.name)
	if not options.dryRun:
		os.system('bsub -q 2nd -o %s.log %s'%(f.name,f.name))

# __main__
per_cat_cards = mk_spin_cat_cards()

if not options.skipWorkspace:
	print 'Using combine tool to make workspace from card', options.card
	if options.noFloatMass:
		os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating --PO=muFloating -m %6.2f -o %s/%s'%(options.card,options.genMass,options.dir,combfile))
		for card in per_cat_cards:
			os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating --PO=muFloating -m %6.2f -o %s/%s'%(card,options.genMass,options.dir,card.replace('.txt','.root')))
	else:
		os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating --PO=muFloating --PO higgsMassRange=122,128 -o %s/%s'%(options.card,options.dir,combfile))
		for card in per_cat_cards:
			os.system('text2workspace.py %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating --PO=muFloating --PO higgsMassRange=122,128 -o %s/%s'%(card,options.dir,card.replace('.txt','.root')))
		
# write scripts
# toys
for cat in range(options.spinCats):
	writeGen(cat,0.0,0,'sm','SMAsimovCat%d'%cat)
	writeGen(cat,0.0,1,'gravgg','GravGGAsimovCat%d'%cat)
	writeGen(cat,0.5,1,'grav50gg50qq','Grav50GG50QQAsimovCat%d'%cat)
	writeGen(cat,1.0,1,'gravqq','GravQQAsimovCat%d'%cat)

	writeFit(cat,0.0,0.,'smtosm','SMtoSMCat%d'%cat,'sm','SMAsimovCat%d'%cat)
	writeFit(cat,0.0,0.,'smtogravgg','SMtoGravGGCat%d'%cat,'gravgg','GravGGAsimovCat%d'%cat)
	writeFit(cat,0.0,0.,'smtograv50gg50qq','SMtoGrav50GG50QQCat%d'%cat,'grav50gg50qq','Grav50GG50QQAsimovCat%d'%cat)
	writeFit(cat,0.0,0.,'smtogravqq','SMtoGravQQCat%d'%cat,'gravqq','GravQQAsimovCat%d'%cat)

	writeFit(cat,0.0,0.,'smtodata','SMtoDataCat%d'%cat,'','',True)
	writeFit(cat,0.0,1.,'gravggtodata','GravGGtoDataCat%d'%cat,'','',True)
	writeFit(cat,0.5,1.,'grav50gg50qqtodata','Grav50GG50QQtoDataCat%d'%cat,'','',True)
	writeFit(cat,1.0,1.,'gravqqtodata','GravQQtoDataCat%d'%cat,'','',True)

writeFitComb(0.0,0.,'smtodata','SMtoData')
writeFitComb(0.0,1.,'gravggtodata','GravGGtoData')
writeFitComb(0.5,1.,'grav50gg50qqtodata','Grav50GG50QQtoData')
writeFitComb(1.0,1.,'gravqqtodata','GravQQtoData')
"""
