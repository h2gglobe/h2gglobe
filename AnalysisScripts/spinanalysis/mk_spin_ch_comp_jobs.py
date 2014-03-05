#!/usr/bin/env python

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-i","--card")
parser.add_option("-d","--dir")
parser.add_option("-k","--kinCats",type="int",default=4)
parser.add_option("-s","--spinCats",type="int",default=5)
parser.add_option("-m","--genMass",type="float",default=125)
parser.add_option("--gridpoints",type="int")
parser.add_option("--skipWorkspace",dest="skipWorkspace",action="store_true",default=False)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
parser.add_option("--noFloatMass",action="store_true",default=False)
(options,args)=parser.parse_args()

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
	card = open(options.card)
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
