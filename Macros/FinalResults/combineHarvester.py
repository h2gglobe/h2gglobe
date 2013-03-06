#!/usr/bin/env python

# e.g
# ./subCombine.py --card2011=datacard_parametric_7TeV_massfac.txt --unblind --run2011 --splitChannels2011="cat0_7TeV cat1_7TeV cat2_7TeV cat3_7TeV" --splitChannels2011="cat4_7TeV"
# ./subCombine.py --card2012=datacard_parametric_8TeV_cutbased_cbvbf.txt --unblind --run2012 --splitChannels2012="cat0 cat1 cat2 cat3" --splitChannels2012="cat4 cat5" --splitChannels2012="cat6 cat7" --splitChannels2012="cat8"
# ./subCombine.py --cardBoth=hgg.combined.txt --unblind --runBoth --splitChannelsBoth="ch1_cat0_7TeV ch1_cat1_7TeV ch1_cat2_7TeV ch1_cat3_7TeV ch2_cat0_8TeV ch2_cat1_8TeV ch2_cat2_8TeV ch2_cat3_8TeV" --splitChannelsBoth="ch1_cat4_7TeV ch2_cat4_8TeV ch2_cat5_8TeV" --splitChannelsBoth="ch2_cat6_8TeV ch2_cat7_8TeV ch2_cat8_8TeV"

# one command to rule them all
# ./subCombine.py --card2011=datacard_parametric_7TeV_massfac.txt --unblind --run2011 --splitChannels2011="cat0 cat1 cat2 cat3" --splitChannels2011="cat4" --card2012=datacard_parametric_8TeV_cutbased_cbvbf.txt --unblind --run2012 --splitChannels2012="cat0 cat1 cat2 cat3" --splitChannels2012="cat4 cat5" --splitChannels2012="cat6 cat7" --splitChannels2012="cat8" --cardBoth=hgg.combined.txt --unblind --runBoth --splitChannelsBoth="ch1_cat0_7TeV ch1_cat1_7TeV ch1_cat2_7TeV ch1_cat3_7TeV ch2_cat0_8TeV ch2_cat1_8TeV ch2_cat2_8TeV ch2_cat3_8TeV" --splitChannelsBoth="ch1_cat4_7TeV ch2_cat4_8TeV ch2_cat5_8TeV" --splitChannelsBoth="ch2_cat6_8TeV ch2_cat7_8TeV ch2_cat8_8TeV" --doSpecials -q 1nh

import os
import numpy
import sys

from optparse import OptionParser

parser = OptionParser()
parser.add_option("","--card2011",dest="card2011",type="string",help="Name of 2011 card")
parser.add_option("","--card2012",dest="card2012",type="string",help="Name of 2012 card")
parser.add_option("","--cardBoth",dest="cardBoth",type="string",help="Name of combined card (default is to combine the two cards given)")
parser.add_option("","--unblind",dest="unblind",default=False,action="store_true",help="Unblind")
parser.add_option("","--parametric",dest="parametric",default=False,action="store_true",help="Run parametric model")
parser.add_option("-q","--queue",dest="queue",type="string",default="1nh",help="queue")
parser.add_option("-M","--methods",dest="methods",default=[],action="append",help="Combine methods - default does several (can pass multiple times)")
parser.add_option("-L","--mHlow",dest="mHlow",type="float",default=110,help="Lowest MH")
parser.add_option("-H","--mHhigh",dest="mHhigh",type="float",default=150,help="Highest MH")
parser.add_option("-S","--mHstep",dest="mHstep",type="float",default=0.5,help="MH step")
parser.add_option("","--run2011",dest="run2011",default=False,action="store_true",help="Run only 2011")
parser.add_option("","--run2012",dest="run2012",default=False,action="store_true",help="Run only 2012")
parser.add_option("","--runBoth",dest="runBoth",default=False,action="store_true",help="Run only 2011+2012")
parser.add_option("","--skipDatacard",dest="skipDatacard",default=False,action="store_true",help="Skip the original datacard (if running sub channels)")
parser.add_option("","--splitChannels2011",dest="splitChannels2011",default=[],action="append",help="Run these channels for 2011 (can pass multiple times)")
parser.add_option("","--splitChannels2012",dest="splitChannels2012",default=[],action="append",help="Run these channels for 2012 (can pass multiple times)")
parser.add_option("","--splitChannelsBoth",dest="splitChannelsBoth",default=[],action="append",help="Run these channels for Both (can pass multiple times)")
parser.add_option("","--doSpecials",dest="doSpecials",type="float",help="Do special plots e.g. BestFitMass, RVRF, MuMH. Can pass a mass here")
parser.add_option("","--dryRun",dest="dryRun",default=False,action="store_true",help="Do not submit")
parser.add_option("-v","--verbose",dest="verbose",default=False,action="store_true")
(options,args)=parser.parse_args()

path = os.getcwd()

# define a method for splitting the card into sub channels if requested
def cardSplitting(card,subCatJobs):
  categories = os.popen('grep bin %s | grep -v max | grep -v combine | grep -v CMS-HGG | grep -v shapes | grep -v observation | grep -v Combination | head -1'%card).readlines()[0].strip('bin')
  allCats = categories.split()
  splitCards=[]
  for job in subCatJobs:
    mycats = job.split(" ")
    veto = ""
    for cat in allCats:
      if cat in mycats: continue
      else: veto += "|ch1_"+cat
    veto=veto[1:]
    newcardname = '%s_%s.txt'%(card.replace('.txt',''),job.replace(' ',''))
    if options.verbose: print 'combineCards.py --xc="%s" %s > %s'%(veto,card,newcardname)
    os.system('combineCards.py --xc="%s" %s > %s'%(veto,card,newcardname))
    splitCards.append([job.replace(' ',''),newcardname])
  return splitCards

def writePreamble(file,path,card,mh,type=''):
  file.write('#!/bin/bash\n')
  file.write('cd %s\n'%path)
  file.write('touch %s.run\n'%file.name)
  file.write('echo ---------------------------------------------------\n')
  file.write('echo   Running %s at mass %5.1f '%(card,mh))
  if type!='': file.write('with options -- %s \n'%(type))
  else: file.write('\n')
  file.write('echo ---------------------------------------------------\n')
  file.write('eval `scramv1 runtime -sh`\n')

def writeAsymptotic(file,path,card,folder,mh):
  if 'Asymptotic' in options.methods:
    m = '%5.1f'%mh
    m = m.replace('.0','')
    file.write('if ( \n')
    if options.unblind: file.write('\t combine %s/%s -M Asymptotic -m %s -n %s\n'%(path,card,m,folder))  
    else: file.write('\t combine %s/%s -M Asymptotic -m %s -n %s --run=expected\n'%(path,card,m,folder))  
    file.write(') then\n')
    file.write('\t mv higgsCombine%s.Asymptotic.mH%s.root %s/%s/Asymptotic\n'%(folder,m,path,folder))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')

def writeProfileLikelhood(file,path,card,folder,mh):
  m = '%5.1f'%mh
  m = m.replace('.0','')
  if 'ExpProfileLikelihood' in options.methods:
    file.write('if ( \n')
    file.write('\t combine %s/%s -M ProfileLikelihood -m %s -n %sExp --signif --pval -t -1 --expectSignal=1\n'%(path,card,m,folder))  
    file.write(') then\n')
    file.write('\t mv higgsCombine%sExp.ProfileLikelihood.mH%s.root %s/%s/ExpProfileLikelihood\n'%(folder,m,path,folder))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')
  if 'ExpProfileLikelihoodm125' in options.methods and options.parametric:
    file.write('if ( \n')
    file.write('\t combine %s/%s -M ProfileLikelihood -m %s -n %sExpm125 --signif --pval -t -1 --expectSignal=1 --expectSignalMass=125\n'%(path,card,m,folder))  
    file.write(') then\n')
    file.write('\t mv higgsCombine%sExpm125.ProfileLikelihood.mH%s.root %s/%s/ExpProfileLikelihoodm125\n'%(folder,m,path,folder))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')
  if 'ProfileLikelihood' in options.methods and options.unblind:
    file.write('if ( \n')
    file.write('\t combine %s/%s -M ProfileLikelihood -m %s -n %s --signif --pval\n'%(path,card,m,folder))  
    file.write(') then\n')
    file.write('\t mv higgsCombine%s.ProfileLikelihood.mH%s.root %s/%s/ProfileLikelihood\n'%(folder,m,path,folder))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')

def writeProfileLikelhoodSplitCard(file,path,splitcard,folder,mh):
  m = '%5.1f'%mh
  m = m.replace('.0','')
  if 'ExpProfileLikelihood' in options.methods:
    file.write('if ( \n')
    file.write('\t combine %s/%s -M ProfileLikelihood -m %s -n %sExp%s --signif --pval -t -1 --expectSignal=1\n'%(path,splitcard[1],m,folder,splitcard[0]))  
    file.write(') then\n')
    file.write('\t mv higgsCombine%sExp%s.ProfileLikelihood.mH%s.root %s/%s/ExpProfileLikelihood/%s\n'%(folder,splitcard[0],m,path,folder,splitcard[0]))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')
  if 'ProfileLikelihood' in options.methods and options.unblind:
    file.write('if ( \n')
    file.write('\t combine %s/%s -M ProfileLikelihood -m %s -n %s%s --signif --pval\n'%(path,splitcard[1],m,folder,splitcard[0]))  
    file.write(') then\n')
    file.write('\t mv higgsCombine%s%s.ProfileLikelihood.mH%s.root %s/%s/ProfileLikelihood/%s\n'%(folder,splitcard[0],m,path,folder,splitcard[0]))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')

def writeMaxLikelihoodFit(file,path,card,folder,mh):
  m = '%5.1f'%mh
  m = m.replace('.0','')
  if 'MaxLikelihoodFit' in options.methods and options.unblind:
    file.write('if ( \n')
    #file.write('\t combine %s/%s -M MaxLikelihoodFit --minimizerAlgo=Minuit -m %s -n %s --rMin=-3. --rMax=3.\n'%(path,card,m,folder))
    file.write('\t combine %s/%s -M MaxLikelihoodFit -m %s -n %s --rMin=-3. --rMax=3.\n'%(path,card,m,folder))
    file.write(') then\n')
    file.write('\t mv higgsCombine%s.MaxLikelihoodFit.mH%s.root %s/%s/MaxLikelihoodFit\n'%(folder,m,path,folder))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')

def writeChannelCompatibility(path,card,folder,mh):
  m = '%5.1f'%mh
  m = m.replace('.0','')
  file = open('%s/subChanComp.sh'%(folder),'w')
  if not options.dryRun:
    os.system('rm -f %s.done\n'%file.name)
    os.system('rm -f %s.log\n'%file.name)
    os.system('rm -f %s.fail\n'%file.name)
    os.system('rm -f %s.run\n'%file.name)
  writePreamble(file,path,card,mh,'mu Scan')
  file.write('if ( \n')
  file.write('\t combine %s/%s -M ChannelCompatibilityCheck -m %s -n %s --rMin=-25 --verbose=1 --saveFitResult\n'%(path,card,m,folder))
  file.write(') then\n')
  file.write('\t mv higgsCombine%s.ChannelCompatibilityCheck.mH%s.root %s/%s/Specials\n'%(folder,m,path,folder))
  file.write('\t touch %s.done\n'%file.name)
  file.write('else\n')
  file.write('\t touch %s.fail\n'%file.name)
  file.write('fi\n')
  file.write('rm %s.run\n'%file.name)
  file.close()
  os.system('chmod +x %s'%file.name)
  if not options.dryRun: os.system('bsub -q 1nh -o %s/%s.log %s/%s'%(path,file.name,path,file.name))

def rejigVal(ws,name,cent,down,up):
  ws.var(name).setVal(cent)
  ws.var(name).setMin(down)
  ws.var(name).setMax(up)
  if options.verbose: ws.var(name).Print()

def writeNLLmuScan(path,card,folder,mh):
  m = '%5.1f'%mh
  m = m.replace('.0','')
  file = open('%s/subMuScan.sh'%(folder),'w')
  if not options.dryRun:
    os.system('rm -f %s.done\n'%file.name)
    os.system('rm -f %s.log\n'%file.name)
    os.system('rm -f %s.fail\n'%file.name)
    os.system('rm -f %s.run\n'%file.name)
  writePreamble(file,path,card,mh,'mu Scan')
  print 'Creating workspace for NLL mu scan...'
  os.system('text2workspace.py %s/%s -o %s/floatingMu.root'%(path,card,path))
  file.write('if ( \n')
  file.write('\t combine %s/floatingMu.root -M MultiDimFit -m %s --rMin=-2 --rMax=4 --points=100 -n %sMuSyst\n'%(path,m,folder))
  file.write(') then\n')
  file.write('\t mv higgsCombine%sMuSyst.MultiDimFit.mH%s.root %s/%s/Specials\n'%(folder,m,path,folder))
  file.write('\t touch %s.done\n'%file.name)
  file.write('else\n')
  file.write('\t touch %s.fail\n'%file.name)
  file.write('fi\n')
  file.write('if ( \n')
  file.write('\t combine %s/floatingMu.root -M MultiDimFit -m %s --rMin=-2 --rMax=4 --points=100 -n %sMuStat --fastScan\n'%(path,m,folder))
  file.write(') then\n')
  file.write('\t mv higgsCombine%sMuStat.MultiDimFit.mH%s.root %s/%s/Specials\n'%(folder,m,path,folder))
  file.write('\t touch %s.done\n'%file.name)
  file.write('else\n')
  file.write('\t touch %s.fail\n'%file.name)
  file.write('fi\n')
  file.write('rm %s.run\n'%file.name)
  file.close()
  os.system('chmod +x %s'%file.name)
  if not options.dryRun: os.system('bsub -q 1nh -o %s/%s.log %s/%s'%(path,file.name,path,file.name))

def writeNLLmassScan(path,card,folder,mh):
  ml = int(mh-5)
  mu = int(mh+5)
  m = '%5.1f'%mh
  m = m.replace('.0','')
  file = open('%s/subMHScan.sh'%(folder),'w')
  if not options.dryRun:
    os.system('rm -f %s.done\n'%file.name)
    os.system('rm -f %s.log\n'%file.name)
    os.system('rm -f %s.fail\n'%file.name)
    os.system('rm -f %s.run\n'%file.name)
  writePreamble(file,path,card,mh,'MH Scan')
  print 'Creating workspace for NLL mass scan...'
  os.system('text2workspace.py %s/%s -o %s/floatingMass.root -m %s -P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs --PO higgsMassRange=%d,%d\n'%(path,card,path,m,ml,mu))
  file.write('if ( \n')
  file.write('\t combine %s/floatingMass.root --preFitValue=0. --saveNLL --floatOtherPOIs=1 -M MultiDimFit -m %s --algo=grid --points=100 -n %sMassSyst -P MH\n'%(path,m,folder))
  file.write(') then\n')
  file.write('\t mv higgsCombine%sMassSyst.MultiDimFit.mH%s.root %s/%s/Specials\n'%(folder,m,path,folder))
  file.write('\t touch %s.done\n'%file.name)
  file.write('else\n')
  file.write('\t touch %s.fail\n'%file.name)
  file.write('fi\n')
  file.write('if ( \n')
  file.write('\t combine %s/floatingMass.root --preFitValue=0. --saveNLL --floatOtherPOIs=1 -M MultiDimFit -m %s --algo=grid --points=100 -n %sMassStat -P MH --fastScan\n'%(path,m,folder))
  file.write(') then\n')
  file.write('\t mv higgsCombine%sMassStat.MultiDimFit.mH%s.root %s/%s/Specials\n'%(folder,m,path,folder))
  file.write('\t touch %s.done\n'%file.name)
  file.write('else\n')
  file.write('\t touch %s.fail\n'%file.name)
  file.write('fi\n')
  file.write('rm %s.run\n'%file.name)
  file.close()
  os.system('chmod +x %s'%file.name)
  if not options.dryRun: os.system('bsub -q %s -o %s/%s.log %s/%s'%(options.queue,path,file.name,path,file.name))

def writeNLLmumhScan(path,card,folder,mh,npoints,njobs):
  ml = int(mh-5)
  mu = int(mh+5)
  m = '%5.1f'%mh
  m = m.replace('.0','')
  pointsperjob = int(npoints/njobs)
  print 'Creating workspace for NLL mu vs mass scan...'
  os.system('text2workspace.py %s/%s -o %s/floatingMuMH.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass'%(path,card,path))
  tf = r.TFile('%(path)s/floatingMuMH.root'%locals(),'UPDATE')
  ws = tf.Get('w')
  rejigVal(ws,'MH',mh,ml,mu)
  rejigVal(ws,'r',1.,0.,2.5)
  ws.Write(ws.GetName(),r.TObject.kWriteDelete)
  tf.Close()
  for j in range(njobs):
    file = open('%s/subMuMHScan_%d.sh'%(folder,j),'w')
    if not options.dryRun:
      os.system('rm -f %s.done\n'%file.name)
      os.system('rm -f %s.log\n'%file.name)
      os.system('rm -f %s.fail\n'%file.name)
      os.system('rm -f %s.run\n'%file.name)
    writePreamble(file,path,card,mh,'MuMH Scan job %d'%j)
    firstPoint = j*pointsperjob
    lastPoint = (j*pointsperjob)+pointsperjob-1
    file.write('if ( \n')
    file.write('\tcombine %s/floatingMuMH.root --preFitValue=0. --saveNLL -P r -P MH -M MultiDimFit -m %s --algo=grid --points=%d --firstPoint=%d --lastPoint=%d -n %sMuMHScan_%d\n'%(path,m,npoints,firstPoint,lastPoint,folder,j))
    file.write(') then\n')
    file.write('\t mv higgsCombine%sMuMHScan_%d.MultiDimFit.mH%s.root %s/%s/Specials\n'%(folder,j,m,path,folder))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')
    file.write('rm %s.run\n'%file.name)
    file.close()
    os.system('chmod +x %s'%file.name)
    if not options.dryRun: os.system('bsub -q %s -o %s/%s.log %s/%s'%(options.queue,path,file.name,path,file.name))

def writeNLLrvrfScan(path,card,folder,mh,npoints,njobs):
  m = '%5.1f'%mh
  m = m.replace('.0','')
  pointsperjob = int(npoints/njobs)
  print 'Creating workspace for NLL Rv vs Rf scan...'
  os.system('text2workspace.py %(path)s/%(card)s -o %(path)s/floatingRvRf.root -m %(m)s -P HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs\n'%locals())
  tf = r.TFile('%(path)s/floatingRvRf.root'%locals(),'UPDATE')
  ws = tf.Get('w')
  rejigVal(ws,'RV',1.,-5,5.)
  rejigVal(ws,'RF',1.,-5,5.)
  ws.Write(ws.GetName(),r.TObject.kWriteDelete)
  tf.Close()
  for j in range(njobs):
    file = open('%s/subRvRfScan_%d.sh'%(folder,j),'w')
    if not options.dryRun:
      os.system('rm -f %s.done\n'%file.name)
      os.system('rm -f %s.log\n'%file.name)
      os.system('rm -f %s.fail\n'%file.name)
      os.system('rm -f %s.run\n'%file.name)
    writePreamble(file,path,card,mh,'RvRf Scan job %d'%j)
    firstPoint = j*pointsperjob
    lastPoint = (j*pointsperjob)+pointsperjob-1
    file.write('if ( \n')
    file.write('\tcombine %s/floatingRvRf.root --preFitValue=0. --saveNLL -M MultiDimFit -m %s --algo=grid --points=%d --firstPoint=%d --lastPoint=%d -n %sRvRfScan_%d | tee combine_%s_grid%d.log \n'%(path,m,npoints,firstPoint,lastPoint,folder,j,m,firstPoint))
    file.write(') then\n')
    file.write('\t mv higgsCombine%sRvRfScan_%d.MultiDimFit.mH%s.root %s/%s/Specials\n'%(folder,j,m,path,folder))
    file.write('\t mv  combine_%s_grid%d.log %s/%s/Specials\n'%(m,firstPoint,path,folder))
    file.write('\t touch %s.done\n'%file.name)
    file.write('else\n')
    file.write('\t touch %s.fail\n'%file.name)
    file.write('fi\n')
    file.write('rm %s.run\n'%file.name)
    file.close()
    os.system('chmod +x %s'%file.name)
    if not options.dryRun: os.system('bsub -q %s -o %s/%s.log %s/%s'%(options.queue,path,file.name,path,file.name))
  
# MAIN STARTS HERE
# if a specific year has not been set run everything
if not options.run2011 and not options.run2012 and not options.runBoth:
  options.run2011=True
  options.run2012=True
  options.runBoth=True

# figure out which methods to run
allowed_methods = ['Asymptotic','ProfileLikelihood','ExpProfileLikelihood','MaxLikelihoodFit']
if options.parametric:
  allowed_methods.append('ExpProfileLikelihoodm125')
# if no specific combine methods then do all
if len(options.methods)==0:
  options.methods = allowed_methods
else:
  for option in options.methods:
    if option not in allowed_methods:
      sys.exit('%s is not in allowed_methods: %s'%(option,allowed_methods))

# auto generate the combined card if it hasn't been specified
if options.runBoth and not options.cardBoth:
  options.cardBoth = 'datacard_autocombined.txt'
  print 'Auto generating combined card - ', options.cardBoth
  os.system("combineCards.py %s %s > %s"%(options.card2011,options.card2012,options.cardBoth))

# setup some common loops over the years
config=[]
if options.run2011:
  config.append([options.card2011,'7TeV',options.splitChannels2011])
if options.run2012:
  config.append([options.card2012,'8TeV',options.splitChannels2012])
if options.runBoth:
  config.append([options.cardBoth,'7and8TeV',options.splitChannelsBoth])

# loop over the years (i.e 2011, 2012 and 2011+2012)
for card, folder, splitChannels in config:
  splitCards = cardSplitting(card,splitChannels)
  os.system('mkdir -p %s'%folder)
  for method in options.methods: os.system('mkdir -p %s/%s'%(folder,method))
  for splitcard in splitCards:
    if 'ProfileLikelihood' in options.methods: os.system('mkdir -p %s/ProfileLikelihood/%s'%(folder,splitcard[0]))
    if 'ExpProfileLikelihood' in options.methods: os.system('mkdir -p %s/ExpProfileLikelihood/%s'%(folder,splitcard[0]))

  # loop over each hypothesised mass and make a sub script for each
  for m in numpy.arange(options.mHlow,options.mHhigh+options.mHstep,options.mHstep):
    f = open('%s/%s/sub_m%5.1f.sh'%(path,folder,m),'w')
    if not options.dryRun:
      os.system('rm -f %s.done\n'%f.name)
      os.system('rm -f %s.log\n'%f.name)
      os.system('rm -f %s.fail\n'%f.name)
      os.system('rm -f %s.run\n'%f.name)
    writePreamble(f,path,card,m)
    
    # write out the combine commands
    if not options.skipDatacard:
      writeAsymptotic(f,path,card,folder,m)
      writeProfileLikelhood(f,path,card,folder,m)
      writeMaxLikelihoodFit(f,path,card,folder,m)
    for splitcard in splitCards:
      writeProfileLikelhoodSplitCard(f,path,splitcard,folder,m)

    f.write('rm %s.run\n'%f.name)
    f.close()
    os.system('chmod +x %s'%f.name)

    if not options.dryRun and (not options.skipDatacard or len(splitCards)!=0): os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))
  
  # these scripts should really be run locally
  if options.doSpecials:
    os.system('mkdir -p %s/%s/Specials'%(path,folder))
    m=options.doSpecials
    import ROOT as r
    r.gSystem.Load('libHiggsAnalysisCombinedLimit')
    writeChannelCompatibility(path,card,folder,m)
    writeNLLmuScan(path,card,folder,m)
    writeNLLmassScan(path,card,folder,m)
    writeNLLmumhScan(path,card,folder,m,2000,20)
    writeNLLrvrfScan(path,card,folder,m,2000,20)

