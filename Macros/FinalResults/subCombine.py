#!/usr/bin/env python

# e.g
# ./subCombine.py --card2011=datacard_parametric_7TeV_massfac.txt --unblind --run2011 --splitChannels2011="cat0_7TeV cat1_7TeV cat2_7TeV cat3_7TeV" --splitChannels2011="cat4_7TeV"
# ./subCombine.py --card2012=datacard_parametric_8TeV_cutbased_cbvbf.txt --unblind --run2012 --splitChannels2012="cat0 cat1 cat2 cat3" --splitChannels2012="cat4 cat5" --splitChannels2012="cat6 cat7" --splitChannels2012="cat8"
# ./subCombine.py --cardBoth=hgg.combined.txt --unblind --runBoth --splitChannelsBoth="ch1_cat0_7TeV ch1_cat1_7TeV ch1_cat2_7TeV ch1_cat3_7TeV ch2_cat0_8TeV ch2_cat1_8TeV ch2_cat2_8TeV ch2_cat3_8TeV" --splitChannelsBoth="ch1_cat4_7TeV ch2_cat4_8TeV ch2_cat5_8TeV" --splitChannelsBoth="ch2_cat6_8TeV ch2_cat7_8TeV ch2_cat8_8TeV"

# one command to rule them all
# ./subCombine.py --card2011=datacard_parametric_7TeV_massfac.txt --unblind --run2011 --splitChannels2011="cat0 cat1 cat2 cat3" --splitChannels2011="cat4" --card2012=datacard_parametric_8TeV_cutbased_cbvbf.txt --unblind --run2012 --splitChannels2012="cat0 cat1 cat2 cat3" --splitChannels2012="cat4 cat5" --splitChannels2012="cat6 cat7" --splitChannels2012="cat8" --cardBoth=hgg.combined.txt --unblind --runBoth --splitChannelsBoth="ch1_cat0_7TeV ch1_cat1_7TeV ch1_cat2_7TeV ch1_cat3_7TeV ch2_cat0_8TeV ch2_cat1_8TeV ch2_cat2_8TeV ch2_cat3_8TeV" --splitChannelsBoth="ch1_cat4_7TeV ch2_cat4_8TeV ch2_cat5_8TeV" --splitChannelsBoth="ch2_cat6_8TeV ch2_cat7_8TeV ch2_cat8_8TeV"

import os
import numpy
import sys

from optparse import OptionParser

parser = OptionParser()
parser.add_option("","--card2011",dest="card2011",type="string",help="Name of 2011 card")
parser.add_option("","--card2012",dest="card2012",type="string",help="Name of 2012 card")
parser.add_option("","--cardBoth",dest="cardBoth",type="string",help="Name of combined card (default is to combine the two cards given)")
parser.add_option("","--unblind",dest="unblind",default=False,action="store_true",help="Unblind")
parser.add_option("-q","--queue",dest="queue",type="string",default="8nm",help="queue")
parser.add_option("-M","--method",dest="method",type="string",help="Combine method (default does several)")
parser.add_option("-L","--mHlow",dest="mHlow",type="float",default=110,help="Lowest MH")
parser.add_option("-H","--mHhigh",dest="mHhigh",type="float",default=150,help="Highest MH")
parser.add_option("-S","--mHstep",dest="mHstep",type="float",default=1,help="MH step")
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

# if a specific year has not been set run everything
if not options.run2011 and not options.run2012 and not options.runBoth:
  options.run2011=True
  options.run2012=True
  options.runBoth=True

if options.runBoth and not options.cardBoth:
  options.cardBoth = 'datacard_autocombined.txt'
  print 'Auto generating combined card - ', options.cardBoth
  os.system("combineCards.py %s %s > %s"%(options.card2011,options.card2012,options.cardBoth))

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

config=[]
if options.run2011:
  config.append([options.card2011,'7TeV',options.splitChannels2011])
if options.run2012:
  config.append([options.card2012,'8TeV',options.splitChannels2012])
if options.runBoth:
  config.append([options.cardBoth,'7and8TeV',options.splitChannelsBoth])

for card, folder, splitChannels in config:
  splitCards = cardSplitting(card,splitChannels)
  os.system('mkdir -p %s'%folder)
  if options.method:
    os.system('mkdir -p %s/%s'%(folder,options.method))
  else:
    os.system('mkdir -p %s/Asymptotic'%folder)
    os.system('mkdir -p %s/ProfileLikelihood'%folder)
    os.system('mkdir -p %s/ExpProfileLikelihood'%folder)
    os.system('mkdir -p %s/MaxLikelihoodFit'%folder)
    os.system('mkdir -p %s/ChannelCompatibility'%folder)
    for splitcard in splitCards:
      os.system('mkdir -p %s/ProfileLikelihood/%s'%(folder,splitcard[0]))
      os.system('mkdir -p %s/ExpProfileLikelihood/%s'%(folder,splitcard[0]))
 
  for m in numpy.arange(options.mHlow,options.mHhigh+options.mHstep,options.mHstep):
    f = open('%s/%s/sub_m%5.1f.sh'%(path,folder,m),'w')
    if not options.dryRun:
      os.system('rm -f %s.done\n'%f.name)
      os.system('rm -f %s.log\n'%f.name)
      os.system('rm -f %s.fail\n'%f.name)
      os.system('rm -f %s.run\n'%f.name)
    f.write('#!/bin/bash\n')
    f.write('cd %s\n'%path)
    f.write('touch %s.run\n'%f.name)
    f.write('echo --------------------------------------------\n')
    f.write('echo   Running %s at mass %5.1f \n'%(card,m)) 
    f.write('echo --------------------------------------------\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('if ( \n')
    if options.unblind:
      if not options.skipDatacard:
        f.write('\t combine %s/%s -M Asymptotic -m %5.1f -n %s\n'%(path,card,m,folder))  
        f.write('\t combine %s/%s -M ProfileLikelihood -m %5.1f -n %s --signif --pval\n'%(path,card,m,folder))  
        f.write('\t combine %s/%s -M ProfileLikelihood -m %5.1f -n %sExp --signif --pval -t -1 --expectSignal=1\n'%(path,card,m,folder))  
        f.write('\t combine %s/%s -M MaxLikelihoodFit -m %5.1f -n %s --rMin=-3. --rMax=3.\n'%(path,card,m,folder))
        f.write('\t combine %s/%s -M ChannelCompatibilityCheck -m %5.1f -n %sChComp --rMin=-25 --verbose=1 --saveFitResult -s -1\n'%(path,card,m,folder))
      for splitcard in splitCards:
        f.write('\t combine %s/%s -M ProfileLikelihood -m %5.1f -n %s%s --signif --pval\n'%(path,splitcard[1],m,folder,splitcard[0]))  
        f.write('\t combine %s/%s -M ProfileLikelihood -m %5.1f -n %sExp%s --signif --pval -t -1 --expectSignal=1\n'%(path,splitcard[1],m,folder,splitcard[0]))  
    else:
      if not options.skipDatacard:
        f.write('\t combine %s/%s -M Asymptotic -m %5.1f -n %s --run=expected\n'%(path,card,m,folder))  
        f.write('\t combine %s/%s -M ProfileLikelihood -m %5.1f -n %sExp --signif --pval -t -1 --expectSignal=1\n'%(path,card,m,folder))  
      for splitcard in splitCards:
        f.write('\t combine %s/%s -M ProfileLikelihood -m %5.1f -n %sExp%s --signif --pval -t -1 --expectSignal=1\n'%(path,splitcard[1],m,folder,splitcard[0]))  
    f.write(') then\n')
    f.write('\t touch %s.done\n'%f.name)
    f.write('else\n')
    f.write('\t touch %s.fail\n'%f.name)
    f.write('fi\n')
    for splitcard in splitCards:
      f.write('mv higgsCombine%sExp*%s.ProfileLikelihood.*root %s/%s/ExpProfileLikelihood/%s\n'%(folder,splitcard[0],path,folder,splitcard[0]))
      f.write('mv higgsCombine%s*%s.ProfileLikelihood.*root %s/%s/ProfileLikelihood/%s\n'%(folder,splitcard[0],path,folder,splitcard[0]))
    f.write('mv higgsCombine%s*.Asymptotic.*root %s/%s/Asymptotic\n'%(folder,path,folder))
    f.write('mv higgsCombine%sExp*.ProfileLikelihood.*root %s/%s/ExpProfileLikelihood\n'%(folder,path,folder))
    if options.unblind:
      f.write('mv higgsCombine%s*.ProfileLikelihood.*root %s/%s/ProfileLikelihood\n'%(folder,path,folder))
      f.write('mv higgsCombine%s*.MaxLikelihoodFit.*root %s/%s/MaxLikelihoodFit\n'%(folder,path,folder))
      f.write('mv higgsCombine%s*.MaxLikelihoodFit.*root %s/%s/MaxLikelihoodFit\n'%(folder,path,folder))
      f.write('mv higgsCombine%s*.ChannelCompatibilityCheck.*root %s/%s/ChannelCompatibility\n'%(folder,path,folder))
    f.write('rm %s.run\n'%f.name)
    f.close()
    os.system('chmod +x %s'%f.name)
    if not options.dryRun: os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))

  if options.doSpecials:
    m=options.doSpecials
    print 'This will run multiDimFit.py which itself uses parallel - make sure they are both available in this directory.'
    print 'You may want to consider changing some of the options in multiDimFit.py'
    # Best Fit Mass
    os.system('mkdir -p %s/%s/BestFitMass/stat_syst'%(path,folder))
    sf = open('%s/%s/sub_BestFitMass.sh'%(path,folder),'w')
    if not options.dryRun:
      os.system('rm -f %s.done\n'%sf.name)
      os.system('rm -f %s.log\n'%sf.name)
      os.system('rm -f %s.fail\n'%sf.name)
      os.system('rm -f %s.run\n'%sf.name)
    sf.write('#!/bin/bash\n')
    sf.write('cd %s\n'%path)
    sf.write('touch %s.run\n'%sf.name)
    sf.write('echo --------------------------------------------\n')
    sf.write('echo   Running %s for best fit mass \n'%(card)) 
    sf.write('echo --------------------------------------------\n')
    sf.write('eval `scramv1 runtime -sh`\n')
    sf.write('if ( \n')
    sf.write('\t ./multiDimFit.py -w %s/%s/BestFitMass/stat_syst -f -d %s/%s -n 100 --model mH\n'%(path,folder,path,card)) 
    sf.write('\t mv higgsCombine*grid*MultiDimFit*.root %s/%s/BestFitMass/stat_syst/\n'%(path,folder))
    sf.write('\t mv combine*grid*.log %s/%s/BestFitMass/stat_syst/\n'%(path,folder))
    sf.write('\t ./multiDimFit.py -w %s/%s/BestFitMass/stat_syst -f -d %s/%s -n 100 --model mH --statOnly\n'%(path,folder,path,card)) 
    sf.write('\t mv higgsCombine*grid*MultiDimFit*.root %s/%s/BestFitMass/stat_only/\n'%(path,folder))
    sf.write('\t mv combine*grid*.log %s/%s/BestFitMass/stat_only/\n'%(path,folder))
    sf.write(') then\n')
    sf.write('\t touch %s.done\n'%sf.name)
    sf.write('else\n')
    sf.write('\t touch %s.fail\n'%sf.name)
    sf.write('fi\n')
    sf.write('rm %s.run\n'%sf.name)
    sf.close()
    os.system('chmod +x %s'%sf.name)
