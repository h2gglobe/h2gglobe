#!/usr/bin/env python
import os
import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-M","--method",dest="methods",default=[],action="append")
parser.add_option("-d","--datacard",dest="datacard")
parser.add_option("-m","--mass",dest="mass",type=float)
parser.add_option("-D","--directory",dest="directory",default="jobs")
parser.add_option("-f","--fitNuis",dest="fitNuis",type="int",default=0)
parser.add_option("-n","--njobs",dest="njobs",type=int)
parser.add_option("-t","--toysperjob",dest="toysperjob",type=int)
parser.add_option("-w","--workspaceName",dest="wsname",default="floatMu.root")
parser.add_option("-s","--skipWorkspace",dest="skipWorkspace",action="store_true",default=False)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
(options,args)=parser.parse_args()

import ROOT as r

os.system('mkdir -p %s'%options.directory)
options.wsname = os.getcwd()+'/'+options.directory+'/'+options.wsname
print options.wsname
if not options.skipWorkspace:
  print 'Using combine tool to make workspace from card', options.datacard
  os.system('text2workspace.py -m %s %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating -o %s'%(options.mass,options.datacard,options.wsname))

massString = ('%3.1f'%options.mass).replace('.0','')

def doChannelCompability():
  print 'Running channel compatibility..'
  f = open('%s/%s/chcomp.sh'%(os.getcwd(),options.directory),'w')
  f.write('#!/bin/bash\n')
  f.write('cd %s/%s\n'%(os.getcwd(),options.directory))
  f.write('rm -f %s.done\n'%f.name)
  f.write('rm -f %s.fail\n'%f.name)
  f.write('touch %s.run\n'%f.name)
  f.write('eval `scramv1 runtime -sh`\n')
 
  if options.fitNuis==0:
    # throw SM asimov
    f.write('if ( combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=0.0,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r,fqq,x --expectSignal=1 --toysFrequentist -n SMAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
    # throw grav_gg asimov
    f.write('  && combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=0.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r,fqq,x --expectSignal=1 --toysFrequentist -n GravGGAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
    # throw gqq asimov
    f.write('  && combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=1.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r,fqq,x --expectSignal=1 --toysFrequentist -n GravQQAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
    # throw g50gg50qq asimov
    f.write('  && combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=0.5,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r,fqq,x --expectSignal=1 --toysFrequentist -n Grav50GG50QQAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
  
  elif options.fitNuis==1:
    # throw SM asimov
    f.write('if ( combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=0.0,x=0 --freezeNuisances fqq,x --redefineSignalPOIs fqq,x --expectSignal=0 --toysFrequentist -n SMAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
    # throw grav_gg asimov
    f.write('  && combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=0.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs fqq,x --expectSignal=0 --toysFrequentist -n GravGGAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
    # throw grav_qq asimov
    f.write('  && combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=1.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs fqq,x --expectSignal=1 --toysFrequentist -n GravQQAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
    # throw grav_50gg50qq asimov
    f.write('  && combine %s -M GenerateOnly -m %s -t -1 --setPhysicsModelParameters fqq=0.5,x=1 --freezeNuisances fqq,x --redefineSignalPOIs fqq,x --expectSignal=0.5 --toysFrequentist -n Grav50GG50QQAsimov --saveToys -s 0\\\n'%(options.wsname,massString)) 
  
  else:
    sys.exit('fitNuis value %d not recognised'%options.fitNuis)
   
  # FIT BACK SM 
  # fit SM asimov 
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineSMAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n SMtoSM --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_gg asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGravGGAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n SMtoGravGG --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_qq asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGravQQAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n SMtoGravQQ --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_50gg50qq asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGrav50GG50QQAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n SMtoGrav50GG50QQ --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit data
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s --setPhysicsModelParameters fqq=0.0,x=0 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n SMtoData --saveFitResult\\\n'%(options.wsname,massString)) 

  # FIT BACK GRAV_GG 
  # fit SM asimov 
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineSMAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravGGtoSM --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_gg asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGravGGAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravGGtoGravGG --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_qq asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGravQQAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravGGtoGravQQ --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_50gg50qq asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGrav50GG50QQAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=0.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravGGtoGrav50GG50QQ --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit data
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s --setPhysicsModelParameters fqq=0.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravGGtoData --saveFitResult\\\n'%(options.wsname,massString)) 

  # FIT BACK GRAV_QQ 
  # fit SM asimov 
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineSMAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=1.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravQQtoSM --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_gg asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGravGGAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=1.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravQQtoGravGG --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_qq asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGravQQAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=1.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravQQtoGravQQ --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit grav_50gg50qq asimov
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s -t -1 --toysFile higgsCombineGrav50GG50QQAsimov.GenerateOnly.mH%s.0.root --setPhysicsModelParameters fqq=1.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravQQtoGrav50GG50QQ --saveFitResult\\\n'%(options.wsname,massString,massString))
  # fit data
  f.write('  && combine %s -M ChannelCompatibilityCheck -m %s --setPhysicsModelParameters fqq=1.0,x=1 --freezeNuisances fqq,x --redefineSignalPOIs r --setPhysicsModelParameterRanges r=-10.,10. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravQQtoData --saveFitResult\\\n'%(options.wsname,massString)) 

  
  f.write(')\n')
  f.write('\tthen touch %s.done\n'%f.name)
  f.write('\telse touch %s.fail\n'%f.name)
  f.write('fi\n')
  f.write('rm -f %s.run\n'%f.name)
  os.system('chmod +x %s'%f.name)
  f.close()
  if not options.dryRun: 
    os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))

def doqqbar():
  print 'Running qqbar...'
  print 'Making jobs...'

  # do best fit nll scan
  bf = open('%s/%s/subBFqqbar.sh'%(os.getcwd(),options.directory),'w')
  bf.write('#!/bin/bash\n')
  bf.write('cd %s/%s\n'%(os.getcwd(),options.directory))
  bf.write('rm -f %s.done\n'%bf.name)
  bf.write('rm -f %s.fail\n'%bf.name)
  bf.write('touch %s.run\n'%bf.name)
  bf.write('eval `scramv1 runtime -sh`\n')
  bf.write('if ( combine %s -M MultiDimFit -m %s --redefineSignalPOIs fqq --setPhysicsModelParameters x=1.,MH=%s --setPhysicsModelParameterRanges fqq=0.,1. --algo grid --points 100 )\n'%(options.wsname,massString,massString)) 
  bf.write('\tthen touch %s.done\n'%bf.name)
  bf.write('\telse touch %s.fail\n'%bf.name)
  bf.write('fi\n')
  bf.write('rm -f %s.run\n'%bf.name)
  os.system('chmod +x %s'%bf.name)
  if not options.dryRun:
    os.system('bsub -q 8nh -o %s.log %s'%(bf.name,bf.name))

  # do likelihood ratio at fqq points
  for j in range(options.njobs):
    f = open('%s/%s/sub%dqqbar.sh'%(os.getcwd(),options.directory,j),'w')
    f.write('#!/bin/bash\n')
    f.write('cd %s/%s\n'%(os.getcwd(),options.directory))
    f.write('rm -f %s.done\n'%f.name)
    f.write('rm -f %s.fail\n'%f.name)
    f.write('touch %s.run\n'%f.name)
    f.write('eval `scramv1 runtime -sh`\n')
    
    f.write('if ( combine %s -m %s -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.0 --freezeNuisances fqq -s 0 -n job%dfqq0.00 '%(options.wsname,massString,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %s -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.25 --freezeNuisances fqq -s 0 -n job%dfqq0.25 '%(options.wsname,massString,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %s -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.50 --freezeNuisances fqq -s 0 -n job%dfqq0.50 '%(options.wsname,massString,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %s -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.75 --freezeNuisances fqq -s 0 -n job%dfqq0.75 '%(options.wsname,massString,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %s -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=1.00 --freezeNuisances fqq -s 0 -n job%dfqq1.00 '%(options.wsname,massString,options.fitNuis,options.toysperjob,j))
    f.write(' )\n')

    f.write('\tthen touch %s.done\n'%f.name)
    f.write('\telse touch %s.fail\n'%f.name)
    f.write('fi\n')
    f.write('rm -f %s.run\n'%f.name)
    os.system('chmod +x %s'%f.name)

    if not options.dryRun:
      os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))

if len(options.methods)==0: options.methods=['ChannelCompatibility','qqbar']
if 'ChannelCompatibility' in options.methods: doChannelCompability()
if 'qqbar' in options.methods: doqqbar()

# fqq scan 
#os.system('text2workspace.py -m %s %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating -o %s'%(massString,options.datacard,options.wsname))
#f.write('if ( combine %s -m %s -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameter fqq=0.025 --freezeNuisances fqq -s 0 -n job%d );\n'%(options.wsname,massString,options.toysperjob,j))

# compatibility
# combineCards.py -xc "ALT" datacard_spinanalysis
#combine combine datacard_spinanalysis_justSM.txt -M GenerateOnly -m 125 -t 4000 --expectSignal=1 -n SMToys --saveToys -s -1 
#combine datacard_spinanalysis_justSM.txt -M ChannelCompatibilityCheck -m 125 -t 4000 --toysFile -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 
#combine datacard_spinanalysis_justGrav.txt -M ChannelCompatibilityCheck -m 125 -t 4000 --toysFile -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4


# CAN DO IT ALL FROM QQBAR workspace "
# e.g combine floatMu_qqbar.root -M MaxLikelihoodFit --setPhysicsModelParameters fqq=0.0,x=0 -m 125 --rMin=-5. --freezeNuisances fqq,x --redefineSignalPOIs r



