#!/usr/bin/env python
import os

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--datacard",dest="datacard")
parser.add_option("-m","--mass",dest="mass",type=float)
parser.add_option("-D","--directory",dest="directory",default="jobs")
parser.add_option("-n","--njobs",dest="njobs",type=int)
parser.add_option("-t","--toysperjob",dest="toysperjob",type=int)
parser.add_option("-w","--workspaceName",dest="wsname",default="floatMu.root")
parser.add_option("-s","--skipWorkspace",dest="skipWorkspace",action="store_true",default=False)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
(options,args)=parser.parse_args()

if not options.skipWorkspace:
  print 'Using combine tool to make workspace from card...'
  os.system('text2workspace.py -m %3.1f %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=muFloating -o %s'%(options.mass,options.datacard,options.wsname))

os.system('mkdir -p %s'%options.directory)
options.wsname = os.getcwd()+'/'+options.wsname
print 'Making jobs...'

for j in range(options.njobs):
  f = open('%s/%s/sub%d.sh'%(os.getcwd(),options.directory,j),'w')
  f.write('#!/bin/bash\n')
  f.write('cd %s/%s\n'%(os.getcwd(),options.directory))
  f.write('rm -f %s.done'%f.name)
  f.write('rm -f %s.fail'%f.name)
  f.write('touch %s.run\n'%f.name)

  f.write('eval `scramv1 runtime -sh`\n')
  f.write('if ( combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys -s 0 -n job%d );\n'%(options.wsname,options.mass,options.toysperjob,j))
  f.write('\tthen touch %s.done\n'%f.name)
  f.write('\telse touch %s.fail\n'%f.name)
  f.write('fi\n')
  f.write('rm -f %s.run'%f.name)
  os.system('chmod +x %s'%f.name)

  if not options.dryRun:
    os.system('bsub -q 1nh -o %s.log %s'%(f.name,f.name))
  
 
# fqq scan 
#os.system('text2workspace.py -m %3.1f %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating -o %s'%(options.mass,options.datacard,options.wsname))
#f.write('if ( combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameter fqq=0.025 --freezeNuisances fqq -s 0 -n job%d );\n'%(options.wsname,options.mass,options.toysperjob,j))

# compatibility
# combineCards.py -xc "ALT" datacard_spinanalysis
#combine datacard_smonly.txt -M GenerateOnly -m 125 -T 1000 --expectSignal=1 -n toyFile --saveToys
#combine datacard_smonly.txt -M ChannelCompatibilityCheck -m 125 -T 1000 --toysFile -g 
#combine datacard_gravonly.txt -M ChannelCompatibilityCheck -m 125 -T 1000 --toysFile -g 
