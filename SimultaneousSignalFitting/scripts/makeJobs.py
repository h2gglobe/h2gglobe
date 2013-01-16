#!/usr/bin/env python

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-i","--inputfile",dest="inputfile",type="str",default="dat/config.dat",help="Input config file")
parser.add_option("-o","--outputdir",dest="outputdir",type="str",help="Output directory location")
parser.add_option("-f","--fork",dest="fork",type="int",default=8,help="Fork NLL calcs over multiple CPUs")
parser.add_option("-r","--recursive",dest="recursive",action="store_true",default=False)
parser.add_option("-q","--queue",dest="queue",type="str",default="1nh")
parser.add_option("--submit",action="store_true",dest="submit",default=False)
(options,args)=parser.parse_args()

import os

os.system('mkdir -p %s'%options.outputdir)
os.system('mkdir -p scripts/jobs')

combinefile = open('scripts/jobs/filestocombine.dat','w')

datfile = open(options.inputfile)
for line in datfile.readlines():
  if line.startswith('#') or line=='\n': continue
  if line.startswith('infile'):
    inWSfile=line.split('=')[1].strip()
    print 'Input file:'
    print '\t', inWSfile
    print 'Output file location:'
    print '\t', options.outputdir
    print 'Output files:'
  
  line_els  = line.split()
  if len(line_els)!=7: continue
  print '\t', line,
  outfile   = os.path.join(options.outputdir,line_els[0])
  proc      = line_els[1]
  cat       = int(line_els[2])
  nGaus     = int(line_els[3])
  dmOrder   = int(line_els[4])
  sigOrder  = int(line_els[5])
  fracOrder = int(line_els[6])
  combinefile.write('%s %s_cat%d\n'%(outfile,proc,cat))
  script = open('%s/scripts/jobs/sub_%s_cat%d.sh'%(os.getcwd(),proc,cat),'w')
  os.system('rm -f %s.done'%(script.name))
  os.system('rm -f %s.fail'%(script.name))
  os.system('rm -f %s.run'%(script.name))
  os.system('rm -f %s.log'%(script.name))
  script.write('#!/bin/bash\n')
  script.write('cd %s\n'%os.getcwd())
  script.write('eval `scramv1 runtime -sh`\n')
  script.write('touch %s.run\n'%(script.name))
  script.write('if ( ./bin/SimultaneousSignalFit -i %s -o %s -p %s -c %d -g %d --dmOrder %d --sigmaOrder %d --fracOrder %d '%(inWSfile,outfile,proc,cat,nGaus,dmOrder,sigOrder,fracOrder))
  if options.fork: script.write('--fork %d '%options.fork)
  if options.recursive: script.write('--recursive ')
  script.write(') then \n')
  script.write('\ttouch %s.done\n'%(script.name))
  script.write('\trm -f %s.run\n'%(script.name))
  script.write('else\n')
  script.write('\ttouch %s.fail\n'%(script.name))
  script.write('fi\n')
  script.close()
  os.system('chmod +x %s'%(script.name))
  #print 'bsub -q %s -o %s.log %s'%(options.queue,script.name,script.name)
  if options.submit: os.system('bsub -q %s -o %s.log %s'%(options.queue,script.name,script.name))

if not options.submit: print 'To submit run with --submit'
