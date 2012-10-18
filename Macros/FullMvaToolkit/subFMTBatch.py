#!/usr/bin/env python

import os
import sys

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p","--path",dest="path")
parser.add_option("-i","--infile",dest="inputFile")
parser.add_option("-o","--outfile",dest="outputFile")
parser.add_option("-d","--datfile",dest="datFile")
parser.add_option("-w","--webfolder",dest="web")
parser.add_option("","--dryRun",dest="dryRun",default=False,action="store_true")
(options,args)=parser.parse_args()

dir = options.path
infile = options.inputFile
outfile = options.outputFile
datfile = options.datFile

f = open('batchFMT.sh','w')
f.write('#!/bin/bash\n')
f.write('cd %s\n'%dir)
f.write('eval `scramv1 runtime -sh`\n')
f.write('\n')
f.write('# Clear plot directory\n')
f.write('rm -rf %s/plots/*\n\n'%(dir))

f.write('# histos from trees first \n')
f.write('if ( %s/bin/RunFullMvaAnalysis -i %s/%s -o %s/CMS-HGG_fmtHistos.root --skipRebin --histosFromTrees -v --useDat %s/%s ); \n'%(dir,dir,infile,dir,dir,datfile))
f.write('\tthen echo \'histo file written to CMS-HGG_fmtHistos.root\'\n')
f.write('else\n\ttouch %s/%s.fail \n\texit \nfi\n\n'%(dir,f.name))

f.write('# now do rebinning and fitting\n')
f.write('if ( %s/bin/RunFullMvaAnalysis -i %s/CMS-HGG_fmtHistos.root -o %s/%s --dumpDatFile %s/dat/mvaanalysis_fmtedges.dat --diagnose -v --useDat %s/%s ); \n'%(dir,dir,dir,outfile,dir,dir,datfile))
f.write('\tthen cp %s/%s %s/CMS-HGG_fmtRebinned.root\n'%(dir,outfile,dir))
f.write('\techo \'rebinned file written to %s/%s and copied to %s/CMS-HGG_fmtRebinned.root\'\n'%(dir,outfile,dir))
f.write('else\n\ttouch %s/%s.fail \n\texit \nfi\n\n'%(dir,f.name))

f.write('# now do background model\n')
f.write('if ( %s/bin/RunFullMvaAnalysis -i %s/CMS-HGG_fmtHistos.root -o %s/%s --useDat %s/dat/mvaanalysis_fmtedges.dat --skipRebin --bkgModel --diagnose -v ); \n'%(dir,dir,dir,outfile,dir))
f.write('\tthen cp %s/%s %s/CMS-HGG_fmtLastSuceess.root\n'%(dir,outfile,dir))
f.write('\techo \'corrected background file written to %s/%s and copied to %s/CMS-HGG_fmtLastSuccess.root\'\n'%(dir,outfile,dir))
f.write('else\n\ttouch %s/%s.fail \n\texit \nfi\n\n'%(dir,f.name))

f.write('# now do interpolation\n')
f.write('if ( %s/bin/RunFullMvaAnalysis -i %s/CMS-HGG_fmtHistos.root -o %s/%s --useDat %s/dat/mvaanalysis_fmtedges.dat --skipRebin --interp --diagnose -v ); \n'%(dir,dir,dir,outfile,dir))
f.write('\tthen cp %s/%s %s/CMS-HGG_fmtLastSuceess.root\n'%(dir,outfile,dir))
f.write('\techo \'interpolated file written to %s/%s and copied to %s/CMS-HGG_fmtLastSuccess.root\'\n'%(dir,outfile,dir))
f.write('else\n\ttouch %s/%s.fail \n\texit \nfi\n\n'%(dir,f.name))

f.write('# now do write datacards\n')
f.write('if ( %s/bin/RunFullMvaAnalysis -i %s/CMS-HGG_fmtHistos.root -o %s/%s --useDat %s/dat/mvaanalysis_fmtedges.dat --skipRebin --datacards --diagnose -v ); \n'%(dir,dir,dir,outfile,dir))
f.write('\tthen echo \'datacards written to %s/mva-datacards-grad\'\n'%(dir))
f.write('else\n\ttouch %s/%s.fail \n\texit \nfi\n\n'%(dir,f.name))

f.write('# now make plots\n')
f.write('if ( %s/bin/RunFullMvaAnalysis -i %s/CMS-HGG_fmtHistos.root -o %s/%s --useDat %s/dat/mvaanalysis_fmtedges.dat --skipRebin --doPlot --diagnose -v --www %s ); \n'%(dir,dir,dir,outfile,dir,options.web))
f.write('\tthen cp %s/%s %s/CMS-HGG_fmtComplete.root\n'%(dir,outfile,dir))
f.write('\techo \'complete file written to %s/%s and copied to %s/CMS-HGG_fmtComplete.root\'\n'%(dir,outfile,dir))
f.write('else\n\ttouch %s/%s.fail \n\texit \nfi\n\n'%(dir,f.name))

f.write('# now run combine\n')
f.write('if ( %s/bin/RunFullMvaAnalysis -i %s/CMS-HGG_fmtHistos.root -o %s/%s --useDat %s/dat/mvaanalysis_fmtedges.dat --skipRebin --runCombine --diagnose -v ); \n'%(dir,dir,dir,outfile,dir))
f.write('\tthen echo \'combine successful out put written to %s\''%(dir))
f.write('else\n\ttouch %s/%s.fail \n\texit \nfi\n\n'%(dir,f.name))

f.write('touch %s/%s.success\n\n'%(dir,f.name))

f.close()
os.system('chmod +x %s'%f.name)
if not options.dryRun: os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))


