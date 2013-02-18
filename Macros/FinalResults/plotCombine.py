#!/usr/bin/env python

# e.g. ./plotCombine.py --unblind --run2011 --lumi2011=5.1 --splitChannels2011="cat0_7TeVcat1_7TeVcat2_7TeVcat3_7TeV:Inclusive:4" --splitChannels2011="cat4_7TeV:Dijet Tag:2"
# e.g. /plotCombine.py --unblind --run2012 --lumi2012=19.6 --splitChannels2012="cat0cat1cat2cat3:Inclusive:4" --splitChannels2012="cat4cat5:Dijet Tag:2" --splitChannels2012="cat6cat7:Lepton Tag:3" --splitChannels2012="cat8:MET Tag:6"
# e.g. ./plotCombine.py --unblind --runBoth --lumi2011=5.1 --lumi2012=19.6 --splitChannelsBoth="ch1_cat0_7TeVch1_cat1_7TeVch1_cat2_7TeVch1_cat3_7TeVch2_cat0_8TeVch2_cat1_8TeVch2_cat2_8TeV ch2_cat3_8TeV:Inclusive:4" --splitChannelsBoth="ch1_cat4_7TeVch2_cat4_8TeVch2_cat5_8TeV:Dijet Tag:2" --splitChannelsBoth="ch2_cat6_8TeVch2_cat7_8TeVch2_cat8_8TeV:Lepton/MET Tag:6"

# one command to rule them all
# e.g. ./plotCombine.py --unblind --run2011 --lumi2011=5.1 --splitChannels2011="cat0_7TeVcat1_7TeVcat2_7TeVcat3_7TeV:Inclusive:4" --splitChannels2011="cat4_7TeV:Dijet Tag:2" --run2012 --lumi2012=19.6 --splitChannels2012="cat0cat1cat2cat3:Inclusive:4" --splitChannels2012="cat4cat5:Dijet Tag:2" --splitChannels2012="cat6cat7:Lepton Tag:3" --splitChannels2012="cat8:MET Tag:6" --runBoth --splitChannelsBoth="ch1_cat0_7TeVch1_cat1_7TeVch1_cat2_7TeVch1_cat3_7TeVch2_cat0_8TeVch2_cat1_8TeVch2_cat2_8TeV ch2_cat3_8TeV:Inclusive:4" --splitChannelsBoth="ch1_cat4_7TeVch2_cat4_8TeVch2_cat5_8TeV:Dijet Tag:2" --splitChannelsBoth="ch2_cat6_8TeVch2_cat7_8TeVch2_cat8_8TeV:Lepton/MET Tag:6"


import os
import numpy

from optparse import OptionParser

parser = OptionParser()
parser.add_option("","--unblind",dest="unblind",default=False,action="store_true",help="Unblind")
parser.add_option("","--run2011",dest="run2011",default=False,action="store_true",help="Run only 2011")
parser.add_option("","--run2012",dest="run2012",default=False,action="store_true",help="Run only 2012")
parser.add_option("","--runBoth",dest="runBoth",default=False,action="store_true",help="Run only 2011+2012")
parser.add_option("","--lumi2011",dest="lumi2011",type="float",default=5.3,help="2011 luminosity")
parser.add_option("","--lumi2012",dest="lumi2012",type="float",default=19.6,help="2012 luminosity")
parser.add_option("","--splitChannels2011",dest="splitChannels2011",default=[],action="append",help="Run these channels for 2011 (can pass multiple times) Split with a colon - dir:leg:col")
parser.add_option("","--splitChannels2012",dest="splitChannels2012",default=[],action="append",help="Run these channels for 2012 (can pass multiple times) Split with a colon - dir:leg:col")
parser.add_option("","--splitChannelsBoth",dest="splitChannelsBoth",default=[],action="append",help="Run these channels for Both (can pass multiple times) Split with a colon - dir:leg:col")
parser.add_option("","--skipHadd",dest="skipHadd",default=False,action="store_true",help="Skip hadding of combine output")
parser.add_option("-M","--method",dest="method",type="string",help="Combine method (default does several)")
(options,args)=parser.parse_args()

os.system('mkdir -p plots')

# if a specific year has not been set run everything
if not options.run2011 and not options.run2012 and not options.runBoth:
  options.run2011=True
  options.run2012=True
  options.runBoth=True

config=[]
if options.run2011:
  config.append(['7TeV','#sqrt{s}=7TeV L=%3.1ffb^{-1}'%(options.lumi2011),options.splitChannels2011])
if options.run2012:
  config.append(['8TeV','#sqrt{s}=8TeV L=%3.1ffb^{-1}'%(options.lumi2012),options.splitChannels2012])
if options.runBoth:
  config.append(['7and8TeV','#splitline{#sqrt{s}=7TeV L=%3.1ffb^{-1}}{#sqrt{s}=8TeV L=%3.1ffb^{-1}}'%(options.lumi2011,options.lumi2012),options.splitChannelsBoth])

for folder, lumistring, splitChannels in config:
  if not options.skipHadd:
    if options.method: os.system('hadd -f %s/%s/%s.root %s/%s/higgsCombine*.root'%(folder,options.method,options.method,folder,options.method))
    else:
      os.system('hadd -f %s/Asymptotic/Asymptotic.root %s/Asymptotic/higgsCombine*.root'%(folder,folder))
      os.system('hadd -f %s/ExpProfileLikelihood/ExpProfileLikelihood.root %s/ExpProfileLikelihood/higgsCombine*.root'%(folder,folder))
      for job in splitChannels:
        os.system('hadd -f %s/ExpProfileLikelihood/%s/ExpProfileLikelihood.root %s/ExpProfileLikelihood/%s/higgsCombine*.root'%(folder,job.split(':')[0].replace(' ',''),folder,job.split(':')[0].replace(' ','')))
      if options.unblind:
        os.system('hadd -f %s/ProfileLikelihood/ProfileLikelihood.root %s/ProfileLikelihood/higgsCombine*.root'%(folder,folder))
        os.system('hadd -f %s/MaxLikelihoodFit/MaxLikelihoodFit.root %s/MaxLikelihoodFit/higgsCombine*.root'%(folder,folder))
        for job in splitChannels:
          os.system('hadd -f %s/ProfileLikelihood/%s/ProfileLikelihood.root %s/ProfileLikelihood/%s/higgsCombine*.root'%(folder,job.split(':')[0].replace(' ',''),folder,job.split(':')[0].replace(' ','')))
    
  if options.unblind:
    # plot pval
    add_string = ""
    for job in splitChannels:
      add_string += ' -f %s/ExpProfileLikelihood/%s/ExpProfileLikelihood.root -c %d -s 2 -w 2 -n -1'%(folder,job.split(':')[0].replace(' ',''),int(job.split(':')[2]))
      if options.unblind:
        add_string += ' -f %s/ProfileLikelihood/%s/ProfileLikelihood.root -c %d -s 1 -w 2 -n %s'%(folder,job.split(':')[0].replace(' ',''),int(job.split(':')[2]),job.split(':')[1])
    os.system('./makeCombinePlots.py -f %s/ExpProfileLikelihood/ExpProfileLikelihood.root -c 1 -s 2 -w 2 -n "Expected 1xSM" -f %s/ProfileLikelihood/ProfileLikelihood.root -c 1 -s 1 -w 2 -n "Observed" --text="%s" --pval %s'%(folder,folder,lumistring,add_string))
    os.system('cp pval.pdf plots/%s_obspval.pdf'%folder)
    os.system('cp pval.png plots/%s_obspval.png'%folder)
    # plot maxlh
    os.system('./makeCombinePlots.py -f %s/MaxLikelihoodFit/MaxLikelihoodFit.root --text="%s" --maxlh'%(folder,lumistring))
    os.system('cp maxlh.pdf plots/%s_maxlh.pdf'%folder)
    os.system('cp maxlh.png plots/%s_maxlh.png'%folder)
    # plot limit
    os.system('./makeCombinePlots.py -f %s/Asymptotic/Asymptotic.root --text="%s" --limit'%(folder,lumistring))
    os.system('cp limit.pdf plots/%s_obslimit.pdf'%folder)
    os.system('cp limit.png plots/%s_obslimit.png'%folder)

  else:
    for job in splitChannels:
      add_string += ' -f %s/ExpProfileLikelihood/%s/ExpProfileLikelihood.root -c %d -s 2 -w 2 -n -1'%(folder,job.split(':')[0].replace(' ',''),int(job.split(':')[2]))
    os.system('./makeCombinePlots.py -f %s/ExpProfileLikelihood/ExpProfileLikelihood.root -c 1 -s 2 -w 2 -n "Expected 1xSM" --text="%s" --pval %s'%(folder,folder,lumistring,add_string))
    os.system('cp pval.pdf plots/%s_exppval.pdf'%folder)
    os.system('cp pval.png plots/%s_exppval.png'%folder)
    os.system('./makeCombinePlots.py -f %s/Asymptotic/Asymptotic.root --text="%s" --limit -e'%(folder,lumistring))
    os.system('cp limit.pdf plots/%s_explimit.pdf'%folder)
    os.system('cp limit.png plots/%s_explimit.png'%folder)

  if options.run2011 and options.run2012 and options.runBoth and folder=='7and8TeV':
    if options.unblind:
      os.system('./makeCombinePlots.py -f 7TeV/ExpProfileLikelihood/ExpProfileLikelihood.root -c 2 -s 2 -w 2 -n -1 -f 8TeV/ExpProfileLikelihood/ExpProfileLikelihood.root -c 4 -s 2 -w 2 -n -1 -f 7and8TeV/ExpProfileLikelihood/ExpProfileLikelihood.root -c 1 -s 2 -w 2 -n "Expected 1xSM" -f 7TeV/ProfileLikelihood/ProfileLikelihood.root -c 2 -s 1 -w 2 -n "7TeV Observed" -f 8TeV/ProfileLikelihood/ProfileLikelihood.root -c 4 -s 1 -w 2 -n "8TeV Observed" -f 7and8TeV/ProfileLikelihood/ProfileLikelihood.root -c 1 -s 1 -w 2 -n "Combined Observed" --text="%s" --pval'%(lumistring))
      os.system('cp pval.pdf plots/obspval.pdf')
      os.system('cp pval.png plots/obspval.png')
    

