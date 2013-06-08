#!/bin/env python

import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-Q","--qqbarCard",dest="qqbarCard",action="store_true",default=False)
parser.add_option("-s","--justSM",dest="justSM",action="store_true",default=False)
parser.add_option("-g","--justGravGG",dest="justGravGG",action="store_true",default=False)
parser.add_option("-q","--justGravQQ",dest="justGravQQ",action="store_true",default=False)
parser.add_option("-c","--nCosThetaCats",dest="cTcats",type="int",default=5)
parser.add_option("-C","--nKinCats",dest="kCats",type="int",default=4)
parser.add_option("-n","--nameOfCard",dest="cardname")
parser.add_option("-S","--sigfile",default="CMS-HGG_interpolated.root")
parser.add_option("-B","--bkgfile",default="CMS-HGG.root")
(options,args)=parser.parse_args()

ncats=options.cTcats*options.kCats
card = open(options.cardname,'w')

card.write('CMS-HGG spin card for for use with combine. RooDataHist+Parametrised Background\n')
if options.qqbarCard: card.write('For qqbar scan\n')
else: card.write('For separation and model dependence plots\n')
card.write('\n------------------------------\nimax *\njmax *\nkmax *\n------------------------------\n')
card.write('shapes data_obs    *        %s cms_hgg_workspace:roohist_data_mass_$CHANNEL\n'%options.bkgfile)

if not options.justGravGG and not options.justGravQQ:
  card.write('shapes ggH         *      %s cms_hgg_workspace:roohist_sig_ggh_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_ggh_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)
  card.write('shapes qqH         *      %s cms_hgg_workspace:roohist_sig_vbf_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_vbf_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)
  card.write('shapes  VH         *      %s cms_hgg_workspace:roohist_sig_wzh_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_wzh_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)
  card.write('shapes ttH         *      %s cms_hgg_workspace:roohist_sig_tth_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_tth_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)

if not options.justSM and not options.justGravQQ:
  card.write('shapes ggH_ALT     *      %s cms_hgg_workspace:roohist_sig_gg_grav_mass_m$MASS_$CHANNEL cms_hgg_workspace:roohist_sig_gg_grav_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)

if options.qqbarCard or options.justGravQQ: 
  card.write('shapes qqbarH_ALT  *      %s cms_hgg_workspace:roohist_sig_qq_grav_mass_m$MASS_$CHANNEL cms_hgg_workspace:roohist_sig_qq_grav_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)

card.write('shapes bkg         *        %s cms_hgg_workspace:pdf_data_pol_model_8TeV_$CHANNEL\n'%options.bkgfile) 
card.write('------------------------------\n')

procs=[]
iprocs=[]
rates=[]

# normal signals
if not options.justGravGG and not options.justGravQQ:
  procs += ['ggH','qqH','VH','ttH']
  iprocs += ['0','-1','-2','-3']
  rates += ['-1','-1','-1','-1']
# ggh grav signal
if not options.justSM and not options.justGravQQ:
  procs.append('ggH_ALT')
  iprocs.append('-4')
  rates.append('-1')
# qqbar grav signal
if options.qqbarCard or options.justGravQQ: 
  procs.append('qqbarH_ALT')
  iprocs.append('-5')
  rates.append('-1')
# background
procs.append('bkg')
iprocs.append('1')
rates.append('1')

lumi_err = '1.044'
scale_ggH = '0.930/1.076'
scale_qqH = '0.972/1.026'
scale_VH  = '0.958/1.042'
scale_ttH = '0.920/1.080'
pdf_ggH   = '0.918/1.076'
pdf_qqH   = '0.992/1.003'
pdf_VH    = '0.982/1.021'
pdf_ttH   = '0.906/1.041'

binned_systs = ['E_res','E_scale','idEff','triggerEff','vtxEff','r9Eff','ptSpin']

card.write('\nbin     ')
for cat in range(ncats): card.write('cat%d '%cat)

card.write('\nobservation    ')
for cat in range(ncats): card.write('-1   ')
card.write('\n------------------------------')

card.write('\nbin     ')
for cat in range (ncats):
  for i,proc in enumerate(procs):
    card.write('cat%d '%cat)
    
card.write('\nprocess  ')
for cat in range (ncats):
  for i,proc in enumerate(procs):
    card.write('%s '%proc)
    
card.write('\nprocess  ')
for cat in range (ncats):
  for i,proc in enumerate(procs):
    card.write('%s '%iprocs[i])
    
card.write('\nrate  ')
for cat in range (ncats):
  for i,proc in enumerate(procs):
    card.write('%s '%rates[i])
    
card.write('\n------------------------------')
card.write('\nlumi    lnN  ')
for cat in range (ncats):
  for i,proc in enumerate(procs):
    if proc=='bkg': card.write(' -')
    else: card.write(' %s'%lumi_err)
   
card.write('\nQCDscale_ggH    lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'ggH' in proc and 'ggH_ALT' not in proc: card.write(' %s'%scale_ggH)
    else: card.write(' -')
   
card.write('\nQCDscale_qqH    lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'qqH' in proc: card.write(' %s'%scale_qqH)
    else: card.write(' -')
   
card.write('\nQCDscale_VH     lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'VH' in proc: card.write(' %s'%scale_VH)
    else: card.write(' -')
   
card.write('\nQCDscale_ttH    lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'ttH' in proc: card.write(' %s'%scale_ttH)
    else: card.write(' -')
   
card.write('\nPDF_ggH    lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'ggH' in proc and 'ggH_ALT' not in proc: card.write(' %s'%pdf_ggH)
    else: card.write(' -')
   
card.write('\nPDF_qqH    lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'qqH' in proc: card.write(' %s'%pdf_qqH)
    else: card.write(' -')

card.write('\nPDF_VH     lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'VH' in proc: card.write(' %s'%pdf_VH)
    else: card.write(' -')
   
card.write('\nPDF_ttH    lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'ttH' in proc: card.write(' %s'%pdf_ttH)
    else: card.write(' -')
   
for syst in binned_systs:
  card.write('\n%s    shape  '%syst)
  for cat in range (ncats):
    for i,proc in enumerate(procs):
      if proc=='bkg': card.write(' 0')
      else: card.write(' 0.333333')
     
card.close()

catMap = {}
for c in range(options.kCats):
  for s in range(options.cTcats):
    catMap[(c*options.cTcats)+s] = (c,s)

if options.justSM or options.justGravGG or options.justGravQQ or options.qqbarCard:
  os.system('combineCards.py %s > %s'%(card.name,(options.cardname).replace('.txt','_temp.txt')))
  old_card = open('%s'%(options.cardname).replace('.txt','_temp.txt'))
  new_card = open('%s'%options.cardname,'w')
  for line in old_card.readlines():
    for cat in range(ncats-1,-1,-1):
        line = line.replace('ch1_cat%d'%cat,'ch1_kinCat%d_spinCat%d'%(catMap[cat][0],catMap[cat][1]))
    new_card.write(line)
  new_card.close()
  os.system('rm -f %s'%old_card.name)



  
