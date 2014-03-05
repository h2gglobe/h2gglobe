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
parser.add_option("--isMultiPdf",default=False,action="store_true")
parser.add_option("--sqrtS",type="int",default=8)
(options,args)=parser.parse_args()

ncats=options.cTcats*options.kCats
card = open(options.cardname,'w')

bkgWS = 'cms_hgg_workspace'
dataWS = 'cms_hgg_workspace'
if options.isMultiPdf: 
	bkgWS = 'multipdf'
	dataWS = 'multipdf'

card.write('CMS-HGG spin card for for use with combine. RooDataHist+Parametrised Background\n')
if options.qqbarCard: card.write('For qqbar scan\n')
else: card.write('For separation and model dependence plots\n')
card.write('\n------------------------------\nimax *\njmax *\nkmax *\n------------------------------\n')
card.write('shapes data_obs    *        %s %s:roohist_data_mass_$CHANNEL\n'%(options.bkgfile,bkgWS))

if not options.justGravGG and not options.justGravQQ:
  card.write('shapes ggH         *      %s cms_hgg_workspace:roohist_sig_ggh_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_ggh_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)
  card.write('shapes qqH         *      %s cms_hgg_workspace:roohist_sig_vbf_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_vbf_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)
  card.write('shapes  WH         *      %s cms_hgg_workspace:roohist_sig_wh_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_wh_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)
  card.write('shapes  ZH         *      %s cms_hgg_workspace:roohist_sig_zh_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_zh_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)
  card.write('shapes ttH         *      %s cms_hgg_workspace:roohist_sig_tth_mass_m$MASS_$CHANNEL      cms_hgg_workspace:roohist_sig_tth_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)

if not options.justSM and not options.justGravQQ:
  card.write('shapes ggH_ALT     *      %s cms_hgg_workspace:roohist_sig_gg_grav_mass_m$MASS_$CHANNEL cms_hgg_workspace:roohist_sig_gg_grav_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)

if options.qqbarCard or options.justGravQQ: 
  card.write('shapes qqbarH_ALT  *      %s cms_hgg_workspace:roohist_sig_qq_grav_mass_m$MASS_$CHANNEL cms_hgg_workspace:roohist_sig_qq_grav_mass_m$MASS_$CHANNEL_$SYSTEMATIC01_sigma\n'%options.sigfile)

if options.isMultiPdf:
	card.write('shapes bkg         *        %s multipdf:CMS_hgg_$CHANNEL_%dTeV_bkgshape\n'%(options.bkgfile,options.sqrtS))
else:
	card.write('shapes bkg         *        %s cms_hgg_workspace:pdf_data_pol_model_%dTeV_$CHANNEL\n'%(options.bkgfile,options.sqrtS)) 
card.write('------------------------------\n')

procs=[]
iprocs=[]
rates=[]

# normal signals
if not options.justGravGG and not options.justGravQQ:
  #procs += ['ggH','qqH','VH','ttH']
  procs += ['ggH','qqH','WH','ZH','ttH']
  iprocs += ['0','-1','-2','-3','-4']
  rates += ['-1','-1','-1','-1','-1']
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

lumi_err = '1.025' if options.sqrtS==8 else '1.022'
br_err = [1.050,0.951]

theorySyst = {}
theorySyst['QCDscale_ggH'] = {}
theorySyst['QCDscale_qqH'] = {}
theorySyst['QCDscale_VH'] = {}
theorySyst['QCDscale_ttH'] = {}
theorySyst['pdf_gg'] = {}
theorySyst['pdf_qqbar'] = {}
# 7 TeV
if options.sqrtS==7:
	# scale
	theorySyst['QCDscale_ggH']['ggH'] = [0.071,-0.078] 
	theorySyst['QCDscale_ggH']['ggH_ALT'] = [0.071,-0.078] 
	theorySyst['QCDscale_qqH']['qqH'] = [0.003,-0.003]
	theorySyst['QCDscale_qqH']['qqbarH_ALT'] = [0.003,-0.003]
	theorySyst['QCDscale_VH']['WH'] = [0.009,-0.009]
	theorySyst['QCDscale_VH']['ZH'] = [0.029,-0.029]
	theorySyst['QCDscale_ttH']['ttH'] = [0.032,-0.093]
	# pdf
	theorySyst['pdf_gg']['ggH'] = [0.076,-0.071] 
	theorySyst['pdf_gg']['ggH_ALT'] = [0.076,-0.071] 
	theorySyst['pdf_qqbar']['qqH'] = [0.025,-0.021]
	theorySyst['pdf_qqbar']['qqbarH_ALT'] = [0.025,-0.021]
	theorySyst['pdf_qqbar']['WH'] = [0.026,-0.026]
	theorySyst['pdf_qqbar']['ZH'] = [0.027,-0.027]
	theorySyst['pdf_gg']['ttH'] = [0.81,-0.81]
# 8 TeV
elif options.sqrtS==8:
	# scale
	theorySyst['QCDscale_ggH']['ggH'] = [0.072,-0.078]
	theorySyst['QCDscale_ggH']['ggH_ALT'] = [0.072,-0.078]
	theorySyst['QCDscale_qqH']['qqH'] = [0.002,-0.002]
	theorySyst['QCDscale_qqH']['qqbarH_ALT'] = [0.002,-0.002]
	theorySyst['QCDscale_VH']['WH'] = [0.010,-0.010]
	theorySyst['QCDscale_VH']['ZH'] = [0.031,-0.031]
	theorySyst['QCDscale_ttH']['ttH'] = [0.038,-0.093]
	# pdf
	theorySyst['pdf_gg']['ggH'] = [0.075,-0.069] 
	theorySyst['pdf_gg']['ggH_ALT'] = [0.075,-0.069] 
	theorySyst['pdf_qqbar']['qqH'] = [0.026,-0.028]
	theorySyst['pdf_qqbar']['qqbarH_ALT'] = [0.026,-0.028]
	theorySyst['pdf_qqbar']['WH'] = [0.023,-0.023]
	theorySyst['pdf_qqbar']['ZH'] = [0.025,-0.025]
	theorySyst['pdf_gg']['ttH'] = [0.081,-0.081]
else:
	sys.exit(str(options.sqrtS)+'TeV is not a valid sqrtS')

"""
scale_ggH = '0.930/1.076'
scale_qqH = '0.972/1.026'
#scale_VH  = '0.958/1.042'
scale_WH  = '0.958/1.042'
scale_ZH  = '0.958/1.042'
scale_ttH = '0.920/1.080'
pdf_ggH   = '0.918/1.076'
pdf_qqH   = '0.992/1.003'
#pdf_VH    = '0.982/1.021'
pdf_WH    = '0.982/1.021'
pdf_ZH    = '0.982/1.021'
pdf_ttH   = '0.906/1.041'
"""

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
card.write('\n')
   
card.write('\nCMS_hgg_BR    lnN  ')
for cat in range (ncats):
  for i,proc in enumerate(procs):
    if proc=='bkg': card.write(' -')
    else: card.write(' %s/%s'%(br_err[1],br_err[0]))
card.write('\n')
   
for systName, systDetails in theorySyst.items():
	card.write('%-35s   lnN   '%systName)
	for c in range(ncats):
		for p in procs:
			if p in systDetails.keys():
				card.write('%5.3f/%5.3f '%(1.+systDetails[p][1],1.+systDetails[p][0]))
			else:
				card.write('- ')
	card.write('\n')
card.write('\n')

"""
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
   
card.write('\nQCDscale_WH     lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'WH' in proc: card.write(' %s'%scale_WH)
    else: card.write(' -')
   
card.write('\nQCDscale_ZH     lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'ZH' in proc: card.write(' %s'%scale_ZH)
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

card.write('\nPDF_WH     lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'WH' in proc: card.write(' %s'%pdf_WH)
    else: card.write(' -')
   
card.write('\nPDF_ZH     lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'ZH' in proc: card.write(' %s'%pdf_ZH)
    else: card.write(' -')
   
card.write('\nPDF_ttH    lnN  ')
for cat in range(ncats):
  for i,proc in enumerate(procs):
    if 'ttH' in proc: card.write(' %s'%pdf_ttH)
    else: card.write(' -')
"""
   
for syst in binned_systs:
  card.write('\n%s    shape  '%syst)
  for cat in range (ncats):
    for i,proc in enumerate(procs):
      if proc=='bkg': card.write(' 0')
      else: card.write(' 0.333333')

if options.isMultiPdf:
	card.write('\n')
	for cat in range(ncats):
		card.write('pdfindex_%d_%dTeV discrete\n'%(cat,options.sqrtS))

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
        line = line.replace('ch1_cat%d'%cat,'kinCat%d_spinCat%d_%dTeV'%(catMap[cat][0],catMap[cat][1],options.sqrtS))
    new_card.write(line)
  new_card.close()
  os.system('rm -f %s'%old_card.name)


  
