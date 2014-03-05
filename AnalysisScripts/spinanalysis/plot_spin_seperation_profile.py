#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--dir")
parser.add_option("-S","--sqrtS",default="comb")
parser.add_option("-u","--unblind",default=False,action="store_true")
parser.add_option("-r","--rebin",type="int",default=1000)
parser.add_option("-x","--xrange",default="-10,10.")
parser.add_option("-b","--batch",default=False,action="store_true")
(options,args) = parser.parse_args()

options.xrange = [float(options.xrange.split(',')[0]),float(options.xrange.split(',')[1])]

import ROOT as r
import sys, os

if options.batch: r.gROOT.SetBatch()
r.gROOT.ProcessLine(".L extractSignificanceStats.C+g")
from ROOT import extractSignificanceStats

cls={}
sigma={}
range={}
if options.sqrtS=='7' or options.sqrtS=='7TeV':
	range[0.00] = [-10.,10.]
	range[0.25] = [-8.,8.]
	range[0.50] = [-6.,6.]
	range[0.75] = [-6.,6.]
	range[1.00] = [-8.,8.]
elif options.sqrtS=='8' or options.sqrtS=='8TeV':
	range[0.00] = [-12.,12.]
	range[0.25] = [-8.,8.]
	range[0.50] = [-6.,6.]
	range[0.75] = [-6.,6.]
	range[1.00] = [-8.,8.]
elif options.sqrtS=='comb':
	range[0.00] = [-15.,15.]
	range[0.25] = [-12.,12.]
	range[0.50] = [-10.,10.]
	range[0.75] = [-10.,10.]
	range[1.00] = [-12.,12.]
else:
	sys.exit('Invalid sqrtS'+options.sqrtS)

for fqq in [0.0,0.25,0.5,0.75,1.0]:
	fname = options.dir+'/testStat_fqq%4.2f.root'%fqq
	if not os.path.exists(fname): fname = options.dir+'/qmu_qqbar%4.2f.root'%fqq
	legname = '2^{+}_{m}(%3.0f%% q#bar{q})'%(fqq*100.)
	altname = '2pm%4.2f'%fqq
	expCLs = r.Double(0)
	expSigSM = r.Double(0)
	expSigALT = r.Double(0)
	obsCLs = r.Double(0)
	obsSigSM = r.Double(0)
	obsSigALT = r.Double(0)
	outname = options.dir+'/'+altname+'_blind'
	if options.unblind: outname = outname.replace('blind','unblind')
	extractSignificanceStats(expCLs,obsCLs,expSigSM,expSigALT,obsSigSM,obsSigALT,options.unblind,legname,outname,fname,options.rebin,range[fqq][0],range[fqq][1],options.sqrtS)
	cls[fqq] = [expCLs,obsCLs]
	sigma[fqq] = [expSigSM,expSigALT,obsSigSM,obsSigALT]

print 'fqq - 1-expCLs - expP(SM<ALT) - expP(ALT>SM) - 1-obsCLs - obsP(SM<Obs) - obsP(ALT>Obs)'
for fqq in [0.0,0.25,0.5,0.75,1.00]:
	print '%4.2f - %5.3f - %4.2f - %4.2f - %5.3f - %4.2f - %4.2f'%(fqq,1-cls[fqq][0],r.RooStats.PValueToSignificance(sigma[fqq][0]),r.RooStats.PValueToSignificance(sigma[fqq][1]),1-cls[fqq][1],r.RooStats.PValueToSignificance(sigma[fqq][2]),r.RooStats.PValueToSignificance(sigma[fqq][3]))
