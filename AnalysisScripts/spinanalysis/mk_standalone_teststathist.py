#!/usr/bin/env python

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-d","--dir")
parser.add_option("--expJobs",type="int",help="Number of jobs expected")
parser.add_option("--isGrid",default=False,action="store_true")
parser.add_option("--doUnpack",default=False,action="store_true")
(options,args)=parser.parse_args()

import os
import ROOT as r
import array
import sys
import fnmatch

def getTestStat(smFileName,altFileName):
	smFile = r.TFile(smFileName)
	altFile = r.TFile(altFileName)
	smTree = smFile.Get('limit')
	altTree = altFile.Get('limit')
	assert(smTree.GetEntries()==altTree.GetEntries())
	res = []
	for e in range(smTree.GetEntries()):
		smTree.GetEntry(e)
		altTree.GetEntry(e)
		res.append(smTree.absNLL - altTree.absNLL)
	smFile.Close()
	altFile.Close()
	return res

def getNLL(file):
	ttree = file.Get('limit')
	ttree.GetEntry(0)
	return ttree.absNLL

thData=[]
thSM=[]
thALT=[]

for i, fqq in enumerate([0.00,0.25,0.50,0.75,1.00]):
	thData.append(r.TH1F('thData_fqqI%d'%i,'',1000,-10,10))
	thSM.append(r.TH1F('thSM_fqqI%d'%i,'',100,-10,10))
	thALT.append(r.TH1F('thALT_fqqI%d'%i,'',100,-10,10))

if options.isGrid:

	if options.doUnpack:
		for fqq,x in [ [0.00,0.] , [0.00,1.], [0.25,1.], [0.50,1.], [0.75,1.], [1.00,1] ]:
			dir = options.dir+'/'+'crab_fqq%4.2f_x%1.0f'%(fqq,x)
			nToys=0
			for root,dirs,files in os.walk(dir+'/res/'):
				for file in fnmatch.filter(files,'outputToy*.tgz'):
					toy = int((file.split('fqq%4.2f_x%1.0f_'%(fqq,x))[1]).split('_')[0])
					nToys+=1
					os.system('mkdir -p %s/crab_fqq%4.2f_x%1.0f/res/outputToy%d'%(options.dir,fqq,x,toy))
					os.system('tar -xvzf %s/crab_fqq%4.2f_x%1.0f/res/%s'%(options.dir,fqq,x,file))
					os.system('mv outputToy_fqq%4.2f_x%1.0f/* %s/crab_fqq%4.2f_x%1.0f/res/outputToy%d'%(fqq,x,options.dir,fqq,x,toy))
			print 'Found', nToys, 'in', dir

	for i, fqq in enumerate([0.00,0.25,0.50,0.75,1.00]):
		outf = r.TFile('%s/testStat_fqq%4.2f.root'%(options.dir,fqq),'RECREATE')
		tree = r.TTree('q','q')
		q = -999.
		type = -999
		qhold = array.array('f',[q])
		thold = array.array('i',[type])
		tree.Branch('q',qhold,'q/F')
		tree.Branch('type',thold,'type/I')

		thData[i].SetLineColor(r.kBlack)
		thData[i].SetFillColor(r.kBlack)
		thSM[i].SetLineWidth(2)
		thSM[i].SetLineColor(r.kBlue)
		thALT[i].SetLineWidth(2)
		thALT[i].SetLineColor(r.kRed)

		# data
		smDataFitF = r.TFile('%s/higgsCombineDataFit_fqq0.00_x0.MultiDimFit.mH120.root'%(options.dir))
		gravDataFitF = r.TFile('%s/higgsCombineDataFit_fqq%4.2f_x1.MultiDimFit.mH120.root'%(options.dir,fqq))
		smDataFitABSNLL = getNLL(smDataFitF)
		gravDataFitABSNLL = getNLL(gravDataFitF)

		# data
		thold[0]=0
		qhold[0] = smDataFitABSNLL - gravDataFitABSNLL
		tree.Fill()
		thData[i].Fill(-2.*qhold[0])

		print 'fqq - ', fqq
		#print '\tData: ', thold[0], qhold[0]

		smToys=0
		gravToys=0
		for jobn in range(options.expJobs):
			
			# SM hypo
			'higgsCombineGen_fqq0.25_x1_Fit_fqq0.00_x0.MultiDimFit.mH120.123456.root'
			smToySMFitName = '%s/crab_fqq0.00_x0/res/outputToy%d/higgsCombineGen_fqq0.00_x0_Fit_fqq0.00_x0.MultiDimFit.mH120.123456.root'%(options.dir,jobn)
			smToyGravFitName = '%s/crab_fqq0.00_x0/res/outputToy%d/higgsCombineGen_fqq0.00_x0_Fit_fqq%4.2f_x1.MultiDimFit.mH120.123456.root'%(options.dir,jobn,fqq)
				
			if os.path.exists(smToySMFitName) and os.path.exists(smToyGravFitName):
				qvals = getTestStat(smToySMFitName,smToyGravFitName)
				for qv in qvals:
					thold[0] = -1
					qhold[0] = qv
					tree.Fill()
					thSM[i].Fill(-2.*qhold[0])
					smToys+=1
					#print '\tSMtoy: ', thold[0], qhold[0]

			# GRAV hypo
			gravToySMFitName = '%s/crab_fqq%4.2f_x1/res/outputToy%d/higgsCombineGen_fqq%4.2f_x1_Fit_fqq0.00_x0.MultiDimFit.mH120.123456.root'%(options.dir,fqq,jobn,fqq)
			gravToyGravFitName = '%s/crab_fqq%4.2f_x1/res/outputToy%d/higgsCombineGen_fqq%4.2f_x1_Fit_fqq%4.2f_x1.MultiDimFit.mH120.123456.root'%(options.dir,fqq,jobn,fqq,fqq)
			if os.path.exists(gravToySMFitName) and os.path.exists(gravToyGravFitName):
				qvals = getTestStat(gravToySMFitName,gravToyGravFitName)
				for qv in qvals:
					thold[0] = 1
					qhold[0] = qv
					tree.Fill()
					thALT[i].Fill(-2.*qhold[0])
					gravToys +=1
					#print '\tALTtoy:', thold[0], qhold[0]
		outf.cd()
		tree.Write()
		print 'Found %d smToys and %d gravToys'%(smToys,gravToys)
		print 'Written to file', outf.GetName()
		outf.Close()
	
	"""
	print thData, thSM, thALT

	canv = r.TCanvas('canv','canv',1200,900)
	canv.Divide(3,2)
	thDummy = r.TH1F("d","",1,-10,10)
	for i, fqq in enumerate([0.00,0.25,0.50,0.75,1.00]):
		canv.cd(i+1)
		thDummy.SetTitle('fqq%4.2f'%fqq)
		thDummy.Draw("AXISG")
		thALT[i].Draw("HISTsame")
		thSM[i].Draw("HISTsame")
		thData[i].Draw("same")
		outf.cd()
		thALT[i].Write()
		thSM[i].Write()
		thData[i].Write()
		
	canv.Update()
	canv.Modified()
	canv.Print("%s/testStats.pdf"%options.dir)
	"""
	
else:

	job_folders=[]
	for root,dirs,files in os.walk(options.dir):
		if 'grid' in root: continue
		if 'job' in os.path.basename(root):
			print root
			print os.path.basename(root)
			print os.path.dirname(root)
			job_folders.append(int(root.split('job')[1]))

	print 'Found %d job folders'%len(job_folders)


	for i, fqq in enumerate([0.00,0.25,0.50,0.75,1.00]):
		outf = r.TFile('%s/testStat_fqq%4.2f.root'%(options.dir,fqq),'RECREATE')
		tree = r.TTree('q','q')
		q = -999.
		type = -999
		qhold = array.array('f',[q])
		thold = array.array('i',[type])
		tree.Branch('q',qhold,'q/F')
		tree.Branch('type',thold,'type/I')

		thData[i].SetLineColor(r.kBlack)
		thData[i].SetFillColor(r.kBlack)
		thSM[i].SetLineWidth(2)
		thSM[i].SetLineColor(r.kBlue)
		thALT[i].SetLineWidth(2)
		thALT[i].SetLineColor(r.kRed)

		# data
		smDataFitF = r.TFile('%s/data_fits/higgsCombineDataFit_fqq0.00_x0.MultiDimFit.mH120.root'%(options.dir))
		gravDataFitF = r.TFile('%s/data_fits/higgsCombineDataFit_fqq%4.2f_x1.MultiDimFit.mH120.root'%(options.dir,fqq))
		smDataFitABSNLL = getNLL(smDataFitF)
		gravDataFitABSNLL = getNLL(gravDataFitF)

		# data
		thold[0]=0
		qhold[0] = smDataFitABSNLL - gravDataFitABSNLL
		tree.Fill()
		thData[i].Fill(-2.*qhold[0])

		print 'fqq - ', fqq
		#print '\tData: ', thold[0], qhold[0]

		for jobn in job_folders:
			
			# SM hypo
			smToySMFitName = '%s/job%d/higgsCombineJob%d_fqq0.00_x0_Fit_fqq0.00_x0.MultiDimFit.mH120.123456.root'%(options.dir,jobn,jobn)
			smToyGravFitName = '%s/job%d/higgsCombineJob%d_fqq0.00_x0_Fit_fqq%4.2f_x1.MultiDimFit.mH120.123456.root'%(options.dir,jobn,jobn,fqq)
				
			if os.path.exists(smToySMFitName) and os.path.exists(smToyGravFitName):
				qvals = getTestStat(smToySMFitName,smToyGravFitName)
				for qv in qvals:
					thold[0] = -1
					qhold[0] = qv
					tree.Fill()
					thSM[i].Fill(-2.*qhold[0])
					#print '\tSMtoy: ', thold[0], qhold[0]

			# GRAV hypo
			gravToySMFitName = '%s/job%d/higgsCombineJob%d_fqq%4.2f_x1_Fit_fqq0.00_x0.MultiDimFit.mH120.123456.root'%(options.dir,jobn,jobn,fqq)
			gravToyGravFitName = '%s/job%d/higgsCombineJob%d_fqq%4.2f_x1_Fit_fqq%4.2f_x1.MultiDimFit.mH120.123456.root'%(options.dir,jobn,jobn,fqq,fqq)
			if os.path.exists(gravToySMFitName) and os.path.exists(gravToyGravFitName):
				qvals = getTestStat(gravToySMFitName,gravToyGravFitName)
				for qv in qvals:
					thold[0] = 1
					qhold[0] = qv
					tree.Fill()
					thALT[i].Fill(-2.*qhold[0])
					#print '\tALTtoy:', thold[0], qhold[0]
		outf.cd()
		tree.Write()
		outf.Close()

	"""
	print thData, thSM, thALT

	canv = r.TCanvas('canv','canv',1200,900)
	canv.Divide(3,2)
	thDummy = r.TH1F("d","",1,-10,10)
	for i, fqq in enumerate([0.00,0.25,0.50,0.75,1.00]):
		canv.cd(i+1)
		thDummy.SetTitle('fqq%4.2f'%fqq)
		thDummy.Draw("AXISG")
		thALT[i].Draw("HISTsame")
		thSM[i].Draw("HISTsame")
		thData[i].Draw("same")
		outf.cd()
		thALT[i].Write()
		thSM[i].Write()
		thData[i].Write()
		
	canv.Update()
	canv.Modified()
	canv.Print("%s/testStats.pdf"%options.dir)
	"""
