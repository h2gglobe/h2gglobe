import ROOT as r
import sys

class BdtToyMaker:


	def __init__(self,fileName,pdfName):

		#r.RooMsgService.instance().setGlobalKillBelow(r.WARNING)
		r.gROOT.SetStyle("Plain")
		r.gROOT.SetBatch(True)

		self.hasFit_ = False

		self.fileName_= fileName
		self.thefile_ = r.TFile(fileName)
		self.ws_      = self.thefile_.Get("cms_hgg_workspace")
		self.var_     = self.ws_.var("CMS_hgg_mass")
		self.data_    = self.ws_.data("data_mass_cat0")

		data_th1_     = self.thefile_.Get("th1f_data_mass_cat0")

		#self.nEvents_ = data_th1_.Integral()# cannot use this since BDT has to me more inclusive than mass for keys pdf to workv

		self.pdfName_ = pdfName
		self.pdf_     = self.ws_.pdf("pdf_%s_cat0"%self.pdfName_)

		self.rand_    = r.TRandom3(0)

	def fitData(self):

		# Fit the data full range
		self.pdf_.getVal()
		#self.var_.setRange("window124L",100,124.0*0.98)
		#self.var_.setRange("window124H",124.0*1.02,180)
#		self.pdf_.fitTo(self.data_,r.RooFit.Range("window124L,window124H"))
		self.pdf_.fitTo(self.data_)
	 	self.hasFit_=True	

	def genData(self):

		print "BdtToyMaker -- Generating Toy Dataset"	
		if not self.hasFit_: self.fitData()
		self.toyNevents_ = self.rand_.Poisson(self.nEvents_)
		self.gendata_    = self.pdf_.generate(r.RooArgSet(self.var_),self.toyNevents_)
		self.genbdtdata_ = self.pdfbdt_.generate(r.RooArgSet(self.bdtvar_),self.toyNevents_)
		
	def getN(self,mH,size):
		
	#	if self.gendata_:
		  self.var_.setRange("rngeForMe",mH*(1.-size),mH*(1.+size))
		  self.var_.setRange("FULL",100,180)
		  int_ = self.pdf_.createIntegral(r.RooArgSet(self.var_),"rngeForMe")	
		  intfull_ = self.pdf_.createIntegral(r.RooArgSet(self.var_),"FULL")	
		  return self.data_.sumEntries()*int_.getVal()/intfull_.getVal()
	#	else: sys.exit("No Data generated yet! first run powerLaw.genData()")

	def createHistPdf(self,HIST):
		self.bdthist_ = HIST.Clone()
		self.bdtvar_=r.RooRealVar("diphotonMVA","diphotonMVA",HIST.GetBinLowEdge(1),HIST.GetBinLowEdge(HIST.GetNbinsX()+1))
		self.bdtdatahist_=r.RooDataHist("roohist_bdt","roohist_bdt",r.RooArgList(self.bdtvar_),self.bdthist_)
		self.pdfbdt_=r.RooHistPdf("pdf_bdt","pdf_bdt",r.RooArgSet(self.bdtvar_),self.bdtdatahist_)

	def createKeysPdf(self,tree):
		#self.nEvents_ = tree.GetEntries()
		self.bdtvar_=r.RooRealVar("diphotonMVA","diphotonMVA",-1.0,1.0)	# go a littel below actual cut to let keysPdf fit properly 
		self.bdtdatahist_=r.RooDataSet("roohist_bdt","roohist_bdt",tree,r.RooArgSet(self.bdtvar_),"diphotonMVA>-1.0")	
		self.nEvents_=self.bdtdatahist_.sumEntries()
		self.pdfbdt_=r.RooKeysPdf("pdf_bdt","pdf_bdt",self.bdtvar_,self.bdtdatahist_)

	def loadKeysPdf(self,fi):
		#self.nEvents_ = tree.GetEntries()
		self.ws_=fi.Get("bdtworkspace")
		self.bdtvar_=self.ws_.var("diphotonMVA")
		self.bdtdatahist_=self.ws_.data("roohist_bdt")	
		self.pdfbdt_=self.ws_.pdf("pdf_bdt")
		#self.bdtvar_=r.RooRealVar("diphotonMVA","diphotonMVA",-0.1,1.0)	# go a littel below actual cut to let keysPdf fit properly 
#i#		self.bdtdatahist_=r.RooDataSet("roohist_bdt","roohist_bdt",tree,r.RooArgSet(self.bdtvar_),"diphotonMVA>-0.1")	
		self.nEvents_=self.bdtdatahist_.sumEntries()
#		self.pdfbdt_=r.RooKeysPdf("pdf_bdt","pdf_bdt",self.bdtvar_,self.bdtdatahist_)

	def saveBdtWorkspace(self,fi):
		saveFile = r.TFile("bdtws.root","RECREATE")		
		self.ws_=r.RooWorkspace("bdtworkspace")
		getattr(self.ws_,'import')(self.bdtdatahist_)
		getattr(self.ws_,'import')(self.pdfbdt_)
		saveFile.cd()
		self.ws_.Write()
		saveFile.Close()
		print "Saved BDT Pdf in workspace to ", fi

	def returnWindowToyData(self,mH,size):# iterate through the 2-variable dataset and return the values

		returnList = []

		mHL = mH*(1.-size)
		mHH = mH*(1.+size)

  		for  i in range(self.toyNevents_):
			val_m = (self.gendata_.get(i)).getRealValue("CMS_hgg_mass");
			val_b = (self.genbdtdata_.get(i)).getRealValue("diphotonMVA");
			if val_m > mHL and val_m < mHH and val_b >=0.05: returnList.append((val_b,((val_m-mH)/mH)))

		return returnList
			
		

	def plotRealData(self,nbins):


		can = r.TCanvas("c","c",1100,600)
		can.Divide(2,1)

		can.cd(1);frame1 = self.var_.frame();self.data_.plotOn(frame1,r.RooFit.Binning(int(nbins)))
		self.pdf_.plotOn(frame1); frame1.Draw()

		can.cd(2);frame2 = self.bdtvar_.frame();self.bdtdatahist_.plotOn(frame2)
		self.pdfbdt_.plotOn(frame2); frame2.Draw(); frame2.GetXaxis().SetRangeUser(0.05,1.0)
		can.SaveAs("datafitFullSpec_%s.pdf"%(self.fileName_))	

	def plotGenData(self,nbins):
	     if self.gendata_:
		can = r.TCanvas("c","c",1100,600)
		can.Divide(2,1)

		can.cd(1);frame1 = self.var_.frame();self.gendata_.plotOn(frame1,r.RooFit.Binning(int(nbins)))
		self.pdf_.plotOn(frame1); frame1.Draw()

		can.cd(2);frame2 = self.bdtvar_.frame();self.genbdtdata_.plotOn(frame2)
		self.pdfbdt_.plotOn(frame2); frame2.Draw()

		can.SaveAs("gendataFromFit_%s.pdf"%(self.fileName_))	
	     else: sys.exit("No Data generated yet! first run powerLaw.genData()")


	
