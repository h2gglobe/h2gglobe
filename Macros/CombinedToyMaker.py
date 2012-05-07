import ROOT as r
import os,sys

r.gROOT.SetStyle("Plain")
r.gROOT.SetBatch(True)

r.gROOT.ProcessLine(".L makeToyWS.C+g")
from ROOT import makeToyWS

class CombinedToyMaker:

  def __init__(self,mitFileName):
    
    r.gROOT.SetStyle("Plain")
    r.gROOT.SetBatch(True)

    self.mitFileName_ = mitFileName
    self.mitFile_     = r.TFile(mitFileName)
    self.mitWS_       = self.mitFile_.Get("cms_hgg_workspace")
    self.mitvar_      = self.mitWS_.var("CMS_hgg_mass")
    self.mitpdfs_     =[]
    self.mitsighists_ =[]
    self.mitsigpdfs_  =[]
    # get background pdf for each mit mass cat
    for cat in range(5): 
      self.mitpdfs_.append(self.mitWS_.pdf("pdf_data_pol_model_cat%i"%cat))
    # get signal histpdf for each mit mass cat
    for cat in range(5):
      temp = self.mitWS_.data("roohist_sig_ggh_mass_m124_cat%i"%cat)
      temp.add(self.mitWS_.data("roohist_sig_vbf_mass_m124_cat%i"%cat))
      temp.add(self.mitWS_.data("roohist_sig_wzh_mass_m124_cat%i"%cat))
      temp.add(self.mitWS_.data("roohist_sig_tth_mass_m124_cat%i"%cat))
      self.mitsigpdfs_.append(r.RooHistPdf("sigmasspdf_cat%i"%cat,"sigmasspdf_cat%i"%cat,r.RooArgSet(self.mitvar_),temp))

    self.hasFit_      = False
    self.rand_        = r.TRandom3(0)

  def createPdfs(self,infile,outfile,expSig=0):
    makeToyWS(infile,outfile)
    self.loadPdfs(outfile,expSig)

  def loadPdfs(self,filename,expSig=0):
    fi = r.TFile(filename)
    self.keyws_           = fi.Get("fits_workspace")
    self.bdtvar_          = self.keyws_.var("bdtoutput")
    self.bdtmgg_          = self.keyws_.var("CMS_hgg_mass")
    self.bdtvbf_          = self.keyws_.var("vbf")
    self.bdtdata_         = self.keyws_.data("data_bdt_novbf")
    self.bdtdatavbf_      = self.keyws_.data("data_bdt_vbf")
    self.bdtdatacut_      = self.keyws_.data("data_bdt_cut_all")
    self.icpdf_           = self.keyws_.pdf("data_pow_model")
    self.bdtpdf_          = self.keyws_.pdf("data_pdf")
    self.dataNoVBFEvents_ = self.bdtdata_.sumEntries("bdtoutput>=0.05")
    self.dataVBFEvents_   = self.bdtdatavbf_.sumEntries("bdtoutput>=0.05")
    if expSig>0:
      self.bdtweight_     = self.keyws_.var("weight")
      self.bdtsigpdf_     = self.keyws_.pdf("sig_pdf")
      self.bdtsigdata_    = self.keyws_.data("sig_bdt_novbf")
      self.bdtsigdatavbf_ = self.keyws_.data("sig_bdt_vbf")
      self.sigNoVBFEvents_= self.bdtsigdata_.sumEntries("bdtoutput>=0.05")
      self.sigVBFEvents_  = self.bdtsigdatavbf_.sumEntries("bdtoutput>=0.05")
    
    print "Non VBF events (data):  ", self.dataNoVBFEvents_
    print "VBF events (data):      ",self.dataVBFEvents_
    if expSig>0:
      print "Non VBF events (sig):   ", self.sigNoVBFEvents_
      print "VBF events (sig):       ",self.sigVBFEvents_
    self.hasFit_      = True

  def savePdfWorkspace(self,filename,expSig):
    fi = r.TFile(filename,"RECREATE")
    fi.cd()
    self.ws_          = r.RooWorkspace("toysworkspace")
    getattr(self.ws_,'import')(self.bdtdata_)
    getattr(self.ws_,'import')(self.bdtdatavbf_)
    getattr(self.ws_,'import')(self.bdtpdf_)
    if expSig>0:
      getattr(self.ws_,'import')(self.bdtsigdata_)
      getattr(self.ws_,'import')(self.bdtsigdatavbf_)
      getattr(self.ws_,'import')(self.bdtsigpdf_)
    self.ws_.Write()
    print "Saved PDF and DataSets to workspace in ", fi.GetName()
    fi.Close()

  def saveToyWorkspace(self,filename):
    fi = r.TFile(filename,"RECREATE")
    fi.cd()
    self.toyws_       = r.RooWorkspace("gentoyworkspace")
    getattr(self.toyws_,'import')(self.genbdtdata_)
    for i, genmdat in enumerate(self.genmassdata_):
      getattr(self.toyws_,'import')(genmdat)
      getattr(self.toyws_,'import')(self.mitpdfs_[i])

    getattr(self.toyws_,'import')(self.allgenmassdata_)
    self.toyws_.Write()
    fi.Close()

  def getEventsPerCat(self,gendata):
    events=[]
    events.append(gendata.sumEntries("bdtoutput>=0.89"))
    events.append(gendata.sumEntries("bdtoutput>=0.74 && bdtoutput<0.89"))
    events.append(gendata.sumEntries("bdtoutput>=0.545 && bdtoutput<0.74"))
    events.append(gendata.sumEntries("bdtoutput>=0.05 && bdtoutput<0.545"))
    for i, ev in enumerate(events):
      print 'Cat %d has %5d events'%(i,ev)
    return events

  def genData(self,outwsname,expSig):
    if not self.hasFit_: sys.exit('Diphoton output not yet fitted. Bailing out.')

    # first gen around bdt keys pdf
    self.toyNoVBFEvents_ = self.rand_.Poisson(int(self.dataNoVBFEvents_))
    self.bdtvar_.setRange(0.05,1.)
    self.genbdtdata_  = self.bdtpdf_.generate(r.RooArgSet(self.bdtvar_),int(self.toyNoVBFEvents_))
    self.genbdtdata_.SetName("gen_diphoton_output")
   
    eventsPerCat = self.getEventsPerCat(self.genbdtdata_)
    if len(eventsPerCat)!=len(self.mitpdfs_)-1: sys.exit('Different numbers of categories. Bailing out.')
    
    # given nevs generate around each non VBF cat
    can = r.TCanvas()
    self.genmassdata_ = []
    for i,ev in enumerate(eventsPerCat):
      self.genmassdata_.append(self.mitpdfs_[i].generate(r.RooArgSet(self.bdtmgg_),int(ev)))
    
    # combine data in each non VBF cat for IC analysis
    for i,dataset in enumerate(self.genmassdata_):
      dataset.SetName("gen_mass_data_cat%d"%i)
      if i==0: self.allgenmassdatanovbf_=r.RooDataSet("gen_combmassdata_novbf","gen_combmassdata_novbf",dataset,r.RooArgSet(self.bdtmgg_))
      else: self.allgenmassdatanovbf_.append(dataset)
    
    # throw toy in VBF cat
    self.toyVBFEvents_ = self.rand_.Poisson(int(self.dataVBFEvents_))
    self.genmassdata_.append(self.mitpdfs_[4].generate(r.RooArgSet(self.bdtmgg_),int(self.toyVBFEvents_)))

    # combine data in VBF cat for IC analysis
    self.allgenmassdata_ = r.RooDataSet("gen_combmassdata_all","gen_combmassdata_all",self.allgenmassdatanovbf_,r.RooArgSet(self.bdtmgg_))
    self.allgenmassdata_.append(self.genmassdata_[4])

    print 'Background toy thrown'
    for i, massdat in enumerate(self.genmassdata_):
      print 'Cat %d toy has  %5d events'%(i,massdat.sumEntries())
    print 'No vbf toy has %5d events'%self.allgenmassdatanovbf_.sumEntries()
    print 'Comb toy has   %5d events'%self.allgenmassdata_.numEntries()

    # then throw signal around bdt correlated with total mass and split into mass categories
    self.gensigmassdata_ = []
    if expSig>0:
      self.toyNoVBFEventsSig_ = self.rand_.Poisson(expSig*int(self.sigNoVBFEvents_))
      self.gensigdatanovbf_ = self.bdtsigpdf_.generate(r.RooArgSet(self.bdtvar_,self.bdtmgg_),self.toyNoVBFEventsSig_)
      self.gensigdatanovbf_.SetName("gen_sig_output_novbf")
      self.gensigmassdata_.append(r.RooDataSet("gen_sig_mass_cat0","gen_sig_mass_cat0",self.gensigdatanovbf_,r.RooArgSet(self.bdtvar_,self.bdtmgg_),"bdtoutput>=0.89"))
      self.gensigmassdata_.append(r.RooDataSet("gen_sig_mass_cat1","gen_sig_mass_cat0",self.gensigdatanovbf_,r.RooArgSet(self.bdtvar_,self.bdtmgg_),"bdtoutput>=0.74 && bdtoutput<0.89"))
      self.gensigmassdata_.append(r.RooDataSet("gen_sig_mass_cat2","gen_sig_mass_cat0",self.gensigdatanovbf_,r.RooArgSet(self.bdtvar_,self.bdtmgg_),"bdtoutput>=0.545 && bdtoutput<0.74"))
      self.gensigmassdata_.append(r.RooDataSet("gen_sig_mass_cat3","gen_sig_mass_cat0",self.gensigdatanovbf_,r.RooArgSet(self.bdtvar_,self.bdtmgg_),"bdtoutput>=0.05 && bdtoutput<0.545"))

    # throw signal toy in VBF cat
    if expSig>0:
      self.toyVBFEventsSig_ = self.rand_.Poisson(expSig*int(self.sigVBFEvents_))
      self.gensigmassdata_.append(self.mitsigpdfs_[4].generate(r.RooArgSet(self.bdtmgg_),int(self.toyVBFEventsSig_)))

      # combine signal in no VBF cats for IC analysis
      for i,dataset in enumerate(self.gensigmassdata_):
        dataset.SetName("gen_sig_mass_data_cat%d"%i)
        if i==4: continue
        elif i==0: self.allgensigmassdatanovbf_=r.RooDataSet("gen_sig_combmassdata_novbf","gen_sig_combmassdata_novbf",dataset,r.RooArgSet(self.bdtmgg_))
        else: self.allgensigmassdatanovbf_.append(dataset)

      # combine signal in VBF cats for IC analysis
      self.allgensigmassdata_=r.RooDataSet("gen_sig_combmassdata","gen_sig_combmassdata",self.allgensigmassdatanovbf_,r.RooArgSet(self.bdtmgg_))
      self.allgensigmassdata_.append(self.gensigmassdata_[4])
      
    # combine signal and background
      self.genbdtdata_.append(self.gensigdatanovbf_)
      for i, dataset in enumerate(self.genmassdata_):
        dataset.append(self.gensigmassdata_[i])
      self.allgenmassdatanovbf_.append(self.allgensigmassdatanovbf_)
      self.allgenmassdata_.append(self.allgensigmassdata_)
    
    # name gen datasets and make gen datahists to match workspace and save workspace
    self.outwsFile_ = r.TFile(outwsname,"RECREATE")
    self.outws_ = r.RooWorkspace("cms_hgg_workspace")
    self.genmassdatahist_=[]
    for i,dataset in enumerate(self.genmassdata_):
      dataset.SetName("data_mass_cat%i"%i)
      self.genmassdatahist_.append(r.RooDataHist("roohist_data_mass_cat%i"%i,"roohist_data_mass_cat%i"%i,r.RooArgSet(self.bdtmgg_),dataset))
      getattr(self.outws_,'import')(dataset)
      getattr(self.outws_,'import')(self.genmassdatahist_[i])
    self.outwsFile_.cd()
    self.outws_.Write()
    print 'Mass fac toy workspace written to file ', self.outwsFile_.GetName()
    self.outwsFile_.Close()

  def getToyVBFevents(self):
    print 'Toy VBF events', self.toyVBFEvents_
    return int(self.toyVBFEvents_)

  def returnWindowToyData(self,mH,size):
    returnList = []

    mHL = mH*(1.-size)
    mHH = mH*(1.+size)

    # CHEATING NICK TEST !!!!!!!!!!!!!!!!!!!
    
    
    for i in range(self.toyNoVBFEvents_):
      val_m = (self.allgenmassdatanovbf_.get(i)).getRealValue("CMS_hgg_mass");
      val_b = (self.genbdtdata_.get(i)).getRealValue("bdtoutput");
      if val_m > mHL and val_m < mHH and val_b >=0.05: returnList.append((val_b,((val_m-mH)/mH)))

    for i in range(self.toyVBFEvents_):
      val_m = (self.genmassdata_[4].get(i)).getRealValue("CMS_hgg_mass");
      val_b = 1.01
      if val_m > mHL and val_m < mHH and val_b >=0.05: returnList.append((val_b,((val_m-mH)/mH)))

    return returnList

  def plotData(self,nbinsMass,nbinsBDT):
    
    if not os.path.isdir("StuffForToys"):
      os.makedirs("StuffForToys")
    can = r.TCanvas()
    frame2 = self.bdtmgg_.frame()
    self.bdtdatacut_.plotOn(frame2,r.RooFit.Binning(nbinsMass))
    self.icpdf_.plotOn(frame2)
    frame2.SetTitle("Mass in data for bdtoutput>=0.05")
    frame2.Draw()
    can.SaveAs("StuffForToys/data_mass.pdf")

    frame1 = self.bdtvar_.frame()
    frame1.SetTitle("BDT output in data")
    self.bdtdata_.plotOn(frame1,r.RooFit.Binning(nbinsBDT))
    self.bdtdatavbf_.plotOn(frame1,r.RooFit.Binning(nbinsBDT),r.RooFit.MarkerColor(4),r.RooFit.LineColor(4))
    self.bdtdata_.plotOn(frame1,r.RooFit.Binning(nbinsBDT))
    self.bdtpdf_.plotOn(frame1)
    frame1.Draw()
    can.SaveAs("StuffForToys/data_bdt.pdf")

  def plotToy(self,nbinsMass,nbinsBDT):
    
    if not os.path.isdir("StuffForToys"):
      os.makedirs("StuffForToys")
    
    self.bdtvar_.setRange(0.05,1.)
    can = r.TCanvas()
    frame1 = self.bdtvar_.frame()
    self.genbdtdata_.plotOn(frame1,r.RooFit.Binning(nbinsBDT))
    self.bdtpdf_.plotOn(frame1)
    frame1.SetTitle("BDT output toy no VBF")
    frame1.Draw()
    can.SaveAs("StuffForToys/toy_bdt_novbf.pdf")

    # mass
    for i, dataset in enumerate(self.genmassdata_):
      frame1 = self.bdtmgg_.frame()
      frame1.SetTitle("Mass toy for category %i"%i)
      dataset.plotOn(frame1,r.RooFit.Binning(nbinsMass))
      self.mitpdfs_[i].plotOn(frame1)
      frame1.Draw()
      can.SaveAs("StuffForToys/toy_mass_cat%i.pdf"%i)

    frame3 = self.bdtmgg_.frame()
    frame3.SetTitle("Mass toy for all non VBF categories combined")
    self.allgenmassdatanovbf_.plotOn(frame3,r.RooFit.Binning(nbinsMass))
    self.icpdf_.plotOn(frame3,r.RooFit.Normalization(self.toyNoVBFEvents_,r.RooAbsReal.NumEvent),r.RooFit.LineColor(r.kMagenta))
    frame3.Draw()
    can.SaveAs("StuffForToys/toy_mass_comb_noVBF.pdf")

    frame4 = self.bdtmgg_.frame()
    frame4.SetTitle("Mass toy for all categories combined")
    self.allgenmassdata_.plotOn(frame4,r.RooFit.Binning(nbinsMass))
    self.icpdf_.plotOn(frame4,r.RooFit.Normalization(self.toyNoVBFEvents_+self.toyVBFEvents_,r.RooAbsReal.NumEvent),r.RooFit.LineColor(r.kMagenta))
    frame4.Draw()
    can.SaveAs("StuffForToys/toy_mass_comb.pdf")
      
    #self.icr1_       = r.RooRealVar("r1","r1",-8.,-10.,0.)
    #self.icr2_       = r.RooRealVar("r2","r2",-0.05,-10.,0.)
    #self.icf1_       = r.RooRealVar("f1","f1",0.01,0.,1.)
    #self.icpdfmass_  = r.RooGenericPdf("data_pow_model","data_pow_model","(1-@3)*TMath::Power(@0,@1)+@3*TMath::Power(@0,@2)",r.RooArgList(self.mitvar_,self.icr1_,self.icr2_,self.icf1_));

  """
  def createSigPdf(self,tree):
    if not self.hasFit_:
      sys.exit('There is no background PDF made or loaded yet. Can\'t create signal. Bailing out')

    # also need to get event weight for signal
    self.weight_ = r.RooRealVar("weight","weight",0,1)
    self.bdtsigdata_    = r.RooDataSet("sigdata_bdt","sigdata_bdt",r.RooArgSet(self.bdtvar_,self.var_,self.bdtvbf_,self.weight_),r.RooFit.Import(tree),r.RooFit.Cut("vbf==0"),"weight")
    self.bdtsigdatavbf_ = r.RooDataSet("sigdata_bdtVBF","sigdata_bdtVBF",r.RooArgSet(self.bdtvar_,self.var_,self.bdtvbf_,self.weight_),r.RooFit.Import(tree),r.RooFit.Cut("vbf==1"),"weight")

    self.sigNoVBFEvents_= self.bdtsigdata_.sumEntries("bdtoutput>=0.05")
    self.sigVBFEvents_  = self.bdtsigdatavbf_.sumEntries("bdtoutput>=0.05")
    self.var_.setBins(320)
    self.bdtvar_.setBins(200)
    self.bdtsighist_    = r.RooDataHist("sigbdthist","sigbdthist",r.RooArgSet(self.bdtvar_,self.var_),self.bdtsigdata_)
    self.bdtsigpdf_     = r.RooHistPdf("sigpdf","sigpdf",r.RooArgSet(self.bdtvar_,self.var_),self.bdtsighist_)

    print 'not vbf sig above -0.8 ', self.bdtsigdata_.sumEntries("bdtoutput>-0.8")
    print 'not vbf sig above 0.05 ', self.bdtsigdata_.sumEntries("bdtoutput>=0.05")
    print 'vbf sig above -0.8 ', self.bdtsigdatavbf_.sumEntries("bdtoutput>-0.8")
    print 'vbf sig above 0.05 ', self.bdtsigdatavbf_.sumEntries("bdtoutput>=0.05")
  
  def createKeysPdf(self,tree):
    self.bdtvar_      = r.RooRealVar("bdtoutput","bdtoutput",-1,1)
    self.bdtmgg_      = r.RooRealVar("mgg","mgg",100,180)
    self.bdtvbf_      = r.RooRealVar("vbf","vbf",0,1)
    self.dataforan_   = r.RooDataSet("foran","foran",r.RooArgSet(self.bdtvar_,self.bdtmgg_),r.RooFit.Import(tree),r.RooFit.Cut("bdtoutput>=0.05"))
    print 'all data above 0.05 ', self.dataforan_.sumEntries()
    
    self.bdtdata_     = r.RooDataSet("roodata_bdt","roodata_bdt",r.RooArgSet(self.bdtvar_,self.bdtmgg_,self.bdtvbf_),r.RooFit.Import(tree),r.RooFit.Cut("bdtoutput>-0.8 && vbf==0"))
    print 'not vbf data above -0.8  ', self.bdtdata_.sumEntries()
    print 'not vbf data above 0.05 ', self.bdtdata_.sumEntries("bdtoutput>=0.05")

    self.bdtdatavbf_  = r.RooDataSet("roodata_bdtVBF","roodata_bdtVBF",r.RooArgSet(self.bdtvar_,self.bdtmgg_,self.bdtvbf_),r.RooFit.Import(tree),r.RooFit.Cut("bdtoutput>-0.8 && vbf==1"))
    print 'vbf data above 0.8  ', self.bdtdatavbf_.sumEntries()
    print 'vbf data above 0.05 ', self.bdtdatavbf_.sumEntries("bdtoutput>=0.05")
    
    self.bdtvar_.setRange(-0.8,1.)
    self.bdtpdf_ = r.RooKeysPdf("bdtpdf","bdtpdf",self.bdtvar_,self.bdtdata_)

    self.dataNoVBFEvents_=self.bdtdata_.sumEntries("bdtoutput>=0.05")
    self.dataVBFEvents_=self.bdtdatavbf_.sumEntries("bdtoutput>=0.05")
    print self.dataNoVBFEvents_
    print self.dataVBFEvents_
    self.hasFit_      = True

  """

