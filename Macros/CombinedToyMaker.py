import ROOT as r
import sys

class CombinedToyMaker:

  def __init__(self,mitFileName):
    
    r.gROOT.SetStyle("Plain")
    r.gROOT.SetBatch(True)

    self.mitFileName_ = mitFileName
    self.mitFile_     = r.TFile(mitFileName)
    self.mitWS_       = self.mitFile_.Get("cms_hgg_workspace")
    self.var_         = self.mitWS_.var("CMS_hgg_mass")
    self.mitpdfs_     =[]
    for i in range(5): self.mitpdfs_.append(self.mitWS_.pdf("pdf_data_pol_model_cat%i"%i))
    self.hasFit_      = False
    self.rand_        = r.TRandom3(0)

    #self.icr1_       = r.RooRealVar("r1","r1",-8.,-10.,0.)
    #self.icr2_       = r.RooRealVar("r2","r2",-0.05,-10.,0.)
    #self.icf1_       = r.RooRealVar("f1","f1",0.01,0.,1.)
    #self.icpdfmass_  = r.RooGenericPdf("data_pow_model","data_pow_model","(1-@3)*TMath::Power(@0,@1)+@3*TMath::Power(@0,@2)",r.RooArgList(self.mitvar_,self.icr1_,self.icr2_,self.icf1_));

  def createKeysPdf(self,tree):
    self.bdtvar_      = r.RooRealVar("bdtoutput","bdtoutput",-1,1)
    self.bdtmgg_      = r.RooRealVar("mgg","mgg",100,180)
    self.bdtvbf_      = r.RooRealVar("vbf","vbf",0,1)
    self.dataforan_   = r.RooDataSet("foran","foran",r.RooArgSet(self.bdtvar_,self.bdtmgg_),r.RooFit.Import(tree),r.RooFit.Cut("bdtoutput>=0.05"))
    print 'all data above 0.05 ', self.dataforan_.sumEntries()
    
    self.bdtdata_     = r.RooDataSet("roodata_bdt","roodata_bdt",r.RooArgSet(self.bdtvar_,self.bdtmgg_,self.bdtvbf_),r.RooFit.Import(tree),r.RooFit.Cut("bdtoutput>-0.8 && vbf==0"))
    print 'not vbf data above 0.8  ', self.bdtdata_.sumEntries()
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

  def loadKeysPdf(self,filename):
    fi = r.TFile(filename)
    self.keyws_       = fi.Get("toysworkspace")
    self.bdtvar_      = self.keyws_.var("bdtoutput")
    self.bdtmgg_      = self.keyws_.var("mgg")
    self.bdtvbf_      = self.keyws_.var("vbf")
    self.bdtdata_     = self.keyws_.data("roodata_bdt")
    self.bdtdatavbf_  = self.keyws_.data("roodata_bdtVBF")
    self.bdtpdf_      = self.keyws_.pdf("bdtpdf")
    self.dataforan_   = self.keyws_.data("foran")
    self.dataNoVBFEvents_=self.bdtdata_.sumEntries("bdtoutput>=0.05")
    self.dataVBFEvents_=self.bdtdatavbf_.sumEntries("bdtoutput>=0.05")
    print self.dataNoVBFEvents_
    print self.dataVBFEvents_
    self.hasFit_      = True

  def savePdfWorkspace(self,filename):
    fi = r.TFile(filename,"RECREATE")
    fi.cd()
    self.ws_          = r.RooWorkspace("toysworkspace")
    getattr(self.ws_,'import')(self.bdtdata_)
    getattr(self.ws_,'import')(self.bdtdatavbf_)
    getattr(self.ws_,'import')(self.bdtpdf_)
    getattr(self.ws_,'import')(self.dataforan_)
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
      print ev, ' events in cat ', i
    return events

  def genData(self,outwsname):
    if not self.hasFit_: sys.exit('Diphoton output not yet fitted. Bailing out.')

    # first gen around bdt keys pdf
    self.toyNoVBFNEvents_ = self.rand_.Poisson(int(self.dataNoVBFEvents_))
    self.bdtvar_.setRange(0.05,1.)
    self.genbdtdata_  = self.bdtpdf_.generate(r.RooArgSet(self.bdtvar_),int(self.toyNoVBFNEvents_))
    self.genbdtdata_.SetName("gen_diphoton_output")
   
    eventsPerCat = self.getEventsPerCat(self.genbdtdata_)
    if len(eventsPerCat)!=len(self.mitpdfs_)-1: sys.exit('Different numbers of categories. Bailing out.')
    
    # given nevs generate around each non VBF cat
    can = r.TCanvas()
    self.genmassdata_ = []
    for i,ev in enumerate(eventsPerCat):
      self.genmassdata_.append(self.mitpdfs_[i].generate(r.RooArgSet(self.var_),int(ev)))
    
    # combine data in each non VBF cat for IC analysis
    for i,dataset in enumerate(self.genmassdata_):
      dataset.SetName("gen_mass_data_cat%d"%i)
      if i==0: self.allgenmassdata_=r.RooDataSet("gen_combmassdata","gen_combmassdata",dataset,r.RooArgSet(self.var_))
      else: self.allgenmassdata_.append(dataset)

    # throw toy in VBF cat
    self.toyVBFEvents_ = self.rand_.Poisson(int(self.dataVBFEvents_))
    self.genmassdata_.append(self.mitpdfs_[4].generate(r.RooArgSet(self.var_),int(self.toyVBFEvents_)))

    # name gen datasets to match workspace and save workspace
    self.outwsFile_ = r.TFile(outwsname,"RECREATE")
    self.outws_ = r.RooWorkspace("cms_hgg_workspace")
    for i,dataset in enumerate(self.genmassdata_):
      frameM = self.var_.frame()
      dataset.SetName("roohist_data_mass_cat%i"%i)
      getattr(self.outws_,'import')(dataset)
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

    for i in range(self.toyNoVBFNEvents_):
      val_m = (self.allgenmassdata_.get(i)).getRealValue("CMS_hgg_mass");
      val_b = (self.genbdtdata_.get(i)).getRealValue("bdtoutput");
      if val_m > mHL and val_m < mHH and val_b >=0.05: returnList.append((val_b,((val_m-mH)/mH)))

    for i in range(self.toyVBFEvents_):
      val_m = (self.genmassdata_[4].get(i)).getRealValue("CMS_hgg_mass");
      val_b = 1.01
      if val_m > mHL and val_m < mHH and val_b >=0.05: returnList.append((val_b,((val_m-mH)/mH)))

    return returnList

  def plotData(self,nbinsMass,nbinsBDT):
    
    can = r.TCanvas()
    frame2 = self.bdtmgg_.frame()
    self.dataforan_.plotOn(frame2,r.RooFit.Binning(nbinsMass))
    frame2.SetTitle("Mass in data for bdtoutput>=0.05")
    frame2.Draw()
    can.SaveAs("StuffForToys/data_mass.pdf")
   
    self.bdtvar_.setRange(-0.8,1.)
    frame1 = self.bdtvar_.frame()
    frame1.SetTitle("BDT output in data")
    self.dataforan_.plotOn(frame1,r.RooFit.Binning(nbinsBDT))
    self.bdtdatavbf_.plotOn(frame1,r.RooFit.Binning(nbinsBDT),r.RooFit.MarkerColor(4),r.RooFit.LineColor(4))
    self.bdtdata_.plotOn(frame1,r.RooFit.Binning(nbinsBDT),r.RooFit.MarkerColor(2),r.RooFit.LineColor(2))
    self.bdtpdf_.plotOn(frame1)
    frame1.Draw()
    can.SaveAs("StuffForToys/data_bdt.pdf")

  def plotToy(self,nbinsMass,nbinsBDT):
    can = r.TCanvas()
    frame1 = self.bdtvar_.frame()
    self.genbdtdata_.plotOn(frame1,r.RooFit.Binning(nbinsBDT))
    self.bdtpdf_.plotOn(frame1)
    frame1.SetTitle("BDT output toy no VBF")
    frame1.Draw()
    can.SaveAs("StuffForToys/toy_bdt_novbf.pdf")

    # mass
    for i, dataset in enumerate(self.genmassdata_):
      frame1 = self.var_.frame()
      frame1.SetTitle("Mass toy for category %i"%i)
      dataset.plotOn(frame1,r.RooFit.Binning(nbinsMass))
      self.mitpdfs_[i].plotOn(frame1)
      frame1.Draw()
      can.SaveAs("StuffForToys/toy_mass_cat%i.pdf"%i)

    frame3 = self.var_.frame()
    frame3.SetTitle("Mass toy for all non VBF categories combined")
    self.allgenmassdata_.plotOn(frame3,r.RooFit.Binning(nbinsMass))
    frame3.Draw()
    can.SaveAs("StuffForToys/toy_mass_comb_noVBF.pdf")

      
