#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "TFile.h"
#include "TMath.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TAxis.h"
#include "TStopwatch.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooBernstein.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooRandom.h"

using namespace std;
using namespace RooFit;

double getTotalEvents(map<string,double> events){
  
  double total=0.;
  for (map<string,double>::iterator it=events.begin(); it!=events.end(); it++){
    total += it->second;
  }
  return total;
}

void Plot(RooRealVar *mass, RooCategory *category, RooDataSet *data, RooSimultaneous *pdf, int nBDTCats, int nSpinCats, bool sm, string name){
  
  TCanvas *canv = new TCanvas();
  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      RooPlot *plot = mass->frame();
      data->plotOn(plot,Cut(Form("category==category::cat%d_spin%d",c,s)));
      pdf->plotOn(plot,Slice(*category,Form("cat%d_spin%d",c,s)),Components(Form("data_pol_model_cat%d_spin%d",c,s)),ProjWData(*category,*data),LineStyle(kDashed),LineColor(kRed));
      if (sm) pdf->plotOn(plot,Slice(*category,Form("cat%d_spin%d",c,s)),Components(Form("sigModel_SM_cat%d_spin%d",c,s)),ProjWData(*category,*data),LineColor(kRed));
      else pdf->plotOn(plot,Slice(*category,Form("cat%d_spin%d",c,s)),Components(Form("sigModel_GRAV_cat%d_spin%d",c,s)),ProjWData(*category,*data),LineColor(kBlue));
      pdf->plotOn(plot,Slice(*category,Form("cat%d_spin%d",c,s)),ProjWData(*category,*data),LineColor(kRed));
      plot->GetXaxis()->SetTitle("m_{#gamma#gamma}");
      plot->SetTitle(Form("cat%d_spin%d",c,s));
      plot->Draw();
      canv->Print(Form("%s_cat%d_spin%d.pdf",name.c_str(),c,s));
      canv->Print(Form("%s_cat%d_spin%d.png",name.c_str(),c,s));
    }
  }
  delete canv;
  
}

int main(int argc, char* argv[]){
 
  int nToys;
  string filename;
  int nBDTCats;
  int nSpinCats;
  bool globePDFs=false;
  if (argc!=5 && argc!=6) {
    cout << "usage ./bin/diySeparation <ntoys> <filename> <nBDTcats> <nSpinCats> --globePDFs" << endl;
    exit(1);
  }
  else {
    nToys=atoi(argv[1]);
    filename=string(argv[2]);
    nBDTCats=atoi(argv[3]);
    nSpinCats=atoi(argv[4]);
    for (int i=0; i<argc; i++) {
      if (string(argv[i])=="--globePDFs") globePDFs=true;
    }
  }

  // set plotting style
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  TFile *inFile = TFile::Open(filename.c_str());
  RooWorkspace *work = (RooWorkspace*)inFile->Get("HggSpinStudies");

  TFile *outFile = new TFile("LLOut.root","RECREATE");
  TTree *tree_ = new TTree("limit","limit");

  // define vars to save in output tree
  double q_data_=-999.;
  double q_smtoy_=-999.;
  double q_gravtoy_=-999.;
  double muSMSM_=-999.;
  double muSMGRAV_=-999.;
  double muGRAVSM_=-999.;
  double muGRAVGRAV_=-999.;
  double muSMData_=-999.;
  double muGRAVData_=-999.;
  double muSM_perCTbin[nSpinCats];
  double muGRAV_perCTbin[nSpinCats];

  tree_->Branch("q_data",&q_data_);
  tree_->Branch("q_smtoy",&q_smtoy_);
  tree_->Branch("q_gravtoy",&q_gravtoy_);
  tree_->Branch("muSMSM",&muSMSM_);
  tree_->Branch("muSMGRAV",&muSMGRAV_);
  tree_->Branch("muGRAVSM",&muGRAVSM_);
  tree_->Branch("muGRAVGRAV",&muGRAVGRAV_);
  tree_->Branch("muSMData",&muSMData_);
  tree_->Branch("muGRAVData",&muGRAVData_);
  tree_->Branch("nSpinCats",&nSpinCats);
  tree_->Branch("muSM_perCTbin",muSM_perCTbin,"muSM_perCTbin[nSpinCats]/D");
  tree_->Branch("muGRAV_perCTbin",muGRAV_perCTbin,"muGRAV_perCTbin[nSpinCats]/D");

  // get poi from workspace
  RooRealVar *mass = (RooRealVar*)work->var("mass");
  if (globePDFs) mass->setBins(160);
  else mass->setBins(100);
 
  // map for storing all data
  map<string,RooDataSet*> data;

  // map for storing expEvents
  map<string,double> expEventsSM;
  map<string,double> expEventsGRAV;
  map<string,double> expEventsALL;

  // map for storing yields
  map<string,RooFormulaVar*> sigYieldSM;
  map<string,RooFormulaVar*> sigYieldGRAV;
  map<string,RooRealVar*> bkgYield;

  // map for storing all pdfs;
  map<string,RooAddPdf*> sigSM;
  map<string,RooAddPdf*> sigGRAV;
  map<string,RooAbsPdf*> bkgMod;
  map<string,RooAbsPdf*> sbModSM;
  map<string,RooAbsPdf*> sbModGRAV;

  // define global category
  RooCategory *category = new RooCategory("category","category");
  
  // define per cosTheta cat categories
  vector<RooCategory*> miniCategories;
  
  // vectors for storing maps of pdfs and data in each cos theta cat
  vector<map<string,RooAbsPdf*> > sbModSMvect;
  vector<map<string,RooAbsPdf*> > sbModGRAVvect;
  vector<RooSimultaneous*> simPdfSMvect;
  vector<RooSimultaneous*> simPdfGRAVvect;
  vector<RooDataSet*> combDatavect;

  // mu's
  RooRealVar *muSM = new RooRealVar("muSM","muSM",0.,-10.,10.);
  RooRealVar *muGRAV = new RooRealVar("muGRAV","muGRAV",0.,-10.,10.);

  RooRandom::randomGenerator()->SetSeed(0);

  // first find expected events to calculate scale factor
  for (int s=0; s<nSpinCats; s++){
    for (int c=0; c<nBDTCats; c++){
      string catname = Form("cat%d_spin%d",c,s);
      
      expEventsSM.insert(pair<string,double>(catname,((RooDataSet*)work->data(Form("mcSigSMHiggs_bdt%d_cTh%d",c,s)))->sumEntries()));
      expEventsGRAV.insert(pair<string,double>(catname,((RooDataSet*)work->data(Form("mcSigGraviton_bdt%d_cTh%d",c,s)))->sumEntries()));
      expEventsALL.insert(pair<string,double>(catname,((RooDataSet*)work->data(Form("data_mass_cat%d_spin%d",c,s)))->sumEntries()));

      cout << "BDTCat:" << c << " SpinCat:" << s << " Data: " << expEventsALL[catname] << " Spin0: " << expEventsSM[catname] << " Spin2: " << expEventsGRAV[catname] << endl;
    }
  }

  double allSM = getTotalEvents(expEventsSM);
  double allGRAV = getTotalEvents(expEventsGRAV);
  double allData = getTotalEvents(expEventsALL);
  double scaleFactor = allSM/allGRAV;
  cout << "TOTAL: " << "Data: " << allData << " Spin0: " << allSM << " Spin2: " << allGRAV << " ScaleFac: " << scaleFactor << endl;

  for (int s=0; s<nSpinCats; s++){
    
    // define categories and maps for indepedent cos theta bin fits
    miniCategories.push_back(new RooCategory(Form("mincat_spin%d",s),Form("mincat_spin%d",s)));
    map<string,RooAbsPdf*> sbModSMthisCTheta;
    map<string,RooAbsPdf*> sbModGRAVthisCTheta;
    map<string,RooDataSet*> datathisCTheta;
    
    for (int c=0; c<nBDTCats; c++){
      string catname = Form("cat%d_spin%d",c,s);
      category->defineType(catname.c_str());
      miniCategories[s]->defineType(catname.c_str());

      // add to pdf maps
      sigSM.insert(pair<string,RooAddPdf*>(catname,(RooAddPdf*)work->pdf(Form("sigModel_SM_cat%d_spin%d",c,s))));
      sigGRAV.insert(pair<string,RooAddPdf*>(catname,(RooAddPdf*)work->pdf(Form("sigModel_GRAV_cat%d_spin%d",c,s))));
      if (globePDFs) bkgMod.insert(pair<string,RooAbsPdf*>(catname,(RooBernstein*)work->pdf(Form("data_pol_model_cat%d_spin%d",c,s))));
      else bkgMod.insert(pair<string,RooAbsPdf*>(catname,(RooGenericPdf*)work->pdf(Form("data_pow_model_cat%d_spin%d",c,s))));

      // add to yield maps
      sigYieldSM.insert(pair<string,RooFormulaVar*>(catname,new RooFormulaVar(Form("sigYield_SM_cat%d_spin%d",c,s),Form("sigYield_SM_cat%d_spin%d",c,s),Form("@0*%8.4f",expEventsSM[catname]),RooArgList(*muSM))));
      sigYieldGRAV.insert(pair<string,RooFormulaVar*>(catname,new RooFormulaVar(Form("sigYield_GRAV_cat%d_spin%d",c,s),Form("sigYield_GRAV_cat%d_spin%d",c,s),Form("@0*%8.4f",expEventsGRAV[catname]*scaleFactor),RooArgList(*muGRAV))));
      bkgYield.insert(pair<string,RooRealVar*>(catname,new RooRealVar(Form("bkgYield_cat%d_spin%d",c,s),Form("bkgYield_cat%d_spin%d",c,s),expEventsALL[catname],200.,1.e6)));

      // create signal+background model for SM and GRAV
      RooAddPdf *sbTempSM = new RooAddPdf(Form("sbModSM_cat%d_spin%d",c,s),Form("sbModSM_cat%d_spin%d",c,s),RooArgList(*sigSM[catname],*bkgMod[catname]),RooArgList(*sigYieldSM[catname],*bkgYield[catname]));
      RooAddPdf *sbTempGRAV = new RooAddPdf(Form("sbModGRAV_cat%d_spin%d",c,s),Form("sbModGRAV_cat%d_spin%d",c,s),RooArgList(*sigGRAV[catname],*bkgMod[catname]),RooArgList(*sigYieldGRAV[catname],*bkgYield[catname])); 
     
      // fill global pdf map
      sbModSM.insert(pair<string,RooAbsPdf*>(catname,sbTempSM));
      sbModGRAV.insert(pair<string,RooAbsPdf*>(catname,sbTempGRAV));

      // fill pdf map per cos theta cat
      sbModSMthisCTheta.insert(pair<string,RooAbsPdf*>(catname,sbTempSM));
      sbModGRAVthisCTheta.insert(pair<string,RooAbsPdf*>(catname,sbTempGRAV));

      // get data for this cat
      RooDataSet *tempData = (RooDataSet*)work->data(Form("data_mass_cat%d_spin%d",c,s));

      // fill data global and per cos theta cats data maps
      data.insert(pair<string,RooDataSet*>(catname,tempData));
      datathisCTheta.insert(pair<string,RooDataSet*>(catname,tempData));

    }

    // put the per cos theta cats pdfs in vectors
    sbModSMvect.push_back(sbModSMthisCTheta);
    sbModGRAVvect.push_back(sbModGRAVthisCTheta);
    
    // make per cos theta cat simulataneousPdfs and combDataSets and put them in vectors
    simPdfSMvect.push_back(new RooSimultaneous(Form("simPdfSM_spin%d",s),Form("simPdfSM_spin%d",s),sbModSMvect[s],*miniCategories[s]));
    simPdfGRAVvect.push_back(new RooSimultaneous(Form("simPdfGRAV_spin%d",s),Form("simPdfGRAV_spin%d",s),sbModGRAVvect[s],*miniCategories[s]));

    combDatavect.push_back(new RooDataSet(Form("combData_spin%d",s),Form("combData_spin%d",s),RooArgList(*mass),Index(*miniCategories[s]),Import(datathisCTheta)));
  }

  // make global simPdf and combData
  RooSimultaneous *simPdfSM = new RooSimultaneous("simPdfSM","simPdfSM",sbModSM,*category);
  RooSimultaneous *simPdfGRAV = new RooSimultaneous("simPdfGRAV","simPdfGRAV",sbModGRAV,*category);

  RooDataSet *combData = new RooDataSet("combData","combData",RooArgList(*mass),Index(*category),Import(data));

  // if nToys==0 fit data
  if (nToys==0) {
    RooFitResult *fitResDataSM = simPdfSM->fitTo(*combData,Save(true));
    muSMData_ = muSM->getVal();

    RooFitResult *fitResDataGRAV = simPdfGRAV->fitTo(*combData,Save(true));
    muGRAVData_ = muGRAV->getVal();

    fitResDataSM->floatParsFinal().Print("v");
    fitResDataGRAV->floatParsFinal().Print("v");
    
    q_data_ = 2.*(fitResDataSM->minNll()-fitResDataGRAV->minNll());
    
    delete fitResDataSM;
    delete fitResDataGRAV;

    Plot(mass,category,combData,simPdfSM,nBDTCats,nSpinCats,true,"smpdf_data");
    Plot(mass,category,combData,simPdfGRAV,nBDTCats,nSpinCats,false,"gravpdf_data");

    cout << "Expected events..." << endl;
    for (int c=0; c<nBDTCats; c++){
      for (int s=0; s<nSpinCats; s++){
        string catname = Form("cat%d_spin%d",c,s);
        cout << "c " << c << " s " << s << endl;
        cout << "\tSM:   " << expEventsSM[catname] << endl;
        cout << "\tGRAV: " << expEventsGRAV[catname] << endl;
        cout << "\tDATA: " << expEventsALL[catname] << endl;
      }
    }
    
    cout << "Global fit to data...." << endl;
    cout << "\tmuSM =   " << muSM->getVal() << endl;
    cout << "\tmuGRAV = " << muGRAV->getVal() << endl;

    for (int s=0; s<nSpinCats; s++){
      RooFitResult *fitResDataSMCT = simPdfSMvect[s]->fitTo(*combDatavect[s],Save(true));
      muSM_perCTbin[s] = muSM->getVal();

      RooFitResult *fitResDataGRAVCT = simPdfGRAVvect[s]->fitTo(*combDatavect[s],Save(true));
      muGRAV_perCTbin[s] = muGRAV->getVal();

      fitResDataSMCT->floatParsFinal().Print("v");
      fitResDataGRAVCT->floatParsFinal().Print("v");
      delete fitResDataSMCT;
      delete fitResDataGRAVCT;
    }
    tree_->Fill();
  }

  for (int t=0; t<nToys; t++){
    cout << "---------------------------" << endl;
    cout << "------Running toy " << t << " -------" << endl;
    cout << "---------------------------" << endl;
    TStopwatch sw;
    sw.Start();
   
    // global maps for toys
    map<string,RooDataHist*> toySM;
    map<string,RooDataHist*> toyGRAV;

    for (int s=0; s<nSpinCats; s++){
     
      // per cos theta bin cats for toys
      map<string,RooDataHist*> toySMthisCTheta;
      map<string,RooDataHist*> toyGRAVthisCTheta;

      for (int c=0; c<nBDTCats; c++){
        
        string catname = Form("cat%d_spin%d",c,s);

        // throw toys
        muSM->setVal(1.);
        muGRAV->setVal(1.);
        bkgYield[catname]->setVal(expEventsALL[catname]);
        double sigSMToyEvents = RooRandom::randomGenerator()->PoissonD(expEventsSM[catname]);
        double sigGRAVToyEvents = RooRandom::randomGenerator()->PoissonD(scaleFactor*expEventsGRAV[catname]);
        int bkgToyEvents = RooRandom::randomGenerator()->Poisson(expEventsALL[catname]);

        RooDataSet *bkgToy = (RooDataSet*)bkgMod[catname]->generate(*mass,bkgToyEvents);
        RooDataSet *smToy = (RooDataSet*)sigSM[catname]->generate(*mass,sigSMToyEvents);
        RooDataSet *gravToy = (RooDataSet*)sigGRAV[catname]->generate(*mass,sigGRAVToyEvents);
        
        smToy->append(*bkgToy);
        gravToy->append(*bkgToy);

        RooDataHist *tempToySM = new RooDataHist(Form("sm_toy%d_cat%d_spin%d",t,c,s),Form("sm_toy%d_cat%d_spin%d",t,c,s),RooArgSet(*mass),*smToy);
        RooDataHist *tempToyGRAV = new RooDataHist(Form("sm_toy%d_cat%d_spin%d",t,c,s),Form("sm_toy%d_cat%d_spin%d",t,c,s),RooArgSet(*mass),*gravToy);

        // add toy to global map
        toySM.insert(pair<string,RooDataHist*>(catname,tempToySM));
        toyGRAV.insert(pair<string,RooDataHist*>(catname,tempToyGRAV));

        // add toy to per cos theta cat map
        toySMthisCTheta.insert(pair<string,RooDataHist*>(catname,tempToySM));
        toyGRAVthisCTheta.insert(pair<string,RooDataHist*>(catname,tempToyGRAV));
      }
      
      // make per cos theta cat data
      RooDataHist *combDataSMthisCTheta = new RooDataHist(Form("combDataSM_spin%d_toy%d",s,t),Form("combDataSM_spin%d_toy%d",s,t),RooArgList(*mass),Index(*miniCategories[s]),Import(toySMthisCTheta));
      //RooDataHist *combDataGRAVthisCTheta = new RooDataHist(Form("combDataGRAV_spin%d_toy%d",s,t),Form("combDataGRAV_spin%d_toy%d",s,t),RooArgList(*mass),Index(*miniCategories[s]),Import(toyGRAVthisCTheta));

      // for per cos theta cat only want to fit to the same SM toy
      simPdfSMvect[s]->fitTo(*combDataSMthisCTheta);
      muSM_perCTbin[s] = muSM->getVal();
      
      simPdfGRAVvect[s]->fitTo(*combDataSMthisCTheta);
      muGRAV_perCTbin[s] = muGRAV->getVal();

    }
    
    //make global data
    RooDataHist *combDataSM = new RooDataHist(Form("combDataSM_toy%d",t),Form("combDataSM_toy%d",t),RooArgList(*mass),Index(*category),Import(toySM));
    RooDataHist *combDataGRAV = new RooDataHist(Form("combDataGRAV_toy%d",t),Form("combDataGRAV_toy%d",t),RooArgList(*mass),Index(*category),Import(toyGRAV));

    cout << "---------------------------" << endl;
    cout << "------Fitting toy " << t << " -------" << endl;
    cout << "---------------------------" << endl;
   
    muSM->setVal(1.);
    RooFitResult *fitResSMSM = simPdfSM->fitTo(*combDataSM,Save(true));
    muSMSM_ = muSM->getVal();

    muSM->setVal(1.);
    RooFitResult *fitResSMGRAV = simPdfSM->fitTo(*combDataGRAV,Save(true));
    muSMGRAV_ = muSM->getVal();
    
    muGRAV->setVal(1.);
    RooFitResult *fitResGRAVSM = simPdfGRAV->fitTo(*combDataSM,Save(true));
    muGRAVSM_ = muGRAV->getVal();
    
    muGRAV->setVal(1.);
    RooFitResult *fitResGRAVGRAV = simPdfGRAV->fitTo(*combDataGRAV,Save(true));
    muGRAVGRAV_ = muGRAV->getVal();
   
    cout << "Fits done. Getting NLL...." << endl;
    
    q_smtoy_ = 2.*(fitResSMSM->minNll()-fitResGRAVSM->minNll());
    q_gravtoy_ = 2.*(fitResSMGRAV->minNll()-fitResGRAVGRAV->minNll());

    cout << "TestStat: SM - " << q_smtoy_ << "  --  GRAV - " << q_gravtoy_ << endl;
    cout << "Global: muSMSM - " << muSMSM_ << " - muSMGRAV - " << muSMGRAV_ << " - muGRAVSM - " << muGRAVSM_ << " - muGRAVGRAV - " << muGRAVGRAV_ << endl;
    for (int s=0; s<nSpinCats; s++){
      cout << "SpinCat " << s << ": muSM - " << muSM_perCTbin[s] << " - muGRAV - " << muGRAV_perCTbin[s] << endl;
    }

    tree_->Fill();
    
    delete fitResSMSM;
    delete fitResSMGRAV;
    delete fitResGRAVSM;
    delete fitResGRAVGRAV;
    
    for (map<string,RooDataHist*>::iterator it=toySM.begin(); it!=toySM.end(); it++) delete it->second;
    for (map<string,RooDataHist*>::iterator it=toyGRAV.begin(); it!=toyGRAV.end(); it++) delete it->second;
    delete combDataSM;  
    delete combDataGRAV; 

    sw.Stop();
    cout << "Throwing and fitting toy took:" << endl;
    cout << "\t"; sw.Print();
    cout << q_smtoy_ << " -- " << q_gravtoy_ << endl;
  }
  
  outFile->cd();
  tree_->Write();
  outFile->Close();
  inFile->Close();

  for (map<string,RooAddPdf*>::iterator it=sigSM.begin(); it!=sigSM.end(); it++) delete it->second; 
  for (map<string,RooAddPdf*>::iterator it=sigGRAV.begin(); it!=sigGRAV.end(); it++) delete it->second;
  for (map<string,RooAbsPdf*>::iterator it=bkgMod.begin(); it!=bkgMod.end(); it++) delete it->second;
  
  for (map<string,RooFormulaVar*>::iterator it=sigYieldSM.begin(); it!=sigYieldSM.end(); it++) delete it->second;
  for (map<string,RooFormulaVar*>::iterator it=sigYieldGRAV.begin(); it!=sigYieldGRAV.end(); it++) delete it->second;
  for (map<string,RooRealVar*>::iterator it=bkgYield.begin(); it!=bkgYield.end(); it++) delete it->second;

  for (map<string,RooAbsPdf*>::iterator it=sbModSM.begin(); it!=sbModSM.end(); it++) delete it->second;
  for (map<string,RooAbsPdf*>::iterator it=sbModGRAV.begin(); it!=sbModGRAV.end(); it++) delete it->second;


}
