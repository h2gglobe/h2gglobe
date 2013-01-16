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
  if (argc!=5) {
    cout << "usage ./bin/diySeparation <ntoys> <filename> <nBDTcats> <nSpinCats>" << endl;
    exit(1);
  }
  else {
    nToys=atoi(argv[1]);
    filename=string(argv[2]);
    nBDTCats=atoi(argv[3]);
    nSpinCats=atoi(argv[4]);
  }

  // set plotting style
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  TFile *inFile = TFile::Open(filename.c_str());
  RooWorkspace *work = (RooWorkspace*)inFile->Get("HggSpinStudies");

  TFile *outFile = new TFile("LLOut.root","RECREATE");
  TTree *tree_ = new TTree("limit","limit");

  double minllSMSM_;
  double minllSMGRAV_;
  double minllGRAVSM_;
  double minllGRAVGRAV_;
  double q_smtoy_temp_;
  double q_gravtoy_temp_;
  double q_smtoy_;
  double q_gravtoy_;
  double muSMSM_;
  double muSMGRAV_;
  double muGRAVSM_;
  double muGRAVGRAV_;

  //tree_->Branch("minllSMSM",&minllSMSM_);
  //tree_->Branch("minllSMGRAV",&minllSMGRAV_);
  //tree_->Branch("minllGRAVSM",&minllGRAVSM_);
  //tree_->Branch("minllGRAVGRAV",&minllGRAVGRAV_);
  tree_->Branch("q_smtoy",&q_smtoy_);
  tree_->Branch("q_gravtoy",&q_gravtoy_);
  tree_->Branch("muSMSM",&muSMSM_);
  tree_->Branch("muSMGRAV",&muSMGRAV_);
  tree_->Branch("muGRAVSM",&muGRAVSM_);
  tree_->Branch("muGRAVGRAV",&muGRAVGRAV_);

  RooRealVar *mass = (RooRealVar*)work->var("mass");
  mass->setBins(160);
  
  map<string,RooDataSet*> data;

  map<string,RooAddPdf*> sigSM;
  map<string,RooAddPdf*> sigGRAV;
  map<string,RooBernstein*> bkgMod;
  
  map<string,double> expEventsSM;
  map<string,double> expEventsGRAV;
  map<string,double> expEventsALL;

  map<string,RooFormulaVar*> sigYieldSM;
  map<string,RooFormulaVar*> sigYieldGRAV;
  map<string,RooRealVar*> bkgYield;

  map<string,RooAbsPdf*> sbModSM;
  map<string,RooAbsPdf*> sbModGRAV;

  RooRealVar *muSM = new RooRealVar("muSM","muSM",0.,-1.,4.);
  RooRealVar *muGRAV = new RooRealVar("muGRAV","muGRAV",0.,-1.,4.);

  RooCategory *category = new RooCategory("category","category");

  RooRandom::randomGenerator()->SetSeed(0);

  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      string catname = Form("cat%d_spin%d",c,s);
      category->defineType(catname.c_str());
      
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

  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      string catname = Form("cat%d_spin%d",c,s);
      category->defineType(catname.c_str());

      sigSM.insert(pair<string,RooAddPdf*>(catname,(RooAddPdf*)work->pdf(Form("sigModel_SM_cat%d_spin%d",c,s))));
      sigGRAV.insert(pair<string,RooAddPdf*>(catname,(RooAddPdf*)work->pdf(Form("sigModel_GRAV_cat%d_spin%d",c,s))));
      bkgMod.insert(pair<string,RooBernstein*>(catname,(RooBernstein*)work->pdf(Form("data_pol_model_cat%d_spin%d",c,s))));

      sigYieldSM.insert(pair<string,RooFormulaVar*>(catname,new RooFormulaVar(Form("sigYield_SM_cat%d_spin%d",c,s),Form("sigYield_SM_cat%d_spin%d",c,s),Form("@0*%8.4f",expEventsSM[catname]),RooArgList(*muSM))));
      sigYieldGRAV.insert(pair<string,RooFormulaVar*>(catname,new RooFormulaVar(Form("sigYield_GRAV_cat%d_spin%d",c,s),Form("sigYield_GRAV_cat%d_spin%d",c,s),Form("@0*%8.4f",expEventsGRAV[catname]*scaleFactor),RooArgList(*muGRAV))));
      bkgYield.insert(pair<string,RooRealVar*>(catname,new RooRealVar(Form("bkgYield_cat%d_spin%d",c,s),Form("bkgYield_cat%d_spin%d",c,s),expEventsALL[catname],200.,1.e6)));

      sbModSM.insert(pair<string,RooAbsPdf*>(catname, new RooAddPdf(Form("sbModSM_cat%d_spin%d",c,s),Form("sbModSM_cat%d_spin%d",c,s),RooArgList(*sigSM[catname],*bkgMod[catname]),RooArgList(*sigYieldSM[catname],*bkgYield[catname]))));
      sbModGRAV.insert(pair<string,RooAbsPdf*>(catname, new RooAddPdf(Form("sbModGRAV_cat%d_spin%d",c,s),Form("sbModGRAV_cat%d_spin%d",c,s),RooArgList(*sigGRAV[catname],*bkgMod[catname]),RooArgList(*sigYieldGRAV[catname],*bkgYield[catname]))));

      data.insert(pair<string,RooDataSet*>(catname, (RooDataSet*)work->data(Form("data_mass_cat%d_spin%d",c,s))));

    }
  }
 
  RooSimultaneous *simPdfSM = new RooSimultaneous("simPdfSM","simPdfSM",sbModSM,*category);
  RooSimultaneous *simPdfGRAV = new RooSimultaneous("simPdfGRAV","simPdfGRAV",sbModGRAV,*category);

  RooDataSet *combData = new RooDataSet("combData","combData",RooArgList(*mass),Index(*category),Import(data));

  if (nToys==0) {
    RooFitResult *fitResDataSM = simPdfSM->fitTo(*combData,Save(true));
    RooFitResult *fitResDataGRAV = simPdfGRAV->fitTo(*combData,Save(true));
    /*
    while (fitResDataSM->covQual()!=3){
      fitResDataSM = simPdfSM->fitTo(*combData,Save(true),NumCPU(8));
    }
    while (fitResDataGRAV->covQual()!=3){
      fitResDataGRAV = simPdfGRAV->fitTo(*combData,Save(true),NumCPU(8));
    }
    */
    fitResDataSM->floatParsFinal().Print("v");
    fitResDataGRAV->floatParsFinal().Print("v");

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
    
    cout << "Fit to data...." << endl;
    cout << "\tmuSM =   " << muSM->getVal() << endl;
    cout << "\tmuGRAV = " << muGRAV->getVal() << endl;
     
  }

  for (int t=0; t<nToys; t++){
    cout << "---------------------------" << endl;
    cout << "------Running toy " << t << " -------" << endl;
    cout << "---------------------------" << endl;
    TStopwatch sw;
    sw.Start();
   
    map<string,RooDataHist*> toySM;
    map<string,RooDataHist*> toyGRAV;

    for (int c=0; c<nBDTCats; c++){
      for (int s=0; s<nSpinCats; s++){
        
        string catname = Form("cat%d_spin%d",c,s);

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

        toySM.insert(pair<string,RooDataHist*>(catname, new RooDataHist(Form("sm_toy%d_cat%d_spin%d",t,c,s),Form("sm_toy%d_cat%d_spin%d",t,c,s),RooArgSet(*mass),*smToy)));
        toyGRAV.insert(pair<string,RooDataHist*>(catname, new RooDataHist(Form("sm_toy%d_cat%d_spin%d",t,c,s),Form("sm_toy%d_cat%d_spin%d",t,c,s),RooArgSet(*mass),*gravToy)));
      }
    }

    RooDataHist *combDataSM = new RooDataHist(Form("combDataSM_toy%d",t),Form("combDataSM_toy%d",t),RooArgList(*mass),Index(*category),Import(toySM));
    RooDataHist *combDataGRAV = new RooDataHist(Form("combDataGRAV_toy%d",t),Form("combDataGRAV_toy%d",t),RooArgList(*mass),Index(*category),Import(toyGRAV));

    cout << "---------------------------" << endl;
    cout << "------Fitting toy " << t << " -------" << endl;
    cout << "---------------------------" << endl;
   
    muSM->setVal(1.);
    RooFitResult *fitResSMSM = simPdfSM->fitTo(*combDataSM,Save(true));
    muSMSM_ = muSM->getVal();
    minllSMSM_ = (simPdfSM->createNLL(*combDataSM))->getVal();

    muSM->setVal(1.);
    RooFitResult *fitResSMGRAV = simPdfSM->fitTo(*combDataGRAV,Save(true));
    muSMGRAV_ = muSM->getVal();
    minllSMGRAV_ = (simPdfSM->createNLL(*combDataGRAV))->getVal();
    
    muGRAV->setVal(1.);
    RooFitResult *fitResGRAVSM = simPdfGRAV->fitTo(*combDataSM,Save(true));
    muGRAVSM_ = muGRAV->getVal();
    minllGRAVSM_ = (simPdfGRAV->createNLL(*combDataSM))->getVal();
    
    muGRAV->setVal(1.);
    RooFitResult *fitResGRAVGRAV = simPdfGRAV->fitTo(*combDataGRAV,Save(true));
    muGRAVGRAV_ = muGRAV->getVal();
    minllGRAVGRAV_ = (simPdfGRAV->createNLL(*combDataGRAV))->getVal();
   
    cout << "Fits done. Getting NLL...." << endl;
    
    q_smtoy_temp_ = 2.*(minllSMSM_-minllGRAVSM_);
    q_gravtoy_temp_ = 2.*(minllSMGRAV_-minllGRAVGRAV_);
    
    q_smtoy_ = 2.*(fitResSMSM->minNll()-fitResGRAVSM->minNll());
    q_gravtoy_ = 2.*(fitResSMGRAV->minNll()-fitResGRAVGRAV->minNll());

    cout << "TestStat: SM - " << q_smtoy_temp_ << "  --  GRAV - " << q_gravtoy_temp_ << endl;
    cout << "TestStat: SM - " << q_smtoy_ << "  --  GRAV - " << q_gravtoy_ << endl;

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
  for (map<string,RooBernstein*>::iterator it=bkgMod.begin(); it!=bkgMod.end(); it++) delete it->second;
  
  for (map<string,RooFormulaVar*>::iterator it=sigYieldSM.begin(); it!=sigYieldSM.end(); it++) delete it->second;
  for (map<string,RooFormulaVar*>::iterator it=sigYieldGRAV.begin(); it!=sigYieldGRAV.end(); it++) delete it->second;
  for (map<string,RooRealVar*>::iterator it=bkgYield.begin(); it!=bkgYield.end(); it++) delete it->second;

  for (map<string,RooAbsPdf*>::iterator it=sbModSM.begin(); it!=sbModSM.end(); it++) delete it->second;
  for (map<string,RooAbsPdf*>::iterator it=sbModGRAV.begin(); it!=sbModGRAV.end(); it++) delete it->second;


}
