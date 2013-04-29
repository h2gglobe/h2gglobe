#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "boost/lexical_cast.hpp"

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

#define MAX_REPEAT 2

using namespace std;
using namespace RooFit;

double getTotalEvents(map<string,double> events){

  double total=0.;
  for (map<string,double>::iterator it=events.begin(); it!=events.end(); it++){
    total += it->second;
  }
  return total;
}

void Plot(RooRealVar *mass, RooCategory *category, RooDataSet *data, RooSimultaneous *pdf, int nBDTCats, int nSpinCats, bool sm, string name, bool correlateCosThetaCategories=false){

  TCanvas *canv = new TCanvas();
  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      RooPlot *plot = mass->frame();
      data->plotOn(plot,Cut(Form("category==category::cat%d_spin%d",c,s)));
      if(!correlateCosThetaCategories)
        pdf->plotOn(plot,Slice(*category,Form("cat%d_spin%d",c,s)),Components(Form("data_pol_model_cat%d_spin%d",c,s)),ProjWData(*category,*data),LineStyle(kDashed),LineColor(kRed));
      else
        pdf->plotOn(plot,Slice(*category,Form("cat%d_spin%d",c,s)),Components(Form("data_pol_model_cat%d",c)),ProjWData(*category,*data),LineStyle(kDashed),LineColor(kRed));
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

  int nToys=0;
  string filename;
  int nBDTCats=0;
  int nSpinCats=0;
  bool globePDFs=false;
  std::vector<double> cosThetaBoundaries;
  string boundaries="";
  bool correlateCosThetaCategories=false;
  if (argc!=3) {
    cout << "usage ./bin/diySeparation <datfilename> <ntoys>" << endl;
    exit(1);
  }
  else{
    nToys=atoi(argv[2]);
    ifstream datfile;
    datfile.open(argv[1]);
    if (datfile.fail()){
      cout << "datfile " << argv[1] << " doesn't exist" << endl;
      exit(1);
    }
    while (datfile.good()){
      string line;
      getline(datfile,line);
      if (line.find("wsfile=")!=string::npos) filename = line.substr(line.find("=")+1,string::npos);
      if (line.find("globePDFs=")!=string::npos) globePDFs = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
      if (line.find("nBDTCats=")!=string::npos) nBDTCats = boost::lexical_cast<int>(line.substr(line.find("=")+1,string::npos));
      if (line.find("nSpinCats=")!=string::npos) nSpinCats = boost::lexical_cast<int>(line.substr(line.find("=")+1,string::npos));
      if (line.find("catBoundaries=")!=string::npos) boundaries = line.substr(line.find("=")+1,string::npos);
      if (line.find("correlateCosThetaCategories=")!=string::npos) correlateCosThetaCategories = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
    }
    datfile.close();

    if(boundaries != "")
    {
      size_t position = boundaries.find(" ");
      do
      {
        position = boundaries.find(" ");
        string testnum = boundaries.substr(0,position);
        boundaries = boundaries.substr(position+1);

        std::stringstream os(testnum);

        double temp;
        os >> temp;
        if(!os.fail())
        {
          cout << "Interpreted boundary '" << testnum << "' as: " << temp << endl;
          cosThetaBoundaries.push_back(temp);
          for(UInt_t i = cosThetaBoundaries.size()-1; i > 0; i--)
          {
            if(cosThetaBoundaries[i] > cosThetaBoundaries[i-1])
              std::swap(cosThetaBoundaries[i], cosThetaBoundaries[i-1]);
            else
              break;
          }
        }
        else
          cout << "Was not able understand boundary '" << testnum << "', ignoring it." << endl;
      }while(position!=string::npos);
      nSpinCats = cosThetaBoundaries.size()-1;
    }
    else
    {
      if(nSpinCats == 0)
        nSpinCats = 5;

      cout << "Using " << nSpinCats << " cos(Theta*) categories." << endl;

      //Boundaries are defined in (1-cosTheta) space
      for(int i = 0; i <= nSpinCats; i++)
      {
        cosThetaBoundaries.push_back((nSpinCats-i)*1./nSpinCats);
      }
    }
  }

  // set plotting style
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  RooMsgService::instance().setGlobalKillBelow(RooFit::MsgLevel(RooFit::ERROR));
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooFit::PrintLevel(-200000);
  RooMsgService::instance().setStreamStatus(1,false);

  TFile *inFile = TFile::Open(filename.c_str());
  RooWorkspace *work = (RooWorkspace*)inFile->FindObjectAny("HggSpinStudies");

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
      if(correlateCosThetaCategories)
      {
        if (globePDFs) bkgMod.insert(pair<string,RooAbsPdf*>(catname,(RooBernstein*)work->pdf(Form("data_pol_model_cat%d",c))));
        else bkgMod.insert(pair<string,RooAbsPdf*>(catname,(RooGenericPdf*)work->pdf(Form("data_pow_model_cat%d",c))));
      }
      else
      {
        if (globePDFs) bkgMod.insert(pair<string,RooAbsPdf*>(catname,(RooBernstein*)work->pdf(Form("data_pol_model_cat%d_spin%d",c,s))));
        else bkgMod.insert(pair<string,RooAbsPdf*>(catname,(RooGenericPdf*)work->pdf(Form("data_pow_model_cat%d_spin%d",c,s))));
      }

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

    Plot(mass,category,combData,simPdfSM,nBDTCats,nSpinCats,true,"smpdf_data", correlateCosThetaCategories);
    Plot(mass,category,combData,simPdfGRAV,nBDTCats,nSpinCats,false,"gravpdf_data", correlateCosThetaCategories);

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

  map<string,RooAbsPdf*> sigSM_gen;
  map<string,RooAbsPdf*> sigGRAV_gen;
  map<string,RooAbsPdf*> bkgMod_gen;
  for(int c = 0; c < nBDTCats; c++)
  {
    for(int s = 0; s < nSpinCats; s++)
    {
      string catname = Form("cat%d_spin%d",c,s);

      sigSM_gen.insert(pair<string,RooAbsPdf*>(catname,(RooAbsPdf*)sigSM[catname]->cloneTree()));
      sigGRAV_gen.insert(pair<string,RooAbsPdf*>(catname,(RooAbsPdf*)sigGRAV[catname]->cloneTree()));
      bkgMod_gen.insert(pair<string,RooAbsPdf*>(catname,(RooAbsPdf*)bkgMod[catname]->cloneTree()));
    }
  }

  TStopwatch total;
  total.Start();
  bool repeat = false;
  UInt_t count = 0;

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

        RooDataSet *bkgToy = (RooDataSet*)bkgMod_gen[catname]->generate(*mass,bkgToyEvents);
        RooDataSet *smToy = (RooDataSet*)sigSM_gen[catname]->generate(*mass,sigSMToyEvents);
        RooDataSet *gravToy = (RooDataSet*)sigGRAV_gen[catname]->generate(*mass,sigGRAVToyEvents);

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
      repeat = false;
      count = 0;
      do{
        RooFitResult *fitRes = simPdfSMvect[s]->fitTo(*combDataSMthisCTheta, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"), RooFit::PrintLevel(-1));
        if(fitRes->status() == -1)
          repeat = true;
        count++;
      }while(repeat && count<MAX_REPEAT);
      muSM_perCTbin[s] = muSM->getVal();

      repeat = false;
      count = 0;
      do{
        RooFitResult *fitRes = simPdfGRAVvect[s]->fitTo(*combDataSMthisCTheta, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"), RooFit::PrintLevel(-1));
        if(fitRes->status() == -1)
          repeat = true;
        count++;
      }while(repeat && count<MAX_REPEAT);
      muGRAV_perCTbin[s] = muGRAV->getVal();

    }

    //make global data
    RooDataHist *combDataSM = new RooDataHist(Form("combDataSM_toy%d",t),Form("combDataSM_toy%d",t),RooArgList(*mass),Index(*category),Import(toySM));
    RooDataHist *combDataGRAV = new RooDataHist(Form("combDataGRAV_toy%d",t),Form("combDataGRAV_toy%d",t),RooArgList(*mass),Index(*category),Import(toyGRAV));

    cout << "---------------------------" << endl;
    cout << "------Fitting toy " << t << " -------" << endl;
    cout << "---------------------------" << endl;

    muSM->setVal(1.);
    repeat = false;
    count = 0;
    RooFitResult *fitResSMSM;
    do{
      fitResSMSM = simPdfSM->fitTo(*combDataSM,Save(true), RooFit::Minimizer("Minuit2","Migrad"), RooFit::PrintLevel(-1));
      if(fitResSMSM->status() == -1)
        repeat = true;
      count++;
    }while(repeat && count<MAX_REPEAT);
    muSMSM_ = muSM->getVal();

    muSM->setVal(1.);
    repeat = false;
    count = 0;
    RooFitResult *fitResSMGRAV;
    do{
      fitResSMGRAV = simPdfSM->fitTo(*combDataGRAV,Save(true), RooFit::Minimizer("Minuit2","Migrad"), RooFit::PrintLevel(-1));
      if(fitResSMGRAV->status() == -1)
        repeat = true;
      count++;
    }while(repeat && count<MAX_REPEAT);
    muSMGRAV_ = muSM->getVal();

    muGRAV->setVal(1.);
    repeat = false;
    count = 0;
    RooFitResult *fitResGRAVSM;
    do{
      fitResGRAVSM = simPdfGRAV->fitTo(*combDataSM,Save(true), RooFit::Minimizer("Minuit2","Migrad"), RooFit::PrintLevel(-1));
      if(fitResGRAVSM->status() == -1)
        repeat = true;
      count++;
    }while(repeat && count<MAX_REPEAT);
    muGRAVSM_ = muGRAV->getVal();

    muGRAV->setVal(1.);
    repeat = false;
    count = 0;
    RooFitResult *fitResGRAVGRAV;
    do{
      fitResGRAVGRAV = simPdfGRAV->fitTo(*combDataGRAV,Save(true), RooFit::Minimizer("Minuit2","Migrad"), RooFit::PrintLevel(-1));
      if(fitResGRAVGRAV->status() == -1)
        repeat = true;
      count++;
    }while(repeat && count<MAX_REPEAT);
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
  total.Stop();
  cout << endl << endl << "Total time:" << endl;
  cout << "\t"; total.Print();

  outFile->cd();
  tree_->Write();
  outFile->Close();
  inFile->Close();

  for (map<string,RooAddPdf*>::iterator it=sigSM.begin(); it!=sigSM.end(); it++) delete it->second;
  for (map<string,RooAddPdf*>::iterator it=sigGRAV.begin(); it!=sigGRAV.end(); it++) delete it->second;
  if(!correlateCosThetaCategories) //If we are correlating the backgrounds in the cosTheta categories, it seems RooAbsPdf doesn't handle very well the multiple dependence, generating seg-faults
    for (map<string,RooAbsPdf*>::iterator it=bkgMod.begin(); it!=bkgMod.end(); it++) delete it->second;

  for (map<string,RooFormulaVar*>::iterator it=sigYieldSM.begin(); it!=sigYieldSM.end(); it++) delete it->second;
  for (map<string,RooFormulaVar*>::iterator it=sigYieldGRAV.begin(); it!=sigYieldGRAV.end(); it++) delete it->second;
  for (map<string,RooRealVar*>::iterator it=bkgYield.begin(); it!=bkgYield.end(); it++) delete it->second;

  for (map<string,RooAbsPdf*>::iterator it=sbModSM.begin(); it!=sbModSM.end(); it++) delete it->second;
  for (map<string,RooAbsPdf*>::iterator it=sbModGRAV.begin(); it!=sbModGRAV.end(); it++) delete it->second;


}
