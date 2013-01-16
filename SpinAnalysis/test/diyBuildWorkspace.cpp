#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"
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

pair<int,int> getMassFacCategory(float diphoton_bdt, double costheta_cs){
  
  // fabrice style
  int mycat=-1;
  int costheta_csBin=-1; 
  if (TMath::Abs(costheta_cs)>=0. && TMath::Abs(costheta_cs)<0.2) {
    costheta_csBin=0;
    if (diphoton_bdt >= 0.1 && diphoton_bdt < 0.7) mycat=0;
    else if (diphoton_bdt >= 0.7) mycat=1;
    else mycat=-1;
  }
  else if (TMath::Abs(costheta_cs)>=0.2 && TMath::Abs(costheta_cs)<0.4) {
    costheta_csBin=1;
    if (diphoton_bdt >= -0.05 && diphoton_bdt < 0.65) mycat=0;
    else if (diphoton_bdt >= 0.65) mycat=1;
    else mycat=-1;
  }
  else if (TMath::Abs(costheta_cs)>=0.4 && TMath::Abs(costheta_cs)<0.6) {
    costheta_csBin=2;
    if (diphoton_bdt >= -0.2 && diphoton_bdt < 0.6) mycat=0;
    else if (diphoton_bdt >= 0.6) mycat=1;
    else mycat=-1;
  }
  else if (TMath::Abs(costheta_cs)>=0.6 && TMath::Abs(costheta_cs)<0.8) {
    costheta_csBin=3;
    if (diphoton_bdt >= -0.4 && diphoton_bdt < 0.3) mycat=0;
    else if (diphoton_bdt >= 0.3) mycat=1;
    else mycat=-1;
  }
  else if (TMath::Abs(costheta_cs)>=0.8 && TMath::Abs(costheta_cs)<=1.0) {
    costheta_csBin=4;
    if (diphoton_bdt >= 0.1 && diphoton_bdt < 0.85) mycat=0;
    else if (diphoton_bdt >= 0.85) mycat=1;
    else mycat=-1;
  }
  else {
    costheta_csBin=-1;
    mycat=-1;
  }

  return pair<int,int>(mycat,costheta_csBin);
}

pair<int,int> getCutBasedCategory(double lead_eta, double sublead_eta, float lead_r9, float sublead_r9, double costheta_cs){
  
  // marco style
  int mycat=-1;
  int costheta_csBin=-1; 
  if (TMath::Abs(costheta_cs)>=0. && TMath::Abs(costheta_cs)<0.2) costheta_csBin=0;
  else if (TMath::Abs(costheta_cs)>=0.2 && TMath::Abs(costheta_cs)<0.4) costheta_csBin=1;
  else if (TMath::Abs(costheta_cs)>=0.4 && TMath::Abs(costheta_cs)<0.6) costheta_csBin=2;
  else if (TMath::Abs(costheta_cs)>=0.6 && TMath::Abs(costheta_cs)<0.8) costheta_csBin=3;
  else if (TMath::Abs(costheta_cs)>=0.8 && TMath::Abs(costheta_cs)<=1.0) costheta_csBin=4;
  else costheta_csBin=-1;

  if (TMath::Abs(lead_eta)<=1.4442 && TMath::Abs(sublead_eta)<=1.4442){
    if (lead_r9>=0.94 && sublead_r9>=0.94) mycat=0;
    else mycat=1;
  }
  else {
    if (lead_r9>=0.94 && sublead_r9>=0.94) mycat=2;
    else mycat=3;
  }

  return pair<int,int>(mycat,costheta_csBin);
}

string getCatName(int cat, int scat){
  return Form("cat%d_spin%d",cat,scat); 
}

string getCatName(pair<int,int> cat){
  return Form("cat%d_spin%d",cat.first,cat.second); 
}

void fillDataSet(RooRealVar *mass, map<string,RooDataSet*> &dataMap, TTree *tree, bool isMassFac){
  
  int category;
  float evweight;
  double costheta_cs;
  double lead_calo_eta;
  double sublead_calo_eta;
  float lead_r9;
  float sublead_r9;
  float diphoton_bdt;
  double higgs_mass;

  tree->SetBranchAddress("category",&category);
  tree->SetBranchAddress("evweight",&evweight);
  tree->SetBranchAddress("costheta_cs",&costheta_cs);
  tree->SetBranchAddress("lead_calo_eta",&lead_calo_eta);
  tree->SetBranchAddress("sublead_calo_eta",&sublead_calo_eta);
  tree->SetBranchAddress("lead_r9",&lead_r9);
  tree->SetBranchAddress("sublead_r9",&sublead_r9);
  tree->SetBranchAddress("diphoton_bdt",&diphoton_bdt);
  tree->SetBranchAddress("higgs_mass",&higgs_mass);

  cout << "Filling from " << tree->GetName() << endl;
  for (int e=0; e<tree->GetEntries(); e++){
    tree->GetEntry(e);
    if (e%10000==0) cout << "Entry " << e << " of " << tree->GetEntries() << endl;
    pair<int,int> mycat;
    if (isMassFac) mycat = getMassFacCategory(diphoton_bdt,costheta_cs);
    else mycat = getCutBasedCategory(lead_calo_eta,sublead_calo_eta,lead_r9,sublead_r9,costheta_cs);
    if (mycat.first==-1 || mycat.second==-1) continue; 
    mass->setVal(higgs_mass);
    string catname = getCatName(mycat);
    dataMap[catname]->add(*mass,evweight);
  }
}

void Plot(string name, RooRealVar *mass, RooDataSet* data, RooAbsPdf* pdf){
  
  TCanvas *canv = new TCanvas();
  RooPlot *plot = mass->frame();
  data->plotOn(plot);
  pdf->plotOn(plot);
  plot->Draw();
  canv->Print(Form("%s.pdf",name.c_str()));
  canv->Print(Form("%s.png",name.c_str()));
  delete plot;
  delete canv;
}


void makeBackgroundModel(string name, RooRealVar *mass, map<string,RooDataSet*> data, RooWorkspace *work, int nBDTCategories, int nCosThetaCategories, bool verbose=false){

  string path="plots/"+name; 
  system(Form("mkdir -p %s",path.c_str()));

  for(int a=0; a<nBDTCategories; a++){
    for(int b=0; b<nCosThetaCategories; b++){
      
      string catname = getCatName(a,b);

      RooRealVar *p0 = new RooRealVar(Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);
      RooRealVar *p1 = new RooRealVar(Form("CMS_hgg_quartic1_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);
      RooRealVar *p2 = new RooRealVar(Form("CMS_hgg_quartic2_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);
      RooRealVar *p3 = new RooRealVar(Form("CMS_hgg_quartic3_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);
      
      RooFormulaVar *f0 = new RooFormulaVar(Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p0));
      RooFormulaVar *f1 = new RooFormulaVar(Form("CMS_hgg_modquartic1_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p1));
      RooFormulaVar *f2 = new RooFormulaVar(Form("CMS_hgg_modquartic2_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p2));
      RooFormulaVar *f3 = new RooFormulaVar(Form("CMS_hgg_modquartic3_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p3));

      RooBernstein *bern = new RooBernstein(Form("data_pol_model_cat%d_spin%d",a,b),Form("data_pol_model_cat%d_spin%d",a,b),*mass,RooArgList(*f0,*f1,*f2,*f3));
      RooRealVar *norm = new RooRealVar(Form("%s_norm",bern->GetName()),Form("%s_norm",bern->GetName()),10.,1.e6);
      RooExtendPdf *ext = new RooExtendPdf(Form("pdf_data_pol_model_cat%d_spin%d",a,b),Form("pdf_data_pol_model_cat%d_spin%d",a,b),*bern,*norm);

      cout << data[catname]->GetName() << " : " << data[catname]->sumEntries() << endl;
      
      bern->fitTo(*data[catname],PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));

      RooFitResult *fitRes;
      verbose ? 
        fitRes = ext->fitTo(*data[catname],Save(true)) :
        fitRes = ext->fitTo(*data[catname],Save(true),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1)) 
      ;
      fitRes->floatParsFinal().Print("v");
      delete fitRes;

      work->import(*data[catname]);
      work->import(*ext);

      Plot(Form("%s/bkg_cat%d_spin%d",path.c_str(),a,b),mass,data[catname],ext);
    }
  }
}

void makeSignalModel(string name, RooRealVar *mass, map<string,RooDataSet*> sigMC, RooWorkspace *work, int nBDTCategories, int nCosThetaCategories, bool isSM, bool verbose=false){
 
  string path="plots/"+name; 
  system(Form("mkdir -p %s",path.c_str()));
  for(int a=0; a<nBDTCategories; a++){
    for(int b=0; b<nCosThetaCategories; b++){
     
      string catname = getCatName(a,b);

      string postfix="";
      if (isSM) postfix="SM";
      else postfix="GRAV";

      RooAddPdf *sigPdf;
      vector<RooGaussian*> pdfGaussians;
      vector<RooRealVar*> means;
      vector<RooRealVar*> sigmas;
      vector<RooRealVar*> fracs;
      for(int i=0; i<4;i++){
        string nameMean     = Form("mean%d_%s_cat%d_spin%d"         ,i,postfix.c_str(),a,b);
        string nameSigma    = Form("sigma%d_%s_cat%d_spin%d"        ,i,postfix.c_str(),a,b);
        string nameGauss    = Form("gauss%d_%s_cat%d_spin%d"        ,i,postfix.c_str(),a,b);
        means.push_back(new RooRealVar(nameMean.c_str(),nameMean.c_str(),125.,115.,135.));
        sigmas.push_back(new RooRealVar(nameSigma   .c_str(),nameSigma   .c_str(),  2., 0.5,  20.));
        pdfGaussians.push_back(new RooGaussian(nameGauss.c_str(),nameGauss.c_str(),*mass,*means[i],*sigmas[i]));
        if (i<3){
          string nameFrac     = Form("frac%d_%s_cat%d_spin%d"         ,i,postfix.c_str(),a,b);
          fracs.push_back(new RooRealVar(nameFrac   .c_str(),nameFrac   .c_str(),  0.01, 0.,  1.));
        }
      }

      sigPdf = new RooAddPdf(Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b),Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b), RooArgList(*pdfGaussians[0],*pdfGaussians[1],*pdfGaussians[2],*pdfGaussians[3]),RooArgList(*fracs[0],*fracs[1],*fracs[2]));
     
      cout << sigMC[catname]->GetName() << " : " << sigMC[catname]->sumEntries() << endl;

      RooFitResult *fitRes;
      verbose ? 
        fitRes = sigPdf->fitTo(*sigMC[catname],SumW2Error(true),Save(true)) :
        fitRes = sigPdf->fitTo(*sigMC[catname],SumW2Error(true),Save(true),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1))
      ;
      fitRes->floatParsFinal().Print("v");
      delete fitRes;
      
      for (vector<RooRealVar*>::iterator it=means.begin(); it!=means.end(); it++) (*it)->setConstant(true);
      for (vector<RooRealVar*>::iterator it=sigmas.begin(); it!=sigmas.end(); it++) (*it)->setConstant(true);
      for (vector<RooRealVar*>::iterator it=fracs.begin(); it!=fracs.end(); it++) (*it)->setConstant(true);

      RooRealVar *norm;
      norm = new RooRealVar(Form("%s_norm",sigPdf->GetName()),Form("%s_norm",sigPdf->GetName()),sigMC[catname]->sumEntries());

      RooExtendPdf *ext = new RooExtendPdf(Form("pdf_%s",sigPdf->GetName()),Form("pdf_%s",sigPdf->GetName()),*sigPdf,*norm);
      work->import(*sigMC[catname]);
      work->import(*ext);
      
      if (isSM) Plot(Form("%s/sig_spin0_cat%d_spin%d",path.c_str(),a,b),mass,sigMC[catname],sigPdf);
      else Plot(Form("%s/sig_spin2_cat%d_spin%d",path.c_str(),a,b),mass,sigMC[catname],sigPdf);
    }
  }
}

int main(int argc, char* argv[]){

  string filename;
  bool isMassFac=false;
  if (argc!=2 && argc!=3) {
    cout << "usage ./bin/diyBuildWorkspace <filename> --isMassFac" << endl;
    exit(1);
  }
  if (argc==3) isMassFac=true;
  filename = string(argv[1]);

  TFile *inFile = TFile::Open(filename.c_str());

  string outfilename;
  if (isMassFac) outfilename="HggSpin_MassFac.root";
  else outfilename="HggSpin_CutBased.root";

  TFile *outFile = new TFile(outfilename.c_str(),"RECREATE");
  RooWorkspace *work = new RooWorkspace("HggSpinStudies");

  RooRealVar *mass = new RooRealVar("mass","mass",100,180);
  RooRealVar *weight = new RooRealVar("weight","weight",0,100);

  TTree *dataTree = (TTree*)inFile->Get("spin_trees/Data");
  TTree *spin0Tree = (TTree*)inFile->Get("spin_trees/spin0plus_m125_8TeV");
  TTree *spin2Tree = (TTree*)inFile->Get("spin_trees/grav2pm_m125_8TeV");

  map<string,RooDataSet*> dataDSet;
  map<string,RooDataSet*> spin0DSet;
  map<string,RooDataSet*> spin2DSet;

  int nBDTCats;
  int nSpinCats;
  if (isMassFac) nBDTCats=2;
  else nBDTCats=4;
  nSpinCats=5;

  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      string catname = getCatName(c,s);
      dataDSet.insert(pair<string,RooDataSet*>(catname,new RooDataSet(Form("data_mass_%s",catname.c_str()),Form("data_mass_%s",catname.c_str()),RooArgSet(*mass,*weight),weight->GetName())));
      spin0DSet.insert(pair<string,RooDataSet*>(catname,new RooDataSet(Form("mcSigSMHiggs_bdt%d_cTh%d",c,s),Form("mcSigSMHiggs_bdt%d_cTh%d",c,s),RooArgSet(*mass,*weight),weight->GetName())));
      spin2DSet.insert(pair<string,RooDataSet*>(catname,new RooDataSet(Form("mcSigGraviton_bdt%d_cTh%d",c,s),Form("mcSigGraviton_bdt%d_cTh%d",c,s),RooArgSet(*mass,*weight),weight->GetName())));
    }
  }

  fillDataSet(mass,dataDSet,dataTree,isMassFac);
  fillDataSet(mass,spin0DSet,spin0Tree,isMassFac);
  fillDataSet(mass,spin2DSet,spin2Tree,isMassFac);
  
  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      string catname = getCatName(c,s);
      cout << "c" << c << " s" << s << " Data: " << dataDSet[catname]->sumEntries() << " Spin0: " << spin0DSet[catname]->sumEntries() << " Spin2: " << spin2DSet[catname]->sumEntries() << endl;
      //work->import(*dataDSet[catname]);
      //work->import(*spin0DSet[catname]);
      //work->import(*spin2DSet[catname]);
    }
  }

  string dirname;
  if (isMassFac) dirname="massFac";
  else dirname="cutBased";
  makeBackgroundModel(dirname,mass,dataDSet,work,nBDTCats,nSpinCats);
  makeSignalModel(dirname,mass,spin0DSet,work,nBDTCats,nSpinCats,true);
  makeSignalModel(dirname,mass,spin2DSet,work,nBDTCats,nSpinCats,false);

  for (map<string,RooDataSet*>::iterator it=dataDSet.begin(); it!=dataDSet.end(); it++) delete it->second;
  for (map<string,RooDataSet*>::iterator it=spin0DSet.begin(); it!=spin0DSet.end(); it++) delete it->second;
  for (map<string,RooDataSet*>::iterator it=spin2DSet.begin(); it!=spin2DSet.end(); it++) delete it->second;
  
  inFile->Close();
  outFile->cd();
  work->Write();
  outFile->Close();

  delete inFile;
  delete work;
  delete outFile;
  
  return 0;
}
 
