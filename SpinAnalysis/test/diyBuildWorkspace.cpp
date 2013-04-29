#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <fstream>

#include "boost/lexical_cast.hpp"

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
#include "RooGenericPdf.h"
#include "RooRandom.h"

using namespace std;
using namespace RooFit;

int GetCategory(std::vector<double> & v, float val)
{
    if( val == v[0] ) { return 0; }
    std::vector<double>::iterator bound =  lower_bound( v.begin(), v.end(), val, std::greater<double>  ());
    int cat = ( val >= *bound ? bound - v.begin() - 1 : bound - v.begin() );
    if( cat >= int(v.size()) - 1 ) { cat = -1; }
    return cat;
}

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

pair<int,int> getCutBasedCategory(double lead_eta, double sublead_eta, float lead_r9, float sublead_r9, double costheta_cs, std::vector<double> cosThetaBoundaries){

  // marco style
  int mycat=-1;
  int costheta_csBin=-1;
  /*if (TMath::Abs(costheta_cs)>=0. && TMath::Abs(costheta_cs)<0.2) costheta_csBin=0;
  else if (TMath::Abs(costheta_cs)>=0.2 && TMath::Abs(costheta_cs)<0.4) costheta_csBin=1;
  else if (TMath::Abs(costheta_cs)>=0.4 && TMath::Abs(costheta_cs)<0.6) costheta_csBin=2;
  else if (TMath::Abs(costheta_cs)>=0.6 && TMath::Abs(costheta_cs)<0.8) costheta_csBin=3;
  else if (TMath::Abs(costheta_cs)>=0.8 && TMath::Abs(costheta_cs)<=1.0) costheta_csBin=4;
  else costheta_csBin=-1;*/
  costheta_csBin=GetCategory(cosThetaBoundaries, 1-TMath::Abs(costheta_cs));

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

void fillDataSet(RooRealVar *mass, map<string,RooDataSet*> &dataMap, vector<TTree*> trees, bool isMassFac, double mlow, double mhigh, std::vector<double> cosThetaBoundaries){

  int category;
  float evweight;
  double costheta_cs;
  double lead_calo_eta;
  double sublead_calo_eta;
  float lead_r9;
  float sublead_r9;
  float diphoton_bdt;
  double higgs_mass;

  for (vector<TTree*>::iterator treeIt=trees.begin(); treeIt!=trees.end(); treeIt++){
    TTree *tree = *treeIt;
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
      else mycat = getCutBasedCategory(lead_calo_eta,sublead_calo_eta,lead_r9,sublead_r9,costheta_cs,cosThetaBoundaries);
      if (mycat.first==-1 || mycat.second==-1) continue;
      if (higgs_mass<mlow || higgs_mass>mhigh) continue;
      mass->setVal(higgs_mass);
      string catname = getCatName(mycat);
      dataMap[catname]->add(*mass,evweight);
    }
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


void makeBackgroundModel(string name, RooRealVar *mass, map<string,RooDataSet*> data, RooWorkspace *work, int nBDTCategories, int nCosThetaCategories, bool globePDFs, bool verbose=false, bool correlateCosThetaCategories=false){

  string path="plots/"+name;
  system(Form("mkdir -p %s",path.c_str()));


  for(int a=0; a<nBDTCategories; a++){

    RooRealVar *p0, *p1, *p2, *p3, *pow0;
    RooFormulaVar *f0, *f1, *f2, *f3;
    RooAbsPdf *bkg;

    if(correlateCosThetaCategories)
    {
      if(globePDFs)
      {
        p0 = new RooRealVar(Form("CMS_hgg_quartic0_cat%d",a),Form("CMS_hgg_quartic0_cat%d",a),0.02,-5.0,5.0);
        p1 = new RooRealVar(Form("CMS_hgg_quartic1_cat%d",a),Form("CMS_hgg_quartic0_cat%d",a),0.02,-5.0,5.0);
        p2 = new RooRealVar(Form("CMS_hgg_quartic2_cat%d",a),Form("CMS_hgg_quartic0_cat%d",a),0.02,-5.0,5.0);
        p3 = new RooRealVar(Form("CMS_hgg_quartic3_cat%d",a),Form("CMS_hgg_quartic0_cat%d",a),0.02,-5.0,5.0);

        f0 = new RooFormulaVar(Form("CMS_hgg_modquartic0_cat%d",a),Form("CMS_hgg_modquartic0_cat%d",a),"@0*@0",RooArgList(*p0));
        f1 = new RooFormulaVar(Form("CMS_hgg_modquartic1_cat%d",a),Form("CMS_hgg_modquartic0_cat%d",a),"@0*@0",RooArgList(*p1));
        f2 = new RooFormulaVar(Form("CMS_hgg_modquartic2_cat%d",a),Form("CMS_hgg_modquartic0_cat%d",a),"@0*@0",RooArgList(*p2));
        f3 = new RooFormulaVar(Form("CMS_hgg_modquartic3_cat%d",a),Form("CMS_hgg_modquartic0_cat%d",a),"@0*@0",RooArgList(*p3));

        bkg = new RooBernstein(Form("data_pol_model_cat%d",a),Form("data_pol_model_cat%d",a),*mass,RooArgList(*f0,*f1,*f2,*f3));
      }
      else
      {
        pow0 = new RooRealVar(Form("CMS_hgg_pow0_cat%d",a),Form("CMS_hgg_pow0_cat%d",a),10.,2.,20.);
        bkg = new RooGenericPdf(Form("data_pow_model_cat%d",a),Form("data_pow_model_cat%d",a),"pow(@0,-@1)",RooArgList(*mass,*pow0));
      }
    }

    for(int b=0; b<nCosThetaCategories; b++){

      string catname = getCatName(a,b);

      if(!correlateCosThetaCategories)
      {
        if (globePDFs) {
          p0 = new RooRealVar(Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);
          p1 = new RooRealVar(Form("CMS_hgg_quartic1_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);
          p2 = new RooRealVar(Form("CMS_hgg_quartic2_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);
          p3 = new RooRealVar(Form("CMS_hgg_quartic3_cat%d_spin%d",a,b),Form("CMS_hgg_quartic0_cat%d_spin%d",a,b),0.02,-5.0,5.0);

          f0 = new RooFormulaVar(Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p0));
          f1 = new RooFormulaVar(Form("CMS_hgg_modquartic1_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p1));
          f2 = new RooFormulaVar(Form("CMS_hgg_modquartic2_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p2));
          f3 = new RooFormulaVar(Form("CMS_hgg_modquartic3_cat%d_spin%d",a,b),Form("CMS_hgg_modquartic0_cat%d_spin%d",a,b),"@0*@0",RooArgList(*p3));

          if (b==nCosThetaCategories-1){
            bkg = new RooBernstein(Form("data_pol_model_cat%d_spin%d",a,b),Form("data_pol_model_cat%d_spin%d",a,b),*mass,RooArgList(*f0,*f1,*f2));
          }
          else {
            bkg = new RooBernstein(Form("data_pol_model_cat%d_spin%d",a,b),Form("data_pol_model_cat%d_spin%d",a,b),*mass,RooArgList(*f0,*f1,*f2,*f3));
          }
        }
        else {
          pow0 = new RooRealVar(Form("CMS_hgg_pow0_cat%d_spin%d",a,b),Form("CMS_hgg_pow0_cat%d_spin%d",a,b),10.,2.,20.);
          //RooRealVar *pow1 = new RooRealVar(Form("CMS_hgg_pow1_cat%d_spin%d",a,b),Form("CMS_hgg_pow1_cat%d_spin%d",a,b),10.,2.,20.);
          //RooRealVar *f1 = new RooRealVar(Form("CMS_hgg_f1_cat%d_spin%d",a,b),Form("CMS_hgg_f1_cat%d_spin%d",a,b),0.6,0.5,1.);
          //bkg = new RooGenericPdf(Form("data_pow_model_cat%d_spin%d",a,b),Form("data_pow_model_cat%d_spin%d",a,b),"@3*pow(@0,-@1)+(1.-@3)*pow(@01,-@2)",RooArgList(*mass,*pow0,*pow1,*f1));
          bkg = new RooGenericPdf(Form("data_pow_model_cat%d_spin%d",a,b),Form("data_pow_model_cat%d_spin%d",a,b),"pow(@0,-@1)",RooArgList(*mass,*pow0));
        }
      }

      RooRealVar *norm = new RooRealVar(Form("%s_norm",bkg->GetName()),Form("%s_norm",bkg->GetName()),10.,1.e6);
      RooExtendPdf *ext;
      if (globePDFs) ext = new RooExtendPdf(Form("pdf_data_pol_model_cat%d_spin%d",a,b),Form("pdf_data_pol_model_cat%d_spin%d",a,b),*bkg,*norm);
      else ext = new RooExtendPdf(Form("pdf_data_pow_model_cat%d_spin%d",a,b),Form("pdf_data_pow_model_cat%d_spin%d",a,b),*bkg,*norm);

      cout << data[catname]->GetName() << " : " << data[catname]->sumEntries() << endl;

      bkg->fitTo(*data[catname]);//,PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));

      /*
      RooFitResult *fitRes;
      verbose ?
        fitRes = ext->fitTo(*data[catname],Save(true)) :
        fitRes = ext->fitTo(*data[catname],Save(true),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1))
      ;
      fitRes->floatParsFinal().Print("v");
      delete fitRes;
      */

      if (globePDFs) mass->setBins(320);
      else mass->setBins(200);
      RooDataHist *temp = data[catname]->binnedClone(Form("roohist_%s",data[catname]->GetName()));
      work->import(*temp);
      work->import(*data[catname]);
      work->import(*ext, RooFit::RecycleConflictNodes());

      Plot(Form("%s/bkg_cat%d_spin%d",path.c_str(),a,b),mass,data[catname],bkg);
    }
  }
}

void makeSignalModel(string name, RooRealVar *mass, map<string,RooDataSet*> sigMC, RooWorkspace *work, int nBDTCategories, int nCosThetaCategories, bool globePDFs, bool isSM, bool verbose=false){

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
      int nGaus=4;
      if (globePDFs) nGaus=4;
      for (int i=0; i<nGaus; i++){
        string nameMean     = Form("mean%d_%s_cat%d_spin%d"         ,i,postfix.c_str(),a,b);
        string nameSigma    = Form("sigma%d_%s_cat%d_spin%d"        ,i,postfix.c_str(),a,b);
        string nameGauss    = Form("gauss%d_%s_cat%d_spin%d"        ,i,postfix.c_str(),a,b);
        means.push_back(new RooRealVar(nameMean.c_str(),nameMean.c_str(),125.,115.,135.));
        sigmas.push_back(new RooRealVar(nameSigma   .c_str(),nameSigma   .c_str(),  2., 0.5,  20.));
        pdfGaussians.push_back(new RooGaussian(nameGauss.c_str(),nameGauss.c_str(),*mass,*means[i],*sigmas[i]));
        if (i<nGaus-1){
          string nameFrac     = Form("frac%d_%s_cat%d_spin%d"         ,i,postfix.c_str(),a,b);
          fracs.push_back(new RooRealVar(nameFrac   .c_str(),nameFrac   .c_str(),  0.01, 0.,  1.));
        }
      }

      if (globePDFs) sigPdf = new RooAddPdf(Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b),Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b), RooArgList(*pdfGaussians[0],*pdfGaussians[1],*pdfGaussians[2],*pdfGaussians[3]),RooArgList(*fracs[0],*fracs[1],*fracs[2]));
      else sigPdf = new RooAddPdf(Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b),Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b), RooArgList(*pdfGaussians[0],*pdfGaussians[1],*pdfGaussians[2],*pdfGaussians[3]),RooArgList(*fracs[0],*fracs[1],*fracs[2]));
      //else sigPdf = new RooAddPdf(Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b),Form("sigModel_%s_cat%d_spin%d",postfix.c_str(),a,b), RooArgList(*pdfGaussians[0],*pdfGaussians[1]),RooArgList(*fracs[0]));

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

  string filename="";
  string outfilename="";
  bool isMassFac=false;
  bool globePDFs=false;
  bool fullSMproc=false;
  bool useSMpowheg=false;
  bool useSpin2LP=false;
  int nBDTCats=0;
  int nSpinCats=0;
  std::vector<double> cosThetaBoundaries;
  string boundaries="";
  bool correlateCosThetaCategories = false;
  bool useBackgroundMC = false;

  if (argc!=2){
    cout << "usage ./bin/diyBuildWorkspace <datfilename>" << endl;
    exit(1);
  }
  else {
    ifstream datfile;
    datfile.open(argv[1]);
    if (datfile.fail()){
      cout << "datfile " << argv[1] << " doesn't exist" << endl;
      exit(1);
    }
    while (datfile.good()){
      string line;
      getline(datfile,line);
      if (line=="\n" || line.substr(0,1)=="#" || line==" " || line.empty()) continue;
      if (line.find("treefile=")!=string::npos) filename = line.substr(line.find("=")+1,string::npos);
      if (line.find("wsfile=")!=string::npos) outfilename = line.substr(line.find("=")+1,string::npos);
      if (line.find("isMassFac=")!=string::npos) isMassFac = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
      if (line.find("globePDFs=")!=string::npos) globePDFs = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
      if (line.find("fullSMproc=")!=string::npos) fullSMproc = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
      if (line.find("useSMpowheg=")!=string::npos) useSMpowheg = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
      if (line.find("useSpin2LP=")!=string::npos) useSpin2LP = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
      if (line.find("nBDTCats=")!=string::npos) nBDTCats = boost::lexical_cast<int>(line.substr(line.find("=")+1,string::npos));
      if (line.find("nSpinCats=")!=string::npos) nSpinCats = boost::lexical_cast<int>(line.substr(line.find("=")+1,string::npos));
      if (line.find("catBoundaries=")!=string::npos) boundaries = line.substr(line.find("=")+1,string::npos);
      if (line.find("correlateCosThetaCategories=")!=string::npos) correlateCosThetaCategories = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
      if (line.find("useBackgroundMC=")!=string::npos) useBackgroundMC = boost::lexical_cast<bool>(line.substr(line.find("=")+1,string::npos));
    }
    datfile.close();

    if(boundaries != "")
    {
      //cout << "Testing boundaries: " << boundaries << endl;

      size_t position = boundaries.find(" ");
      do
      {
        position = boundaries.find(" ");
        string testnum = boundaries.substr(0,position);
        boundaries = boundaries.substr(position+1);

        //cout << testnum << endl;
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

  cout << "Reading from file - " << filename << endl;
  if (isMassFac) cout << "Running MassFac selection" << endl;
  else cout << "Running CutBased selection" << endl;

  TFile *inFile = TFile::Open(filename.c_str());

  TFile *outFile = new TFile(outfilename.c_str(),"RECREATE");
  RooWorkspace *work = new RooWorkspace("HggSpinStudies");

  double mlow=100.;
  double mhigh=150.;
  if (globePDFs) mhigh=180.;

  RooRealVar *mass = new RooRealVar("mass","mass",mlow,mhigh);
  RooRealVar *weight = new RooRealVar("weight","weight",0,100);

  vector<TTree*> dataTrees;
  vector<TTree*> spin0Trees;
  vector<TTree*> spin2Trees;

  if(!useBackgroundMC)
    dataTrees.push_back((TTree*)inFile->Get("spin_trees/Data"));
  else
  {
    dataTrees.push_back((TTree*)inFile->Get("spin_trees/gjet_20_8TeV_pf"));
    dataTrees.push_back((TTree*)inFile->Get("spin_trees/gjet_40_8TeV_pf"));
    dataTrees.push_back((TTree*)inFile->Get("spin_trees/gjet_20_8TeV_pp"));
    dataTrees.push_back((TTree*)inFile->Get("spin_trees/gjet_40_8TeV_pp"));
    dataTrees.push_back((TTree*)inFile->Get("spin_trees/diphojet_8TeV"));
    dataTrees.push_back((TTree*)inFile->Get("spin_trees/dipho_Box_25_8TeV"));
  }

  if (useSpin2LP) spin2Trees.push_back((TTree*)inFile->Get("spin_trees/spin2lp_m125_8TeV"));
  else spin2Trees.push_back((TTree*)inFile->Get("spin_trees/ggh_grav_m125_8TeV"));
  if (useSMpowheg){
    spin0Trees.push_back((TTree*)inFile->Get("spin_trees/ggh_m125_8TeV"));
    if (fullSMproc){
      spin0Trees.push_back((TTree*)inFile->Get("spin_trees/vbf_m125_8TeV"));
      spin0Trees.push_back((TTree*)inFile->Get("spin_trees/wzh_m125_8TeV"));
      spin0Trees.push_back((TTree*)inFile->Get("spin_trees/tth_m125_8TeV"));
    }
  }
  else {
    spin0Trees.push_back((TTree*)inFile->Get("spin_trees/spin0plus_m125_8TeV"));
  }

  map<string,RooDataSet*> dataDSet;
  map<string,RooDataSet*> spin0DSet;
  map<string,RooDataSet*> spin2DSet;

  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      string catname = getCatName(c,s);
      dataDSet.insert(pair<string,RooDataSet*>(catname,new RooDataSet(Form("data_mass_%s",catname.c_str()),Form("data_mass_%s",catname.c_str()),RooArgSet(*mass,*weight),weight->GetName())));
      spin0DSet.insert(pair<string,RooDataSet*>(catname,new RooDataSet(Form("mcSigSMHiggs_bdt%d_cTh%d",c,s),Form("mcSigSMHiggs_bdt%d_cTh%d",c,s),RooArgSet(*mass,*weight),weight->GetName())));
      spin2DSet.insert(pair<string,RooDataSet*>(catname,new RooDataSet(Form("mcSigGraviton_bdt%d_cTh%d",c,s),Form("mcSigGraviton_bdt%d_cTh%d",c,s),RooArgSet(*mass,*weight),weight->GetName())));
    }
  }
  cout << "BDTCats: " << nBDTCats << endl << "SpinCats: " << nSpinCats << endl;
  fillDataSet(mass,dataDSet,dataTrees,isMassFac,mlow,mhigh,cosThetaBoundaries);
  fillDataSet(mass,spin0DSet,spin0Trees,isMassFac,mlow,mhigh,cosThetaBoundaries);
  fillDataSet(mass,spin2DSet,spin2Trees,isMassFac,mlow,mhigh,cosThetaBoundaries);

  double totalSM=0.;
  double totalGRAV=0.;
  for (int c=0; c<nBDTCats; c++){
    for (int s=0; s<nSpinCats; s++){
      string catname = getCatName(c,s);
      totalSM += spin0DSet[catname]->sumEntries();
      totalGRAV += spin2DSet[catname]->sumEntries();
      cout << "c" << c << " s" << s << " Data: " << dataDSet[catname]->sumEntries() << " Spin0: " << spin0DSet[catname]->sumEntries() << " Spin2: " << spin2DSet[catname]->sumEntries() << endl;
      //work->import(*dataDSet[catname]);
      //work->import(*spin0DSet[catname]);
      //work->import(*spin2DSet[catname]);
    }
  }
  cout << "TOTAL-- sm: " << totalSM << " grav: " << totalGRAV << endl;
  cout << "SCALE FACTOR = " << totalSM/totalGRAV << endl;

  string dirname;
  if (isMassFac) dirname="massFac";
  else dirname="cutBased";
  makeBackgroundModel(dirname,mass,dataDSet,work,nBDTCats,nSpinCats,globePDFs,false,correlateCosThetaCategories);
  makeSignalModel(dirname,mass,spin0DSet,work,nBDTCats,nSpinCats,globePDFs,true);
  makeSignalModel(dirname,mass,spin2DSet,work,nBDTCats,nSpinCats,globePDFs,false);

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

