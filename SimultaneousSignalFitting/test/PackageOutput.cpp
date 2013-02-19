#include <iostream>
#include <vector>
#include <string>

#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"

#include "TFile.h"
#include "TROOT.h"

#include "RooWorkspace.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooHistFunc.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = boost::program_options;

string infilename_;
string outfilename_;
string wsname_;
int mhLow_=110;
int mhHigh_=150;
int ncats_;
string webdir_;
bool web_;
bool makePlots_=false;
vector<int> allMH_;

vector<int> getAllMH(){
  vector<int> result;
  for (int m=mhLow_; m<=mhHigh_; m+=5){
    cout << "Adding mass: " << m << endl;
    result.push_back(m);
  }
  return result;
}

void OptionParser(int argc, char *argv[]){

  po::options_description desc("Allowed options");

  desc.add_options()
    ("help,h",                                                                                "Show help")
    ("infilename,i", po::value<string>(&infilename_)->default_value("dat/filestocombine.dat"),"Input file name")
    ("outfilename,o", po::value<string>(&outfilename_)->default_value("CMS-HGG_sigfit.root"), "Output file name")
    ("wsname,W", po::value<string>(&wsname_)->default_value("wsig_8TeV"),                     "Output workspace name")
    ("mhLow,L", po::value<int>(&mhLow_)->default_value(110),                                  "Low mass point")
    ("mhHigh,H", po::value<int>(&mhHigh_)->default_value(150),                                "High mass point")
    ("ncats,n", po::value<int>(&ncats_)->default_value(9),                                    "Number of categories")
    ("html,w", po::value<string>(&webdir_),                                                   "Make html in this directory")
    ("makePlots,P",                                                                           "Make AN style signal model plots")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")){ cout << desc << endl; exit(1);}
  if (vm.count("html"))             web_=true;
  if (vm.count("makePlots"))        makePlots_=true;
  allMH_ = getAllMH();
}
 
int main (int argc, char *argv[]){
  
  OptionParser(argc,argv);
  
  vector<string> processes;
  processes.push_back("ggh");
  processes.push_back("vbf");
  processes.push_back("wzh");
  processes.push_back("tth");
  
  map<string,pair<string,int> > filestocombine;

  ifstream datfile;
  datfile.open(infilename_.c_str());
  if (datfile.fail()) {
    cerr << "datfile: " << infilename_ << " not found. Bailing out" << endl;
    exit(1);
  }
  while (datfile.good()){
    string line;
    getline(datfile,line);
    if (line=="\n" || line.substr(0,1)=="#" || line==" " || line.empty()) continue;
    string fname = line.substr(0,line.find(" "));
    string info = line.substr(line.find(" ")+1,string::npos);
    string proc = info.substr(0,info.find("_cat"));
    string cat_string = info.substr(info.find("_cat")+4,string::npos);
    cout << fname << " : " << info << " : " << proc << " : " << cat_string << endl;
    int cat = boost::lexical_cast<int>(cat_string);
    filestocombine.insert(pair<string,pair<string,int> >(fname,pair<string,int>(proc,cat)));
  }
  datfile.close();

  TFile *outFile = new TFile(outfilename_.c_str(),"RECREATE");

  RooRealVar *intLumi = new RooRealVar("IntLumi","IntLumi",19620.,0.,10.e5);
  // first loop files and import all pdfs and dataset into one workspace
  RooWorkspace *work = NULL;
  for (map<string,pair<string,int> >::iterator file=filestocombine.begin(); file!=filestocombine.end(); file++){
    string filename = file->first;
    string proc = file->second.first;
    int cat = file->second.second;

    TFile *thisFile = TFile::Open(file->first.c_str());
    
    if (file==filestocombine.begin()) {
      work = (RooWorkspace*)thisFile->Get("wsig");
      work->SetName(wsname_.c_str());
    }
    else {
      RooWorkspace *tempWS = (RooWorkspace*)thisFile->Get("wsig");
      RooAddPdf *resultPdf = (RooAddPdf*)tempWS->pdf(Form("sigpdfrel_%s_cat%d",proc.c_str(),cat));
      work->import(*resultPdf);
      for (unsigned int i=0; i<allMH_.size(); i++){
        int m = allMH_[i];
        RooDataSet *resultData = (RooDataSet*)tempWS->data(Form("sig_%s_mass_m%d_cat%d",proc.c_str(),m,cat));
        work->import(*resultData);
      }
    }
  }
 
  vector<string> expectedObjectsNotFound;

  // now we want to sum up some datasets and their respective pdfs for sig_eff calculations later
  for (unsigned int i=0; i<allMH_.size(); i++){
    int m = allMH_[i];
    RooDataSet *allDataThisMass = NULL;
    for (int cat=0; cat<ncats_; cat++){
      RooDataSet *allDataThisCat = NULL;
      for (vector<string>::iterator proc=processes.begin(); proc!=processes.end(); proc++){
        RooDataSet *tempData = (RooDataSet*)work->data(Form("sig_%s_mass_m%d_cat%d",proc->c_str(),m,cat));
        if (!tempData) {
          cerr << "WARNING -- dataset: " << Form("sig_%s_mass_m%d_cat%d",proc->c_str(),m,cat) << " not found. It will be skipped" << endl;
          expectedObjectsNotFound.push_back(Form("sig_%s_mass_m%d_cat%d",proc->c_str(),m,cat));
          continue;
        }
        if (cat==0) allDataThisMass = (RooDataSet*)tempData->Clone(Form("sig_mass_m%d_AllCats",m));
        else allDataThisMass->append(*tempData);
        if (*proc=="ggh") allDataThisCat = (RooDataSet*)tempData->Clone(Form("sig_mass_m%d_cat%d",m,cat));
        else allDataThisCat->append(*tempData);
      }
      if (!allDataThisCat) {
        cerr << "WARNING -- allData for cat " << cat << " is NULL. Probably because the relevant datasets couldn't be found. Skipping.. " << endl;
        continue;
      }
      work->import(*allDataThisCat);
    }
    if (!allDataThisMass) {
      cerr << "WARNING -- allData for mass " << m << " is NULL. Probably because the relevant datasets couldn't be found. Skipping.. " << endl;
      continue;
    }
    work->import(*allDataThisMass);
  }
  RooArgList *sumPdfs = new RooArgList();
  for (int cat=0; cat<ncats_; cat++){
    RooArgList *sumPdfsThisCat = new RooArgList();
    for (vector<string>::iterator proc=processes.begin(); proc!=processes.end(); proc++){
      //RooExtendPdf *tempPdf = (RooExtendPdf*)work->pdf(Form("sigpdfrel_%s_cat%d",proc->c_str(),cat));
      RooHistFunc *norm = (RooHistFunc*)work->function(Form("hggpdfrel_%s_cat%d_norm",proc->c_str(),cat));
      RooAddPdf *pdf = (RooAddPdf*)work->pdf(Form("hggpdfrel_%s_cat%d",proc->c_str(),cat));
      RooFormulaVar *thisLumNorm = new RooFormulaVar(Form("hggpdfabs_%s_cat%d_norm",proc->c_str(),cat),Form("hggpdfabs_%s_cat%d_norm",proc->c_str(),cat),"@*@1",RooArgList(*norm,*intLumi));
      if (!norm && !pdf) cout << "AHHHH" << endl;
      RooExtendPdf *tempPdf = new RooExtendPdf(Form("sigpdfabs_%s_cat%d",proc->c_str(),cat),Form("sigpdfabs_%s_cat%d",proc->c_str(),cat),*pdf,*thisLumNorm);
      if (!tempPdf) {
        cerr << "WARNING -- pdf: " << Form("sigpdfrel_%s_cat%d",proc->c_str(),cat) << " not found. It will be skipped" << endl;
        expectedObjectsNotFound.push_back(Form("sigpdfrel_%s_cat%d",proc->c_str(),cat));
        continue;
      }
      sumPdfsThisCat->add(*tempPdf);
    }
    if (sumPdfsThisCat->getSize()==0){
        cerr << "WARNING -- sumPdfs for cat " << cat << " is EMPTY. Probably because the relevant pdfs couldn't be found. Skipping.. " << endl;
        continue;
    }
    RooAddPdf *sumPdfsPerCat = new RooAddPdf(Form("sigpdfrelcat%d_allProcs",cat),Form("sigpdfrelcat%d_allProcs",cat),*sumPdfsThisCat);
    sumPdfs->add(*sumPdfsPerCat);
  }
  if (sumPdfs->getSize()==0){
      cerr << "WARNING -- sumAllPdfs is EMPTY. Probably because the relevant pdfs couldn't be found. Skipping.. " << endl;
  }
  else {
    RooAddPdf *sumPdfsAllCats = new RooAddPdf("sigpdfrelAllCats_allProcs","sigpdfrelAllCats_allProcs",*sumPdfs);
    work->import(*sumPdfsAllCats);
  }

  outFile->cd();
  work->Write();
  outFile->Close();
  delete outFile;

  if (expectedObjectsNotFound.size()==0) cout << "All expected objects found and packaged successfully." << endl;
  else {
    cout << "Some expected objects not found:" << endl;
    for (vector<string>::iterator it=expectedObjectsNotFound.begin(); it!=expectedObjectsNotFound.end(); it++) cout << "\t" << *it << endl;
  }
  cout << "Output written to: " << endl;
  cout << "\t" << outfilename_ << endl;

  if (web_){
    string sitename = webdir_.substr(webdir_.find("www")+4,string::npos);
    sitename = "www.cern.ch/mkenzie/"+sitename;
    cout << "Publishing to web directory " << webdir_ << endl;
    system(Form("rm -rf %s",webdir_.c_str()));
    system(Form("cp -r plots %s",webdir_.c_str()));
    system(Form("make_html.py %s -c -l -n --title \"Simultaneous Signal Fitting\"",webdir_.c_str()));
    cout << "\t" << sitename << endl;
  }

  if (makePlots_){
    ifstream tempfile("../Macros/makeParametricSignalModelPlots.C");
    if (!tempfile) {
      cout << "I'm looking for a file: ../Macros/makeParametricSignalModelPlots.C to make the AN plots but I can't find it. Bailing!" << endl;
      exit(1);
    }
    gROOT->ProcessLine(".L ../Macros/makeParametricSignalModelPlots.C+g");
  }

  return 0;
}  

