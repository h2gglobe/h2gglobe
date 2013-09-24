#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <typeinfo>

#include "TFile.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "TKey.h"
#include "TMacro.h"

#include "../interface/InitialFit.h"
#include "../interface/LinearInterp.h"
#include "../interface/FinalModelConstruction.h"
#include "../interface/Packager.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = boost::program_options;

string filename_;
string outfilename_;
string datfilename_;
string systfilename_;
int mhLow_=110;
int mhHigh_=150;
int nCats_;
float constraintValue_;
int constraintValueMass_;
bool spin_=false;
bool splitVH_=false;
bool splitRVWV_=true;
bool doSecondaryModels_=true;
bool recursive_=false;
int verbose_=0;

void OptionParser(int argc, char *argv[]){
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                			"Show help")
    ("infilename,i", po::value<string>(&filename_),                                           			"Input file name")
    ("outfilename,o", po::value<string>(&outfilename_)->default_value("CMS-HGG_sigfit.root"), 			"Output file name")
    ("datfilename,d", po::value<string>(&datfilename_)->default_value("dat/config.dat"),      			"Configuration file")
		("systfilename,s", po::value<string>(&systfilename_)->default_value("dat/photonCatSyst.dat"),		"Systematic model numbers")
    ("mhLow,L", po::value<int>(&mhLow_)->default_value(110),                                  			"Low mass point")
    ("mhHigh,H", po::value<int>(&mhHigh_)->default_value(150),                                			"High mass point")
    ("nCats,n", po::value<int>(&nCats_)->default_value(9),                                    			"Number of total categories")
    ("constraintValue,C", po::value<float>(&constraintValue_)->default_value(0.1),            			"Constraint value")
    ("constraintValueMass,M", po::value<int>(&constraintValueMass_)->default_value(125),            "Constraint value mass")
    ("skipSecondaryModels",                                                                   			"Turn off creation of all additional models")
    ("splitVH",                                                                               			"Split VH into WH and ZH")
    ("recursive",                                                                             			"Recursively calculate gaussian fractions")
    ("verbose,v", po::value<int>(&verbose_)->default_value(0),                                			"Verbosity level: 0 (lowest) - 3 (highest)")
  ;                                                                                             		
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")){ cout << desc << endl; exit(1);}
  if (vm.count("spin"))                     spin_=true;
  if (vm.count("splitVH"))                  splitVH_=true;
  if (vm.count("nosplitRVWV"))              splitRVWV_=false;
  if (vm.count("skipSecondaryModels"))      doSecondaryModels_=false;
  if (vm.count("splitVH"))                  splitVH_=true;
  if (vm.count("recursive"))                recursive_=true;
}

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

int main(int argc, char *argv[]){
 
  OptionParser(argc,argv);

  TStopwatch sw;
  sw.Start();

  TFile *inFile = TFile::Open(filename_.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->SetTitle("m_{#gamma#gamma}");
  mass->setUnit("GeV");
  RooRealVar *intLumi = (RooRealVar*)inWS->var("IntLumi");
  RooRealVar *MH = new RooRealVar("MH","m_{H}",mhLow_,mhHigh_);
  MH->setUnit("GeV");

  RooRealVar *MH_SM = new RooRealVar("MH_SM","m_{H} (SM)",mhLow_,mhHigh_);
  RooRealVar *DeltaM = new RooRealVar("DeltaM","#Delta m_{H}",0.,-10.,10.);
  DeltaM->setUnit("GeV");
  RooAddition *MH_2 = new RooAddition("MH_2","m_{H} (2)",RooArgList(*MH,*DeltaM));
  RooRealVar *higgsDecayWidth = new RooRealVar("HiggsDecayWidth","#Gamma m_{H}",0.,0.,10.);
 
  TFile *outFile = new TFile(outfilename_.c_str(),"RECREATE");
  RooWorkspace *outWS = new RooWorkspace("wsig_8TeV");

  transferMacros(inFile,outFile);

  // run fits for each line in datfile
  ifstream datfile;
  datfile.open(datfilename_.c_str());
  if (datfile.fail()) exit(1);
  while (datfile.good()){
    string line;
    getline(datfile,line);
    if (line=="\n" || line.substr(0,1)=="#" || line==" " || line.empty()) continue;
    vector<string> els;
    split(els,line,boost::is_any_of(" "));
    assert(els.size()==4);
    string proc = els[0];
    int cat = boost::lexical_cast<int>(els[1]);
    int nGaussiansRV = boost::lexical_cast<int>(els[2]);
    int nGaussiansWV = boost::lexical_cast<int>(els[3]);

    cout << "-----------------------------------------------------------------" << endl;
    cout << Form("Running fits for proc:%s - cat:%d with nGausRV:%d nGausWV:%d",proc.c_str(),cat,nGaussiansRV,nGaussiansWV) << endl;
    cout << "-----------------------------------------------------------------" << endl;
    // get datasets for each MH here
    map<int,RooDataSet*> datasetsRV;
    map<int,RooDataSet*> datasetsWV;
    map<int,RooDataSet*> datasets;

    for (int mh=mhLow_; mh<=mhHigh_; mh+=5){
      RooDataSet *dataRV = (RooDataSet*)inWS->data(Form("sig_%s_mass_m%d_rv_cat%d",proc.c_str(),mh,cat));
      RooDataSet *dataWV = (RooDataSet*)inWS->data(Form("sig_%s_mass_m%d_wv_cat%d",proc.c_str(),mh,cat));
      RooDataSet *data = (RooDataSet*)inWS->data(Form("sig_%s_mass_m%d_cat%d",proc.c_str(),mh,cat));
      datasetsRV.insert(pair<int,RooDataSet*>(mh,dataRV));
      datasetsWV.insert(pair<int,RooDataSet*>(mh,dataWV));
      datasets.insert(pair<int,RooDataSet*>(mh,data));
    }

    // these guys do the fitting
    // right vertex
    InitialFit initFitRV(mass,MH,mhLow_,mhHigh_);
    initFitRV.setVerbosity(verbose_);
    initFitRV.buildSumOfGaussians(Form("%s_cat%d",proc.c_str(),cat),nGaussiansRV);
    initFitRV.setDatasets(datasetsRV);
    initFitRV.runFits(1);
    initFitRV.saveParamsToFileAtMH(Form("dat/in/%s_cat%d_rv.dat",proc.c_str(),cat),constraintValueMass_);
    initFitRV.loadPriorConstraints(Form("dat/in/%s_cat%d_rv.dat",proc.c_str(),cat),constraintValue_);
    initFitRV.runFits(1);
    initFitRV.plotFits(Form("plots/%s_cat%d/rv",proc.c_str(),cat));
    map<int,map<string,RooRealVar*> > fitParamsRV = initFitRV.getFitParams();
    
    // wrong vertex
    InitialFit initFitWV(mass,MH,mhLow_,mhHigh_);
    initFitWV.setVerbosity(verbose_);
    initFitWV.buildSumOfGaussians(Form("%s_cat%d",proc.c_str(),cat),nGaussiansWV,recursive_);
    initFitWV.setDatasets(datasetsWV);
    initFitWV.runFits(1);
    initFitWV.saveParamsToFileAtMH(Form("dat/in/%s_cat%d_wv.dat",proc.c_str(),cat),constraintValueMass_);
    initFitWV.loadPriorConstraints(Form("dat/in/%s_cat%d_wv.dat",proc.c_str(),cat),constraintValue_);
    initFitWV.runFits(1);
    initFitRV.plotFits(Form("plots/%s_cat%d/wv",proc.c_str(),cat));
    map<int,map<string,RooRealVar*> > fitParamsWV = initFitWV.getFitParams();

    //these guys do the interpolation
    // right vertex
    LinearInterp linInterpRV(MH,mhLow_,mhHigh_,fitParamsRV,doSecondaryModels_);
    linInterpRV.setVerbosity(verbose_);
    linInterpRV.setSecondaryModelVars(MH_SM,DeltaM,MH_2,higgsDecayWidth);
    linInterpRV.interpolate(nGaussiansRV);
    map<string,RooSpline1D*> splinesRV = linInterpRV.getSplines();

    // wrong vertex
    LinearInterp linInterpWV(MH,mhLow_,mhHigh_,fitParamsWV,doSecondaryModels_);
    linInterpWV.setVerbosity(verbose_);
    linInterpWV.setSecondaryModelVars(MH_SM,DeltaM,MH_2,higgsDecayWidth);
    linInterpWV.interpolate(nGaussiansWV);
    map<string,RooSpline1D*> splinesWV = linInterpWV.getSplines();

    // this guy constructs the final model with systematics, eff*acc etc.
    FinalModelConstruction finalModel(mass,MH,intLumi,mhLow_,mhHigh_,proc,cat,doSecondaryModels_,systfilename_,verbose_,false);
    finalModel.setSecondaryModelVars(MH_SM,DeltaM,MH_2,higgsDecayWidth);
    finalModel.setRVsplines(splinesRV);
    finalModel.setWVsplines(splinesWV);
    finalModel.setRVdatasets(datasetsRV);
    finalModel.setWVdatasets(datasetsWV);
    //finalModel.setSTDdatasets(datasets);
		finalModel.makeSTDdatasets();
    finalModel.buildRvWvPdf("hggpdfsmrel",nGaussiansRV,nGaussiansWV,recursive_);
    finalModel.getNormalization();
    finalModel.plotPdf("plots");
    finalModel.save(outWS);
  }
  
  datfile.close();

  sw.Stop();
  cout << "Whole fitting process took..." << endl;
  cout << "\t";
  sw.Print();

  sw.Start();
  
  cout << "Starting to combine fits..." << endl;
  // this guy packages everything up
  Packager packager(outWS,splitVH_,nCats_,mhLow_,mhHigh_);
  packager.packageOutput();
  sw.Stop();
  cout << "Combination complete." << endl;
  cout << "Whole process took..." << endl;
  cout << "\t";
  sw.Print();

  cout << "Writing to file..." << endl;
  outFile->cd();
  outWS->Write();
  outFile->Close();
  inFile->Close();
  cout << "Done." << endl;

  return 0;
}
