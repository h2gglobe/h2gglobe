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
string plotDir_;
int mhLow_=110;
int mhHigh_=150;
int nCats_;
float constraintValue_;
int constraintValueMass_;
bool spin_=false;
vector<string> procs_;
string procStr_;
bool isCutBased_=false;
bool is2011_=false;
string massesToSkip_;
vector<int> skipMasses_;
bool splitRVWV_=true;
bool doSecondaryModels_=true;
bool runInitialFitsOnly_=false;
bool recursive_=true;
int verbose_=0;

void OptionParser(int argc, char *argv[]){
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                			"Show help")
    ("infilename,i", po::value<string>(&filename_),                                           			"Input file name")
    ("outfilename,o", po::value<string>(&outfilename_)->default_value("CMS-HGG_sigfit.root"), 			"Output file name")
    ("datfilename,d", po::value<string>(&datfilename_)->default_value("dat/config.dat"),      			"Configuration file")
		("systfilename,s", po::value<string>(&systfilename_)->default_value("dat/photonCatSyst.dat"),		"Systematic model numbers")
    ("plotDir,p",	po::value<string>(&plotDir_)->default_value("plots"),															"Put plots in this directory")
		("mhLow,L", po::value<int>(&mhLow_)->default_value(110),                                  			"Low mass point")
    ("mhHigh,H", po::value<int>(&mhHigh_)->default_value(150),                                			"High mass point")
    ("nCats,n", po::value<int>(&nCats_)->default_value(9),                                    			"Number of total categories")
    ("constraintValue,C", po::value<float>(&constraintValue_)->default_value(0.1),            			"Constraint value")
    ("constraintValueMass,M", po::value<int>(&constraintValueMass_)->default_value(125),            "Constraint value mass")
    ("skipSecondaryModels",                                                                   			"Turn off creation of all additional models")
		("procs", po::value<string>(&procStr_)->default_value("ggh,vbf,wh,zh,tth"),											"Processes (comma sep)")
    ("isCutBased",                                                                               		"Is this the cut based analysis")
    ("is2011",                                                                               				"Is this the 7TeV analysis")
		("skipMasses", po::value<string>(&massesToSkip_)->default_value(""),														"Skip these mass points - used eg for the 7TeV where there's no mc at 145")
		("runInitialFitsOnly",																																					"Just fit gaussians - no interpolation, no systematics - useful for testing nGaussians")
    ("nonRecursive",                                                                             		"Do not recursively calculate gaussian fractions")
    ("verbose,v", po::value<int>(&verbose_)->default_value(0),                                			"Verbosity level: 0 (lowest) - 3 (highest)")
  ;                                                                                             		
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")){ cout << desc << endl; exit(1);}
  if (vm.count("spin"))                     spin_=true;
  if (vm.count("isCutBased"))               isCutBased_=true;
  if (vm.count("is2011"))               		is2011_=true;
  if (vm.count("runInitialFitsOnly"))       runInitialFitsOnly_=true;
  if (vm.count("nosplitRVWV"))              splitRVWV_=false;
  if (vm.count("skipSecondaryModels"))      doSecondaryModels_=false;
  if (vm.count("recursive"))                recursive_=false;
	if (vm.count("skipMasses")) {
		cout << "Masses to skip... " << endl;
		vector<string> els;
		split(els,massesToSkip_,boost::is_any_of(","));
		if (els.size()>0 && massesToSkip_!="") {
			for (vector<string>::iterator it=els.begin(); it!=els.end(); it++) {
				skipMasses_.push_back(boost::lexical_cast<int>(*it));
			}
		}
		cout << "\t";
		for (vector<int>::iterator it=skipMasses_.begin(); it!=skipMasses_.end(); it++) cout << *it << " ";
		cout << endl;
	}
	split(procs_,procStr_,boost::is_any_of(","));
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

bool skipMass(int mh){
	for (vector<int>::iterator it=skipMasses_.begin(); it!=skipMasses_.end(); it++) {
		if (*it==mh) return true;
	}
	return false;
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
  RooWorkspace *outWS;
	if (is2011_) outWS = new RooWorkspace("wsig_7TeV");
	else outWS = new RooWorkspace("wsig_8TeV");

  transferMacros(inFile,outFile);

	system(Form("mkdir -p %s/initialFits",plotDir_.c_str()));
	system("mkdir -p dat/in");

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
			if (skipMass(mh)) continue;
      RooDataSet *dataRV = (RooDataSet*)inWS->data(Form("sig_%s_mass_m%d_rv_cat%d",proc.c_str(),mh,cat));
      RooDataSet *dataWV = (RooDataSet*)inWS->data(Form("sig_%s_mass_m%d_wv_cat%d",proc.c_str(),mh,cat));
      RooDataSet *data = (RooDataSet*)inWS->data(Form("sig_%s_mass_m%d_cat%d",proc.c_str(),mh,cat));
      datasetsRV.insert(pair<int,RooDataSet*>(mh,dataRV));
      datasetsWV.insert(pair<int,RooDataSet*>(mh,dataWV));
      datasets.insert(pair<int,RooDataSet*>(mh,data));
    }

    // these guys do the fitting
    // right vertex
    InitialFit initFitRV(mass,MH,mhLow_,mhHigh_,skipMasses_);
    initFitRV.setVerbosity(verbose_);
    initFitRV.buildSumOfGaussians(Form("%s_cat%d",proc.c_str(),cat),nGaussiansRV,recursive_);
    initFitRV.setDatasets(datasetsRV);
    initFitRV.runFits(1);
		if (!runInitialFitsOnly_) {
			initFitRV.saveParamsToFileAtMH(Form("dat/in/%s_cat%d_rv.dat",proc.c_str(),cat),constraintValueMass_);
			initFitRV.loadPriorConstraints(Form("dat/in/%s_cat%d_rv.dat",proc.c_str(),cat),constraintValue_);
			initFitRV.runFits(1);
		}
    initFitRV.plotFits(Form("%s/initialFits/%s_cat%d_rv",plotDir_.c_str(),proc.c_str(),cat));
    map<int,map<string,RooRealVar*> > fitParamsRV = initFitRV.getFitParams();
    
    // wrong vertex
    InitialFit initFitWV(mass,MH,mhLow_,mhHigh_,skipMasses_);
    initFitWV.setVerbosity(verbose_);
    initFitWV.buildSumOfGaussians(Form("%s_cat%d",proc.c_str(),cat),nGaussiansWV,recursive_);
    initFitWV.setDatasets(datasetsWV);
    initFitWV.runFits(1);
		if (!runInitialFitsOnly_) {
			initFitWV.saveParamsToFileAtMH(Form("dat/in/%s_cat%d_wv.dat",proc.c_str(),cat),constraintValueMass_);
			initFitWV.loadPriorConstraints(Form("dat/in/%s_cat%d_wv.dat",proc.c_str(),cat),constraintValue_);
			initFitWV.runFits(1);
		}
    initFitWV.plotFits(Form("%s/initialFits/%s_cat%d_wv",plotDir_.c_str(),proc.c_str(),cat));
    map<int,map<string,RooRealVar*> > fitParamsWV = initFitWV.getFitParams();

		if (!runInitialFitsOnly_) {
			//these guys do the interpolation
			// right vertex
			LinearInterp linInterpRV(MH,mhLow_,mhHigh_,fitParamsRV,doSecondaryModels_,skipMasses_);
			linInterpRV.setVerbosity(verbose_);
			linInterpRV.setSecondaryModelVars(MH_SM,DeltaM,MH_2,higgsDecayWidth);
			linInterpRV.interpolate(nGaussiansRV);
			map<string,RooSpline1D*> splinesRV = linInterpRV.getSplines();

			// wrong vertex
			LinearInterp linInterpWV(MH,mhLow_,mhHigh_,fitParamsWV,doSecondaryModels_,skipMasses_);
			linInterpWV.setVerbosity(verbose_);
			linInterpWV.setSecondaryModelVars(MH_SM,DeltaM,MH_2,higgsDecayWidth);
			linInterpWV.interpolate(nGaussiansWV);
			map<string,RooSpline1D*> splinesWV = linInterpWV.getSplines();

			// this guy constructs the final model with systematics, eff*acc etc.
			FinalModelConstruction finalModel(mass,MH,intLumi,mhLow_,mhHigh_,proc,cat,doSecondaryModels_,systfilename_,skipMasses_,verbose_,isCutBased_,is2011_);
			finalModel.setSecondaryModelVars(MH_SM,DeltaM,MH_2,higgsDecayWidth);
			finalModel.setRVsplines(splinesRV);
			finalModel.setWVsplines(splinesWV);
			finalModel.setRVdatasets(datasetsRV);
			finalModel.setWVdatasets(datasetsWV);
			//finalModel.setSTDdatasets(datasets);
			finalModel.makeSTDdatasets();
			if (is2011_) {
				finalModel.buildRvWvPdf("hggpdfsmrel_7TeV",nGaussiansRV,nGaussiansWV,recursive_);
			} else {
				finalModel.buildRvWvPdf("hggpdfsmrel_8TeV",nGaussiansRV,nGaussiansWV,recursive_);
			}
			finalModel.getNormalization();
			finalModel.plotPdf(plotDir_);
			finalModel.save(outWS);
		}
  }
  
  datfile.close();

  sw.Stop();
  cout << "Whole fitting process took..." << endl;
  cout << "\t";
  sw.Print();

 	if (!runInitialFitsOnly_) { 
		sw.Start();
		cout << "Starting to combine fits..." << endl;
		// this guy packages everything up
		Packager packager(outWS,procs_,nCats_,mhLow_,mhHigh_,skipMasses_,is2011_,plotDir_);
		packager.packageOutput();
		sw.Stop();
		cout << "Combination complete." << endl;
		cout << "Whole process took..." << endl;
		cout << "\t";
		sw.Print();
	}

  cout << "Writing to file..." << endl;
  outFile->cd();
  outWS->Write();
  outFile->Close();
  inFile->Close();
  cout << "Done." << endl;

  return 0;
}
