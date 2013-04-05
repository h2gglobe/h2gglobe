#include <TFile.h>
#include <TTree.h>
#include <TRint.h>

#include <dlfcn.h>

#include "GenericAnalysis.h"
#include "parser.h"

#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <fstream>

using namespace std;

/** the name of the tree within the input file */
const string inputTreeName = "opttree";

void buildBkgModel(RooContainer* rooContainer,  map<string, string> values) {
  
  std::string postfix = (getUint(values,"dataIs2011")?"":"_8TeV");
  unsigned int nCategories_ = getUint(values, "nCategories");
  std::vector<int> bkgPolOrderByCat = getVint(values, "bkgPolOrderByCat");

  // sanity check
  assert(bkgPolOrderByCat.size() == nCategories_);

  rooContainer->AddRealVar("CMS_hgg_pol6_0" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_1" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_2" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_3" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_4" + postfix, -0.01,- 1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_5" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_0" + postfix, "@0*@0", "CMS_hgg_pol6_0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_1" + postfix, "@0*@0", "CMS_hgg_pol6_1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_2" + postfix, "@0*@0", "CMS_hgg_pol6_2" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_3" + postfix, "@0*@0", "CMS_hgg_pol6_3" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_4" + postfix, "@0*@0", "CMS_hgg_pol6_4" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_5" + postfix, "@0*@0", "CMS_hgg_pol6_4" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_pol5_0" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_1" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_2" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_3" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_4" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_0" + postfix, "@0*@0", "CMS_hgg_pol5_0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_1" + postfix, "@0*@0", "CMS_hgg_pol5_1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_2" + postfix, "@0*@0", "CMS_hgg_pol5_2" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_3" + postfix, "@0*@0", "CMS_hgg_pol5_3" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_4" + postfix, "@0*@0", "CMS_hgg_pol5_4" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_quartic0"+postfix,-0.1,-1.0,1.0);
  rooContainer->AddRealVar("CMS_hgg_quartic1"+postfix,-0.1,-1.0,1.0);
  rooContainer->AddRealVar("CMS_hgg_quartic2"+postfix,-0.1,-1.0,1.0);
  rooContainer->AddRealVar("CMS_hgg_quartic3"+postfix,-0.01,-1.0,1.0);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic0" + postfix, "@0*@0", "CMS_hgg_quartic0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic1" + postfix, "@0*@0", "CMS_hgg_quartic1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic2" + postfix, "@0*@0", "CMS_hgg_quartic2" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic3" + postfix, "@0*@0", "CMS_hgg_quartic3" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_quad0" + postfix, -0.1, -1.5, 1.5);
  rooContainer->AddRealVar("CMS_hgg_quad1" + postfix, -0.01, -1.5, 1.5);
  rooContainer->AddFormulaVar("CMS_hgg_modquad0" + postfix, "@0*@0", "CMS_hgg_quad0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquad1" + postfix, "@0*@0", "CMS_hgg_quad1" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_cubic0" + postfix, -0.1, -1.5, 1.5);
  rooContainer->AddRealVar("CMS_hgg_cubic1" + postfix, -0.1, -1.5, 1.5);
  rooContainer->AddRealVar("CMS_hgg_cubic2" + postfix, -0.01, -1.5, 1.5);
  rooContainer->AddFormulaVar("CMS_hgg_modcubic0" + postfix, "@0*@0", "CMS_hgg_cubic0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modcubic1" + postfix, "@0*@0", "CMS_hgg_cubic1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modcubic2" + postfix, "@0*@0", "CMS_hgg_cubic2" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_lin0" + postfix, -0.01, -1.5, 1.5);
  rooContainer->AddFormulaVar("CMS_hgg_modlin0" + postfix, "@0*@0", "CMS_hgg_lin0" + postfix);
  
  // prefix for models parameters
  std::map<int,std::string> parnames;
  parnames[1] = "modlin";
  parnames[2] = "modquad";
  parnames[3] = "modcubic";
  parnames[4] = "modquartic";
  parnames[5] = "modpol5_";
  parnames[6] = "modpol6_";
    
  // map order to categories flags + parameters names
  std::map<int, std::pair<std::vector<int>, std::vector<std::string> > > catmodels;
  // fill the map
  for(unsigned int icat=0; icat<nCategories_; ++icat) {
    // get the poly order for this category
    int catmodel = bkgPolOrderByCat[icat];
    std::vector<int> & catflags = catmodels[catmodel].first;
    std::vector<std::string> & catpars = catmodels[catmodel].second;
    // if this is the first time we find this order, build the parameters
    if( catflags.empty() ) {
      assert( catpars.empty() );
      // by default no category has the new model
      catflags.resize(nCategories_, 0);
      std::string & parname = parnames[catmodel];
      for(int iorder = 0; iorder<catmodel; ++iorder) {
	catpars.push_back( Form( "CMS_hgg_%s%d%s", parname.c_str(), iorder, +postfix.c_str() ) );
      }
    } else {
      assert( catflags.size() == nCategories_ && catpars.size() == catmodel );
    }
    // chose category order
    catflags[icat] = 1;
  }
  
  // now loop over the models and allocate the pdfs
  for(std::map<int, std::pair<std::vector<int>, std::vector<std::string> > >::iterator modit = catmodels.begin();
      modit!=catmodels.end(); ++modit ) {
    std::vector<int> & catflags = modit->second.first;
    std::vector<std::string> & catpars = modit->second.second;
    
    rooContainer->AddSpecificCategoryPdf(&catflags[0],"data_pol_model" + postfix, "0", "CMS_hll_mass", catpars, 70+catpars.size()); 
    // >= 71 means RooBernstein of order >= 1
  }
}

void initRooContainer(RooContainer* rooContainer, map<string, string> values) { 

  rooContainer->SetNCategories(getUint(values, "nCategories"));
  rooContainer->AddGlobalSystematic("lumi", 1.045, 1.00);
  
  // Create observables for shape-analysis with ranges
  rooContainer->AddObservable("CMS_hll_mass", getFloat(values, "massMin"), getFloat(values, "massMax"));
  rooContainer->AddConstant("IntLumi", getFloat(values, "intlumi"));

  // SM Model
  //rooContainer->AddConstant("XSBR_tth_155", 0.00004370);
  //rooContainer->AddConstant("XSBR_ggh_150", 0.01428);
  //rooContainer->AddConstant("XSBR_vbf_150", 0.001308);
  //rooContainer->AddConstant("XSBR_wzh_150", 0.000641);
  //rooContainer->AddConstant("XSBR_tth_150", 0.000066);
  //rooContainer->AddConstant("XSBR_ggh_145", 0.018820);
  //rooContainer->AddConstant("XSBR_vbf_145", 0.001676);
  //rooContainer->AddConstant("XSBR_wzh_145", 0.000891);
  //rooContainer->AddConstant("XSBR_tth_145", 0.000090);
  //rooContainer->AddConstant("XSBR_ggh_140", 0.0234109);
  //rooContainer->AddConstant("XSBR_vbf_140", 0.00203036);
  //rooContainer->AddConstant("XSBR_wzh_140", 0.001163597);
  //rooContainer->AddConstant("XSBR_tth_140", 0.000117189);
  //rooContainer->AddConstant("XSBR_ggh_135", 0.0278604);
  //rooContainer->AddConstant("XSBR_vbf_135", 0.002343);
  //rooContainer->AddConstant("XSBR_wzh_135", 0.001457559);
  //rooContainer->AddConstant("XSBR_tth_135", 0.000145053);
  //rooContainer->AddConstant("XSBR_ggh_130", 0.0319112);
  //rooContainer->AddConstant("XSBR_vbf_130", 0.00260804);
  //rooContainer->AddConstant("XSBR_wzh_130", 0.001759636);
  //rooContainer->AddConstant("XSBR_tth_130", 0.000173070);
  //rooContainer->AddConstant("XSBR_ggh_125", 0.0350599);
  //rooContainer->AddConstant("XSBR_vbf_125", 0.00277319);
  //rooContainer->AddConstant("XSBR_wzh_125", 0.002035123);
  //rooContainer->AddConstant("XSBR_tth_125", 0.000197718);
  //rooContainer->AddConstant("XSBR_ggh_120", 0.0374175);
  //rooContainer->AddConstant("XSBR_vbf_120", 0.00285525);
  //rooContainer->AddConstant("XSBR_wzh_120", 0.002285775);
  //rooContainer->AddConstant("XSBR_tth_120", 0.00021951);
  //rooContainer->AddConstant("XSBR_ggh_123", 0.0360696);
  //rooContainer->AddConstant("XSBR_vbf_123", 0.00281352);
  //rooContainer->AddConstant("XSBR_wzh_123", 0.00213681);
  //rooContainer->AddConstant("XSBR_tth_123", 0.00020663);
  //rooContainer->AddConstant("XSBR_ggh_121", 0.0369736);
  //rooContainer->AddConstant("XSBR_vbf_121", 0.00284082);
  //rooContainer->AddConstant("XSBR_wzh_121", 0.00223491);
  //rooContainer->AddConstant("XSBR_tth_121", 0.00021510);
  //rooContainer->AddConstant("XSBR_ggh_115", 0.0386169);
  //rooContainer->AddConstant("XSBR_vbf_115", 0.00283716);
  //rooContainer->AddConstant("XSBR_wzh_115", 0.002482089);
  //rooContainer->AddConstant("XSBR_tth_115", 0.000235578);
  //rooContainer->AddConstant("XSBR_ggh_110", 0.0390848);
  //rooContainer->AddConstant("XSBR_vbf_110", 0.00275406);
  //rooContainer->AddConstant("XSBR_wzh_110", 0.002654575);
  //rooContainer->AddConstant("XSBR_tth_110", 0.000247629);
  //rooContainer->AddConstant("XSBR_ggh_105", 0.0387684);
  //rooContainer->AddConstant("XSBR_vbf_105", 0.00262016);
  //rooContainer->AddConstant("XSBR_wzh_105", 0.002781962);
  //rooContainer->AddConstant("XSBR_tth_105", 0.000255074);
  //
  //// FF model 
  //rooContainer->AddConstant("ff_XSBR_vbf_150", 0.00259659);
  //rooContainer->AddConstant("ff_XSBR_wzh_150", 0.00127278);
  //rooContainer->AddConstant("ff_XSBR_vbf_145", 0.00387544);
  //rooContainer->AddConstant("ff_XSBR_wzh_145", 0.00205969);
  //rooContainer->AddConstant("ff_XSBR_vbf_140", 0.00565976);
  //rooContainer->AddConstant("ff_XSBR_wzh_140", 0.003243602);
  //rooContainer->AddConstant("ff_XSBR_vbf_135", 0.00825);
  //rooContainer->AddConstant("ff_XSBR_wzh_135", 0.00513225);
  //rooContainer->AddConstant("ff_XSBR_vbf_130", 0.0122324);
  //rooContainer->AddConstant("ff_XSBR_wzh_130", 0.00825316);
  //rooContainer->AddConstant("ff_XSBR_vbf_125", 0.0186494);
  //rooContainer->AddConstant("ff_XSBR_wzh_125", 0.01368598);
  //rooContainer->AddConstant("ff_XSBR_vbf_123", 0.022212);
  //rooContainer->AddConstant("ff_XSBR_wzh_123", 0.0168696);
  //rooContainer->AddConstant("ff_XSBR_vbf_121", 0.0266484);
  //rooContainer->AddConstant("ff_XSBR_wzh_121", 0.0209646);
  //rooContainer->AddConstant("ff_XSBR_vbf_120", 0.0293139);
  //rooContainer->AddConstant("ff_XSBR_wzh_120", 0.02346729);
  //rooContainer->AddConstant("ff_XSBR_vbf_115", 0.0482184);
  //rooContainer->AddConstant("ff_XSBR_wzh_115", 0.04218386);
  //rooContainer->AddConstant("ff_XSBR_vbf_110", 0.083181);
  //rooContainer->AddConstant("ff_XSBR_wzh_110", 0.08017625);
  //rooContainer->AddConstant("ff_XSBR_vbf_105", 0.151616);
  //rooContainer->AddConstant("ff_XSBR_wzh_105", 0.1609787);
  
  // -----------------------------------------------------
  // Configurable background model
  // if no configuration was given, set some defaults

  //int nCategories_ = getUint(values, "nCategories");
  unsigned int nDataBins = getUint(values, "nDataBins");
  std::vector<int>   bkgPolOrderByCat = getVint(values, "bkgPolOrderByCat");
  std::vector<float> sigPointsToBook  = getVfloat(values, "signalPoints");

  // build the model
  buildBkgModel(rooContainer, values);
  
  rooContainer->CreateDataSet("CMS_hll_mass", "data_mass", nDataBins);
  rooContainer->CreateDataSet("CMS_hll_mass", "bkg_mass" , nDataBins);            
  
  // Create Signal DataSets:
  for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
    int sig = sigPointsToBook[isig];
    rooContainer->CreateDataSet("CMS_hll_mass", Form("sig_ggh_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hll_mass", Form("sig_vbf_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hll_mass", Form("sig_rsg_mass_m%d", sig), nDataBins);    
  }
  
  // Make more datasets representing Systematic Shifts of various quantities
  //for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
  //  int sig = sigPointsToBook[isig];
  //  rooContainer->MakeSystematics("CMS_hll_mass", Form("sig_ggh_mass_m%d", sig), -1);    
  //  rooContainer->MakeSystematics("CMS_hll_mass", Form("sig_vbf_mass_m%d", sig), -1);    
  //  rooContainer->MakeSystematics("CMS_hll_mass", Form("sig_wzh_mass_m%d", sig), -1);    
  //  rooContainer->MakeSystematics("CMS_hll_mass", Form("sig_tth_mass_m%d", sig), -1);    
  //}
}

void fitToData(RooContainer* rooContainer, std::map<std::string, std::string> values) {

  std::string postfix = (getUint(values,"dataIs2011")?"":"_8TeV");
  rooContainer->FitToData("data_pol_model" + postfix, "data_mass");  // Fit to full range of dataset
}

//----------------------------------------------------------------------

void initConfig(std::map<std::string, std::string> values, RooContainer* rooContainer) {

  bool doBlinding = bool(getUint(values, "doBlinding"));
  rooContainer->BlindData(doBlinding);
  initRooContainer(rooContainer, values);
  //bookHistogram(histoContainer, values);
}

//----------------------------------------------------------------------

GenericAnalysis *openAnalysisCode(const string &fname) {
  
  // open the shared library
  void* soHandle = dlopen(fname.c_str(), RTLD_LAZY);

  if (!soHandle) {
    cerr << "could not open analysis shared object file '" << fname << "': " << dlerror() << endl;
    exit(1);
  }
  
  // load the instantiation function (the function which creates
  // an instance of the needed class)
  typedef GenericAnalysis *(*InstantiationFunctionType)();
  
  dlerror();
  const char *instantiationFunctionName = "makeInstance";
  
  InstantiationFunctionType instantiationFunction = (InstantiationFunctionType) dlsym(soHandle, instantiationFunctionName);
  const char *errMessage = dlerror();
  if (errMessage) {
    cerr << "could not find symbol '" << instantiationFunctionName << "' in analysis shared object file '" << fname << "'. Exiting." << endl;
    dlclose(soHandle);
    exit(1);
  }
  
  GenericAnalysis *analysis = instantiationFunction();
  return analysis;
}

//----------------------------------------------------------------------

void usage() {
  cerr << endl
       << "usage:   workspacer config.dat MyAnalysis.so input.root output.root" << endl
       << endl
    ;
  exit(1);
}

//----------------------------------------------------------------------

int main(int argc, char **argv) {
  if (argc != 4 + 1)
    usage();
  
  string configFname = argv[1];
  string analysisCodeFname = argv[2];
  string inputFname = argv[3];
  string outputFname = argv[4];

  RooContainer *rooContainer = new RooContainer();
  
  // read the configuration file and book histograms 
  map<std::string, std::string> values = parseConfigFile(configFname);
  initConfig(values, rooContainer);
  
  //--------------------
  // open the input file
  //--------------------
  TFile *fin = TFile::Open(inputFname.c_str());
  
  if (fin == NULL || !fin->IsOpen()) {
    cerr << "could not open input file '" << inputFname << "'" << endl;
    exit(1);
  }

  gROOT->cd();

  TTree *tree = (TTree*)(fin->Get(inputTreeName.c_str()));
  assert(tree != NULL);
  
  //--------------------
  // open the analysis shared object
  //--------------------
  GenericAnalysis *analysis = openAnalysisCode(analysisCodeFname);
  
  // must also activate branches
  analysis->setBranchAddresses(tree, values);
  
  // loop on the events
  unsigned numEvents = tree->GetEntries();
  for (unsigned i = 0; i < numEvents; ++i) {
    // read ith event
    tree->GetEntry(i);
    
    // call user function
    analysis->analyze(rooContainer);
    
  } // end of loop over events

  //--------------------
  
  fitToData(rooContainer, values);

  TFile *fout = TFile::Open(outputFname.c_str(), "RECREATE");

  // write out histograms to output file
  fout->cd();
  rooContainer->Save();
  
  // close output file
  fout->Close();


}
