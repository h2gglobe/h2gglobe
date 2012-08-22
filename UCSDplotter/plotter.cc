#include <iostream>
#include <cassert>
#include <string>
#include <fstream>

#include <TFile.h>
#include <TTree.h>

#include <dlfcn.h>

#include "GenericAnalysis.h"

#include <cstdlib>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/algorithm/string/regex.hpp>


#include <vector>

#include <TRint.h>

#include "Plotting.h"


using namespace std;

//----------------------------------------------------------------------
// some constants / parameters
//----------------------------------------------------------------------

/** the name of the tree within the input file */
//const string inputTreeName = "opttree";
const string inputTreeName = "ntuple";

bool interactive = true;

//----------------------------------------------------------------------
std::string 
getString(const map<string, string> &values, const std::string &paramName)
{
  map<string, string>::const_iterator it = values.find(paramName);
  if (it == values.end())
    throw std::runtime_error("no parameter named '" + paramName + "' found");

  return it->second;
}

//----------------------------------------------------------------------

/** for optional parameters */
std::string 
getString(const map<string, string> &values, const std::string &paramName, const std::string &defaultValue)
{
  map<string, string>::const_iterator it = values.find(paramName);
  if (it == values.end())
    return defaultValue;
  else
    return it->second;
}

//----------------------------------------------------------------------
unsigned getUint(const map<string, string> &values, const std::string &paramName)
{
  try
    {
      return boost::lexical_cast<unsigned>(getString(values, paramName));
    }
  catch (const boost::bad_lexical_cast &ex)
    {
      throw std::runtime_error("value of parameter '" + paramName + "' is not a valid unsigned integer");
    }
}

//----------------------------------------------------------------------
unsigned getFloat(const map<string, string> &values, const std::string &paramName)
{
  try
    {
      return boost::lexical_cast<float>(getString(values, paramName));
    }
  catch (const boost::bad_lexical_cast &ex)
    {
      throw std::runtime_error("value of parameter '" + paramName + "' is not a valid float");
    }
}

//----------------------------------------------------------------------

#include <boost/algorithm/string.hpp>

std::vector<int> getVint(const map<string, string> &values, const std::string &paramName)
{
  try
    {
      vector<string> stringValues;
      string buf = getString(values, paramName);
      boost::split(stringValues, buf, boost::is_any_of(","));

      vector<int> retval;

      BOOST_FOREACH(string str, stringValues)
      {
	retval.push_back(boost::lexical_cast<int>(str));
      }

      return retval;
    }
  catch (const boost::bad_lexical_cast &ex)
    {
      throw std::runtime_error("values of parameter '" + paramName + "' contain invalid integers");
    }
}


//----------------------------------------------------------------------

void bookHistogram(HistoContainer *container, map<string, string> values)
{
  int htyp = getUint(values,"htyp");
  switch (htyp)
    {
    case 0: // 1D histos
      container->Add(getString(values,"name"),
                     getString(values,"xaxis", ""),
                     getString(values,"yaxis", ""),
                     getUint(values,"ncat"),
                     getUint(values,"xbins"),
                     getFloat(values,"xmin"),
                     getFloat(values,"xmax")
                     );
      break;


    case 1: // 2D histos
      container->Add(getString(values,"name"),
                     getString(values,"xaxis", ""),
                     getString(values,"yaxis", ""),
                     getUint(values,"ncat"),

                     getUint(values,"xbins"),
                     getFloat(values,"xmin"),
                     getFloat(values,"xmax"),

                     getUint(values,"ybins"),
                     getFloat(values,"ymin"),
                     getFloat(values,"ymax")
                     );
      break;

      
    default:
      throw std::runtime_error("unknown htyp " + boost::lexical_cast<string>(htyp));
    }
}

void initRooContainer(RooContainer* rooContainer, map<string, string> values) {

  rooContainer->SetNCategories(getUint(values, "nCategories"));
  rooContainer->nsigmas = getUint(values, "nSystSteps");
  rooContainer->sigmaRange = getUint(values, "systRange");
  
  // FIXME
  //l.rooContainer->MakeSystematicStudy(sys,sys_t); // * many times
  rooContainer->AddGlobalSystematic("lumi", 1.045, 1.00);
  
  // Create observables for shape-analysis with ranges
  rooContainer->AddObservable("CMS_hgg_mass", getFloat(values, "massMin"), getFloat(values, "massMax"));
  rooContainer->AddConstant("IntLumi", getFloat(values, "intlumi"));

  // SM Model
  rooContainer->AddConstant("XSBR_tth_155", 0.00004370);
  rooContainer->AddConstant("XSBR_ggh_150", 0.01428);
  rooContainer->AddConstant("XSBR_vbf_150", 0.001308);
  rooContainer->AddConstant("XSBR_wzh_150", 0.000641);
  rooContainer->AddConstant("XSBR_tth_150", 0.000066);
  rooContainer->AddConstant("XSBR_ggh_145", 0.018820);
  rooContainer->AddConstant("XSBR_vbf_145", 0.001676);
  rooContainer->AddConstant("XSBR_wzh_145", 0.000891);
  rooContainer->AddConstant("XSBR_tth_145", 0.000090);
  rooContainer->AddConstant("XSBR_ggh_140", 0.0234109);
  rooContainer->AddConstant("XSBR_vbf_140", 0.00203036);
  rooContainer->AddConstant("XSBR_wzh_140", 0.001163597);
  rooContainer->AddConstant("XSBR_tth_140", 0.000117189);
  rooContainer->AddConstant("XSBR_ggh_135", 0.0278604);
  rooContainer->AddConstant("XSBR_vbf_135", 0.002343);
  rooContainer->AddConstant("XSBR_wzh_135", 0.001457559);
  rooContainer->AddConstant("XSBR_tth_135", 0.000145053);
  rooContainer->AddConstant("XSBR_ggh_130", 0.0319112);
  rooContainer->AddConstant("XSBR_vbf_130", 0.00260804);
  rooContainer->AddConstant("XSBR_wzh_130", 0.001759636);
  rooContainer->AddConstant("XSBR_tth_130", 0.000173070);
  rooContainer->AddConstant("XSBR_ggh_125", 0.0350599);
  rooContainer->AddConstant("XSBR_vbf_125", 0.00277319);
  rooContainer->AddConstant("XSBR_wzh_125", 0.002035123);
  rooContainer->AddConstant("XSBR_tth_125", 0.000197718);
  rooContainer->AddConstant("XSBR_ggh_120", 0.0374175);
  rooContainer->AddConstant("XSBR_vbf_120", 0.00285525);
  rooContainer->AddConstant("XSBR_wzh_120", 0.002285775);
  rooContainer->AddConstant("XSBR_tth_120", 0.00021951);
  rooContainer->AddConstant("XSBR_ggh_123", 0.0360696);
  rooContainer->AddConstant("XSBR_vbf_123", 0.00281352);
  rooContainer->AddConstant("XSBR_wzh_123", 0.00213681);
  rooContainer->AddConstant("XSBR_tth_123", 0.00020663);
  rooContainer->AddConstant("XSBR_ggh_121", 0.0369736);
  rooContainer->AddConstant("XSBR_vbf_121", 0.00284082);
  rooContainer->AddConstant("XSBR_wzh_121", 0.00223491);
  rooContainer->AddConstant("XSBR_tth_121", 0.00021510);
  rooContainer->AddConstant("XSBR_ggh_115", 0.0386169);
  rooContainer->AddConstant("XSBR_vbf_115", 0.00283716);
  rooContainer->AddConstant("XSBR_wzh_115", 0.002482089);
  rooContainer->AddConstant("XSBR_tth_115", 0.000235578);
  rooContainer->AddConstant("XSBR_ggh_110", 0.0390848);
  rooContainer->AddConstant("XSBR_vbf_110", 0.00275406);
  rooContainer->AddConstant("XSBR_wzh_110", 0.002654575);
  rooContainer->AddConstant("XSBR_tth_110", 0.000247629);
  rooContainer->AddConstant("XSBR_ggh_105", 0.0387684);
  rooContainer->AddConstant("XSBR_vbf_105", 0.00262016);
  rooContainer->AddConstant("XSBR_wzh_105", 0.002781962);
  rooContainer->AddConstant("XSBR_tth_105", 0.000255074);
  
  // FF model 
  rooContainer->AddConstant("ff_XSBR_vbf_150", 0.00259659);
  rooContainer->AddConstant("ff_XSBR_wzh_150", 0.00127278);
  rooContainer->AddConstant("ff_XSBR_vbf_145", 0.00387544);
  rooContainer->AddConstant("ff_XSBR_wzh_145", 0.00205969);
  rooContainer->AddConstant("ff_XSBR_vbf_140", 0.00565976);
  rooContainer->AddConstant("ff_XSBR_wzh_140", 0.003243602);
  rooContainer->AddConstant("ff_XSBR_vbf_135", 0.00825);
  rooContainer->AddConstant("ff_XSBR_wzh_135", 0.00513225);
  rooContainer->AddConstant("ff_XSBR_vbf_130", 0.0122324);
  rooContainer->AddConstant("ff_XSBR_wzh_130", 0.00825316);
  rooContainer->AddConstant("ff_XSBR_vbf_125", 0.0186494);
  rooContainer->AddConstant("ff_XSBR_wzh_125", 0.01368598);
  rooContainer->AddConstant("ff_XSBR_vbf_123", 0.022212);
  rooContainer->AddConstant("ff_XSBR_wzh_123", 0.0168696);
  rooContainer->AddConstant("ff_XSBR_vbf_121", 0.0266484);
  rooContainer->AddConstant("ff_XSBR_wzh_121", 0.0209646);
  rooContainer->AddConstant("ff_XSBR_vbf_120", 0.0293139);
  rooContainer->AddConstant("ff_XSBR_wzh_120", 0.02346729);
  rooContainer->AddConstant("ff_XSBR_vbf_115", 0.0482184);
  rooContainer->AddConstant("ff_XSBR_wzh_115", 0.04218386);
  rooContainer->AddConstant("ff_XSBR_vbf_110", 0.083181);
  rooContainer->AddConstant("ff_XSBR_wzh_110", 0.08017625);
  rooContainer->AddConstant("ff_XSBR_vbf_105", 0.151616);
  rooContainer->AddConstant("ff_XSBR_wzh_105", 0.1609787);
  
  
  // -----------------------------------------------------
  // Make some data sets from the observables to fill in the event loop         
  // Binning is for histograms (will also produce unbinned data sets)
  unsigned int nDataBins = getUint(values, "nDataBins");
  rooContainer->CreateDataSet("CMS_hgg_mass", "data_mass", nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
  rooContainer->CreateDataSet("CMS_hgg_mass", "bkg_mass" , nDataBins);            
  

  std::vector<int> sigPointsToBook = getVint(values, "signalPoints");
  // Create Signal DataSets:
  for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
    int sig = sigPointsToBook[isig];
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_ggh_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_vbf_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_wzh_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_tth_mass_m%d", sig), nDataBins);   
    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_ggh_mass_m%d_rv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_vbf_mass_m%d_rv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_wzh_mass_m%d_rv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_tth_mass_m%d_rv", sig), nDataBins);    
    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_ggh_mass_m%d_wv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_vbf_mass_m%d_wv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_wzh_mass_m%d_wv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_tth_mass_m%d_wv", sig), nDataBins);    
  }
  
  // Make more datasets representing Systematic Shifts of various quantities
  for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
    int sig = sigPointsToBook[isig];
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_ggh_mass_m%d", sig), -1);    
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_vbf_mass_m%d", sig), -1);    
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_wzh_mass_m%d", sig), -1);    
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_tth_mass_m%d", sig), -1);    
  }
}

void fitToData(RooContainer* rooContainer, map<string, string> values) {

  std::string postfix = (getUint(values,"dataIs2011")?"":"_8TeV");
  rooContainer->FitToData("data_pol_model" + postfix, "data_mass");  // Fit to full range of dataset
}

void buildBkgModel(RooContainer* rooContainer,  map<string, string> values) {

  std::string postfix = (getUint(values,"dataIs2011")?"":"_8TeV");
  int nCategories_ = getUint(values, "nCategories");
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
  for(int icat=0; icat<nCategories_; ++icat) {
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
  /// for(size_t imodel=0; imodel<catmodels.size(); ++imodel ) {
  for(std::map<int, std::pair<std::vector<int>, std::vector<std::string> > >::iterator modit = catmodels.begin();
      modit!=catmodels.end(); ++modit ) {
    std::vector<int> & catflags = modit->second.first;
    std::vector<std::string> & catpars = modit->second.second;
    
    rooContainer->AddSpecificCategoryPdf(&catflags[0],"data_pol_model" + postfix, "0", "CMS_hgg_mass", catpars, 70+catpars.size()); 
    // >= 71 means RooBernstein of order >= 1
  }
}

//----------------------------------------------------------------------
void parseConfigFile(const std::string &configFname, HistoContainer *histoContainer, RooContainer* rooContainer)
{
  ifstream infile(configFname.c_str());
  string line;

  int lineNum = 0;

  while ( infile.good() )
  {
    try
      {
        getline (infile,line);
        ++lineNum;
        // cout << line << endl;

        // remove everything after the first #
        size_t pos = line.find('#');
        if (pos != string::npos)
          line.erase(pos);

        // remove leading a trailing white space
        boost::algorithm::trim(line);

        if (line == "")
          continue;

        // split line into name=value pairs, separated by whitespace
        vector<string> parts;
    
        boost::algorithm::split_regex(parts, line,
                                      boost::regex( "\\s+"));

    
        map<string, string> values;

        BOOST_FOREACH(std::string part, parts)
          {
            // loop over all parts of this line
            pos = part.find('=');
            assert(pos != string::npos);

            string key = part.substr(0,pos);
            string value = part.substr(pos+1);
            values[key] = value;
          }

        //histogramDefinitions.push_back(ConfigLineData(lineNum, values));

	// FIXME ADD A CONFIGURABLE SWITCH;
	initRooContainer(rooContainer, values); 
        bookHistogram(histoContainer, values);
      } // try
    catch (const std::exception &ex)
      {
        cerr << "exception caught while reading line " << lineNum << " of file " << configFname << ": " << ex.what() << endl;
        exit(1);
      }

  } // loop over lines
}

//----------------------------------------------------------------------


GenericAnalysis *openAnalysisCode(const string &fname)
{
  // open the shared library
  void* soHandle = dlopen(fname.c_str(), RTLD_LAZY);

  if (! soHandle)
    {
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
    if (errMessage)
    {
      cerr << "could not find symbol '" << instantiationFunctionName << "' in analysis shared object file '" << fname << "'. Exiting." << endl;
      dlclose(soHandle);
      exit(1);
    }

    GenericAnalysis *analysis = instantiationFunction();
    return analysis;

}


//----------------------------------------------------------------------

void usage()
{
  cerr << endl
       << "usage:   plotter plotvariables.dat myanalysis.so input.root output.root" << endl
       << endl
    ;
  exit(1);
}

//----------------------------------------------------------------------

int main(int argc, char **argv)
{
  if (argc != 4 + 1)
    usage();

  string configFname = argv[1];
  string analysisCodeFname = argv[2];
  string inputFname = argv[3];
  string outputFname = argv[4];

  HistoContainer *histoContainer = new HistoContainer();
  RooContainer *rooContainer = new RooContainer();

  // read the configuration file and book histograms 
  parseConfigFile(configFname, histoContainer, rooContainer);

  //--------------------

  TRint *myapp = NULL;
  if (interactive)
  {
    int dummyArgc = 1;
    char* dummyArgv[2];
    dummyArgv[0] = argv[0];
    dummyArgv[1] = NULL;

    myapp = new TRint("rint", &dummyArgc, &dummyArgv[0]);
  }

  //--------------------
  // open the input file
  //--------------------

  TFile *fin = TFile::Open(inputFname.c_str());

  if (fin == NULL || !fin->IsOpen())
  {
    cerr << "could not open input file '" << inputFname << "'" << endl;
    exit(1);
  }

  TTree *tree = (TTree*)(fin->Get(inputTreeName.c_str()));
  assert(tree != NULL);

  //--------------------
  // open the analysis shared object
  //--------------------
  GenericAnalysis *analysis = openAnalysisCode(analysisCodeFname);

  //--------------------

  // must also activate branches
  analysis->setBranchAddresses(tree);

  // loop on the events
  unsigned numEvents = tree->GetEntries();
  for (unsigned i = 0; i < numEvents; ++i)
  {
    // read ith event
    tree->GetEntry(i);

    // call user function
    analysis->analyze(histoContainer, rooContainer);

  } // end of loop over events

  //--------------------
  
  // FIXME ADD LATER
  //fitToData(rooContainer, values);

  TFile *fout = TFile::Open(outputFname.c_str(), "RECREATE");

  // write out histograms to output file
  fout->cd();
  histoContainer->Save();
  rooContainer->Save();

  // do plotting
  Plotting plotter(histoContainer);

  // see http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/Macros/Normalization_8TeV.cc?revision=1.10&view=markup
  // for itype numbers for each mass

  // mh = 125 GeV
  plotter.addSignalItype(-37);
  plotter.addSignalItype(-38);
  plotter.addSignalItype(-39);
  plotter.addSignalItype(-40);


  plotter.plotAll();

  if (interactive)
  {
    myapp->Run(kTRUE);            // run ROOT interactively
  }

  // close output file
  fout->Close();

}
