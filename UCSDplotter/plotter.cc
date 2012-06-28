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


using namespace std;

//----------------------------------------------------------------------
// some constants
//----------------------------------------------------------------------

const string inputTreeName = "opttree";

vector<map<string, string> > histogramDefinitions;

//----------------------------------------------------------------------
void parseConfigFile(const std::string &configFname)
{
  ifstream infile(configFname.c_str());
  string line;

  while ( infile.good() )
  {
    getline (infile,line);
    cout << line << endl;

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

    histogramDefinitions.push_back(values);

  } // loop over lines
}

//----------------------------------------------------------------------

void bookHistograms(HistoContainer *container)
{
  for (vector<map<string, string> >::iterator it = histogramDefinitions.begin();
       it != histogramDefinitions.end();
       ++it)
  {
    map<string, string> &histoDef = *it;

    int htyp = boost::lexical_cast<int>(histoDef["htyp"]);
    switch (htyp)
    {
    case 0: // 1D histos
      container->Add(histoDef["name"],
                     histoDef["xaxis"],
                     histoDef["yaxis"],
                     boost::lexical_cast<unsigned>(histoDef["ncat"]),
                     boost::lexical_cast<unsigned>(histoDef["xbins"]),
                     boost::lexical_cast<unsigned>(histoDef["xmin"]),
                     boost::lexical_cast<unsigned>(histoDef["xmax"])
                     );
      break;


    case 1: // 2D histos
      container->Add(histoDef["name"],
                     histoDef["xaxis"],
                     histoDef["yaxis"],
                     boost::lexical_cast<unsigned>(histoDef["ncat"]),

                     boost::lexical_cast<unsigned>(histoDef["xbins"]),
                     boost::lexical_cast<unsigned>(histoDef["xmin"]),
                     boost::lexical_cast<unsigned>(histoDef["xmax"]),

                     boost::lexical_cast<unsigned>(histoDef["ybins"]),
                     boost::lexical_cast<unsigned>(histoDef["ymin"]),
                     boost::lexical_cast<unsigned>(histoDef["ymax"])
                     );
      break;

      
    default:
      cerr << "unknown htyp " << htyp << endl;
      exit(1);
    }
  }
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

int main(int argc, char **argv)
{
  assert(argc == 4 + 1);

  string configFname = argv[1];
  string analysisCodeFname = argv[2];
  string inputFname = argv[3];
  string outputFname = argv[4];

  parseConfigFile(configFname);


  HistoContainer *histoContainer = new HistoContainer();

  bookHistograms(histoContainer); // use HistoContainer

  //--------------------
  // open the input file
  //--------------------

  TFile *fin = TFile::Open(inputFname.c_str());

  assert(fin != NULL);
  assert(fin->IsOpen());

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
    analysis->analyze(histoContainer);

  } // end of loop over events

  //--------------------

  TFile *fout = TFile::Open(outputFname.c_str(), "RECREATE");

  // write out histograms to output file
  fout->cd();
  histoContainer->Save();

  // close output file
  fout->Close();
}
