#include <iostream>
#include <cassert>
#include <string>
#include <fstream>

#include <TFile.h>
#include <TTree.h>

#include <dlfcn.h>

#include "GenericAnalysis.h"

#include <cstdlib>

using namespace std;

//----------------------------------------------------------------------
// some constants
//----------------------------------------------------------------------

const string inputTreeName = "opttree";


//----------------------------------------------------------------------
void parseConfigFile(const std::string &configFname)
{
  ifstream infile(configFname);
  string line;

  while ( infile.good() )
  {
    getline (infile,line);
    cout << line << endl;

    // remove everything after the first #
    
   }



}

//----------------------------------------------------------------------





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

  // bookHistograms(); // use HistoContainer


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
