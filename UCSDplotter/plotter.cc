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
// some constants / parameters
//----------------------------------------------------------------------

/** the name of the tree within the input file */
const string inputTreeName = "opttree";

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

//----------------------------------------------------------------------
void parseConfigFile(const std::string &configFname, HistoContainer *histoContainer)
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

  // read the configuration file and book histograms 
  parseConfigFile(configFname, histoContainer);

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
