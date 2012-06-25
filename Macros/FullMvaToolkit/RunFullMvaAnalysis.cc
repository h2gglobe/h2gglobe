#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include <typeinfo>
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TFile.h"
#include "TMacro.h"
#include "TObjString.h"

#include "interface/FMTPlots.h"
#include "interface/FMTBase.h"
#include "interface/FMTRebin.h"
#include "interface/FMTFit.h"
#include "interface/FMTSetup.h"
#include "interface/FMTSigInterp.h"

using namespace std;
namespace po = boost::program_options;

string getFileName(int argc, char* argv[]){
  string filename="0";
  for (int i=0; i<argc; i++){
    if (string(argv[i])=="--filename" || string(argv[i])=="-i") filename=string(argv[i+1]);
  }
  return filename;
}

int main(int argc, char* argv[]){

  string filename = getFileName(argc,argv);

	FMTSetup *runner = new FMTSetup(filename);
	runner->OptionParser(argc,argv);
	runner->CheckRunOptions();

	runner->printRunOptions("before.txt");
	runner->runRebinning();
	runner->printRunOptions("after.txt");
	runner->runFitting();
	runner->createCorrBkgModel();
	runner->interpolateBDT();
	runner->writeDataCards();
  runner->makePlots();
	runner->publishToWeb();
	runner->runCombine();

	delete runner;
	
  return 0;
}
