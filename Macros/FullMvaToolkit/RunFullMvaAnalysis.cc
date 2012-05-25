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
  if (filename!="0") system(Form("cp %s %s_beforeFMT.root",filename.c_str(),filename.c_str()));

	FMTSetup *runner = new FMTSetup();
	runner->OptionParser(argc,argv);
	runner->CheckRunOptions();

	runner->printRunOptions("before.txt");
	runner->runRebinning();
	runner->printRunOptions("after.txt");
	runner->runFitting();
	runner->makeNormPlot();
	runner->createCorrBkgModel();
	runner->interpolateBDT();
	runner->writeDataCards();
	runner->publishToWeb();
	runner->runCombine();

  cout << "Original file " << filename << " updated." << endl;
  cout << "Old file copied to " << filename << "_beforeFMT.root" << endl;
	
	delete runner;
	
  return 0;
}
