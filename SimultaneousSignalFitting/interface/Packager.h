#ifndef Packager_h 
#define Packager_h

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

#include "../../Macros/Normalization_8TeV.h"

class Packager {

  public:

    Packager(RooWorkspace *ws, bool splitVH, int nCats, int mhLow, int mhHigh, bool is2011=false, string outDir="plots");
    ~Packager();

    void packageOutput();
		void makePlots();
		void makePlot(RooRealVar *mass, RooRealVar *MH, RooAddPdf *pdf, std::map<int,RooDataSet*> data, std::string name);

  private:
    RooWorkspace *outWS;
    bool splitVH_;
    int nCats_;
    int mhLow_;
    int mhHigh_;
		bool is2011_;
		string outDir_;
		int sqrts_;
    std::vector<std::string> procs;
    Normalization_8TeV *normalization;

};
#endif
