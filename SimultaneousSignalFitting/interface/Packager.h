#ifndef Packager_h 
#define Packager_h

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "RooWorkspace.h"

#include "../../Macros/Normalization_8TeV.h"

class Packager {

  public:

    Packager(RooWorkspace *ws, bool splitVH, int nCats, int mhLow, int mhHigh, bool is2011=false);
    ~Packager();

    void packageOutput();

  private:
    RooWorkspace *outWS;
    bool splitVH_;
    int nCats_;
    int mhLow_;
    int mhHigh_;
		bool is2011_;
		int sqrts_;
    std::vector<std::string> procs;
    Normalization_8TeV *normalization;

};
#endif
