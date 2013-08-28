#ifndef Packager_h 
#define Packager_h

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "RooWorkspace.h"

#include "../interface/Normalization_8TeV.h"

class Packager {

  public:

    Packager(RooWorkspace *ws, bool splitVH, int nCats, int mhLow, int mhHigh);
    ~Packager();

    void packageOutput();

  private:
    RooWorkspace *outWS;
    bool splitVH_;
    int nCats_;
    int mhLow_;
    int mhHigh_;
    std::vector<std::string> procs;
    Normalization_8TeV *normalization;

};
#endif
