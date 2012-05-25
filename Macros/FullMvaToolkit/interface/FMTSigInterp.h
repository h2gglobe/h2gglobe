#ifndef FMTSigInterp_h
#define FMTSigInterp_h

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "FMTBase.h"

using namespace std;

class FMTSigInterp : public FMTBase {

  public:
    FMTSigInterp(string, bool, bool, int, int, double, double, double, int, double, double, int, int, int, double, double, bool, bool, vector<string>, bool, vector<map<int,vector<double> > >, bool blind=false,bool verbose=false);
    ~FMTSigInterp();

    TH1F *Interpolate(double,TH1F*,double,TH1F*,double);

    TH1F* linearBin(TH1F*);
    void plotSidebands(TH1F*,vector<TH1F*>,vector<TH1F*>,string);
    void plotOutput(TH1F*, TH1F*, TH1F*, string);
    void plotFrac(TH1F*,TH1F*,TH1F*,string);

    void runInterpolation();

  private:
    TFile *tFile;
    bool diagnose_;
    bool blind_;
};

#endif
