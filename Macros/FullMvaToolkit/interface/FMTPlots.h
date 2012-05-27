#ifndef FMTPlots_h
#define FMTPlots_h

#include <iostream>

#include "TH1F.h"

#include "FMTBase.h"

using namespace std;

class FMTPlots : public FMTBase {

  public:
    FMTPlots(string, int, int, double, double, double, int, double, double, int, int, int, double, double, bool, bool, vector<string>, bool, vector<map<int,vector<double> > >, bool blind=false,bool verbose=false);
    ~FMTPlots();
  
    TH1F *linearBin(TH1F*);
    void plotAll(double);
    void makePlots(TH1F*,TH1F*,TH1F*,TH1F*,vector<TH1F*>,vector<TH1F*>,vector<TH1F*>,vector<TH1F*>,vector<pair<TH1F*,TH1F*> >,TH1F*,TH1F*,double);

    void plotOutput(TH1F*,TH1F*,TH1F*,TH1F*,TH1F*,double);
    void plotSidebands(TH1F*, vector<TH1F*>, vector<TH1F*>, TH1F*, vector<TH1F*>, vector<TH1F*>,double);
    void plotSystFracs(TH1F*, vector<pair<TH1F*,TH1F*> >,double);
    void plotInterpolation(TH1F*, TH1F*, TH1F*, double);
    void makeNormPlot();

  private:
    
    TFile *tFile;
    bool blind_;
    bool verbose_;
};

#endif
