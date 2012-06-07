#ifndef Normalization_8TeV_h
#define Normalization_8TeV_h

#include <vector>
#include <map>
#include <iostream>

#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TLegend.h"

using namespace std;

class Normalization_8TeV {

  public:
	Normalization_8TeV();
	
	double GetBR(double);
	double GetBR(int);
	double GetXsection(double,TString);
	double GetXsection(double);
	double GetXsection(int);
	double GetNorm(double,TH1F*,double, TH1F*,double);
	double GetMass(int);
	TString GetProcess(int);
	void CheckNorm(double,double,double,TString);
	void FillSignalTypes();
	void PlotExpected(double ,double);	
 private:
	std::map<double,double> BranchingRatioMap;
	std::map<double,double> XSectionMap_ggh;
	std::map<double,double> XSectionMap_vbf;
	std::map<double,double> XSectionMap_wh;
	std::map<double,double> XSectionMap_zh;
	std::map<double,double> XSectionMap_wzh;
	std::map<double,double> XSectionMap_tth;
	std::map<int,std::pair<TString,double > > SignalTypeMap;

};
#endif
