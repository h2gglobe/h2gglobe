#ifndef FMTFit_h
#define FMTFit_h

#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"

#include "TFile.h"

#include "FMTBase.h"

using namespace std;
using namespace RooFit;

class FMTFit : public FMTBase{

	public:
		FMTFit(TFile *tFile, int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, bool includeVBF, bool includeLEP, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int,vector<double> > > AllBinEdges, bool verbose=false);
		~FMTFit();

		pair<double,double> FitPow(double);
		void redoFit(double);
		void Plot(double);

		bool getblind();
		bool getplot();

		void setblind(bool);
		void setplot(bool);

	private:
		RooRealVar *r1, *r2, *f1;
		RooAbsPdf *fit;
		RooRealVar *nBkgInSigReg;
		RooWorkspace *outWS;
		RooRealVar *mass_var;
		RooDataSet *data;

		bool blind_;
		bool plot_;

};

#endif
