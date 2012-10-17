#ifndef FMTTree_h
#define FMTTree_h

#include "TH1F.h"

#include "RooRealVar.h"
#include "RooDataset.h"

#include "FMTFit.h"
#include "FMTBase.h"

using namespace std;

class FMTTree : public FMTBase {
	
	public:
    FMTTree(double,bool,int, int, double, double, double, int, double, double, int, int, int, double, double, int, bool, int, bool, int, vector<string>, bool, vector<map<int,vector<double> > >, bool verbose=false);
    ~FMTTree();

		map<string,TTree*> getSignalTrees(TFile *tFile, string dir="");
		map<string,TTree*> getDataTrees(TFile *tFile, string dir="");
		map<string,TTree*> getBackgroundTrees(TFile *tFile, string dir="");

		void initVariables();
		void setBranchVariables(TTree*);
		int icCat(int);

		void FillHist(string, int, double);
		void FillSigHist(string, double);
		void FillSystHist(string, double);

		string getProc(string);
		int getMH(string);
		void FillMassDatasets();

		void run(TFile*);

   private:
    float mass_;
    float bdtoutput_;
    float weight_;
    int category_;
    vector<pair<float,float> > massSyst_;
    vector<pair<float,float> > bdtouputSyst_;
    vector<pair<float,float> > weightSyst_;
    vector<pair<int,int> > categorySyst_;
		
		map<string,TH1F*> th1fs_;
		
		TFile *inFile_;
		TFile *outFile_;
		RooWorkspace *ws_;
		RooRealVar *massVar_;
		RooDataSet *dataSet_;
};

#endif
