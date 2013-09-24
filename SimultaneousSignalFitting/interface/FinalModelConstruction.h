#ifndef FinalModelConstruction_h 
#define FinalModelConstruction_h

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "TGraph.h"

#include "RooAbsReal.h"
#include "RooAddition.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"
#include "../interface/Normalization_8TeV.h"

class FinalModelConstruction {

  public:
    
    FinalModelConstruction(RooRealVar *massVar, RooRealVar *MHvar, RooRealVar *intL, int mhLow, int mhHigh, std::string proc, int cat, bool doSecMods, std::string systematicsFileName, int verbosity, bool isCB=false);
    ~FinalModelConstruction();

		void loadSignalSystematics(std::string filename);
		void printSignalSystematics();

    void setSecondaryModelVars(RooRealVar *mh_sm, RooRealVar *deltam, RooAddition *mh_2, RooRealVar *width);

    void buildRvWvPdf(std::string name, int nGrv, int nGwv, bool recusive);
    void buildStdPdf(std::string name, int nGaussians, bool recursive);
    std::vector<RooAddPdf*> buildPdf(std::string name, int nGaussians, bool recursive, std::map<std::string,RooSpline1D*> splines, string add="");
    void getRvFractionFunc(std::string name);
    void setupSystematics();
    void getNormalization();

		RooAbsReal *getMeanWithPhotonSyst(RooAbsReal *dm, string name);
		RooAbsReal *getSigmaWithPhotonSyst(RooAbsReal *sig_fit, string name);
		RooAbsReal *getRateWithPhotonSyst(string name);
    
		void setRVsplines(std::map<std::string,RooSpline1D*> splines);
    void setWVsplines(std::map<std::string,RooSpline1D*> splines);
    void setSTDsplines(std::map<std::string,RooSpline1D*> splines);

    void setRVdatasets(std::map<int,RooDataSet*> data);
    void setWVdatasets(std::map<int,RooDataSet*> data);
    void setSTDdatasets(std::map<int,RooDataSet*> data);
		void makeSTDdatasets();

    void plotPdf(std::string outDir);

    void save(RooWorkspace *work);

  private:
    
    RooRealVar *mass;
    RooRealVar *MH;
    RooRealVar *intLumi;
    int mhLow_;
    int mhHigh_;
    std::string proc_;
    int cat_;
    int nIncCats_;
    bool doSecondaryModels;
    bool secondaryModelVarsSet;
    bool isCutBased;
    std::vector<int> allMH_;
    std::vector<int> getAllMH();
    int verbosity_;
    Normalization_8TeV *norm;

    std::map<std::string,RooSpline1D*> stdSplines;
    std::map<std::string,RooSpline1D*> rvSplines;
    std::map<std::string,RooSpline1D*> wvSplines;
    std::map<int,RooDataSet*> rvDatasets;
    std::map<int,RooDataSet*> wvDatasets;
    std::map<int,RooDataSet*> stdDatasets;
    
    RooRealVar *vertexNuisance;
    RooSpline1D *rvFracFunc;
    RooRealVar *globalScale;
    RooRealVar *categoryScale;
    RooConstVar *categorySmear;
    RooRealVar *categoryResolution;
    vector<double> constSmearVals;

    RooAddPdf *finalPdf;
    RooAbsReal *finalNorm;
    RooAbsReal *finalNormThisLum;
    RooExtendPdf *extendPdfRel;
    RooExtendPdf *extendPdf;
    // secondary models
    RooAddPdf *finalPdf_SM;
    RooAddPdf *finalPdf_2;
    RooAddPdf *finalPdf_NW;
    RooAbsReal *finalNorm_SM;
    RooAbsReal *finalNorm_2;
    RooAbsReal *finalNorm_NW;

    bool systematicsSet_;
    bool rvFractionSet_;

    RooSpline1D *graphToSpline(string name, TGraph *graph);
    RooSpline1D *graphToSpline(string name, TGraph *graph, RooAbsReal *var);

    std::map<std::string,RooSpline1D*> xsSplines;
    RooSpline1D *brSpline;
    // secondary models
    std::map<std::string,RooSpline1D*> xsSplines_SM;
    RooSpline1D *brSpline_SM;
    std::map<std::string,RooSpline1D*> xsSplines_2;
    RooSpline1D *brSpline_2;
    std::map<std::string,RooSpline1D*> xsSplines_NW;
    RooSpline1D *brSpline_NW;

    RooRealVar *MH_SM;
    RooRealVar *DeltaM;
    RooAddition *MH_2;
    RooRealVar *higgsDecayWidth;

		// photon systematic stuff
		std::vector<std::string> photonCats;
		std::map<std::string,std::map<int,std::map<std::string,double> > > meanSysts;
		std::map<std::string,std::map<int,std::map<std::string,double> > > sigmaSysts;
		std::map<std::string,std::map<int,std::map<std::string,double> > > rateSysts;

		std::map<string,RooRealVar*> photonSystematics;
		std::map<string,RooRealVar*> photonSystematicConsts;

		// utility funcs
		void stripSpace(std::string &line);
		void printVec(std::vector<std::string> vec);
		void printSystMap(std::map<std::string,std::map<int,std::map<std::string,double> > > &theMap);
		void addToSystMap(std::map<std::string,std::map<int,std::map<std::string,double> > > &theMap, string proc, int diphotonCat, string phoSystName, double var);

};

#endif

