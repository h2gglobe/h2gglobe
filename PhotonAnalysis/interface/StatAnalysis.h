#ifndef __STATANALYSIS__
#define __STATANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class StatAnalysis : public PhotonAnalysis 
{
public:
	
	StatAnalysis();
	virtual ~StatAnalysis();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	void Init(LoopAll&);
	void Term(LoopAll&);
	
	void GetBranches(TTree *, std::set<TBranch *>& );
	
	virtual bool SelectEvents(LoopAll&, int);
	virtual void Analysis(LoopAll&, Int_t);
	
	// Options
	float leadEtCut;
	float subleadEtCut;
	std::string efficiencyFile;
	
	EnergySmearer::energySmearingParameters eSmearPars, eSmearDataPars;
	EfficiencySmearer::efficiencySmearingParameters effSmearPars;
	DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

	double GetDifferentialKfactor(double, int);
	bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doR9Syst, doVtxEffSyst, doTriggerEffSyst;
	float systRange;
	int   nSystSteps;   
	int   nEtaCategories, nR9Categories, nPtCategories;

	TString kfacHist;

	TH1D *thm110,*thm120,*thm130,*thm140;
	int nMasses;

protected:
	void PreselectPhotons(LoopAll& l, int jentry);

	std::vector<BaseSmearer *> photonSmearers_;
	std::vector<BaseSmearer *> systPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> diPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> systDiPhotonSmearers_;
	
	EnergySmearer *eScaleSmearer, *eScaleDataSmearer, *eResolSmearer ;
	EfficiencySmearer *idEffSmearer, *r9Smearer;
	DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer;
	
	std::string name_;
	
	int nCategories_;
	int diPhoCounter_;
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	// RooStuff
	RooContainer *rooContainer;

	// Saclay Text
	ofstream SaclayText;
	//vector<double> weights;
	TFile *kfacFile;
	
};

#endif
