#ifndef __STATANALYSIS__
#define __STATANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"

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
	EnergySmearer::energySmearingParameters eSmearPars;
	bool doEscaleSyst, doEresolSyst;
	float systRange;
	int nSystSteps;   
	int nCategories;

protected:
	void PreselectPhotons(LoopAll& l, int jentry);

	std::vector<BaseSmearer *> photonSmearers_;
	std::vector<BaseSmearer *> systPhotonSmearers_;
	std::vector<BaseSmearer *> diPhotonSmearers_;
	
	EnergySmearer *eScaleSmearer, *eResolSmearer;
	
	std::string name_;
	
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	// RooStuff
	RooContainer *rooContainer;

	//vector<double> weights;
	
};

#endif
