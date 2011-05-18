#ifndef __STATANALYSIS__
#define __STATANALYSIS__

#include "BaseAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

// ------------------------------------------------------------------------------------
class StatAnalysis : public BaseAnalysis 
{
public:
	
	StatAnalysis();
	virtual ~StatAnalysis();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	void Init(LoopAll&);
	void Term(LoopAll&);
	
	void ReducedOutputTree(LoopAll &l, TTree *);
	void GetBranches(TTree *, std::set<TBranch *>& );
	
	void FillReductionVariables(LoopAll& l, int jentry);   
	bool SelectEventsReduction(LoopAll&, int);
	
	virtual bool SkimEvents(LoopAll&, int);
	virtual bool SelectEvents(LoopAll&, int);
	virtual void Analysis(LoopAll&, Int_t);
	
	// Public parameters to be read from config file
	VertexAlgoParameters vtxAlgoParams;	 
	std::vector<std::string> vtxVarNames;
	
	// Preselection indexes
	float presel_scet1, presel_scet2, presel_maxeta;
	float presel_ecaliso_eb, presel_ecaliso_ee, presel_sieie_eb, presel_sieie_ee, presel_hoe;

	std::vector<int> pho_acc;
	std::vector<int> pho_presel;
	std::vector<int> pho_presel_lead;
	std::vector<float> pho_sc_et;
	
	// Other options
	bool runStatAnalysis;
	
protected:
	void PreselectPhotons(LoopAll& l, int jentry);
	
	std::string name_;
	
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	// RooStuff
	RooContainer *rooContainer;
	
};

#endif
