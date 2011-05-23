#ifndef __PHOTONANALYSIS__
#define __PHOTONANALYSIS__

#include "BaseAnalysis.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/PhotonInfo.h"

#include "TriggerSelection.h"

// ------------------------------------------------------------------------------------
class PhotonAnalysis : public BaseAnalysis 
{
public:
	
	PhotonAnalysis();
	virtual ~PhotonAnalysis();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	virtual void Init(LoopAll&);
	virtual void Term(LoopAll&);
	
	virtual void ReducedOutputTree(LoopAll &l, TTree *);
	virtual void GetBranches(TTree *, std::set<TBranch *>& );
	
	virtual void FillReductionVariables(LoopAll& l, int jentry);   
	virtual bool SelectEventsReduction(LoopAll&, int);

	virtual bool SkimEvents(LoopAll&, int);
	virtual bool SelectEvents(LoopAll&, int);
	virtual void Analysis(LoopAll&, Int_t);
	
	// Public parameters to be read from config file
	VertexAlgoParameters vtxAlgoParams;	 
	std::vector<std::string> vtxVarNames;
	
	bool doTriggerSelection; 
	std::vector<TriggerSelection> triggerSelections;
	
	// Preselection indexes
	float presel_scet1, presel_scet2, presel_maxeta;
	float presel_ecaliso_eb, presel_ecaliso_ee, presel_sieie_eb, presel_sieie_ee, presel_hoe;

	std::vector<int> pho_acc;
	std::vector<int> pho_presel;
	std::vector<int> pho_presel_lead;
	std::vector<float> pho_sc_et;
	// Other options
	bool runStatAnalysis;
        TString puHist;//name of pileup reweighting histogram

	
protected:
	void PreselectPhotons(LoopAll& l, int jentry);
	void StatAnalysis(LoopAll &l, int jentry);
	
	std::string name_;
	
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	vector<double> weights;
	int trigCounter_;
	
};

#endif
