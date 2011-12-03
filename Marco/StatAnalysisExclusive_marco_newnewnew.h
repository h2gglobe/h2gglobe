#ifndef __STATANALYSISEXCLUSIVE__
#define __STATANALYSISEXCLUSIVE__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class StatAnalysisExclusive : public PhotonAnalysis 
{
public:
	
	StatAnalysisExclusive();
	virtual ~StatAnalysisExclusive();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	void Init(LoopAll&);
	void Term(LoopAll&);
	
	void GetBranches(TTree *, std::set<TBranch *>& );
	
	virtual bool SelectEvents(LoopAll&, int);
	virtual void Analysis(LoopAll&, Int_t);
	
	// Options
  bool includeVBF;
  bool includeVHad;
  bool includeVHlep;
	bool reRunCiC;

  float leadEtCut;
  float leadEtVBFCut;
  float leadEtVHadCut;
  float leadEtVHlepCut;
  float subleadEtCut;
  float subleadEtVBFCut;
  float subleadEtVHadCut;
  float subleadEtVHlepCut;
  int nVBFEtaCategories;
	int nVHadEtaCategories;
  std::string efficiencyFile;
	
	// EnergySmearer::energySmearingParameters eSmearPars; // gone to PhotonAnalysis GF
	EfficiencySmearer::efficiencySmearingParameters effSmearPars;
	DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

	double GetDifferentialKfactor(double, int);

	void FillSignalLabelMap();
	std::string GetSignalLabel(int) ;

  // for N-1
  float  myVBFLeadJPt;
  float  myVBFSubJPt;
  float  myVBFdEta;
  float  myVBFZep;
  float  myVBFdPhi;
  float  myVBF_Mjj;

  float  myVBF_Mgg;

  float  myVHadLeadJPt;
  float  myVHadSubJPt;
  float  myVHaddEta;
  float  myVHadZep;
  float  myVHaddPhi;
  float  myVHad_Mjj;

  float  myVHad_Mgg;

  float  myAllLeadJPt;
  float  myAllSubJPt;
  float  myAllLeadJEta;
  float  myAllSubJEta;
  float  myAlldEta;
  float  myAllZep;
  float  myAlldPhi;
  float  myAll_Mjj;

  float  myAll_Mgg;
  float  myAllPtHiggs;
  float  myAll_phoet1;
  float  myAll_phoet2;
  float  myAll_phoetom1;
  float  myAll_phoetom2;

  float  myInclusive_Mgg;
  float  myInclusivePtHiggs;

  std::pair<int, int> Select2HighestPtJets(LoopAll&, TLorentzVector& leadpho, TLorentzVector& subleadpho, float jtLMinPt, float jtTMinPt);
  int RescaleJetEnergy(LoopAll&);

	bool  doMCSmearing;
	bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
	bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, doKFactorSmear;
	float systRange;
	int   nSystSteps;   
	int   nEtaCategories, nR9Categories, nPtCategories;
	float massMin, massMax;
	int nDataBins;	
	//float smearing_sigma_EBHighR9       ;
	//float smearing_sigma_EBLowR9        ;
	//float smearing_sigma_EEHighR9       ;
	//float smearing_sigma_EELowR9        ;
	//float smearing_sigma_error_EBHighR9 ;
	//float smearing_sigma_error_EBLowR9  ;
	//float smearing_sigma_error_EEHighR9 ;
	//float smearing_sigma_error_EELowR9  ;
	
	std::string kfacHist;

	TH1D *thm110,*thm120,*thm130,*thm140;
	int nMasses;

protected:
	std::vector<BaseSmearer *> photonSmearers_;
	std::vector<BaseSmearer *> systPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> diPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> systDiPhotonSmearers_;
	std::vector<BaseGenLevelSmearer *> genLevelSmearers_;
	std::vector<BaseGenLevelSmearer *> systGenLevelSmearers_;
	
	EnergySmearer /* *eScaleSmearer,*/ *eResolSmearer ; // moved to PhotonAnalysis GF 
	EfficiencySmearer *idEffSmearer, *r9Smearer;
	DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer;
	KFactorSmearer * kFactorSmearer;
	
	std::string name_;
	std::map<int,std::string> signalLabels;
	float nevents, sumwei, sumaccept, sumsmear, sumev; 
	
	int nCategories_;
	int nPhotonCategories_;
	int diPhoCounter_;
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	// RooStuff
	RooContainer *rooContainer;

	ofstream eventListText;
	ofstream eventListTextVBF;
	//vector<double> weights;
	TFile *kfacFile;
	
};

#endif
