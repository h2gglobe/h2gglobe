#ifndef __MASSFACTORIZEDANALYSIS__
#define __MASSFACTORIZEDANALYSIS__

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
class MassFactorizedMvaAnalysis : public PhotonAnalysis 
{
public:
	
	MassFactorizedMvaAnalysis();
	virtual ~MassFactorizedMvaAnalysis();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	void Init(LoopAll&);
	void Term(LoopAll&);
	
	void GetBranches(TTree *, std::set<TBranch *>& );
	
	virtual bool SelectEvents(LoopAll&, int);
	virtual void ResetAnalysis();
	virtual void Analysis(LoopAll&, Int_t);
	
	// Options
  	bool includeVBF;
	bool reRunCiCForData;
	float leadEtCut;
	float subleadEtCut;
	std::string efficiencyFile;
	
	// EnergySmearer::energySmearingParameters eSmearPars; // gone to PhotonAnalysis GF
	EfficiencySmearer::efficiencySmearingParameters effSmearPars;
	DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

	double GetDifferentialKfactor(double, int);

	void FillSignalLabelMap();
	int GetBDTBoundaryCategory(float,bool,bool);
	std::string GetSignalLabel(int) ;

	bool  doMCSmearing;
	bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst, doPhotonMvaIdSyst;
	bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, doKFactorSmear, doPhotonMvaIdSmear;
	bool doRegressionSmear, doRegressionSyst;
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

	std::string bdtTrainingPhilosophy;
	std::string photonLevelMvaUCSD  ;
	std::string eventLevelMvaUCSD   ;                    
	std::string photonLevelMvaMIT_EB;
	std::string photonLevelMvaMIT_EE;
	std::string eventLevelMvaMIT    ;

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
	
	EnergySmearer /* *eScaleSmearer,*/ *eResolSmearer, *eRegressionSmearer ; // moved to PhotonAnalysis GF 
	EfficiencySmearer *idEffSmearer, *r9Smearer;
	DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer,*photonMvaIdSmearer ;
	KFactorSmearer * kFactorSmearer;
	
	std::string name_;
	std::map<int,std::string> signalLabels;
	float nevents, sumwei, sumaccept, sumsmear, sumev; 
	
	int nPhotonCategories_;
	int diPhoCounter_;
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	// RooStuff
	RooContainer *rooContainer;

	ofstream eventListText;
	//vector<double> weights;
	TFile *kfacFile;
	
};

#endif
