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
#include "KFactorSmearer.h"
#include "InterferenceSmearer.h"
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
    virtual void ResetAnalysis();
    virtual void Analysis(LoopAll&, Int_t);
    
    
    std::string efficiencyFile;

    // mva removed cp march 8
    //bool useMVA;
    //std::string phoIDMVAtype;
    //float phoIDMVAloose;
    //float phoIDMVAtight;
    //int nDiphoEventClasses;

    // EnergySmearer::energySmearingParameters eSmearPars; // gone to PhotonAnalysis GF
    EfficiencySmearer::efficiencySmearingParameters effSmearPars;
    DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

    double GetDifferentialKfactor(double, int);

    void FillSignalLabelMap();
    std::string GetSignalLabel(int) ;

    bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
    bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, doKFactorSmear, doInterferenceSmear;
    float systRange;
    int   nSystSteps;   
    //int   nEtaCategories, nR9Categories, nPtCategories;
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
    bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
		      bool & isCorrectVertex,
		      bool isSyst=false, 
		      float syst_shift=0., bool skipSelection=false,
		      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 
    bool VHmuevent, VHelevent, VBFevent, VHhadevent;
    double genLevWeight; 

    std::vector<float> smeared_pho_energy;
    std::vector<float> smeared_pho_r9;
    std::vector<float> smeared_pho_weight;

    void  computeExclusiveCategory(LoopAll & l, int & category, std::pair<int,int> diphoton_index, float pt);	

    void fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs, float lead_r9, float sublead_r9, 
			  int category, float evweight , LoopAll &);
	
    EnergySmearer /* *eScaleSmearer,*/ *eResolSmearer ; // moved to PhotonAnalysis GF 
    EfficiencySmearer *idEffSmearer, *r9Smearer;
    DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer;
    KFactorSmearer * kFactorSmearer;
    InterferenceSmearer * interferenceSmearer;
    
    std::string name_;
    std::map<int,std::string> signalLabels;
    float nevents, sumwei, sumaccept, sumsmear, sumev; 
    
    int nInclusiveCategories_;
    int nCategories_;
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


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
