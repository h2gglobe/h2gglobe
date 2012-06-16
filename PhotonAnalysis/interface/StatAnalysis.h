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
    virtual bool Analysis(LoopAll&, Int_t);
    
    std::string efficiencyFile;

    EfficiencySmearer::efficiencySmearingParameters effSmearPars;
    DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

    double GetDifferentialKfactor(double, int);

    void FillSignalLabelMap(LoopAll &l);
    std::string GetSignalLabel(int) ;

    bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
    bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, 
	doKFactorSmear, doInterferenceSmear;
    float systRange;
    int   nSystSteps;   
    //int   nEtaCategories, nR9Categories, nPtCategories;
    std::vector<int> cicCutLevels;
    float massMin, massMax;
    int nDataBins;  
    bool dataIs2011;
    bool scaleClusterShapes, scaleR9Only;
    bool dumpAscii, dumpMcAscii;
    float phoidMvaCut;
    std::vector<double> zeePtBinLowEdge, zeePtWeight;

    std::string kfacHist;

    TH1D *thm110,*thm120,*thm130,*thm140;

 protected:
    // Factorized Analysis method + loop over systematics
    virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, 
			      int & category, int & diphoton_id,
			      bool & isCorrectVertex,float &kinematic_bdtout,
			      bool isSyst=false, 
			      float syst_shift=0., bool skipSelection=false,
			      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 

    virtual void FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA, int category, float weight, 
				  bool isCorrectVertex);
    virtual void AccumulateSyst(int cur_type, float mass, float diphotonMVA, int category, float weight,
				std::vector<double> & mass_errors,
				std::vector<double> & mva_errors,
				std::vector<int>    & categories,
				std::vector<double> & weights);
    virtual void FillRooContainerSyst(LoopAll& l, const std::string & name,int cur_type,
				      std::vector<double> & mass_errors, std::vector<double> & mva_errors,
				      std::vector<int>    & categories, std::vector<double> & weights);
    
    bool VHmuevent, VHelevent, VBFevent, VHhadevent, VHmetevent;  //met at analysis step
    double genLevWeight; 

    std::vector<float> smeared_pho_energy;
    std::vector<float> smeared_pho_r9;
    std::vector<float> smeared_pho_weight;

    void  computeExclusiveCategory(LoopAll & l, int & category, std::pair<int,int> diphoton_index, float pt);	

    void fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs, 
			  float lead_r9, float sublead_r9, int diphoton_index, 
			  int category, bool rightvtx, float evweight , LoopAll &);

    void fillSignalEfficiencyPlots(float weight, LoopAll & l );

    void rescaleClusterVariables(LoopAll &l);

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
