#ifndef __MASSFACTORIZEDANALYSIS__
#define __MASSFACTORIZEDANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "StatAnalysis.h"
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
class MassFactorizedMvaAnalysis : public StatAnalysis 
{
 public:
    
    MassFactorizedMvaAnalysis();
    virtual ~MassFactorizedMvaAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    virtual void ResetAnalysis();
    //// virtual void Analysis(LoopAll&, Int_t); 

    virtual int GetBDTBoundaryCategory(float,bool,bool);

    void fillZeeControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
	                     const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
			     float phoid_mvaout_lead, float phoid_mvaout_sublead, 
			     float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
			     int diphoton_id, int category, int selectioncategory, float evweight, LoopAll & l );

    void fillEscaleTree(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
	                     const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
			     float phoid_mvaout_lead, float phoid_mvaout_sublead, 
			     float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
			     int diphoton_id, int category, int selectioncategory, float evweight, LoopAll & l );

    void fillTrainTree(LoopAll&, Int_t, Int_t, Int_t, float, TLorentzVector &, TLorentzVector &, float,float,float,const char*,float photonID_1=-50.,float photonID_2=-50.);

    void fillOptTree(LoopAll& l, const TLorentzVector & lead_p4, const TLorentzVector & sublead_p4, Float_t vtxProb,
		     std::pair<int, int> diphoton_index, Int_t diphoton_id, Float_t phoid_mvaout_lead, Float_t phoid_mvaout_sublead,
		     Float_t weight, Float_t mass, Float_t sigmaMrv, Float_t sigmaMwv,
		     const TLorentzVector & Higgs, Float_t diphobdt_output, Int_t category, bool VBFevent, Float_t myVBF_Mjj, Float_t myVBFLeadJPt, 
		     Float_t myVBFSubJPt, Int_t nVBFDijetJetCategories);

    bool doPhotonMvaIdSyst;
    bool doPhotonMvaIdSmear;
    bool doRegressionSmear, doRegressionSyst;
    bool doPdfWeightSmear, doPdfWeightSyst;
    bool fillEscaleTrees;
    bool fillOptree;
    bool forceStdPlotsOnZee;

    std::string bdtTrainingPhilosophy;
    std::string photonLevelMvaUCSD  ;
    std::string eventLevelMvaUCSD   ;                    
    std::string photonLevelMvaMIT_EB;
    std::string photonLevelMvaMIT_EE;
    std::string eventLevelMvaMIT    ;
    std::string photonLevelNewIDMVA_EB;
    std::string photonLevelNewIDMVA_EE;

    std::vector<float> bdtCategoryBoundaries;
    float photonIDMVAShift_EB, photonIDMVAShift_EE;
    
 protected:

    virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
		      bool & isCorrectVertex, float &kinematic_bdtout,
		      bool isSyst=false, 
		      float syst_shift=0., bool skipSelection=false,
		      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 

    EnergySmearer  *eRegressionSmearer ; 
    DiPhoEfficiencySmearer *photonMvaIdSmearer ;
    
    std::string name_;
    std::map<int,std::string> signalLabels;
    
    
    //HggVertexAnalyzer vtxAna_;
    //HggVertexFromConversions vtxConv_;
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
