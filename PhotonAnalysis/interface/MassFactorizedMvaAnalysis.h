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
   
    virtual int GetBDTBoundaryCategory(float,bool,bool);

    virtual void ResetAnalysis();
    //// virtual void Analysis(LoopAll&, Int_t); 

    void fillZeeControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
	                     const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
			     float phoid_mvaout_lead, float phoid_mvaout_sublead, 
			     float diphobdt_output_up, float diphobdt_output_down,
			     float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
			     int diphoton_id, int category, int selectioncategory, float evweight, LoopAll & l );

    void fillEscaleTree(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
	                     const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
			     float phoid_mvaout_lead, float phoid_mvaout_sublead, 
			     float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
			     int diphoton_id, int category, int selectioncategory, float evweight, LoopAll & l );

    void fillTrainTree(LoopAll&, Int_t, Int_t, Int_t, float, TLorentzVector &, TLorentzVector &, float,float,float,const char*,float photonID_1=-50.,float photonID_2=-50.);

    bool PreselDiphoFillDiphoMVA(LoopAll &l, float *pho_energy_array=0, bool smear=false, float syst_shift=0);
    float GetDiphoMva(LoopAll & l, int diphotonId, bool smear=false, float syst_shift=0);
    void ComputeDiphoMvaInputs(LoopAll &l, float &phoid_mvaout_lead, float &phoid_mvaout_sublead, float &vtxProb, int diphoton_id);

    bool doPhotonMvaIdSyst;
    bool doPhotonMvaIdSmear;
    bool doRegressionSmear, doRegressionSyst;
    bool doPdfWeightSmear, doPdfWeightSyst;
    bool fillEscaleTrees;
    bool fillOptree;
    bool forceStdPlotsOnZee;

    bool applyIdmvaCorrection;
    bool applySigmaECorrection;

    TFile* idmvascaleFile;
    TFile* sigmaescaleFile;
    
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
