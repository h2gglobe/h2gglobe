#ifndef __VertexOptimizationAnalysis__
#define __VertexOptimizationAnalysis__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis/interface/PhotonAnalysis.h"
#include "PhotonAnalysis/interface/StatAnalysis.h"
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
class VertexOptimizationAnalysis : public StatAnalysis 
{
 public:
    
    VertexOptimizationAnalysis();
    virtual ~VertexOptimizationAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    void ReducedOutputTree(LoopAll &l, TTree * outputTree);

    virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, 
			      float & mass, float & evweight, int & category, int & diphoton_id,
			      bool & isCorrectVertex, float &kinematic_bdtout,
			      bool isSyst=false, 
			      float syst_shift=0., bool skipSelection=false,
			      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 
    
    TString minBiasRefName;
    int storeNVert;

private:
    std::vector<std::string> vtxVarNames_;
    std::vector<float> vtxVars_;
    
    TFile * uFile_;
    TTree * uTree_;
    TH1 * hMinBiasSpecturm_, *hHiggsSpecturm_, *hMinBiasRef_;
    bool isClosestToGen_, passCiC_;
    int   nPU_, nVert_, itype_;
    float evWeight_, ksprob_;
    TLorentzVector *pho1_;
    TLorentzVector *pho2_;
    TLorentzVector *dipho_;
    
    TTree *evTree_;
    Float_t dZTrue_, zRMS_, zTrue_;
    Int_t category_, nConv_;
    Float_t mTrue_, mTrueVtx_;
    vector<float> MVA_;
    vector<float> dZ_;
    vector<float> diphoM_;
    vector<float> diphoPt_;
    
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
