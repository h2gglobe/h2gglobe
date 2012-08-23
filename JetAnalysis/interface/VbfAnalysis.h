#ifndef __VbfAnalysis__
#define __VbfAnalysis__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis/interface/MassFactorizedMvaAnalysis.h"
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
class VbfAnalysis : public MassFactorizedMvaAnalysis 
{
 public:
    
    VbfAnalysis();
    virtual ~VbfAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    void FillReductionVariables(LoopAll& l, int jentry);   
    void ReducedOutputTree(LoopAll &l, TTree * outputTree);
    bool SelectEventsReduction(LoopAll&, int);

    bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, 
		      int & diphoton_id,
		      bool & isCorrectVertex,
		      float &kinematic_bdtout,
		      bool isSyst=false, 
		      float syst_shift=0., bool skipSelection=false,
		      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 

    
    bool recomputeJetId, expoMatching, dumpFlatTree;
    
    
private:
    TFile * outputFile_;
    TTree * flatTree_;
    int tree_entry;
    int   tree_nPU,_nVertices;
    float tree_evWeight;

    float tree_diphomva;
    float tree_pho1pt;
    float tree_pho2pt;
    float tree_diphopt;
    float tree_diphoM;
    float tree_diphoEta;
    float tree_dijetEta;
    float tree_jet1isMatched,tree_jet2isMatched;
    float tree_jet1genPt,tree_jet2genPt;
    float tree_jet1genDr,tree_jet2genDr;
    float tree_jet1pt, tree_jet2pt, tree_jet1eta, tree_jet2eta, tree_zepp, tree_mj1j2, tree_dphi, tree_dphiJJ, tree_dphiJJ2, tree_deltaEta3;    
    bool  tree_jet1PileupID,tree_jet2PileupID ;
    bool  tree_isSignal;
    int   tree_mctype;
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
