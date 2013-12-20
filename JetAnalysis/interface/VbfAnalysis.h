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

class TreeVariables {

public:
    
    TreeVariables();

    int entry;
    int   nPU;
    float weight;
    float sampleweight;
    
    int   jet1, jet2, jet3;
    float diphomva;
    float pho1pt;
    float pho2pt;
    float diphopt;
    float diphoM;
    float diphoEta;
    float dijetEta;
    float jet1isMatched,jet2isMatched;
    float jet1genPt,jet2genPt;
    float jet1genDr,jet2genDr;
    float jet1Pt, jet2Pt, jet1Eta, jet2Eta, zepp, mj1j2, dphi, dphiJJ, dphiJJ2, deltaEta3;
    bool  jet1PileupID,jet2PileupID ;
    bool  isSignal;
    int   mctype;

    float pho1energy;
    float pho2energy;
    float pho1Phi;
    float pho2Phi;
    float pho1Eta;
    float pho2Eta;
    float pho1scEta;
    float pho2scEta;
    float pho1r9;
    float pho2r9;
    float pho1idMva;
    float pho2idMva;
    float pho1sE;
    float pho2sE;
    float pho1sEsmear;
    float pho2sEsmear;
    float diphosM;
    float diphosMwv;
    float diphovtxProb;
    bool pho1Matched, pho2Matched, corrVeretx;
    int pho1CiC, pho2CiC;
    float diphoMVA, vbfMVA, combinedMVA;
    
};

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

    
    bool recomputeJetId, expoMatching, dumpFlatTree, requireTwoJets;
    
    
private:
    TFile * outputFile_;
    TTree * flatTree_;
    int _nVertices;
    TreeVariables tree_, default_;

};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
