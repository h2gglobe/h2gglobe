#ifndef __JetIdAnalysis__
#define __JetIdAnalysis__

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
class JetIdAnalysis : public StatAnalysis 
{
 public:
    
    JetIdAnalysis();
    virtual ~JetIdAnalysis();
    
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

    void fillJetIdPlots(LoopAll & l, int ijet, int icat, float wei, const char * label, bool kinOnly=true);
    
    bool recomputeJetId, expoMatching, dumpFlatTree;
    
    
private:
    std::vector<std::string> vtxVarNames_;
    std::vector<float> vtxVars_;

    TFile * outputFile_;
    TTree * flatTree_;
    TH1 * hMinBiasSpecturm_, *hHiggsSpecturm_, *hMinBiasRef_;
    bool isClosestToGen_, passCiC_;
    int   nPU_, nVert_;
    float evWeight_, ksprob_;
    TLorentzVector *pho1_;
    TLorentzVector *pho2_;
    TLorentzVector *dipho_;
    
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
