#ifndef __JetHandler__
#define __JetHandler__

#include "FWCore/Utilities/interface/GCC11Compatibility.h"

#ifndef noexcept
#define noexcept
#endif

#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "../interface/JetResponseChange.h"

#include <string>

class LoopAll;

// ------------------------------------------------------------------------------------
class JetHandler
{
public:
    
    JetHandler(const std::string & cfg, LoopAll &l);
    
    void computeBetas(int ijet, int ivx=-1);
    void fillFromJet(int ijet, int ivtx=-1);
    void computeMva(int ijet, int ivtx=-1);
    void computeWp(int ijet, int ivtx=-1);
    
    void recomputeJec(int ijet, bool correct=false);
    void emulateJetResponse(int ijet);
    void setJetResponseStep(float x) { lumiStep_ = x; };
    void applyJecUncertainty(int ijet, float shift);
    void applyJerUncertainty(int ijet, float shift);

    void bookFlatTree(TTree * tree);
    
    virtual ~JetHandler();
    
private:
    LoopAll & l_;
    PileupJetIdAlgo * cutbased, * simple, * full;
    PileupJetIdentifier internalId_;
    
    FactorizedJetCorrector *jecCorData_, *jecCorMc_;
    JetCorrectionUncertainty *jecUnc_;
    std::vector<JetCorrectorParameters> jetCorParsData_, jetCorParsMc_;

    std::vector<double> jerEtaBins_, jerResSf_, jerResSfErr_;
    
    TTree * flatTree;
    boost::shared_ptr<edm::ParameterSet> myPset;

    int nLumiBins_;
    float lumiStep_;
    JetResponseChange *jetResponseChange;
 
};

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
