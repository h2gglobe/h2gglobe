#ifndef __JetHandler__
#define __JetHandler__

#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include <string>

class LoopAll;

// ------------------------------------------------------------------------------------
class JetHandler
{
public:
    
    JetHandler(const std::string & cfg, LoopAll &l);
    
    void compute_betas(int ijet, int ivx=-1);
    void compute_mva(int ijet, int ivtx=-1);
    void  compute_wp(int ijet, int ivtx=-1);
    
    void recompute_jec();
    
    virtual ~JetHandler();
    
private:
    LoopAll & l_;
    PileupJetIdAlgo * cutbased, * simple, * full;

    FactorizedJetCorrector *jecCorData_, *jecCorMc_;
    std::vector<JetCorrectorParameters> jetCorParsData_, jetCorParsMc_;
    
    
};

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
