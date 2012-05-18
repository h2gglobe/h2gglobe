#ifndef __JetHandler__
#define __JetHandler__

#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include <string>

class LoopAll;

// ------------------------------------------------------------------------------------
class JetHandler
{
public:
    
    JetHandler(const std::string & cfg, const std::string & pset, LoopAll &l);
    
    void compute_beta(int ijet, int ivx=-1);
    void compute_mva(int ijet, int ivtx=-1);
    void  compute_wp(int ijet, int ivtx=-1);
    
    void recompute_jec();
    
    virtual ~JetHandler();
    
private:
    LoopAll & l_;
    PileupJetIdAlgo * cutbased, * simple, * full;
    
};

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
