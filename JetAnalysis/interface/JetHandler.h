#ifndef __JetHandler__
#define __JetHandler__

#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include <string>

// ------------------------------------------------------------------------------------
class JetHandler
{
public:
    
    JetHandler(const std::string & cfg, const std::string & pset);
    virtual ~JetHandler();
    
private:
    PileupJetIdAlgo * jetIdAlgo;
  
};

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
