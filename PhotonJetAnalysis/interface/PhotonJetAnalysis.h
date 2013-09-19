#ifndef __PHOTONJETANALYSIS__
#define __PHOTONJETANALYSIS__

#include "PhotonAnalysis/interface/PhotonAnalysis.h"

// ------------------------------------------------------------------------------------
class PhotonJetAnalysis : public PhotonAnalysis 
{
 public:

    bool SkimEvents(LoopAll&, int);
    bool SelectEventsReduction(LoopAll&, int);

    PhotonJetAnalysis();
    virtual ~PhotonJetAnalysis();

};

#endif
