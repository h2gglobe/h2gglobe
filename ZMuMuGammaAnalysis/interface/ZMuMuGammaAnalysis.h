#ifndef __ZMUMUGAMMAANALYSIS__
#define __ZMUMUGAMMAANALYSIS__

#include "PhotonAnalysis/interface/PhotonAnalysis.h"

// ------------------------------------------------------------------------------------
class ZMuMuGammaAnalysis : public PhotonAnalysis 
{
 public:

  bool SkimEvents(LoopAll&, int);
  bool SelectEventsReduction(LoopAll&, int);
  void ReducedOutputTree(LoopAll &l, TTree *);

    ZMuMuGammaAnalysis();
    virtual ~ZMuMuGammaAnalysis();

};

#endif
