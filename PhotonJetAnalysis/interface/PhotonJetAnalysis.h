#ifndef __PHOTONJETANALYSIS__
#define __PHOTONJETANALYSIS__

#include "BaseAnalysis.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/PhotonInfo.h"
#include "PhotonAnalysis/interface/PhotonAnalysis.h"

#include "TriggerSelection.h"
#include "EnergySmearer.h"

#include "MassResolution.h"

#include "TMVA/Reader.h"
#include "PhotonFix.h"
#include <stdio.h>
// #include "HiggsToGammaGamma/interface/GBRForest.h"
//#include "../../../../HiggsToGammaGamma/interface/GBRForest.h"
//#include "HiggsAnalysis/HiggsToGammaGamma/interface/GBRForest.h"

// ------------------------------------------------------------------------------------
class PhotonJetAnalysis : public PhotonAnalysis 
{
 public:

    bool SkimEvents(LoopAll&, int);
    bool SelectEventsReduction(LoopAll&, int);
    void ReducedOutputTree(LoopAll &l, TTree *);
    
};

#endif
