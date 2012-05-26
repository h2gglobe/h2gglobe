#ifndef __EMPTYANALYSIS__
#define __EMPTYANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "StatAnalysis.h"
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
class EmptyAnalysis : public StatAnalysis 
{
 public:
    
    EmptyAnalysis();
    virtual ~EmptyAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    virtual void FillReductionVariables(LoopAll& l, int jentry);   
    virtual bool SelectEventsReduction(LoopAll&, int);
    virtual void ReducedOutputTree(LoopAll &l, TTree * outputTree); 

    virtual bool SkimEvents(LoopAll&, int);
    virtual bool SelectEvents(LoopAll&, int);
    virtual bool Analysis(LoopAll&, Int_t);

};

#endif


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
