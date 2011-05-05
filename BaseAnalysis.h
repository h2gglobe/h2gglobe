#ifndef __BASEANALYSIS__
#define __BASEANALYSIS__

#include "LoopAll.h"

// ------------------------------------------------------------------------------------
class BaseAnalysis {
public:
  BaseAnalysis();
  virtual ~BaseAnalysis();
  
  virtual const std::string & name() const = 0; 
  
  operator const std::string & () const { return this->name(); };
  
  virtual void ReducedOutputTree(TTree *) = 0;
  virtual void GetBranches(TTree *, std::set<TBranch *>& ) = 0;
  
  virtual void Init(LoopAll&) = 0;
  virtual void Term(LoopAll&) = 0;

  virtual void FillReductionVariables(LoopAll&, int) = 0;

  virtual bool SelectEventsReduction(LoopAll&, int) = 0;
  virtual bool SelectEvents(LoopAll&, int) = 0;
  
  virtual void Analysis(LoopAll&, Int_t) = 0;
  
};

bool operator == (BaseAnalysis * a, const std::string & b);

#endif
