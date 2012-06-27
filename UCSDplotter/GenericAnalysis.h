#ifndef _GenericAnalysis_h
#define _GenericAnalysis_h


#include "HistoContainer.h"

class GenericAnalysis
{
public:
  /** set branch addresses for the given tree
      and set active / inactivate the branches */
  virtual void setBranchAddresses(TTree *tree) = 0;

  virtual void analyze(HistoContainer *container) = 0;

};

// GenericAnalysis *makeInstance();
#endif
