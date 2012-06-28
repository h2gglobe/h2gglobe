#ifndef _JimAnalysis_h
#define _JimAnalysis_h

#include "GenericAnalysis.h"
#include "HistoContainer.h"

#include "TTree.h"

/** base class from which user analyses must inherit */
class JimAnalysis : public GenericAnalysis {
 public:
  /** set branch addresses for the given tree
      and set active / inactivate the branches */
  void setBranchAddresses(TTree *tree);

  /** this is called for each event in the input tree
      and is typically the place to fill histograms
      using container->Fill(..) and container->Fill2D(..)
  */
  void analyze(HistoContainer *container) = 0;


  // List of opttree branches

  Int_t run;
  Int_t lumis;
  Int_t event;
  Int_t itype;
  Float_t sigmaMrvoM;
  Float_t sigmaMwvoM;
  Float_t ptoM1;
  Float_t ptoM2;
  Float_t vtxprob;
  Float_t eta1;
  Float_t eta2;
  Float_t cosphi;
  Float_t idmva1;
  Float_t idmva2;
  Float_t rho;
  Float_t xsec_weight;
  Float_t pu_weight;
  Float_t full_weight;
  Int_t pu_n;
  Float_t mass;
  Int_t vbfcat;
  Float_t MET;
  Float_t MET_phi;
  Float_t bdt;
  Int_t bdt_cat;


};

/*

  a shared object with a user analysis must also contain the following
  function:

  GenericAnalysis *makeInstance();

  This must create an instance of the user analysis class and return it.

*/
#endif
