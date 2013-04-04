#ifndef _MyAnalysis_h
#define _MyAnalysis_h

#include "GenericAnalysis.h"
#include "../RooContainer.h"
#include "parser.h"

#include "TTree.h"

/** base class from which user analyses must inherit */
class MyAnalysis : public GenericAnalysis {
 public:
  /** set branch addresses for the given tree
      and set active / inactivate the branches */
  void setBranchAddresses(TTree *tree, std::map<std::string, std::string>& values);

  void defineLabels();
  std::string GetSignalLabel(int id);

  /** this is called for each event in the input tree
      and is typically the place to fill histograms
      using container->Fill(..) and container->Fill2D(..)
  */
  void analyze(RooContainer* rooContainer);

  void FillRooContainer(RooContainer* rooContainer);
  void FillRooContainerSyst(RooContainer* rooContainer);

  // List of opttree branches

  Float_t run;
  Float_t lumis;
  Double_t event;
  Float_t itype;
  Float_t et1;
  Float_t et2;
  Int_t muid1;
  Int_t muid2;
  Float_t eta1;
  Float_t eta2;
  Float_t mass;
  Int_t cat;
  Float_t full_weight;

  std::string* name1;

  std::map<int,std::string> signalLabels;
};

/*

  a shared object with a user analysis must also contain the following
  function:

  GenericAnalysis *makeInstance();

  This must create an instance of the user analysis class and return it.

*/
#endif
