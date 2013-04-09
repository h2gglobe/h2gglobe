#ifndef _NewCategoriesMVA_h
#define _NewCategoriesMVA_h

#include "GenericAnalysis.h"
#include "../RooContainer.h"
#include "parser.h"

#include "TTree.h"

/** base class from which user analyses must inherit */
class NewCategoriesMVA : public GenericAnalysis {
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

  Float_t itype;
  Float_t full_weight;
  Float_t mass;
  Float_t bdt;
  Float_t bdt_cat;
  
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
