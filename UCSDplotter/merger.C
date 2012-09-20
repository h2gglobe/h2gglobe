#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"

#include<iostream>
#include<vector>
#include<string>

int merger(char* filename = "all.root", char* outfilename = "opttree.root", char* treename = "opttree") {

  int nTrees = 0;
  TTree *tree[100];
  TList *list = new TList;
  TString firstName;

  TFile* f = new TFile(filename);
  
  // Create an iterator on the list of keys
  TIter nextTopLevelKey(f->GetListOfKeys());
  TKey *keyTopLevel;
  while (keyTopLevel = (TKey*)nextTopLevelKey()) {
    
    TString name(keyTopLevel->GetName());
    TString className(keyTopLevel->GetClassName());
    
    if (className.CompareTo("TTree") == 0) {
      
      if ((name.CompareTo("lumi") != 0) && 
	  (name.CompareTo("plotvariables") != 0) && 
	  (name.CompareTo("inputfiles") != 0)) {
	std::cout << "Adding " << keyTopLevel->GetName() << std::endl;
	if (nTrees == 0) {
	  firstName = name;
	  //std::cout << firstName << std::endl;
	}
	tree[nTrees] = (TTree*) f->Get(name);
	list->Add(tree[nTrees]);
	nTrees++;
      }
    }
  }

  TFile* out = new TFile(outfilename, "recreate");
  out->cd();
  TTree *opttree = TTree::MergeTrees(list);
  std::cout << list->GetName() << std::endl;
  opttree->SetName(treename);
  f->Close();
  
  opttree->Write();
  out->Close();
  

   out = new TFile(outfilename, "update");
  TIter nextKey(out->GetListOfKeys());
  TKey* key;
  
  while (key = (TKey*)nextKey()) {
    TString name(key->GetName());
    TString className(key->GetClassName());
    if (className.CompareTo("TTree") == 0) {
      if (name.CompareTo(treename) != 0) {
	std::cout << "Cleaning " << key->GetName() << std::endl;
	key->Delete();
      } else {
	key->SetTitle("Opttree");
      }
    }
  }

  out->Close();
  
  return 0;
}
