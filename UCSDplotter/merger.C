#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"

#include<iostream>
#include<vector>
#include<string>


std::vector<std::string> checkNames(TFile* f, const char* s) {

  std::vector<std::string> results;
  char a[100];
  for (int i=0; i<100; i++) {
    sprintf(a, "%s;%d", s, i);  
    std::cout << a << std::endl;
    if (f->GetObjectUnchecked(a) != 0) {
      std::cout << a << std::endl;
      results.push_back(a);
    }
  }

  return results;
}

int merger(const char* filename = "all.root") {

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

  TFile* out = new TFile("opttree.root", "recreate");
  TTree *opttree = TTree::MergeTrees(list);
  opttree->SetName("opttree");
  f->Close();
  
  opttree->Write();
  out->Close();
  
  //out = new TFile("opttree.root", "update");
  //std::vector<std::string> temp2 = checkNames(f, firstName);
  //for (unsigned int y=0; y<temp2.size()-1; y++) 
  //  f->Delete(temp2[y].c_str());
  //out->Write();
  //out->Close();
  
  return 0;
}
