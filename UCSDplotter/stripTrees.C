#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"

#include<iostream>
#include<vector>
#include<string>

int stripTrees(char* filename = "all.root") {
  
  char outfilename[100];
  std::string filenameString(filename);
  std::cout << filenameString.size() << std::endl; 
  std::string filenameStripped = filenameString.substr(filenameString.size()-5, filenameString.size());
  //  filenameStripped.erase(filenameStripped
  sprintf(outfilename, "%s_histograms.root", filenameString.c_str());

  std::vector<TTree*> trees;
  std::vector<TH1F*> th1fs;
  std::vector<TH2F*> th2fs;
  std::vector<TProfile*> tprofiles;

  TFile* f = new TFile(filename);

  // Create an iterator on the list of keys
  TIter nextTopLevelKey = NULL;
  nextTopLevelKey=(f->GetListOfKeys());
  
  TKey *keyTopLevel;
  while (keyTopLevel = (TKey*)nextTopLevelKey()) {
    
    TString name(keyTopLevel->GetName());
    TString className(keyTopLevel->GetClassName());
    
    if ((className.CompareTo("TTree") != 0) || 
	((name.CompareTo("plotvariables") == 0) || (name.CompareTo("inputfiles") == 0))) {
      std::cout << "Adding " << keyTopLevel->GetName() << std::endl;
      
      if (className.CompareTo("TTree") == 0)
	trees.push_back((TTree*) f->Get(name));
      else if (className.CompareTo("TH1F") == 0)
	th1fs.push_back((TH1F*) f->Get(name));
      else if (className.CompareTo("TH2F") == 0)
	th2fs.push_back((TH2F*) f->Get(name));
      else if (className.CompareTo("TProfile") == 0)
	tprofiles.push_back((TProfile*) f->Get(name));
    }
  }

  TFile* out = new TFile(outfilename, "recreate");
  out->cd();
  for (unsigned int i=0; i<trees.size(); i++) 
    trees[i]->Write();

  for (unsigned int i=0; i<th1fs.size(); i++) 
    th1fs[i]->Write();

  for (unsigned int i=0; i<th2fs.size(); i++) 
    th2fs[i]->Write();

  for (unsigned int i=0; i<tprofiles.size(); i++) 
    tprofiles[i]->Write();

  out->Close();
  
  return 0;
}
