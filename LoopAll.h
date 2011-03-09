#ifndef LoopAll_h
#define LoopAll_h

#include "CommonParameters.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TCanvas.h"

#include <TH1.h> 
#include <TH2.h> 
#include <TProfile.h> 

#include <TClonesArray.h>
#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

// binary numbers (HLT)
#include <bitset>
#include <map>

#if defined (__MAKECINT__ )
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::vector<int> >+;
#endif 

class Util;

#include "HistoContainer.h"
#include "../interface/Limits.h"

class LoopAll {
 public :
  TTree          *fChain;   
  Int_t           fCurrent; 


#include "branchdef/branchdef.h"
#include "branchdef/treedef.h"

  std::vector<HistoContainer*> histoContainer;
  

  LoopAll(TTree *tree=0);
  virtual ~LoopAll();
  virtual Int_t  Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void   Init(Int_t typerunpass, TTree *tree);
  virtual void   InitReal(Util* ut, Int_t typerunpass);
  virtual void   TermReal(Int_t typerunpass);
  virtual void   Loop();
  virtual Bool_t Notify();
  virtual void   Show(Long64_t entry = -1);
  
  Util *UtilInstance;
  Int_t typerun;
  Int_t currentindexfiles;
  
  Float_t counters[2];
  Float_t countersred[2];
   
  TFile *hfile;
  TFile * outputFile;
  TTree * outputTree;
  TTree * outputTreePar;
  Int_t makeOutputTree;
  Int_t outputEvents;

#include "GeneralFunctions_h.h"
#include "PhotonAnalysis/PhotonAnalysisFunctions_h.h"
 
 void Loop(Util*);
 void setUtilInstance(Util*);
 void myWritePlot(Util*);
 
 int FillAndReduce(Util*, int);
};

#endif

#ifdef LoopAll_cxx
LoopAll::LoopAll(TTree *tree) {
  //histoContainer = new HistoContainer();
}

LoopAll::~LoopAll() {
  if (!fChain) 
    return;
  delete fChain->GetCurrentFile();
}

Int_t LoopAll::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t LoopAll::LoadTree(Long64_t entry) {

  // Set the environment to read one entry
  if (!fChain) return -5;
  Int_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->IsA() != TChain::Class()) 
    return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

Bool_t LoopAll::Notify() {
  return kTRUE;
}

void LoopAll::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t LoopAll::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}


#endif

