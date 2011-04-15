#ifndef LoopAll_h
#define LoopAll_h

#include "CommonParameters.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TCanvas.h"

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include <bitset>
#include <map>

class Util;

#include "HistoContainer.h"
#include "CounterContainer.h"
#include "SampleContainer.h"
#include "Cut.h"
#include "branchdef/Limits.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/PhotonInfo.h"

class LoopAll {
 public :
  TTree          *fChain;   
  
#include "branchdef/branchdef.h"
#include "branchdef/treedef.h"
#include "GeneralFunctions_h.h"
#include "PhotonAnalysis/PhotonAnalysisFunctions_h.h"

  std::vector<HistoContainer> histoContainer;
  std::vector<CounterContainer> counterContainer;
  std::vector<SampleContainer> sampleContainer;
  std::vector<Cut> cutContainer;
  RooContainer *rooContainer;

  LoopAll(TTree *tree=0);
  virtual ~LoopAll();
  //virtual Int_t  Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void   Init(Int_t typerunpass, TTree *tree);
  virtual void   InitReal(Int_t typerunpass);
  virtual void   TermReal(Int_t typerunpass);
  virtual void   Loop(Double_t a);
  virtual Bool_t Notify();
  virtual void   Show(Long64_t entry = -1);
  virtual void   InitHistos();
  virtual void   BookHisto(int,int,int,int,int,int
			  ,float,float,float,float
			  ,char*);
 
  Util *utilInstance;
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

  VertexAlgoParameters vtxAlgoParams;	 
  std::vector<std::string> vtxVarNames;

  void Loop(Int_t);
  void setUtilInstance(Util*);
  void myWritePlot();
  void myWriteFits();
  void myWriteCounters();
  
  int FillAndReduce(int);
  
  //void BookHistos();
  void AddCut(char*,int,int,int,float*,float*);
  void InitCounters();
  void AddCounter(int,char*,char*,char*,char*);
  
  int ApplyCut(int, float, int);
  int ApplyCut(std::string, float, int);
  
  void FillHist(std::string, float);
  void FillHist2D(std::string, float, float);

  void FillHist(std::string, int, float); 
  void FillHist2D(std::string, int, float, float);
  
  void FillCounter(std::string, int);
  void FillCounter(std::string);
};

#endif
