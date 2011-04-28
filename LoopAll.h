#ifndef LoopAll_h
#define LoopAll_h

#include "CommonParameters.h"

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"

#include <sstream>

#include "HistoContainer.h"
#include "CounterContainer.h"
#include "SampleContainer.h"
#include "Cut.h"
#include "branchdef/Limits.h"
//#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/HggVertexFromConversions.h"
#include "VertexAnalysis/interface/VertexAlgoParameters.h"
#include "VertexAnalysis/interface/PhotonInfo.h"
#include "VertexAnalysis/interface/VertexAlgoParameters.h"

#define DEBUG 0

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
  //RooContainer *rooContainer;
  
  LoopAll(TTree *tree=0);
  virtual ~LoopAll();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void   Init(Int_t typerunpass, TTree *tree);
  virtual void   InitReal(Int_t typerunpass);
  virtual void   TermReal(Int_t typerunpass);
  //virtual void   Loop(Double_t a);
  virtual Bool_t Notify();
  virtual void   Show(Long64_t entry = -1);
  virtual void   InitHistos();
  virtual void   BookHisto(int,int,int,int,int,int
			  ,float,float,float,float
			  ,char*);
 
  //Util *utilInstance;
  void LoopAndFillHistos(TString treename="event");
  void WriteHist();  
  void WriteFits();  
  void WriteCounters();  
  void SetTypeRun(int, const char* n);
  //void SetOutputNames(const char* n, const char* n2="");
  void AddFile(std::string,int);
  void ReadInput(int t=0);
  void DefineSamples(const char*,int,int,int,int, long long
		    ,float,float,float,float,float);

  void Term(); 

  std::vector<std::string> files;
  std::vector<int> itype;
  //int lumireal[MAXFILES];
  int nfiles;
  float intlumi;

  std::vector<TTree*> Trees;
  std::vector<TFile*> Files;
  std::vector<TTree*> TreesPar;

  TTree * outputTree;
  TTree * outputTreePar;

  TFile * outputFile;
  TString outputFileName;
  TString histFileName;
  Int_t makeOutputTree;

  char inputFilesName[1024];
  char countersName[1024];

  Int_t typerun;
  Int_t typeread;

  std::map<int,int> type2HistVal;

  Int_t        current; //current file
  Int_t        tot_events;
  Int_t        sel_events;

  Int_t        outputParTot_Events, outputParSel_Events, outputParType, outputParVersion, outputParReductions, outputParRed_Events[20];
  std::vector<std::string>* outputParParameters;
  std::string* outputParJobMaker;

  Int_t currentindexfiles;
  
  Float_t counters[2];
  Float_t countersred[2];
   
  TFile *hfile;
  Int_t outputEvents;

  VertexAlgoParameters vtxAlgoParams;	 
  std::vector<std::string> vtxVarNames;
  HggVertexAnalyzer vtxAna;
  HggVertexFromConversions vtxConv;

  void Loop(Int_t);
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
  
  std::set<TBranch *> branchesToRead; 
};

#endif
