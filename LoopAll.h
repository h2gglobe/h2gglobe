#ifndef LoopAll_h
#define LoopAll_h

#include "CommonParameters.h"

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TClonesArray.h"

#include <sstream>
#include <set>
#include <string>
#include <list>

#include <THStack.h>
#include <TLegend.h>
#include <TAxis.h>

class BaseAnalysis;

#include "HistoContainer.h"
#include "CounterContainer.h"
#include "SampleContainer.h"
#include "Cut.h"
#include "branchdef/Limits.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/HggVertexFromConversions.h"
#include "VertexAnalysis/interface/VertexAlgoParameters.h"
#include "VertexAnalysis/interface/PhotonInfo.h"
#include "VertexAnalysis/interface/VertexAlgoParameters.h"


#define BRANCH_DICT(NAME) branchDict[# NAME] = branch_info_t(& b ## _ ## NAME, & LoopAll::SetBranchAddress ## _ ## NAME, & LoopAll::Branch ## _ ## NAME )

#define DEBUG 1

class LoopAll {
 public :
  TTree          *fChain;

#include "branchdef/branchdef.h"
#include "branchdef/treedef.h"
#include "branchdef/setbranchaddress.h"
#include "branchdef/treebranch.h"
#include "GeneralFunctions_h.h"

  std::vector<HistoContainer> histoContainer;
  std::vector<CounterContainer> counterContainer;
  std::vector<SampleContainer> sampleContainer;
  std::vector<Cut> cutContainer;
  RooContainer *rooContainer;
  
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

  void LoopAndFillHistos(TString treename="event");
  void MergeContainers();
  //void WriteHist();  
  //void WriteFits();  
  //void WriteCounters();  

  /** @param type is the type of the analysis:
             0 for looper.py (filling histograms)
             1 for reduce.py (reduction step)
      @param n is the name of the output file 
  */
  void SetTypeRun(int type, const char* n);

  //void SetOutputNames(const char* n, const char* n2="");
  void StoreProcessedLumis(TTree * tree);
  void AddFile(std::string,int);
  void ReadInput(int t=0);

  /** adds a new entry to sampleContainer and returns a reference
      to this entry. */
  SampleContainer & DefineSamples(const char *filesshortnam,
                            int type, int histtoindfromfiles, int histoplotit,
                            int nred, long long ntot, float intlumi,
                            float lumi, float xsec, float kfactor,
			    float scale, bool addnevents=false);
  void Term(); 

  std::vector<std::string> files;
  std::vector<int> itype;
  //int lumireal[MAXFILES];
  int nfiles;
  float intlumi;
  float intlumi_;

  std::vector<TTree*> Trees;
  std::vector<TTree*> LumiTrees;
  std::vector<TFile*> Files;
  std::vector<TTree*> TreesPar;

  TTree * outputTree;
  TTree * outputTreeLumi;
  TTree * outputTreePar;
  TTree * plotvartree;
  TTree * inputfiletree;

  TH1D  * pileup;

  TFile * outputFile;
  TString outputFileName;
  TString histFileName;
  Int_t makeOutputTree;

  char inputFilesName[1024];
  char countersName[1024];

  enum runtypes { kFill=0, kReduce=1, kFillReduce=2 };
  Int_t typerun;
  Int_t typeread;

  std::map<int,int> type2HistVal;

  Int_t        current; //current file
  Int_t        current_sample_index; //current file
  Int_t        tot_events;
  Int_t        sel_events;
  std::vector<std::string>  globalCountersNames; 
  std::vector<int> globalCounters; 
  std::vector<int> fileGlobalCounters;
  
  Int_t        outputParTot_Events, outputParSel_Events, outputParType, outputParVersion, outputParReductions, outputParRed_Events[20];
  std::vector<std::string>* outputParParameters;
  std::string* outputParJobMaker;

  Int_t currentindexfiles;
  
  std::vector<Float_t> counters;
  std::vector<Float_t> countersred;
   
  TFile *hfile;
  Int_t outputEvents;

  void Loop(Int_t);
  void WriteHist();
  void WriteFits();
  void WriteCounters();
  
  int FillAndReduce(int);
  
  //void BookHistos();
  void AddCut(char*,int,int,int,float*,float*);
  void InitCounters();

  void AddCounter(int, const char*, const char*, const char*, const char*);

  int ApplyCut(int, float, int);
  int ApplyCut(std::string, float, int);

  void myPrintCounters();
  void myPrintCountersNew();

  virtual void   BookHisto(int, int, int, int, int, int,
                           float, float, float, float, const char*,
			   const char* xaxis="", const char* yaxis="");
 

  void WritePI();
  void AddCut2(char*, int, int, int, float*, float*, int, int, int, float, float, char*, char*);

  int ApplyCut(int icut, int icat);
  int ApplyCut(TString cutname, int icat);

  float GetCutValue(TString cutname, int icat, int highcut=0);

  //int ApplyCut(int icut, float var, int icat);
  //int ApplyCut(TString cutname, float var, int icat);

  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCut(int icut, int * passcategory); //returns the number of categories
  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCut(TString cutname, int * passcategory); //returns the number of categories
  
  int ApplyCutsFill(int icat, int cutset, int & ncutsapplied, int & ncutspassed,  int & ncutsfailed, float histweight=1.0, float countweight=1.0);
  int ApplyCuts(int icat, int cutset, int & ncutsapplied, int & ncutspassed,  int & ncutsfailed);
  int ApplyCutsFill(int icat, int cutset, float histweight=1.0, float countweight=1.0);
  int ApplyCuts(int icat, int cutset);
  
  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCuts(int cutset, int * passcategory, int * ncutsapplied, int * ncutspassed,  int * ncutsfailed); //returns the number of categories
  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCuts(int cutset, int * passcategory); //returns the number of categories
  
  int SetCutVariables(int i, float * variables);
  int SetCutVariables(TString cutname, float * variables);  

  void FillHist(std::string, float);
  void FillHist2D(std::string, float, float);

  void FillHist(std::string, int, float, float wt = 1.0); 
  void FillHist2D(std::string, int, float, float, float wt = 1.0);

  void FillCounter(std::string name, float weight=1., int cat=0);

  BaseAnalysis* AddAnalysis(BaseAnalysis*); 
  template<class T> T * GetAnalysis( const std::string & name ) {
	  std::vector<BaseAnalysis*>::iterator it=find( analyses.begin(), analyses.end(), name);
	  if( it != analyses.end() ) { return dynamic_cast<T*>( *it ); }
	  return 0;
  }
  
  /// void SkimBranch(const std::string & name)   { skimBranchNames.insert(name);  };
  void InputBranch(const std::string & name, int typ)  { inputBranchNames.insert(std::pair<std::string,int>(name,typ)); };
  void OutputBranch(const std::string & name) { if( find(outputBranchNames.begin(), outputBranchNames.end(), name)==outputBranchNames.end() ) { outputBranchNames.push_back(name); } };

  void GetBranches(std::map<std::string,int> & names, std::set<TBranch *> & branches);
  void SetBranchAddresses(std::map<std::string,int> & names);
  void Branches(std::list<std::string> & names);
  void GetEntry(std::set<TBranch *> & branches, int jentry);

  bool CheckLumiSelection( int run, int lumi );
  bool CheckEventList( int run, int lumi, int event );

#ifndef __CINT__
  typedef void (LoopAll::*branch_io_t) (TTree *);
  struct branch_info_t {
	  branch_info_t(TBranch ** b=0, branch_io_t r=0, branch_io_t w=0 ) :
		  branch(b), read(r), write(w) 
		  {};
	  TBranch ** branch;
	  branch_io_t read, write;
  };
  typedef std::map<std::string, branch_info_t> dict_t;
  dict_t branchDict;
#endif

  //// std::set<std::string> skimBranchNames;
  //// std::set<TBranch *> skimBranches; 
  std::map<std::string,int> inputBranchNames;
  std::set<TBranch *> inputBranches; 
  std::list<std::string> outputBranchNames;
  
  /** list of the analyses to be performed */
  std::vector<BaseAnalysis*> analyses;
};

#endif
