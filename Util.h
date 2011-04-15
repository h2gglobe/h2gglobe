#ifndef Util_h
#define Util_h

#include "CommonParameters.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TArrow.h>
#include <iostream>
#include <map>

using namespace std;

class LoopAll;

class Util {
  
 public:

  Util();
 
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

  void CallAddCut(char*,int,int,int,float*,float*);
  void CallInitHistos();
  void CallBookHisto(int,int,int,int,int,int
		    ,float,float,float,float
		    ,char *);  
  void CallInitCounters();
  void CallAddCounter(int,char*,char*,char*,char*);

  std::vector<std::string> files;
  std::vector<int> itype;
  //int lumireal[MAXFILES];
  int nfiles;
  float intlumi;

  std::vector<TTree*> Trees;
  std::vector<TFile*> Files;
  std::vector<TTree*> TreesPar;

  LoopAll * loops;

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

  map<int,int> type2HistVal;

  Int_t        current; //current file
  Int_t        tot_events;
  Int_t        sel_events;

  Int_t        outputParTot_Events, outputParSel_Events, outputParType, outputParVersion, outputParReductions, outputParRed_Events[20];
  std::vector<std::string>* outputParParameters;
  std::string* outputParJobMaker;
};
#endif

