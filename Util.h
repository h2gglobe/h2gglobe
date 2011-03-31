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
  void WriteCounters();  
  void SetTypeRun(int, const char* n);
  void SetOutputNames(const char* n, const char* n2="");
  void AddFile(char*,int);
  void ReadInput(int t=0);
  

  char* files[MAXFILES];
  int itype[MAXFILES];
  //int lumireal[MAXFILES];
  int nfiles;
  float intlumi;

  TFile *Files[MAXFILES];
  TTree *Trees[MAXFILES];
  TTree *TreesPar[MAXFILES];

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

};
#endif

