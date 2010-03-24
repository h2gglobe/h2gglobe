#include "Util.h"
#include "LoopAll.h"

#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>

#define DEBUG 0

Util::Util() {
  loops = new LoopAll();
  loops->setUtilInstance(this);
#include "branchdef/newclonesarray.h"
}

void Util::SetTypeRun(int t, const char* name) {

  typerun = t;

  if (t == 1) {
    makeOutputTree = 1;
    outputFileName = TString(name);
    histFileName  = TString("");
  } else {
    makeOutputTree = 0;
    outputFileName = TString("");
    histFileName  = TString(name);
  }
  
  if (makeOutputTree) {
    cout << "CREATE " << outputFileName<<endl;
    outputFile=new TFile(outputFileName,"recreate");
    outputFile->cd();
    if(outputFile) 
      if (DEBUG)
	cout<<"output file defined"<<endl;
  }
  
  int hasoutputfile=0;
  
  if(typerun == 0 || typerun == 2) {
    hasoutputfile=0;
  } else if(typerun == 1) {

    hasoutputfile=1;
    outputTree = new TTree("event","Reduced tree"); 
    outputTreePar = new TTree("global_variables","Parameters");
    
    loops->PhotonAnalysisReducedOutputTree();
  } 
}

void Util::SetOutputNames(const char* name, const char* name2) {

  if (name2 != "")
    makeOutputTree = 1;
  histFileName  = TString(name);
  outputFileName = TString(name2);
}

void Util::AddFile(char* name) {
  files[nfiles] = name;
  nfiles++;
}

void Util::LoopAndFillHistos(TString treename) {

  int i=0;
  
  if (DEBUG)
    cout<<"LoopAndFillHistos: calling InitReal "<<endl;

  loops->InitReal(typerun);

  while (i<nfiles) {
    
    cout<<"LoopAndFillHistos: opening " << i << " " << files[i]<<endl;

    Files[i] = TFile::Open(files[i]);
    tot_events=1;
    sel_events=1;
    if(typerun == 1) { //this is a reduce job

      if(Files[i])
	TreesPar[i]=(TTree*) Files[i]->Get("global_variables");

      if(TreesPar[i]) {
	TBranch        *b_tot_events;
	TBranch        *b_sel_events;
	TreesPar[i]->SetBranchAddress("tot_events",&tot_events, &b_tot_events);
	TreesPar[i]->SetBranchAddress("sel_events",&sel_events, &b_sel_events);
	b_tot_events->GetEntry(0);
	b_sel_events->GetEntry(0);
      } else {
	cout<<"REDUCE JOB, no global_variables tree ... the C-step job must have crashed ... SKIP THE FILE"<<endl;
	tot_events=0;
	sel_events=0;
      }
    }

    if(tot_events!=0) {

      if(Files[i])
	Trees[i]=(TTree*) Files[i]->Get(treename);

      loops->Init(typerun, Trees[i]);
    }

    loops->Loop(this);

    if(tot_events != 0) {
      Trees[i]->Delete("");
    }
    
    if(Files[i])
      Files[i]->Close();
    
    i++;
  }
 
  loops->TermReal(typerun);
}

void Util::WriteHist() {
    loops->myWritePlot(this);
}

