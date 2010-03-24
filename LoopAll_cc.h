#define DEBUG  0

#include "TH2.h"
#include "TGraph.h"
#include "TStyle.h"
#include <iostream>
#include <math.h>
#include "Util.h"
#include "TRandom.h"

#include <TLorentzVector.h>
#include <TVector3.h>

#include "stdlib.h"

using namespace std;

#include "Tools.h"

#include "PhotonAnalysis/PhotonAnalysisReducedOutputTree.h"

void LoopAll::InitReal(Int_t typerunpass) {

  // Set branch addresses
  typerun=typerunpass;
  hfile = (TFile*)gROOT->FindObject(UtilInstance->histFileName); 
  if (hfile) 
    hfile->Close();

  InitRealPhotonAnalysis(typerun);

  if (UtilInstance->makeOutputTree) 
    UtilInstance->outputFile->cd();
}

void LoopAll::TermReal(Int_t typerunpass) {

  typerun=typerunpass;

  TermRealPhotonAnalysis(typerun);
  
  if (UtilInstance->makeOutputTree) 
    UtilInstance->outputFile->cd();
}


void LoopAll::Init(Int_t typerunpass, TTree *tree) {
  
  // Set branch addresses
  typerun=typerunpass;
  
  fChain = tree;
  if (tree == 0) 
    return;
 
  fCurrent = -1;
  fChain->SetMakeClass(1);

  mySetBranchAddressRedPhotonAnalysis();

  Notify();
}


void LoopAll::setUtilInstance(Util *ut) {
  UtilInstance = ut;
}

void LoopAll::Loop(Util* ut) {
  
  makeOutputTree = ut->makeOutputTree;
  if (makeOutputTree) {
    outputTree = ut->outputTree;
    outputFile = ut->outputFile;
    if(outputFile) 
      outputFile->cd();
  }

  Int_t nentries = 0;
  if(fChain) 
    nentries = Int_t(fChain->GetEntriesFast());

  Int_t nbytes = 0, nb = 0;

  outputEvents=0;

  int hasoutputfile=0;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    
    if(jentry%10000==0) 
      cout << "Entry: "<<jentry<<endl;
    
    if(DEBUG) 
      cout<<"call LoadTree"<<endl;
    
    Int_t ientry = LoadTree(jentry);
  
    if (ientry < 0) 
      break;
   
    if(DEBUG) 
      cout<<"Loopall_cc RRRRRRRRRRRR calling GetEntry"<<endl;
      
    if(ut->typerun == 1) {
      nb = fChain->GetEntry(jentry);
      
      if(DEBUG) 
	cout<<"Loopall_cc RRRRRRRRRRR called GetEntry"<<endl;
    } else {
      //HERE NEED TO DO A GLOBAL GETENTRY
      nb=0;
    }

    nbytes += nb;

    if (Cut(ientry) < 0) 
	continue;

    if(DEBUG) 
      cout<<"Call FillandReduce "<<endl;
      
    hasoutputfile = FillAndReduce(ut, jentry);
      
    if(DEBUG) 
      cout<<"Called FillandReduce "<<endl;
  }

  if(hasoutputfile) {
    if(outputFile) {
      outputFile->cd();
      if (DEBUG)
	cout<<"LoopAll_cc writing outputTree"<<endl;
      outputTree->Write(0,TObject::kWriteDelete);
    }

    
    if(outputFile) {
      outputFile->cd();
      if(ut->TreesPar[0]) {
	std::vector<std::string> *parameters = new std::vector<std::string>;
	Int_t tot_events, sel_events, type, version, reductions;
	Int_t red_events[20];

	ut->TreesPar[0]->SetBranchAddress("tot_events", &tot_events);
	ut->TreesPar[0]->SetBranchAddress("sel_events", &sel_events);
	ut->TreesPar[0]->SetBranchAddress("type", &type);
	ut->TreesPar[0]->SetBranchAddress("version", &version);
	ut->TreesPar[0]->SetBranchAddress("parameters", &parameters);
	if (ut->TreesPar[0]->FindBranch("reductions")) {
	  ut->TreesPar[0]->SetBranchAddress("reductions", &reductions);
	  ut->TreesPar[0]->SetBranchAddress("red_events", &red_events);
	}
	
	TTree* newtree = new TTree("global_variables", "Global Parameters");
	newtree->Branch("tot_events", &tot_events, "tot_events/I");
	newtree->Branch("sel_events", &sel_events, "sel_events/I");
	newtree->Branch("type", &type, "type/I");
	newtree->Branch("version", &version, "version/I");
	newtree->Branch("parameters", "std::vector<string>", &parameters);
	newtree->Branch("reductions", &reductions, "reductions/I");
	newtree->Branch("red_events", &red_events, "red_events[reductions]/I");
	
	ut->TreesPar[0]->GetEntry(0);

	if (!ut->TreesPar[0]->FindBranch("reductions")) {
	  reductions = 1;
	  for (int i=0; i<20; i++) {
	    red_events[i] = -1;
	  }
	  red_events[0] = (int)countersred[1];
	} else {
	  red_events[reductions] = (int)countersred[1];
	  reductions++;
	}

	newtree->Fill();
	newtree->Write();
      } else {
	std::cerr << "Cannot write Parameter tree." << std::endl;
      }
    }
  }
  
  int oldnentries=nentries;
  if (nentries == ut->sel_events) {
    nentries = ut->tot_events;
  }

  if(countersred[1] || oldnentries==0) {
    printf("red: %d_%d \n",(int)countersred[0], (int) countersred[1]);
  } else { 
    printf("norm: %d \n",(int)counters[0]);
  }
}

void LoopAll::myWritePlot(Util * ut) {

  hfile = new TFile(UtilInstance->histFileName, "RECREATE", "Globe ROOT file with histograms");

  hfile->cd();
  hfile->cd();
  histoContainer->Save();
  if (UtilInstance->makeOutputTree) 
    UtilInstance->outputFile->cd();
}

int LoopAll::FillAndReduce(Util * ut, int jentry) {

  int hasoutputfile = 0;

  if(ut->typerun == 1) {
    hasoutputfile = 1;
    if(DEBUG) 
      cout<<"call myReduce"<<endl;
    myReducePhotonAnalysis(ut, jentry);
    if(DEBUG) 
      cout<<"called myReduce"<<endl;
  } else if (ut->typerun == 0) {
    hasoutputfile = 0;
    if(DEBUG) 
      cout<<"call myFillHist"<<endl;
    myFillHistPhotonAnalysis(ut, jentry);
    if(DEBUG) 
      cout<<"called myFillHist"<<endl;
  } else if (ut->typerun == 2) {
    hasoutputfile = 0;
    if(DEBUG) 
      cout<<"call myFillHistRed"<<endl;
    myFillHistPhotonAnalysisRed(ut, jentry);
    if(DEBUG) 
      cout<<"called myFillHistRed"<<endl;
  }
  
  return hasoutputfile;
}















