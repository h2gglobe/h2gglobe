#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/MethodBase.h"
#endif

void photonid_bdt(char* outputFileName = "TMVA_photonid") {
  
  TMVA::Tools::Instance();
  char a[100];
  sprintf(a, "%s.root", outputFileName);
  TFile* outputFile = TFile::Open(a, "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory(outputFileName, outputFile, "!V:!Silent");

  factory->AddVariable("sieie", 'F');
  //factory->AddVariable("good_iso", 'F');
  //factory->AddVariable("bad_iso", 'F');
  factory->AddVariable("goodpf_iso", 'F');
  factory->AddVariable("badpf_iso", 'F');
  factory->AddVariable("drtotk", 'F');
  factory->AddVariable("hoe", 'F');
  //factory->AddVariable("tkiso", 'F');
  factory->AddVariable("tkisopf", 'F');
  factory->AddVariable("r9", 'F');
  factory->AddVariable("ptom", 'F');
  factory->AddVariable("eta", 'F');
  factory->AddVariable("nvtx", 'F');
  factory->AddVariable("scwidtheta", 'F');
  factory->AddVariable("scwidthphi", 'F');
  //factory->AddSpectator("isLeading", 'F');

  TFile* inputSignal = TFile::Open("signal_ID.root");
  TFile* inputBkg = TFile::Open("bg_ID.root");
  
  TTree *signal     = (TTree*)inputSignal->Get("events");
  TTree *background = (TTree*)inputBkg->Get("events");
  
  factory->AddSignalTree    (signal);
  factory->AddBackgroundTree(background);

  factory->SetBackgroundWeightExpression("w");
  factory->SetSignalWeightExpression("w");
  
  //factory->PrepareTrainingAndTestTree("itype > 109000 && itype < 136000", "",
  factory->PrepareTrainingAndTestTree("", "",
				      "nTrain_Signal=0:nTrain_Background=0::nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V" );

  //factory->BookMethod(TMVA::Types::kBDT, "Adaptive", "!H:!V:NTrees=200:BoostType=AdaBoost");
  factory->BookMethod(TMVA::Types::kBDT, "Gradient", "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=10:BoostType=Grad:nCuts=20");
  //factory->BookMethod(TMVA::Types::kBDT, "Gradient2000", "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=9:BoostType=AdaBoost:nCuts=2000");

  //TString vars = "sieie:goodpf_iso:badpf_iso:drtotk:hoe:tkisopf:pt:r9:eta:";
  /*
  TMVA::MethodCategory* mcat2 = 0;
  TMVA::MethodBase* fiCat2 = factory->BookMethod(TMVA::Types::kCategory, "Category_Gradient", "!V");
  mcat2 = dynamic_cast<TMVA::MethodCategory*>(fiCat2);
  mcat2->AddMethod("isLeading < 0",
		   vars,
		   TMVA::Types::kBDT,
		   "Category_BDTCat_1",
		   "!H:!V:NTrees=400:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
  mcat2->AddMethod("isLeading > 0",
		   vars,
		   TMVA::Types::kBDT,
		   "Category_BDTCat_2",
		   "!H:!V:NTrees=400:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
  */
  /*
    mcat2->AddMethod("cat == 0",
    vars,
    TMVA::Types::kBDT,
    "Category_BDTCat_1",
    "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
    mcat2->AddMethod("cat == 1",
    vars,
		  TMVA::Types::kBDT,
		  "Category_BDTCat_2",
		  "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
  mcat2->AddMethod("cat == 2",
		   vars,
		  TMVA::Types::kBDT,
		  "Category_BDTCat_3",
		  "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
  mcat2->AddMethod("cat == 3",
		  vars,
		  TMVA::Types::kBDT,
		  "Category_BDTCat_4",
		  "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
  mcat2->AddMethod("cat == 4",
		   vars,
		   TMVA::Types::kBDT,
		  "Category_BDTCat_5",
		  "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
  mcat2->AddMethod("cat == 5",
		   vars,
		   TMVA::Types::kBDT,
		  "Category_BDTCat_6",
		  "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=9:BoostType=Grad:nCuts=9"); 
  */
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  outputFile->Close();

  delete factory;
}
