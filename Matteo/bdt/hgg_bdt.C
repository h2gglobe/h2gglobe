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
#endif

void hgg_bdt(char* outputFileName = "TMVA_photonid") {
  
  TMVA::Tools::Instance();
  char a[100];
  sprintf(a, "%s.root", outputFileName);
  TFile* outputFile = TFile::Open(a, "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory(outputFileName, outputFile, "!V:!Silent");

  factory->AddVariable("subleadptomass", 'F');
  factory->AddVariable("diphoptom", 'F');
  factory->AddVariable("sumptom", 'F');
  factory->AddVariable("subleadmva2", 'F');
  factory->AddVariable("leadmva2", 'F');
  factory->AddVariable("leadeta", 'F');
  factory->AddVariable("subleadeta", 'F');
  factory->AddVariable("leadr9", 'F');
  factory->AddVariable("subleadr9", 'F');
  factory->AddVariable("dmom", 'F');
  factory->AddVariable("pvtx", 'F');
  factory->AddVariable("nvtx", 'F');
  factory->AddVariable("sigma_mz", 'F');

  //factory->AddSpectator("diphocat2r92eta", 'I');

  TFile* inputSignal = TFile::Open("signal.root");
  TFile* inputBkg = TFile::Open("bg.root");
  
  TTree *signal     = (TTree*)inputSignal->Get("events");
  TTree *background = (TTree*)inputBkg->Get("events");
  
  factory->AddSignalTree    (signal);
  factory->AddBackgroundTree(background);

  factory->SetBackgroundWeightExpression("w");
  factory->SetSignalWeightExpression("w");
  
  //factory->PrepareTrainingAndTestTree("subleadmva > 0.5 && leadmva > 0.5", "subleadmva > 0.5 && leadmva > 0.5",
  //factory->PrepareTrainingAndTestTree("itype != -14 && itype != -24", "itype != -14 && itype != -24",
  factory->PrepareTrainingAndTestTree("", "",
				      "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V");
                                      //"nTrain_Signal=0:nTrain_Background=0:nTest_Signal=1:nTest_Background=1:!V");

  factory->BookMethod(TMVA::Types::kBDT, "Gradient", "!H:!V:NTrees=200:NNodesMax=20:MaxDepth=20:BoostType=Grad:nCuts=20");



  //factory->BookMethod(TMVA::Types::kBDT, "Gradient2000", "!H:!V:NTrees=400:NNodesMax=20:MaxDepth=20:BoostType=Grad:nCuts=2000");
  //factory->BookMethod(TMVA::Types::kBDT, "AdaBoost", "!H:!V:NTrees=400");
  
  //TString vars = "subleadptomass:diphoptom:sumptom:subleadmva:leadmva:leadeta:subleadeta:leadr9:subleadr9:";
  //TString vars = "subleadptomass:diphoptom:sumptom:subleadmva:leadmva:";
  /*
  TMVA::MethodCategory* mcat2 = 0;
  TMVA::MethodBase* fiCat2 = factory->BookMethod(TMVA::Types::kCategory, "Category_Gradient", "!V");
  mcat2 = dynamic_cast<TMVA::MethodCategory*>(fiCat2);
  mcat2->AddMethod("diphocat2r92eta < 0.5",
		   vars,
		  TMVA::Types::kBDT,
		  "Category_BDTCat_1",
		  "!H:!V:NTrees=100:MaxDepth=15:BoostType=Grad:nCuts=20"); 

  mcat2->AddMethod("diphocat2r92eta > 0.5 && diphocat2r92eta < 1.5",
		   vars,
		  TMVA::Types::kBDT,
		  "Category_BDTCat_2",
		  "!H:!V:NTrees=100:MaxDepth=15:BoostType=Grad:nCuts=20"); 

  mcat2->AddMethod("diphocat2r92eta > 1.5 && diphocat2r92eta < 2.5",
		   vars,
		  TMVA::Types::kBDT,
		  "Category_BDTCat_3",
		  "!H:!V:NTrees=100:MaxDepth=15:BoostType=Grad:nCuts=20"); 

  mcat2->AddMethod("diphocat2r92eta > 2.5 && diphocat2r92eta < 3.5",
		  vars,
		  TMVA::Types::kBDT,
		  "Category_BDTCat_4",
		  "!H:!V:NTrees=100:MaxDepth=15:BoostType=Grad:nCuts=20"); 
  */

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  outputFile->Close();

  delete factory;
}
