///==== include ====

#include "TFile.h"
#include "TChain.h"
#include "TMinuit.h"
#include <sstream>
#include <iostream>

#include "TMVA/Factory.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace std;
 
// --------- MAIN -------------------

int main(int argc, char** argv)
{ 

  // options
  bool  useDiphotonPt = true;
  bool  usePhotonsPt  = true;
  
  // Chain
  TChain* tree = new TChain("flatTree");
  //tree->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/vbfTest-17-02-2012/vbfAnalysisTree_vbfanalysis.root");
  // tree->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/vbfTest-18-07-2012/vbfAnalysisTree_vbfanalysis.root");
  tree->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/vbfTest-25-07-2012/vbfAnalysisTree_vbfanalysis.root");
  
  
  // Declaration of leaf types
  int entry ;
  float pho1pt;
  float pho2pt;
  float diphopt;
  float diphoM;
  float diphoEta;
  float dijetEta;
  float jet1pt;
  float jet2pt;
  float jet1eta;
  float jet2eta;
  float mj1j2;
  float zepp;
  float dphi;
  bool  isSignal;
  int   mctype;
    
  tree->SetBranchAddress("entry",   &entry);
  tree->SetBranchAddress("pho1pt",  &pho1pt);
  tree->SetBranchAddress("pho2pt",  &pho2pt);
  tree->SetBranchAddress("diphopt", &diphopt);
  tree->SetBranchAddress("diphoM",  &diphoM);
  tree->SetBranchAddress("diphoEta",&diphoEta);
  tree->SetBranchAddress("dijetEta",&dijetEta);  
  tree->SetBranchAddress("jet1pt",  &jet1pt);
  tree->SetBranchAddress("jet2pt",  &jet2pt);
  tree->SetBranchAddress("jet1eta", &jet1eta);
  tree->SetBranchAddress("jet2eta", &jet2eta);
  tree->SetBranchAddress("mj1j2",   &mj1j2);
  tree->SetBranchAddress("zepp",    &zepp);
  tree->SetBranchAddress("dphi",    &dphi);
  tree->SetBranchAddress("isSignal",&isSignal);
  tree->SetBranchAddress("mctype",  &mctype);


  // Create a new root output file.
  string outputFileName = argv[1];
  if (useDiphotonPt) outputFileName = outputFileName +"_diphopt";
  if (usePhotonsPt)  outputFileName = outputFileName +"_phopt";
  if (useDiphotonPt && usePhotonsPt) outputFileName = outputFileName;
  
  TFile* outputFile = TFile::Open((outputFileName+".root").c_str(), "RECREATE" );
  TMVA::Factory* factory = new TMVA::Factory(outputFileName.c_str(), outputFile, 
					     "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


  // -- variables

  factory->AddVariable( "jet1pt");
  factory->AddVariable( "jet2pt" );
  factory->AddVariable( "deta:=abs(jet1eta-jet2eta)" );
  factory->AddVariable( "mj1j2");
  factory->AddVariable( "zepp" );
  factory->AddVariable( "dphi" );
  
  if (useDiphotonPt)
    factory->AddVariable( "diphoptOverM:=diphopt/diphoM" );
  
  if (usePhotonsPt){
    factory->AddVariable( "pho1ptOverM:=pho1pt/diphoM" );
    factory->AddVariable( "pho2ptOverM:=pho2pt/diphoM" );
  }
  
  
  // -- spectators
  factory->AddSpectator("diphoM");
    

  //event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  

  // ====== To give different trees for training and testing, do as follows:
  //TFile *tmp = new TFile("tmvatrainingtemp.root","RECREATE");
  //tmp->ls();
   
  TTree *trainingSignal     = tree->CopyTree("isSignal==1 && mctype == -58");// mH = 124
  TTree *testSignal         = tree->CopyTree("isSignal==1 && mctype == -38");// mH = 125
  
  TTree *trainingBackground = tree->CopyTree("isSignal==0 && entry%3 != 0");
  TTree *testBackground     = tree->CopyTree("isSignal==0 && entry%3 == 0");

  factory->AddSignalTree( trainingSignal, signalWeight, "Training" );
  factory->AddSignalTree( testSignal,     signalWeight,  "Test" );
  factory->AddBackgroundTree( trainingBackground, backgroundWeight, "Training" );
  factory->AddBackgroundTree( testBackground,     backgroundWeight,  "Test" );
    
  // ====== register trees ====================================================
  //
  // the following method is the prefered one:
  // you can add an arbitrary number of signal or background trees
  //factory->AddSignalTree    ( tree,     signalWeight     );
  //factory->AddBackgroundTree( tree,     backgroundWeight );
   

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "isSignal == 1 && pho1pt/diphoM > (60./120.) && pho2pt/diphoM> (30./120.)"; 
   TCut mycutb = "isSignal == 0 && pho1pt/diphoM > (60./120.) && pho2pt/diphoM> (30./120.)"; 

   //... is photons PT is given as input to the MVA, use looser cuts 
   if ( usePhotonsPt ) {
     mycuts = "isSignal == 1 && pho1pt/diphoM > (40./120.) && pho2pt/diphoM> (25./120.)"; 
     mycutb = "isSignal == 0 && pho1pt/diphoM > (40./120.) && pho2pt/diphoM> (25./120.)"; 
   }
   
   // tell the factory to use all remaining events in the trees after training for testing:
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "SplitMode=Random:NormMode=NumEvents:!V" );
   
   
   // Boosted Decision Trees: use BDTG ( Gradient Boost )
   factory->BookMethod( TMVA::Types::kBDT, "BDTG",
			"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5:MaxDepth=3" );
   //"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=15:MaxDepth=5" );
   
   
   // book Cuts
   //factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
   //			"H:!V:FitMethod=GA:CutRangeMin[0]=20:CutRangeMax[0]=500:CutRangeMin[1]=20:CutRangeMax[1]=500:VarProp=FSmart:VarProp[4]=FMin:EffSel:Steps=30:Cycles=3:PopSize=500:SC_steps=10:SC_rate=5:SC_factor=0.95" );
		   
   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
}
