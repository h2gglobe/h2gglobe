///==== include ====

#include "TFile.h"
#include "TChain.h"
#include "TMinuit.h"
#include <sstream>
#include <iostream>
#include <cstdlib> 
#include "TMVA/Factory.h"
//#include "TMVAMultiClassGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace std;
 
// --------- MAIN -------------------

int main(int argc, char** argv)
{ 
  
  if (argc < 5) std::cout << "Missing arguments " << std::endl;

  // options
  bool   useDiphotonPt = atoi(argv[2]);
  bool   usePhotonsPt  = atoi(argv[3]);
  string str_mjjCut    = argv[4]; 

  cout << "Using diphopt : " << useDiphotonPt << std::endl;
  cout << "Using phopt   : " << usePhotonsPt << std::endl;
  cout << "preselection M(jj) > " << str_mjjCut << std::endl;

  // Chain
  TChain* tree = new TChain("flatTree");
  tree->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VBF/vbfTest-08-08-2012/vbfAnalysisTree_vbfanalysis.root");
  
  
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
					     "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );


  // -- variables

  factory->AddVariable( "jet1pt");
  factory->AddVariable( "jet2pt" );
  factory->AddVariable( "deta:=abs(jet1eta-jet2eta)" );
  //factory->AddVariable( "jet1eta" );
  //factory->AddVariable( "jet2eta" );
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
  // -- signal sample
  TTree *signal     = tree->CopyTree("isSignal==1 && mctype == -58");// VBF mH = 124
    
  // -- background sample 0 : diphoton jets
  TTree *background0 = tree->CopyTree("isSignal==0 "); // DiPhoton+jets
  
  // -- background sample 1 : GluGlu
  TTree *background1 = tree->CopyTree("isSignal==1 && mctype == -57"); // GluGlu mH = 124
  
    

  // ====== register trees ====================================================
  //
  
  // Apply additional cuts on the signal and background samples (can be different)
  string str_cut  = "pho1pt/diphoM > (60./120.) && pho2pt/diphoM> (30./120.) && mj1j2 > "+str_mjjCut;
  TCut mycut = str_cut.c_str();
  
  //... is photons PT is given as input to the MVA, use looser cuts 
  if ( usePhotonsPt ) {
    str_cut = "pho1pt/diphoM > (40./120.) && pho2pt/diphoM> (25./120.) && mj1j2 > "+str_mjjCut;
    mycut = str_cut.c_str();
  }
   
  factory->AddTree    (signal,"Signal",  signalWeight, mycut );
  factory->AddTree    (background0,"bg0", backgroundWeight, mycut);
  factory->AddTree    (background1,"bg1", backgroundWeight, mycut);
 
   // tell the factory to use all remaining events in the trees after training for testing:
  factory->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );
  
   
   // Boosted Decision Trees: use BDTG ( Gradient Boost )
   factory->BookMethod( TMVA::Types::kBDT, "BDTG",
			"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.50:nCuts=20:NNodesMax=8");
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
