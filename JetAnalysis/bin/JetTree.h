//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct  8 09:36:12 2012 by ROOT version 5.32/00
// from TTree flatTree/flatTree
// found on file: jetid_jetanalysis.root
//////////////////////////////////////////////////////////

#ifndef JetTree_h
#define JetTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         dRMean;
   Float_t         frac01;
   Float_t         frac02;
   Float_t         frac03;
   Float_t         frac04;
   Float_t         frac05;
   Float_t         frac06;
   Float_t         frac07;
   Float_t         nNeutrals;
   Float_t         beta;
   Float_t         betaStar;
   Float_t         dZ;
   Float_t         nCharged;
   Float_t         dR2Mean;
   Float_t         betaStarClassic;
   Float_t         jetPt;
   Float_t         jetEta;
   Float_t         nvtx;
   Int_t           ievent;
   Int_t           ijet;
   Bool_t          isMatched;
   Float_t         jetGenPt;
   Float_t         jetHenDr;
   Float_t         puweight;
   Int_t           npu;
   Int_t           isData;
   Float_t         njets;
   Bool_t          jetLooseID;
   Int_t           simpleId;
   Int_t           fullId;
   Int_t           cutbasedId;
   Float_t         simpleDiscriminant;
   Float_t         fullDiscriminant;
   Float_t         cutbasedDiscriminant;
   Float_t         dimuonMass;
   Float_t         dimuonPt;
   Float_t         dphiZJet;

   // List of branches
   TBranch        *b_dRMean;   //!
   TBranch        *b_frac01;   //!
   TBranch        *b_frac02;   //!
   TBranch        *b_frac03;   //!
   TBranch        *b_frac04;   //!
   TBranch        *b_frac05;   //!
   TBranch        *b_frac06;   //!
   TBranch        *b_frac07;   //!
   TBranch        *b_nNeutrals;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_betaStar;   //!
   TBranch        *b_dZ;   //!
   TBranch        *b_nCharged;   //!
   TBranch        *b_dR2Mean;   //!
   TBranch        *b_betaStarClassic;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_ievent;   //!
   TBranch        *b_ijet;   //!
   TBranch        *b_isMatched;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetHenDr;   //!
   TBranch        *b_puweight;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_jetLooseID;   //!
   TBranch        *b_simpleId;   //!
   TBranch        *b_fullId;   //!
   TBranch        *b_cutbasedId;   //!
   TBranch        *b_simpleDiscriminant;   //!
   TBranch        *b_fullDiscriminant;   //!
   TBranch        *b_cutbasedDiscriminant;   //!
   TBranch        *b_dimuonMass;   //!
   TBranch        *b_dimuonPt;   //!
   TBranch        *b_dphiZJet;   //!

   JetTree(TTree *tree=0);
   virtual ~JetTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif


JetTree::JetTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("jetid_jetanalysis.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("jetid_jetanalysis.root");
      }
      f->GetObject("flatTree",tree);

   }
   Init(tree);
}

JetTree::~JetTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("dRMean", &dRMean, &b_dRMean);
   fChain->SetBranchAddress("frac01", &frac01, &b_frac01);
   fChain->SetBranchAddress("frac02", &frac02, &b_frac02);
   fChain->SetBranchAddress("frac03", &frac03, &b_frac03);
   fChain->SetBranchAddress("frac04", &frac04, &b_frac04);
   fChain->SetBranchAddress("frac05", &frac05, &b_frac05);
   fChain->SetBranchAddress("frac06", &frac06, &b_frac06);
   fChain->SetBranchAddress("frac07", &frac07, &b_frac07);
   fChain->SetBranchAddress("nNeutrals", &nNeutrals, &b_nNeutrals);
   fChain->SetBranchAddress("beta", &beta, &b_beta);
   fChain->SetBranchAddress("betaStar", &betaStar, &b_betaStar);
   fChain->SetBranchAddress("dZ", &dZ, &b_dZ);
   fChain->SetBranchAddress("nCharged", &nCharged, &b_nCharged);
   fChain->SetBranchAddress("dR2Mean", &dR2Mean, &b_dR2Mean);
   fChain->SetBranchAddress("betaStarClassic", &betaStarClassic, &b_betaStarClassic);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("ievent", &ievent, &b_ievent);
   fChain->SetBranchAddress("ijet", &ijet, &b_ijet);
   fChain->SetBranchAddress("isMatched", &isMatched, &b_isMatched);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetHenDr", &jetHenDr, &b_jetHenDr);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("jetLooseID", &jetLooseID, &b_jetLooseID);
   fChain->SetBranchAddress("simpleId", &simpleId, &b_simpleId);
   fChain->SetBranchAddress("fullId", &fullId, &b_fullId);
   fChain->SetBranchAddress("cutbasedId", &cutbasedId, &b_cutbasedId);
   fChain->SetBranchAddress("simpleDiscriminant", &simpleDiscriminant, &b_simpleDiscriminant);
   fChain->SetBranchAddress("fullDiscriminant", &fullDiscriminant, &b_fullDiscriminant);
   fChain->SetBranchAddress("cutbasedDiscriminant", &cutbasedDiscriminant, &b_cutbasedDiscriminant);
   fChain->SetBranchAddress("dimuonMass", &dimuonMass, &b_dimuonMass);
   fChain->SetBranchAddress("dimuonPt", &dimuonPt, &b_dimuonPt);
   fChain->SetBranchAddress("dphiZJet", &dphiZJet, &b_dphiZJet);
   Notify();
}

Bool_t JetTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
