//***** Simple macro to compute systematics on di-jet tagged categories from PU-jetID efficiency.
// Inputs needed: 
// (1) file containing mini-tree with jet related variables for ggh, vbf, wzh, tth
// (2) data/MC scale factors derived from Z->mumu, e.g.: AnalysisScripts/aux/JetIdScaleFactor_ZmumuJets_12fb_hcp.root

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1D.h"
#include "TCanvas.h"

int etaregion(float jeteta){
  int region=-1;
  if (fabs(jeteta)<2.50 )                     region = 0; //TK
  if (fabs(jeteta)>2.50 && fabs(jeteta)<2.75) region = 1; //HEin
  if (fabs(jeteta)>2.75 && fabs(jeteta)<3.0)  region = 2; //HEout
  if (fabs(jeteta)>3.00 )                     region = 3; //HF
  return region;
}


void makeJetIdEffUncertainties(string infilename, int mass=125)
{

  TFile *inFile = TFile::Open(infilename.c_str());
  
  //--- load trees
  vector<TTree*> trees;
  trees.push_back((TTree*)inFile->Get(Form("ggh_m%d_8TeV",mass)));
  trees.push_back((TTree*)inFile->Get(Form("vbf_m%d_8TeV",mass)));
  trees.push_back((TTree*)inFile->Get(Form("wzh_m%d_8TeV",mass)));
  trees.push_back((TTree*)inFile->Get(Form("tth_m%d_8TeV",mass)));
  
  TList *listOfTrees = new TList();
  for (vector<TTree*>::iterator it=trees.begin(); it!=trees.end(); it++) listOfTrees->Add(*it);
  
  TFile output("out.root", "RECREATE");
  
  TTree *tree = TTree::MergeTrees(listOfTrees);
  tree->SetName(Form("AllSignal_m%d",mass)); 
  cout << "Merged trees.... " << endl;
  //tree->Print();

  int sampleType;
  float leadPtOverM, subleadPtOverM;
  float leadJPt, subleadJPt;
  float leadJEta, subleadJEta;
  float vbfMVA, diphotonMVA;
  int category ;
  float weight;
  float Zep;
  float MJJ;
  float deltaPhiJJGamGam;
  float deltaEtaJJ;

  tree->SetBranchAddress("sampleType",&sampleType);
  tree->SetBranchAddress("leadPtOverM",&leadPtOverM);
  tree->SetBranchAddress("subleadPtOverM",&subleadPtOverM);
  tree->SetBranchAddress("leadJPt",&leadJPt);
  tree->SetBranchAddress("subleadJPt",&subleadJPt);
  tree->SetBranchAddress("leadJEta",&leadJEta);
  tree->SetBranchAddress("subleadJEta",&subleadJEta);
  tree->SetBranchAddress("MJJ",&MJJ);
  tree->SetBranchAddress("deltaEtaJJ",&deltaEtaJJ);
  tree->SetBranchAddress("deltaPhiJJGamGam",&deltaPhiJJGamGam);
  tree->SetBranchAddress("Zep",&Zep);
  tree->SetBranchAddress("vbfMVA",&vbfMVA);
  tree->SetBranchAddress("diphotonMVA",&diphotonMVA);
  tree->SetBranchAddress("category",&category);
  tree->SetBranchAddress("weight",&weight);
 
  int ncats=9;

  //--- process names
  vector<string> procnames;
  procnames.push_back("ggH");
  procnames.push_back("qqH");
  procnames.push_back("ttH");
  procnames.push_back("VH");

  
  //--- book histograms
  TH1F *h[4][9];
  TH1F *hs[4][9];
  char histo[100];
  for (int iproc=0; iproc < 4; iproc++){
    for (int icat=0; icat < ncats; icat++){
      // -- histograms
      sprintf(histo,"h_%s_cat%d", (procnames.at(iproc)).c_str(), icat);
      h[iproc][icat]= new TH1F(histo, histo, 100, -1,1);
      // -- histograms after re-weighting for PU-jetID
      sprintf(histo,"hs_%s_cat%d",(procnames.at(iproc)).c_str(), icat);
      hs[iproc][icat]= new TH1F(histo, histo, 100, -1,1);
    }
  }

  //--- load data/MC scale factors
  TFile *f = TFile::Open("JetIdScaleFactor_ZmumuJets_12fb_hcp.root"); // new 
  //TFile *f = TFile::Open("../../../AnalysisScripts/aux/JetIdScaleFactor_ZmumuJets.root"); // scale factors used for ICHEP
  TH1F *hJetIdScaleFactor[4];
  hJetIdScaleFactor[0]=(TH1F*)f->Get("hJetIdScaleFactor_TK");
  hJetIdScaleFactor[1]=(TH1F*)f->Get("hJetIdScaleFactor_HEin");
  hJetIdScaleFactor[2]=(TH1F*)f->Get("hJetIdScaleFactor_HEout");
  hJetIdScaleFactor[3]=(TH1F*)f->Get("hJetIdScaleFactor_HF");


  float w1 = 1.; // weight for leading jet
  float w2 = 1.; // weight for sub-leading jet
  float w  = 1.;
  int proc = 0;
  int cat  = 0;
  int r1, r2;

  //--- loop over entries
  for (int ientry = 0 ; ientry < tree->GetEntries(); ientry++ ){

    tree->GetEntry(ientry);

    if (ientry%100000==0) cout << "Analyzing entry:" << ientry << endl;

    //-- 
    if (!(leadJPt>0 && subleadJPt>0)) continue; 
    
    w1=1;
    w2=1;
    w =1;
    r1=-1;
    r2=-1;

    //--- test effect of PU-jetID uncertainties on cut based di-jet tag
    //         bool tight = false;
    //         bool loose = false;
    //         category=-1;    
    //         if ( leadPtOverM> 60./120. && subleadPtOverM > 30./120. && leadJPt>30. && subleadJPt>30. && fabs(deltaEtaJJ)>3 && Zep<2.5 && MJJ>500 && fabs(deltaPhiJJGamGam)>2.6) {
    //           tight = true;
    //           category = 4;
    //         }
    //         if ( leadPtOverM> 60./120. && subleadPtOverM > 30./120. && leadJPt>30. && subleadJPt>20. && fabs(deltaEtaJJ)>3 && Zep<2.5 && MJJ>250 && fabs(deltaPhiJJGamGam)>2.6 && !tight) {
    //           loose = true;
    //           category = 5;
    //         }
    
    
    //-- select event in category 4 or 5 (di-jet tagged categories)
    if (category==4 || category==5){ 
      //-- find the jet eta pseudorapidity region
      r1=etaregion(leadJEta);
      r2=etaregion(subleadJEta);
      if (r1!=-1) {
	w1 = hJetIdScaleFactor[r1]->GetBinContent(hJetIdScaleFactor[r1]->FindBin(leadJPt));
	if (leadJPt>100) w1 = hJetIdScaleFactor[r1]->GetBinContent(hJetIdScaleFactor[r1]->GetNbinsX());
      }
      if (r2!=-1) {
	w2 = hJetIdScaleFactor[r2]->GetBinContent(hJetIdScaleFactor[r2]->FindBin(subleadJPt));
	if (subleadJPt>100) w2 = hJetIdScaleFactor[r2]->GetBinContent(hJetIdScaleFactor[r2]->GetNbinsX());
      }
      
      // -- to be conservative if the data/MC scale factor is > 1, do not correct. 
      if (w1>1) w1=1.;
      if (w2>1) w2=1.;
      w  = w1*w2;
      
      //-- fill histograms with weights
      proc = abs (sampleType+37);
      h[proc][category]->Fill(0, weight);
      hs[proc][category]->Fill(0, (w*weight));
    }
  }
  
  cout << "process   category   N   Nweighted    1-Nweighted/N "  << endl; 
  float r;
  for (int iproc=0; iproc < 4; iproc++){
    for (int icat=4; icat < 6; icat++){
      r=0;
      if (hs[iproc][icat]-> GetSumOfWeights()>0 && h[iproc][icat]-> GetSumOfWeights()>0)
	r = 1.-hs[iproc][icat]-> GetSumOfWeights()/h[iproc][icat]-> GetSumOfWeights() ;
      cout << (procnames.at(iproc)).c_str() << "        " << icat << "      " 
	   << h[iproc][icat] -> GetSumOfWeights() << "    " 
	   << hs[iproc][icat]-> GetSumOfWeights() << "    "
	   << r << " " 
	   << endl;
    }
  }
  
  

}
