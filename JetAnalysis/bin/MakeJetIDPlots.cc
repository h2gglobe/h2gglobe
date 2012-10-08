#include "./JetTree.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"

#include "TMath.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>


// *** 
bool pass_level(int id, int level) { return ( id & (1 << level) ) != 0 ; }

// *** MAIN
int main(int argc, char** argv)
{ 

  if(argc<8){ 
    std::cout << "Usage : " << argv[0]  
	      << "  <input file>  "  
	      << "  <input tree name>  "  
	      << "  <output file name>  " 
	      << "  <eta region: ALL TK HEin HEout HF>" 
	      << "  <min jet pt>  " 
	      << "  <max jet pt>  "  
	      << "  <data flag>  "
	      << std::endl; 
    std::cout <<"Example: " << argv[0] << "jetid_jetanalysis.root flatTree testjetid ALL 20 100 0" << std::endl;
    exit(0); 
  }
   
  std::string inputFile_         = argv[1];
  std::string inputTree_         = argv[2];
  std::string outputRootFileName_= argv[3];
  std::string etaRange_          = argv[4];
  double minJetPt_               = atof(argv[5]);
  double maxJetPt_               = atof(argv[6])  ;
  int dataFlag_                  = atoi(argv[7]);

  outputRootFileName_+="_"+etaRange_+"_pt"+argv[5]+"to"+argv[6]+".root";


  std::cout << "Input parameters: " << std::endl;
  std::cout << "   input file name : "<< argv[1]  << std::endl;
  std::cout << "   input tree name : "<< argv[2]  << std::endl;
  std::cout << "   ouput file name : "<< argv[3]  << std::endl;
  std::cout << "   eta range       : "<< argv[4]  << std::endl;
  std::cout << "   min jet pt      : "<< argv[5]  << std::endl;
  std::cout << "   max jet pt      : "<< argv[6]  << std::endl;
  std::cout << "   is data         : "<< argv[7]  << std::endl;

  float minJetEta_, maxJetEta_;

  if (etaRange_=="TK") {
    minJetEta_=0;
    maxJetEta_=2.5;
  }
  
  if (etaRange_=="HEin") {
    minJetEta_=2.5;
    maxJetEta_=2.75;
  }
  
  if (etaRange_=="HEout") {
    minJetEta_=2.75;
    maxJetEta_=3.0;
  }
  
  if (etaRange_=="HF") {
    minJetEta_=3.0;
    maxJetEta_=5.0;
  }

  if (etaRange_=="ALL") {
    minJetEta_=0.0;
    maxJetEta_=5.0;
  }

  
  //--- PU reweighting
  //edm::LumiReWeighting lumiWeights_ = edm::LumiReWeighting( mcPuFile_.c_str(), dataPuFile_.c_str(), mcPuHisto_.c_str(), dataPuHisto_.c_str());


  std::cout<< "Output file name: " << outputRootFileName_ << std::endl;


  // --- open tree
  std::cout << "Reading tree..." << std::endl;
  std::string treeName = inputTree_ ;
  TChain* chain = new TChain(treeName.c_str());
  chain->Add(inputFile_.c_str());
  JetTree t((chain));

  std::cout << "Number of entries : " << chain-> GetEntries() << std::endl;
 
  //--- book histograms
  TH1F *hNvtx      = new TH1F("hNvtx","Number of vertices", 100, 0, 100.);
  
  // -- jet eta, pt
  TH1F *hjetPt      = new TH1F("hjetPt","jet pt", 100, 0, 100.);
  TH1F *hjetPt_NoPU = new TH1F("hjetPt_NoPU","jet pt (NoPU)", 100, 0, 100.);
  TH1F *hjetPt_PU   = new TH1F("hjetPt_PU","jet pt (PU)", 100, 0, 100.);

  TH1F *hjetEta      = new TH1F("hjetEta","jet eta", 100, -5, 5.);
  TH1F *hjetEta_NoPU = new TH1F("hjetEta_NoPU","jet eta (NoPU)", 100, -5, 5.);
  TH1F *hjetEta_PU   = new TH1F("hjetEta_PU","jet eta (PU)", 100, -5, 5.);

  // -- deltaR
  TH1F *hdR2Mean      = new TH1F("hdR2Mean","dR2Mean", 100, 0, 1.);
  TH1F *hdR2Mean_NoPU = new TH1F("hdR2Mean_NoPU","dR2Mean (NoPU)", 100, 0, 1.);
  TH1F *hdR2Mean_PU   = new TH1F("hdR2Mean_PU","dR2Mean (PU)", 100, 0, 1.);

  // -- fractions in annula of 0.1
  TH1F *hfrac01      = new TH1F("hfrac01","frac01", 100, 0, 1.);
  TH1F *hfrac01_NoPU = new TH1F("hfrac01_NoPU","frac01 (NoPU)", 100, 0, 1.);
  TH1F *hfrac01_PU   = new TH1F("hfrac01_PU","frac01 (PU)", 100, 0, 1.);

  // -- fractions in annula of 0.2
  TH1F *hfrac02      = new TH1F("hfrac02","frac02", 100, 0, 1.);
  TH1F *hfrac02_NoPU = new TH1F("hfrac02_NoPU","frac02 (NoPU)", 100, 0, 1.);
  TH1F *hfrac02_PU   = new TH1F("hfrac02_PU","frac02 (PU)", 100, 0, 1.);

  // -- fractions in annula of 0.3
  TH1F *hfrac03      = new TH1F("hfrac03","frac03", 100, 0, 1.);
  TH1F *hfrac03_NoPU = new TH1F("hfrac03_NoPU","frac03 (NoPU)", 100, 0, 1.);
  TH1F *hfrac03_PU   = new TH1F("hfrac03_PU","frac03 (PU)", 100, 0, 1.);

  // -- fractions in annula of 0.4
  TH1F *hfrac04      = new TH1F("hfrac04","frac04", 100, 0, 1.);
  TH1F *hfrac04_NoPU = new TH1F("hfrac04_NoPU","frac04 (NoPU)", 100, 0, 1.);
  TH1F *hfrac04_PU   = new TH1F("hfrac04_PU","frac04 (PU)", 100, 0, 1.);

  // -- fractions in annula of 0.5
  TH1F *hfrac05      = new TH1F("hfrac05","frac05", 100, 0, 1.);
  TH1F *hfrac05_NoPU = new TH1F("hfrac05_NoPU","frac05 (NoPU)", 100, 0, 1.);
  TH1F *hfrac05_PU   = new TH1F("hfrac05_PU","frac05 (PU)", 100, 0, 1.);

  // -- beta and betaStar
  TH1F *hbeta      = new TH1F("hbeta","beta", 400, 0, 1.);
  TH1F *hbeta_NoPU = new TH1F("hbeta_NoPU","beta (NoPU)", 400, 0, 1.);
  TH1F *hbeta_PU   = new TH1F("hbeta_PU","beta (PU)", 400, 0, 1.);

  TH1F *hbetaStar      = new TH1F("hbetaStar","betaStar", 400, 0, 1.);
  TH1F *hbetaStar_NoPU = new TH1F("hbetaStar_NoPU","betaStar (NoPU)", 400, 0, 1.);
  TH1F *hbetaStar_PU   = new TH1F("hbetaStar_PU","betaStar (PU)", 400, 0, 1.);


  // -- mva output

  TH1F *hsimpleDiscriminant      = new TH1F("hsimpleDiscriminant","simpleDiscriminant", 200, -1, 1.);
  TH1F *hsimpleDiscriminant_NoPU = new TH1F("hsimpleDiscriminant_NoPU","simpleDiscriminant (NoPU)", 200, -1, 1.);
  TH1F *hsimpleDiscriminant_PU   = new TH1F("hsimpleDiscriminant_PU","simpleDiscriminant (PU)", 200, -1, 1.);

  TH1F *hfullDiscriminant      = new TH1F("hfullDiscriminant","fullDiscriminant", 200, -1, 1.);
  TH1F *hfullDiscriminant_NoPU = new TH1F("hfullDiscriminant_NoPU","fullDiscriminant (NoPU)", 200, -1, 1.);
  TH1F *hfullDiscriminant_PU   = new TH1F("hfullDiscriminant_PU","fullDiscriminant (PU)", 200, -1, 1.);

  TH1F *hcutbasedDiscriminant      = new TH1F("hcutbasedDiscriminant","cutbasedDiscriminant", 200, -1, 1.);
  TH1F *hcutbasedDiscriminant_NoPU = new TH1F("hcutbasedDiscriminant_NoPU","cutbasedDiscriminant (NoPU)", 200, -1, 1.);
  TH1F *hcutbasedDiscriminant_PU   = new TH1F("hcutbasedDiscriminant_PU","cutbasedDiscriminant (PU)", 200, -1, 1.);


  //--- for efficiency studies
  TH2F *hPtRatio_vs_Dphi = new TH2F("hPtRatio_vs_Dphi","hPtRatio_vs_Dphi",100,-3.2,3.2,100,0,10);

  char hname[100];
  std::string suffLevel[3] = {"Tight","Medium","Loose"};
  std::string suffId[3]    = {"simpleId","fullId","cutbasedId"};


  //--- pt ratio 
  TH1F *hPtRatio_matched_noId  = new TH1F("hPtRatio_matched_noId","Pt Ratio (MC matched jets)",100,0,10); 
  TH1F *hPtRatio_noId          = new TH1F("hPtRatio_noId","Pt Ratio",100,0,10); 
  TH1F* hPtRatio_matched[3][3];
  TH1F* hPtRatio[3][3];
   
  for ( int iid = 0; iid < 3 ; iid++ ){
    for ( int ilevel = 0; ilevel < 3 ; ilevel++ ){
      sprintf(hname,"hPtRatio_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hPtRatio_matched[iid][ilevel] = new TH1F(hname,hname,100,0,10); 
      sprintf(hname,"hPtRatio_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hPtRatio[iid][ilevel] = new TH1F(hname,hname,100,0,10); 
    }
  }

  //-- jet pt for "true" jets (= matched jets)
  TH1F *hJetPt_matched_noId = new TH1F("hJetPt_matched_noId","hJetPt_matched",100,0,100);
  
  TH1F *hJetPt_noId         = new TH1F("hJetPt_noId","hJetPt",100,0,100);
  TH1F *hJetPt_noId_bkg     = new TH1F("hJetPt_noId_bkg","hJetPt",100,0,100);
  TH1F *hJetPt_noId_bkgSub  = new TH1F("hJetPt_noId_bkgSub","hJetPt",100,0,100);
  
  TH1F *hJetPt_matched[3][3];
  TH1F *hJetPt[3][3];
  TH1F *hJetPt_bkg[3][3];
  TH1F *hJetPt_bkgSub[3][3];
  for ( int iid = 0; iid < 3 ; iid++ ){
    for ( int ilevel = 0; ilevel < 3 ; ilevel++ ){
      sprintf(hname,"hJetPt_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetPt_matched[iid][ilevel] = new TH1F(hname,hname,100,0,100);
      sprintf(hname,"hJetPt_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetPt[iid][ilevel] = new TH1F(hname,hname,100,0,100);
      sprintf(hname,"hJetPt_bkg_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetPt_bkg[iid][ilevel] = new TH1F(hname,hname,100,0,100);
      sprintf(hname,"hJetPt_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetPt_bkgSub[iid][ilevel] = new TH1F(hname,hname,100,0,100);
    }
  }

  //-- jet eta for "true" jets (= matched jets)
  TH1F *hJetEta_matched_noId = new TH1F("hJetEta_matched_noId","hJetEta_matched",100,-5,5);
  TH1F *hJetEta_noId         = new TH1F("hJetEta_noId","hJetEta",100,-5,5);
  TH1F *hJetEta_noId_bkg     = new TH1F("hJetEta_noId_bkg","hJetEta",100,-5,5);
  TH1F *hJetEta_noId_bkgSub  = new TH1F("hJetEta_noId_bkgSub","hJetEta",100,-5,5);  

  TH1F *hJetEta_matched[3][3];
  TH1F *hJetEta[3][3];
  TH1F *hJetEta_bkg[3][3];
  TH1F *hJetEta_bkgSub[3][3];
  for ( int iid = 0; iid < 3 ; iid++ ){
    for ( int ilevel = 0; ilevel < 3 ; ilevel++ ){
      sprintf(hname,"hJetEta_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetEta_matched[iid][ilevel] = new TH1F(hname,hname,100,-5,5);
      sprintf(hname,"hJetEta_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetEta[iid][ilevel] = new TH1F(hname,hname,100,-5,5);
      sprintf(hname,"hJetEta_bkg_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetEta_bkg[iid][ilevel] = new TH1F(hname,hname,100,-5,5);
      sprintf(hname,"hJetEta_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetEta_bkgSub[iid][ilevel] = new TH1F(hname,hname,100,-5,5);
    }
  }

  //-- nvtx
  TH1F *hNumberOfVertices_matched_noId = new TH1F("hNumberOfVertices_matched_noId","hNumberOfVertices_matched",100,0,100);
  TH1F *hNumberOfVertices_noId         = new TH1F("hNumberOfVertices_noId","hNumberOfVertices",100,0,100);
  TH1F *hNumberOfVertices_noId_bkg     = new TH1F("hNumberOfVertices_noId_bkg","hNumberOfVertices",100,0,100);
  TH1F *hNumberOfVertices_noId_bkgSub  = new TH1F("hNumberOfVertices_noId_bkgSub","hNumberOfVertices",100,0,100);

  TH1F *hNumberOfVertices_matched[3][3];
  TH1F *hNumberOfVertices[3][3];
  TH1F *hNumberOfVertices_bkg[3][3];
  TH1F *hNumberOfVertices_bkgSub[3][3];
  for ( int iid = 0; iid < 3 ; iid++ ){
    for ( int ilevel = 0; ilevel < 3 ; ilevel++ ){
      sprintf(hname,"hNumberOfVertices_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hNumberOfVertices_matched[iid][ilevel] = new TH1F(hname,hname,100,0,100);
      sprintf(hname,"hNumberOfVertices_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hNumberOfVertices[iid][ilevel] = new TH1F(hname,hname,100,0,100);
      sprintf(hname,"hNumberOfVertices_bkg_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hNumberOfVertices_bkg[iid][ilevel] = new TH1F(hname,hname,100,0,100);
      sprintf(hname,"hNumberOfVertices_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hNumberOfVertices_bkgSub[iid][ilevel] = new TH1F(hname,hname,100,0,100);
    }
  }


  //-- dphi - to be used for background subtraction
  TH1F *hDphi_noId         = new TH1F("hDphi_noId","hDphi",300,0,3.2);
  TH1F *hDphi_noId_PU      = new TH1F("hDphi_noId_PU","hDphi",300,0,3.2);
  TH1F *hDphi_noId_noPU      = new TH1F("hDphi_noId_noPU","hDphi",300,0,3.2);
  TH1F *hDphi[3][3];
  for ( int iid = 0; iid < 3 ; iid++ ){
    for ( int ilevel = 0; ilevel < 3 ; ilevel++ ){
      sprintf(hname,"hDphi_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hDphi[iid][ilevel] = new TH1F(hname,hname,300,0,3.2);
    }
  }

  //--- profile of  ptratio vs ptz 
  TProfile *pPtRatio_noId = new TProfile("pPtRatio_noId","Pt Ratio vs ptZ",100,0,100,0,10); 
  TProfile *pPtRatio[3][3];
   
  for ( int iid = 0; iid < 3 ; iid++ ){
    for ( int ilevel = 0; ilevel < 3 ; ilevel++ ){
      sprintf(hname,"pPtRatio_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      pPtRatio[iid][ilevel] = new TProfile(hname,hname,100,0,100,0,10); 
    }
  }





  float w = 1;  

  float dphiCut     = 2.5;
  float dphiCutB    = 1.5;
  float scaleFactor =  (TMath::Pi() - dphiCut)/dphiCutB;

  for (int ientry = 0; ientry  < chain->GetEntries(); ientry++ ){
    
    t.GetEntry(ientry);

    if (t.isData!=dataFlag_) continue;
    
    if (ientry%100000==0) std::cout << "Analyzing entry : " << ientry << "\r" << std::flush;

    if (!t.jetLooseID) continue;
    
    if ( t.jetPt < minJetPt_ ) continue;
    if ( t.jetPt > maxJetPt_ ) continue;

    if ( fabs(t.jetEta) < minJetEta_ ) continue;
    if ( fabs(t.jetEta) > maxJetEta_ ) continue;

    
    if (dataFlag_) w = 1;
    else w=t.puweight;
      
    // -- fill histograms for leading jet
    if ( t.ijet==1 ) {
      hNvtx->Fill(t.nvtx,w);
      float ptratio = 0;
      if ( t.dimuonPt > 0. ) ptratio = t.jetPt/t.dimuonPt;

      if ( ptratio > 0. )
	hPtRatio_vs_Dphi -> Fill(t.dphiZJet,ptratio,w);
      
      // -- fill ptRatio 
      if (ptratio > 0. && fabs(t.dphiZJet)>dphiCut){
	
	// -- fill plost for jets matching to gen jets
	if ( dataFlag_==0 && ptratio > 0. && t.isMatched && t.jetGenPt > 10.){
	  hPtRatio_matched_noId -> Fill(ptratio,w);  
	  for (int ilevel = 0; ilevel < 3; ilevel++){
	    if ( pass_level(t.simpleId,ilevel)  )    hPtRatio_matched[0][ilevel] -> Fill(ptratio,w);   
	    if ( pass_level(t.fullId,ilevel)    )    hPtRatio_matched[1][ilevel] -> Fill(ptratio,w);   
	    if ( pass_level(t.cutbasedId,ilevel))    hPtRatio_matched[2][ilevel] -> Fill(ptratio,w);   
	  }
	}
	
	// -- fill plost for selected jets
	hPtRatio_noId -> Fill(ptratio,w);  
	for (int ilevel = 0; ilevel < 3; ilevel++){
	  if ( pass_level(t.simpleId,ilevel)  )  hPtRatio[0][ilevel]-> Fill(ptratio,w);   
	  if ( pass_level(t.fullId,ilevel)    )  hPtRatio[1][ilevel]-> Fill(ptratio,w);   
	  if ( pass_level(t.cutbasedId,ilevel))  hPtRatio[2][ilevel]-> Fill(ptratio,w);   
	}
      }

 //      // for TEST
//       // -- fill plots for selected jets matching to gen jets
//       if ( dataFlag_==0 && ptratio > 0. && t.isMatched && t.jetGenPt > 10.){
// 	hJetPt_matched_noId   -> Fill(t.jetPt,w);  
// 	hJetEta_matched_noId  -> Fill(t.jetEta,w);
// 	  hNumberOfVertices_matched_noId->Fill(t.nvtx,w);
	  
// 	  for (int ilevel = 0; ilevel < 3; ilevel++){
// 	    if ( pass_level(t.simpleId,ilevel)  )  {
// 	      hJetPt_matched[0][ilevel]   -> Fill(t.jetPt,w);  
// 	      hJetEta_matched[0][ilevel]  -> Fill(t.jetEta,w);
// 	      hNumberOfVertices_matched[0][ilevel]-> Fill(t.nvtx,w);
// 	    }
// 	    if ( pass_level(t.fullId,ilevel)    )  {
// 	      hJetPt_matched[1][ilevel]   -> Fill(t.jetPt,w);  
// 	      hJetEta_matched[1][ilevel]  -> Fill(t.jetEta,w);
// 	      hNumberOfVertices_matched[1][ilevel]-> Fill(t.nvtx,w);
// 	    }
// 	    if ( pass_level(t.cutbasedId,ilevel))  {
// 	      hJetPt_matched[2][ilevel]   -> Fill(t.jetPt,w);  
// 	      hJetEta_matched[2][ilevel]  -> Fill(t.jetEta,w);
// 	      hNumberOfVertices_matched[2][ilevel]-> Fill(t.nvtx,w);
// 	    }
// 	  }	
//       }
      
      // -- fill jet-Z dphi 
      if ( ptratio > 0.5 && ptratio < 1.5 ){
	hDphi_noId->Fill(fabs(t.dphiZJet),w);
	if ( t.isMatched == 0 ) hDphi_noId_PU->Fill(fabs(t.dphiZJet),w);
	else hDphi_noId_noPU->Fill(fabs(t.dphiZJet),w);
      }
      
      // -- fill jet pT , eta, nvtx plots in the control region
      if ( ptratio > 0.5 && ptratio < 1.5 && fabs(t.dphiZJet) > dphiCut ) {
	
	// -- fill plots for selected jets matching to gen jets
	if ( dataFlag_==0 && ptratio > 0. && t.isMatched && t.jetGenPt > 10.){
	  hJetPt_matched_noId   -> Fill(t.jetPt,w);  
	  hJetEta_matched_noId  -> Fill(t.jetEta,w);
	  hNumberOfVertices_matched_noId->Fill(t.nvtx,w);
	
	  for (int ilevel = 0; ilevel < 3; ilevel++){
	    if ( pass_level(t.simpleId,ilevel)  )  {
	      hJetPt_matched[0][ilevel]   -> Fill(t.jetPt,w);  
	      hJetEta_matched[0][ilevel]  -> Fill(t.jetEta,w);
	      hNumberOfVertices_matched[0][ilevel]-> Fill(t.nvtx,w);
	    }
	    if ( pass_level(t.fullId,ilevel)    )  {
	      hJetPt_matched[1][ilevel]   -> Fill(t.jetPt,w);  
	      hJetEta_matched[1][ilevel]  -> Fill(t.jetEta,w);
	      hNumberOfVertices_matched[1][ilevel]-> Fill(t.nvtx,w);
	    }
	    if ( pass_level(t.cutbasedId,ilevel))  {
	      hJetPt_matched[2][ilevel]   -> Fill(t.jetPt,w);  
	      hJetEta_matched[2][ilevel]  -> Fill(t.jetEta,w);
	      hNumberOfVertices_matched[2][ilevel]-> Fill(t.nvtx,w);
	    }
	  }	
	}
	
	//-- fill plots for all selected jets
	pPtRatio_noId -> Fill(t.dimuonPt, ptratio, w);
	hJetPt_noId   -> Fill(t.jetPt,w);  
	hJetEta_noId  -> Fill(t.jetEta,w);
	hNumberOfVertices_noId->Fill(t.nvtx,w);

	for (int ilevel = 0; ilevel < 3; ilevel++){
	  if ( pass_level(t.simpleId,ilevel) )  {
	    pPtRatio[0][ilevel] -> Fill(t.dimuonPt, ptratio, w);
	    hJetPt[0][ilevel]   -> Fill(t.jetPt,w);  
	    hJetEta[0][ilevel]  -> Fill(t.jetEta,w); 
	    hNumberOfVertices[0][ilevel]->Fill(t.nvtx,w);
	  } 
	  if ( pass_level(t.fullId,ilevel)   )  {
	    pPtRatio[1][ilevel] -> Fill(t.dimuonPt, ptratio, w);
	    hJetPt[1][ilevel]   -> Fill(t.jetPt,w);  
	    hJetEta[1][ilevel]  -> Fill(t.jetEta,w);  
	    hNumberOfVertices[1][ilevel]->Fill(t.nvtx,w);
	  }
	  if ( pass_level(t.cutbasedId,ilevel) ) {
	    pPtRatio[2][ilevel] -> Fill(t.dimuonPt, ptratio, w);
	    hJetPt[2][ilevel]   -> Fill(t.jetPt,w);  
	    hJetEta[2][ilevel]  -> Fill(t.jetEta,w);  
	    hNumberOfVertices[2][ilevel]->Fill(t.nvtx,w);
	  }
	}
      }

      // --- fill plots for backgroun estimation
      if ( ptratio > 0.5 && ptratio < 1.5 && fabs(t.dphiZJet) < dphiCutB ) {
	hJetPt_noId_bkg    -> Fill(t.jetPt,w);  
	hJetEta_noId_bkg   -> Fill(t.jetEta,w);  
	hNumberOfVertices_noId_bkg -> Fill(t.nvtx,w);  

	for (int ilevel = 0; ilevel < 3; ilevel++){
	  if ( pass_level(t.simpleId,ilevel) )  {
	    hJetPt_bkg[0][ilevel] -> Fill(t.jetPt,w);  
	    hJetEta_bkg[0][ilevel]-> Fill(t.jetEta,w); 
	    hNumberOfVertices_bkg[0][ilevel]->Fill(t.nvtx,w);
	  } 
	  if ( pass_level(t.fullId,ilevel)   )  {
	    hJetPt_bkg[1][ilevel] -> Fill(t.jetPt,w);  
	    hJetEta_bkg[1][ilevel]-> Fill(t.jetEta,w);  
	    hNumberOfVertices_bkg[1][ilevel]->Fill(t.nvtx,w);
	  }
	  if ( pass_level(t.cutbasedId,ilevel) ) {
	    hJetPt_bkg[2][ilevel] -> Fill(t.jetPt,w);  
	    hJetEta_bkg[2][ilevel]-> Fill(t.jetEta,w);  
	    hNumberOfVertices_bkg[2][ilevel]->Fill(t.nvtx,w);
	  }
	}
      }
      
    }
    
    //-- kinematic and jet id variables 
    hjetPt        -> Fill(t.jetPt,w);
    hjetEta       -> Fill(t.jetEta,w);
    hdR2Mean      -> Fill(t.dR2Mean,w);
    hfrac01       -> Fill(t.frac01,w);
    hfrac02       -> Fill(t.frac02,w);
    hfrac03       -> Fill(t.frac03,w);
    hfrac04       -> Fill(t.frac04,w);
    hfrac05       -> Fill(t.frac05,w);
    hbeta         -> Fill(t.beta,w);
    hbetaStar     -> Fill(t.betaStar,w);
    hsimpleDiscriminant -> Fill(t.simpleDiscriminant,w);
    hfullDiscriminant   -> Fill(t.fullDiscriminant,w);
    hcutbasedDiscriminant -> Fill(t.cutbasedDiscriminant,w);

    if ( t.isMatched ) {
      hjetPt_NoPU        -> Fill(t.jetPt,w);
      hjetEta_NoPU       -> Fill(t.jetEta,w);
      hdR2Mean_NoPU      -> Fill(t.dR2Mean,w);
      hfrac01_NoPU       -> Fill(t.frac01,w);
      hfrac02_NoPU       -> Fill(t.frac02,w);
      hfrac03_NoPU       -> Fill(t.frac03,w);
      hfrac04_NoPU       -> Fill(t.frac04,w);
      hfrac05_NoPU       -> Fill(t.frac05,w);
      hbeta_NoPU         -> Fill(t.beta,w);
      hbetaStar_NoPU     -> Fill(t.betaStar,w);
      hsimpleDiscriminant_NoPU  -> Fill(t.simpleDiscriminant,w);
      hfullDiscriminant_NoPU    -> Fill(t.fullDiscriminant,w);
      hcutbasedDiscriminant_NoPU  -> Fill(t.cutbasedDiscriminant,w);
    }
    else{
      hjetPt_PU        -> Fill(t.jetPt,w);
      hjetEta_PU       -> Fill(t.jetEta,w);
      hdR2Mean_PU      -> Fill(t.dR2Mean,w);
      hfrac01_PU       -> Fill(t.frac01,w);
      hfrac02_PU       -> Fill(t.frac02,w);
      hfrac03_PU       -> Fill(t.frac03,w);
      hfrac04_PU       -> Fill(t.frac04,w);
      hfrac05_PU       -> Fill(t.frac05,w);
      hbeta_PU         -> Fill(t.beta,w);
      hbetaStar_PU     -> Fill(t.betaStar,w);
      hsimpleDiscriminant_PU  -> Fill(t.simpleDiscriminant,w);
      hfullDiscriminant_PU    -> Fill(t.fullDiscriminant,w);
      hcutbasedDiscriminant_PU  -> Fill(t.cutbasedDiscriminant,w);
    }

  }// end loop over entries

  // scale bkg 
  
  //--- compute efficiencies
  TH1F *hEff_vs_PtRatio_matched[3][3]; 
  TH1F *hEff_vs_PtRatio[3][3]; 
  
  TH1F *hEff_vs_JetPt_matched[3][3]; 
  TH1F *hEff_vs_JetPt[3][3]; 
  TH1F *hEff_vs_JetPt_bkgSub[3][3]; 

  TH1F *hEff_vs_JetEta_matched[3][3]; 
  TH1F *hEff_vs_JetEta[3][3]; 
  TH1F *hEff_vs_JetEta_bkgSub[3][3]; 
  
  TH1F *hEff_vs_NumberOfVertices_matched[3][3]; 
  TH1F *hEff_vs_NumberOfVertices[3][3]; 
  TH1F *hEff_vs_NumberOfVertices_bkgSub[3][3]; 

    
  hPtRatio_noId        ->Sumw2();
  hPtRatio_matched_noId->Sumw2();
  hJetPt_noId          ->Sumw2();
  hJetPt_matched_noId  ->Sumw2();
  hJetEta_noId         ->Sumw2();
  hJetEta_matched_noId ->Sumw2();


  // -- bkg subtracted histograms
  hJetPt_noId_bkgSub  -> Sumw2();
  hJetEta_noId_bkgSub -> Sumw2();
  hNumberOfVertices_noId_bkgSub -> Sumw2();

  hJetPt_noId_bkgSub = (TH1F*)hJetPt_noId->Clone("hJetPt_noId_bkgSub");
  hJetPt_noId_bkg    -> Scale(scaleFactor);
  hJetPt_noId_bkgSub -> Add(hJetPt_noId_bkgSub, hJetPt_noId_bkg,1, -1.);

  hJetEta_noId_bkgSub = (TH1F*)hJetEta_noId->Clone("hJetEta_noId_bkgSub");
  hJetEta_noId_bkg    -> Scale(scaleFactor);
  hJetEta_noId_bkgSub -> Add(hJetEta_noId_bkgSub, hJetEta_noId_bkg,1, -1.);
  
  hNumberOfVertices_noId_bkgSub = (TH1F*)hNumberOfVertices_noId->Clone("hNumberOfVertices_noId_bkgSub");
  hNumberOfVertices_noId_bkg    -> Scale(scaleFactor);
  hNumberOfVertices_noId_bkgSub -> Add(hNumberOfVertices_noId_bkgSub, hNumberOfVertices_noId_bkg,1, -1.);

  for ( int iid = 0; iid < 3 ; iid++ ){
    for (int ilevel = 0; ilevel < 3; ilevel++){
  
      // -  matched jets --> "true" efficiency
      if (dataFlag_==0){
	hPtRatio_matched[iid][ilevel]->Sumw2();
	sprintf(hname,"hEff_vs_PtRatio_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
	hEff_vs_PtRatio_matched[iid][ilevel]=(TH1F*)hPtRatio_matched[iid][ilevel]->Clone(hname);
	hEff_vs_PtRatio_matched[iid][ilevel]->Divide(hEff_vs_PtRatio_matched[iid][ilevel],hPtRatio_matched_noId,1,1,"B");
	
	hJetPt_matched[iid][ilevel]->Sumw2();
	sprintf(hname,"hEff_vs_JetPt_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
	hEff_vs_JetPt_matched[iid][ilevel]=(TH1F*)hJetPt_matched[iid][ilevel]->Clone(hname);
	hEff_vs_JetPt_matched[iid][ilevel]->Divide(hEff_vs_JetPt_matched[iid][ilevel],hJetPt_matched_noId,1,1,"B");

	hJetEta_matched[iid][ilevel]->Sumw2();
	sprintf(hname,"hEff_vs_JetEta_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
	hEff_vs_JetEta_matched[iid][ilevel]=(TH1F*)hJetEta_matched[iid][ilevel]->Clone(hname);
	hEff_vs_JetEta_matched[iid][ilevel]->Divide(hEff_vs_JetEta_matched[iid][ilevel],hJetEta_matched_noId,1,1,"B");

	hNumberOfVertices_matched[iid][ilevel]->Sumw2();
	sprintf(hname,"hEff_vs_NumberOfVertices_matched_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
	hEff_vs_NumberOfVertices_matched[iid][ilevel]=(TH1F*)hNumberOfVertices_matched[iid][ilevel]->Clone(hname);
	hEff_vs_NumberOfVertices_matched[iid][ilevel]->Divide(hEff_vs_NumberOfVertices_matched[iid][ilevel],hNumberOfVertices_matched_noId,1,1,"B");
      }

      hPtRatio[iid][ilevel]->Sumw2();
      sprintf(hname,"hEff_vs_PtRatio_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hEff_vs_PtRatio[iid][ilevel]=(TH1F*)hPtRatio[iid][ilevel]->Clone(hname);
      hEff_vs_PtRatio[iid][ilevel]->Divide(hEff_vs_PtRatio[iid][ilevel],hPtRatio_noId,1,1,"B");
      
      hJetPt[iid][ilevel]->Sumw2();
      sprintf(hname,"hEff_vs_JetPt_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hEff_vs_JetPt[iid][ilevel]=(TH1F*)hJetPt[iid][ilevel]->Clone(hname);
      hEff_vs_JetPt[iid][ilevel]->Divide(hEff_vs_JetPt[iid][ilevel],hJetPt_noId,1,1,"B");
      
      hJetEta[iid][ilevel]->Sumw2();
      sprintf(hname,"hEff_vs_JetEta_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hEff_vs_JetEta[iid][ilevel]=(TH1F*)hJetEta[iid][ilevel]->Clone(hname);
      hEff_vs_JetEta[iid][ilevel]->Divide(hEff_vs_JetEta[iid][ilevel],hJetEta_noId,1,1,"B");
      
      hNumberOfVertices[iid][ilevel]->Sumw2();
      sprintf(hname,"hEff_vs_NumberOfVertices_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hEff_vs_NumberOfVertices[iid][ilevel]=(TH1F*)hNumberOfVertices[iid][ilevel]->Clone(hname);
      hEff_vs_NumberOfVertices[iid][ilevel]->Divide(hEff_vs_NumberOfVertices[iid][ilevel],hNumberOfVertices_noId,1,1,"B");
      

      // --- compute efficiencies after background subtraction
      // -- pt
      sprintf(hname,"hJetPt_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetPt_bkgSub[iid][ilevel]= (TH1F*)hJetPt[iid][ilevel]->Clone(hname);
      hJetPt_bkgSub[iid][ilevel]->Sumw2();

      hJetPt_bkg[iid][ilevel]   -> Scale(scaleFactor);
      hJetPt_bkgSub[iid][ilevel]-> Add(hJetPt_bkgSub[iid][ilevel], hJetPt_bkg[iid][ilevel], 1., -1.);
      sprintf(hname,"hEff_vs_JetPt_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hEff_vs_JetPt_bkgSub[iid][ilevel]=(TH1F*)hJetPt_bkgSub[iid][ilevel]->Clone(hname);
      hEff_vs_JetPt_bkgSub[iid][ilevel]->Divide(hEff_vs_JetPt_bkgSub[iid][ilevel],hJetPt_noId_bkgSub,1,1,"B");
      
      // -- eta
      sprintf(hname,"hJetEta_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hJetEta_bkgSub[iid][ilevel]= (TH1F*)hJetEta[iid][ilevel]->Clone(hname);
      hJetEta_bkgSub[iid][ilevel]->Sumw2();

      hJetEta_bkg[iid][ilevel]   -> Scale(scaleFactor);
      hJetEta_bkgSub[iid][ilevel]-> Add(hJetEta_bkgSub[iid][ilevel], hJetEta_bkg[iid][ilevel], 1., -1.);
      sprintf(hname,"hEff_vs_JetEta_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hEff_vs_JetEta_bkgSub[iid][ilevel]=(TH1F*)hJetEta_bkgSub[iid][ilevel]->Clone(hname);
      hEff_vs_JetEta_bkgSub[iid][ilevel]->Divide(hEff_vs_JetEta_bkgSub[iid][ilevel],hJetEta_noId_bkgSub,1,1,"B");
 
      // -- number of vertices
      sprintf(hname,"hNumberOfVertices_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hNumberOfVertices_bkgSub[iid][ilevel]= (TH1F*)hNumberOfVertices[iid][ilevel]->Clone(hname);
      hNumberOfVertices_bkgSub[iid][ilevel]->Sumw2();

      hNumberOfVertices_bkg[iid][ilevel]   -> Scale(scaleFactor);
      hNumberOfVertices_bkgSub[iid][ilevel]-> Add(hNumberOfVertices_bkgSub[iid][ilevel], hNumberOfVertices_bkg[iid][ilevel], 1., -1.);
      sprintf(hname,"hEff_vs_NumberOfVertices_bkgSub_%s_%s", suffId[iid].c_str(), suffLevel[ilevel].c_str());
      hEff_vs_NumberOfVertices_bkgSub[iid][ilevel]=(TH1F*)hNumberOfVertices_bkgSub[iid][ilevel]->Clone(hname);
      hEff_vs_NumberOfVertices_bkgSub[iid][ilevel]->Divide(hEff_vs_NumberOfVertices_bkgSub[iid][ilevel],hNumberOfVertices_noId_bkgSub,1,1,"B");
 




    }
  }

  
  // save the histograms
 
  //std::cout << "Saving histograms on file : " << outputRootFilePath_+outputRootFileName_ << std::endl;
  //  TFile* outputRootFile = new TFile((outputRootFilePath_+outputRootFileName_).c_str(), "RECREATE");

  std::cout << "Saving histograms on file : " << outputRootFileName_ << std::endl;
  TFile* outputRootFile = new TFile((outputRootFileName_).c_str(), "RECREATE");

  outputRootFile -> cd();

  hNvtx ->Write();
  
  hjetPt        -> Write();
  hjetEta       -> Write();
  hdR2Mean      -> Write();
  hfrac01       -> Write();
  hfrac02       -> Write();
  hfrac03       -> Write();
  hfrac04       -> Write();
  hfrac05       -> Write();
  hbeta         -> Write();
  hbetaStar     -> Write();
  hsimpleDiscriminant->Write();
  hfullDiscriminant  ->Write();
  hcutbasedDiscriminant->Write();


  hjetPt_NoPU       -> Write();
  hjetEta_NoPU      -> Write();
  hdR2Mean_NoPU     -> Write();
  hfrac01_NoPU      -> Write();
  hfrac02_NoPU      -> Write();
  hfrac03_NoPU      -> Write();
  hfrac04_NoPU      -> Write();
  hfrac05_NoPU      -> Write();
  hbeta_NoPU        -> Write();
  hbetaStar_NoPU    -> Write();
  hsimpleDiscriminant_NoPU->Write();
  hfullDiscriminant_NoPU  ->Write();
  hcutbasedDiscriminant_NoPU->Write();

  hjetPt_PU        -> Write();
  hjetEta_PU       -> Write();
  hdR2Mean_PU      -> Write();
  hfrac01_PU       -> Write();
  hfrac02_PU       -> Write();
  hfrac03_PU       -> Write();
  hfrac04_PU       -> Write();
  hfrac05_PU       -> Write();
  hbeta_PU         -> Write();
  hbetaStar_PU     -> Write();
  hsimpleDiscriminant_PU->Write();
  hfullDiscriminant_PU  ->Write();
  hcutbasedDiscriminant_PU->Write();

  hDphi_noId->Write();
  hDphi_noId_PU->Write();
  hDphi_noId_noPU->Write();


  TDirectory *efficiency  =   outputRootFile -> mkdir("efficiency","efficiency");
  efficiency->cd();

  hPtRatio_matched_noId -> Write();   
  hJetPt_matched_noId -> Write();   
  hJetEta_matched_noId -> Write();   
  hNumberOfVertices_matched_noId -> Write();   

  hPtRatio_vs_Dphi ->Write();
  hPtRatio_noId -> Write();   
  hJetPt_noId  -> Write();   
  hJetEta_noId -> Write();   
  hNumberOfVertices_noId -> Write();   

  pPtRatio_noId -> Write();
 

  hJetPt_noId_bkg     -> Write();  
  hJetPt_noId_bkgSub  -> Write();  
  hJetEta_noId_bkg    -> Write();  
  hJetEta_noId_bkgSub -> Write();  
  hNumberOfVertices_noId_bkg    -> Write();  
  hNumberOfVertices_noId_bkgSub -> Write();  


  for ( int iid = 0; iid < 3 ; iid++ ){
    for (int ilevel = 0; ilevel < 3; ilevel++){
      if ( dataFlag_==0 ){
	hPtRatio_matched[iid][ilevel]->Write();   
	hJetPt_matched[iid][ilevel]->Write();   
	hJetEta_matched[iid][ilevel]->Write();   
	hNumberOfVertices_matched[iid][ilevel]->Write();   
	
	hEff_vs_PtRatio_matched[iid][ilevel]->Write(); 
	hEff_vs_JetPt_matched[iid][ilevel]->Write(); 
	hEff_vs_JetEta_matched[iid][ilevel]->Write(); 
	hEff_vs_NumberOfVertices_matched[iid][ilevel]->Write(); 
      }
      
      pPtRatio[iid][ilevel]->Write();   
      hPtRatio[iid][ilevel]->Write();   
      hJetPt[iid][ilevel]->Write();   
      hJetEta[iid][ilevel]->Write();   
      hNumberOfVertices[iid][ilevel]->Write();   
      
      hEff_vs_PtRatio[iid][ilevel]-> Write(); 
      hEff_vs_JetPt[iid][ilevel]  -> Write(); 
      hEff_vs_JetEta[iid][ilevel] -> Write(); 
      hEff_vs_NumberOfVertices[iid][ilevel]-> Write(); 

      // --- bgk subtraction plots
      hJetPt_bkg[iid][ilevel] -> Write();   
      hJetEta_bkg[iid][ilevel]-> Write();   
      hNumberOfVertices_bkg[iid][ilevel]-> Write();   

      hJetPt_bkgSub[iid][ilevel]  -> Write();   
      hJetEta_bkgSub[iid][ilevel] -> Write();   
      hNumberOfVertices_bkgSub[iid][ilevel]-> Write();   
      
      hEff_vs_JetPt_bkgSub[iid][ilevel]  -> Write(); 
      hEff_vs_JetEta_bkgSub[iid][ilevel] -> Write(); 
      hEff_vs_NumberOfVertices_bkgSub[iid][ilevel] -> Write(); 

   }
  }

  outputRootFile -> Close();
    
  return 0;


 }

