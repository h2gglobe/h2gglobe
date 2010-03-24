#define LoopAll_cxx
#include "LoopAll.h"
#include "TH2.h"
#include "TStyle.h" 
#include "TCanvas.h"

#include "LoopAll_cc.h"
#include "GeneralFunctions_cc.h"

#include "PhotonAnalysis/PhotonAnalysisFunctions_cc.h"


void LoopAll::Loop() {

  if (fChain == 0) 
    return;

  Int_t nentries = Int_t(fChain->GetEntriesFast());
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) 
      break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
  }
}

