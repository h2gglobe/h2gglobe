#include "../interface/FutureAnalysis.h"
#include "Sorters.h"
#include "../../LoopAll.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
FutureAnalysis::FutureAnalysis()  
{
    name_ = "FutureAnalysis";
    deno = 0;
    num = 0;

}

// ----------------------------------------------------------------------------------------------------
FutureAnalysis::~FutureAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::Term(LoopAll& l) 
{
  std::cout << "num=" << num << "  deno=" << deno << std::endl;
  double reco_eff = (double(num)/double(deno)) ;
  std::cout << "reconstruction eff : " << reco_eff << std::endl;
 
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::Init(LoopAll& l) 
{
  std::cout << "Init -> FutureAnalysis "<< std::endl;
}

void FutureAnalysis::ResetAnalysis(){}

void FutureAnalysis::GetBranches(TTree * outputTree, std::set<TBranch *>& ) 
{}


// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
  int type = l.itype[l.current];
  
  /// Loop over reco photons ///
  for (int ipho=0; ipho<l.pho_n; ++ipho) {
    int ncat = l.PhotonCategory(ipho, 2, 2);
    TLorentzVector* ph = dynamic_cast<TLorentzVector*>(l.pho_p4->At(ipho));
    FillFlatTree(l, type, ipho);
    l.FillHist("all_pho_pt", ncat, ph->Pt(), 1);
    
    /// loop over all gen particles ///
    for (int igen=0; igen<l.gp_n; igen++) {
      if (l.gp_pdgid[igen]==22 && l.gp_status[igen]==1) { /// select only stable/final state photons ///
	TLorentzVector* GenPho_p4 = dynamic_cast<TLorentzVector*>(l.gp_p4->At(igen));
	float dr = ph->DeltaR(*GenPho_p4); /// deltaR between reco photon and gen photon ///
	l.FillHist("deltaR_all", 0, dr, 1); ///  weight=1
	double pT_recoPho = ph->Pt();
	double pT_genPho = GenPho_p4->Pt();
	double AbsDeltaPt = fabs(pT_recoPho-pT_genPho);
	l.FillHist("AbsDeltaPt", 0, fabs(pT_recoPho-pT_genPho), 1); ///cat=0, weight=1
      }
    } 
  }
  /// Loop over gen particles ///
  for (int i=0; i<l.gp_n; i++) {
    TLorentzVector* p4_gen = dynamic_cast<TLorentzVector*>(l.gp_p4->At(i));
    /// Pick gen photons with pt>20 GeV ///
    if ((p4_gen->Pt())>20.0 && l.gp_status[i]==1 && l.gp_pdgid[i]==22) {
      deno++ ;
      double thisEta = fabs(p4_gen->Eta());
      double thisPt  = fabs(p4_gen->Pt());
      l.FillHist("deno_hist_absEta",0, thisEta, 1);
      l.FillHist("deno_hist_Pt",0, thisPt, 1);

      int mt = MatchedWithReco(*p4_gen,l) ;
      if (mt != -1) {
	num++;
	l.FillHist("num_hist_absEta", 0, thisEta, 1);
	l.FillHist("num_hist_Pt", 0, thisPt, 1);
      }
    }
  }
  return true;
}

int FutureAnalysis::MatchedWithReco(TLorentzVector& GenPho, LoopAll& l) {
  float dr_start = 0.1;
  float abs_dPt_start = 5.0;
  int matchedIndx = -1 ;
  for (int m=0; m<l.pho_n; ++m) {
    TLorentzVector* recoPho = dynamic_cast<TLorentzVector*>(l.pho_p4->At(m));
    float dr = GenPho.DeltaR(*recoPho);
    float abs_dPt = fabs(GenPho.Pt()-recoPho->Pt());
    if (dr<dr_start && abs_dPt<abs_dPt_start) {
      dr_start=dr;
      abs_dPt_start=abs_dPt;
      matchedIndx = m;
    }
  }
  return matchedIndx;
}

// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::SelectEvents(LoopAll&, int)
{
  return true;
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
}
   
// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::SelectEventsReduction(LoopAll&, int)
{
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::SkimEvents(LoopAll&, int)
{
  return true;
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
}

void FutureAnalysis::FillFlatTree(LoopAll& l, Int_t type, Int_t pho_indx) {
  l.FillTree("run", l.run, "Future_myTry");
  l.FillTree("lumis", l.lumis, "Future_myTry");
  l.FillTree("event", l.event, "Future_myTry");
  l.FillTree("itype", type, "Future_myTry");

  if (pho_indx != -1) {
    TLorentzVector* p4_photon = (TLorentzVector*) l.pho_p4->At(pho_indx);
    //    std::cout << "eta " << p4_photon->Eta() << std::endl;
    l.FillTree("eta_pho", (float)p4_photon->Eta(), "Future_myTry");
  }
}

// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
