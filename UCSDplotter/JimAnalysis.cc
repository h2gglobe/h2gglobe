#include "JimAnalysis.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
void JimAnalysis::setBranchAddresses(TTree* tree) {
  
  tree->SetBranchAddress("run"         , &run);	       
  tree->SetBranchAddress("lumis"       , &lumis);       
  tree->SetBranchAddress("event"       , &event);       
  tree->SetBranchAddress("itype"       , &itype);       
  tree->SetBranchAddress("sigmaMrvoM"  , &sigmaMrvoM);  
  tree->SetBranchAddress("sigmaMwvoM"  , &sigmaMwvoM);  
  tree->SetBranchAddress("ptoM1"       , &ptoM1);       
  tree->SetBranchAddress("ptoM2"       , &ptoM2);       
  tree->SetBranchAddress("vtxprob"     , &vtxprob);     
  tree->SetBranchAddress("eta1"	       , &eta1);	       
  tree->SetBranchAddress("eta2"	       , &eta2);	       
  tree->SetBranchAddress("cosphi"      , &cosphi);      
  tree->SetBranchAddress("idmva1"      , &idmva1);      
  tree->SetBranchAddress("idmva2"      , &idmva2);      
  tree->SetBranchAddress("rho"	       , &rho);	       
  tree->SetBranchAddress("xsec_weight" , &xsec_weight); 
  tree->SetBranchAddress("pu_weight"   , &pu_weight);   
  tree->SetBranchAddress("full_weight" , &full_weight); 
  tree->SetBranchAddress("pu_n"	       , &pu_n);	       
  tree->SetBranchAddress("mass"	       , &mass);	       
  tree->SetBranchAddress("vbfcat"      , &vbfcat);      
  tree->SetBranchAddress("MET"	       , &MET);	       
  tree->SetBranchAddress("MET_phi"     , &MET_phi);     
  tree->SetBranchAddress("bdt"	       , &bdt);	       
  tree->SetBranchAddress("bdt_cat"     , &bdt_cat);      
}

//----------------------------------------------------------------------

void JimAnalysis::analyze(HistoContainer* container) {

  // cout << "run=" << run << endl;
  // container->Fill("run", 0, run, full_weight);

  // fill mass per BDT category
  if (bdt_cat >= 0)
  {
    container->Fill("reco_mass", int(itype + 0.5), int(bdt_cat + 0.5), mass, full_weight);
  }
  

}

//----------------------------------------------------------------------
extern "C" GenericAnalysis *makeInstance()
{
  return new JimAnalysis();
}
//----------------------------------------------------------------------
