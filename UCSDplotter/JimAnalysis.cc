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
  tree->SetBranchAddress("w" , &full_weight); 
  tree->SetBranchAddress("pu_n"	       , &pu_n);	       
  tree->SetBranchAddress("mass"	       , &mass);	       
  tree->SetBranchAddress("vbfcat"      , &vbfcat);      
  tree->SetBranchAddress("MET"	       , &MET);	       
  tree->SetBranchAddress("MET_phi"     , &MET_phi);     
  tree->SetBranchAddress("bdt"	       , &bdt);	       
  tree->SetBranchAddress("diphocat2r92eta"     , &bdt_cat);      

  defineLabels();
}

//----------------------------------------------------------------------

void JimAnalysis::analyze(HistoContainer* container, RooContainer* rooContainer) {

  // cout << "run=" << run << endl;
  // container->Fill("run", 0, run, full_weight);

  // fill mass per BDT category
  if (bdt_cat >= 0) {
    container->Fill("reco_mass", int(itype + 0.5), int(bdt_cat + 0.5), mass, full_weight);
    FillRooContainer(rooContainer, true);
  }
}

//----------------------------------------------------------------------
extern "C" GenericAnalysis *makeInstance()
{
  return new JimAnalysis();
}
//----------------------------------------------------------------------

void JimAnalysis::FillRooContainer(RooContainer* rooContainer, bool isCorrectVertex) {

  if (itype == 0 ) {
    rooContainer->InputDataPoint("data_mass", bdt_cat, mass);
  } else if (itype < 0 ) {
    //if(doMcOptimization) {
    //rooContainer->InputDataPoint("data_mass", bdt_cat, mass, full_weight);
    //} else if (itype != 3 && itype != 4 ) {
    //rooContainer->InputDataPoint("bkg_mass", bdt_cat, mass, full_weight);
    //}
    // FIXME FOR NEW ITYPES
  } else if (itype > 0) {
    rooContainer->InputDataPoint("sig_" + GetSignalLabel(itype), bdt_cat, mass, full_weight);
    if (isCorrectVertex) 
      rooContainer->InputDataPoint("sig_" + GetSignalLabel(itype) + "_rv", bdt_cat, mass, full_weight);
    else 
      rooContainer->InputDataPoint("sig_" + GetSignalLabel(itype) + "_wv", bdt_cat, mass, full_weight);
  }
}

void JimAnalysis::FillRooContainerSyst(RooContainer* rooContainer, const std::string &name, int cur_type,
				       std::vector<double> & mass_errors, std::vector<double> & mva_errors,
				       std::vector<int>    & categories, std::vector<double> & weights) {    
  if (itype < 0){
    // feed the modified signal model to the RooContainer
    rooContainer->InputSystematicSet("sig_" + GetSignalLabel(cur_type), name, categories, mass_errors, weights);
  }
}
    

std::string JimAnalysis::GetSignalLabel(int id) {

  if (id == 100001) id = -69;
  if (id == 100002) id = -70;
  if (id == 100003) id = -72;
  if (id == 100004) id = -71;

  if (id == 105001) id = -13;
  if (id == 105002) id = -14;
  if (id == 105003) id = -16;
  if (id == 105004) id = -15;

  if (id == 115001) id = -21;
  if (id == 115002) id = -22;
  if (id == 115003) id = -24;
  if (id == 115004) id = -23;

  if (id == 120001) id = -25;
  if (id == 120002) id = -26;
  if (id == 120003) id = -28;
  if (id == 120004) id = -27;

  if (id == 110001) id = -17;
  if (id == 110002) id = -18;
  if (id == 110003) id = -20;
  if (id == 110004) id = -19;

  if (id == 121001) id = -53;
  if (id == 121002) id = -54;
  if (id == 121003) id = -56;
  if (id == 121004) id = -55;

  if (id == 123001) id = -57;
  if (id == 123002) id = -58;
  if (id == 123003) id = -60;
  if (id == 123004) id = -59;

  if (id == 124001) id = -37;
  if (id == 124002) id = -38;
  if (id == 124003) id = -40;
  if (id == 124004) id = -39;

  if (id == 130001) id = -29;
  if (id == 130002) id = -30;
  if (id == 130003) id = -32;
  if (id == 130004) id = -31;

  if (id == 135001) id = -41;
  if (id == 135002) id = -42;
  if (id == 135003) id = -44;
  if (id == 135004) id = -43;

  if (id == 140001) id = -33;
  if (id == 140002) id = -34;
  if (id == 140003) id = -36;
  if (id == 140004) id = -35;

  if (id == 145001) id = -45;
  if (id == 145002) id = -46;
  if (id == 145003) id = -48;
  if (id == 145004) id = -47;

  if (id == 150001) id = -49;
  if (id == 150002) id = -50;
  if (id == 150003) id = -52;
  if (id == 150004) id = -51;

  // Basically A Map of the ID (type) to the signal's name which can be filled Now:
  // For the lazy man, can return a memeber of the map rather than doing it yourself
  std::map<int,std::string>::iterator it = signalLabels.find(id);
  if (it!=signalLabels.end()){
    return it->second;
  } else { 
    std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
    return "NULL";
  }  
}

void JimAnalysis::defineLabels() {

  signalLabels[-57]="ggh_mass_m123";
  signalLabels[-58]="vbf_mass_m123";
  signalLabels[-60]="wzh_mass_m123";
  signalLabels[-59]="tth_mass_m123";
  signalLabels[-53]="ggh_mass_m121";
  signalLabels[-54]="vbf_mass_m121";
  signalLabels[-56]="wzh_mass_m121";
  signalLabels[-55]="tth_mass_m121";
  signalLabels[-65]="ggh_mass_m160";
  signalLabels[-66]="vbf_mass_m160";
  signalLabels[-68]="wzh_mass_m160";
  signalLabels[-67]="tth_mass_m160";
  signalLabels[-61]="ggh_mass_m155";
  signalLabels[-62]="vbf_mass_m155";
  signalLabels[-64]="wzh_mass_m155";
  signalLabels[-63]="tth_mass_m155";
  signalLabels[-49]="ggh_mass_m150";
  signalLabels[-50]="vbf_mass_m150";
  signalLabels[-52]="wzh_mass_m150";
  signalLabels[-51]="tth_mass_m150";
  signalLabels[-45]="ggh_mass_m145";
  signalLabels[-46]="vbf_mass_m145";
  signalLabels[-48]="wzh_mass_m145";
  signalLabels[-47]="tth_mass_m145";
  signalLabels[-33]="ggh_mass_m140";
  signalLabels[-34]="vbf_mass_m140";
  signalLabels[-36]="wzh_mass_m140";
  signalLabels[-35]="tth_mass_m140";
  signalLabels[-41]="ggh_mass_m135";
  signalLabels[-42]="vbf_mass_m135";
  signalLabels[-44]="wzh_mass_m135";
  signalLabels[-43]="tth_mass_m135";
  signalLabels[-29]="ggh_mass_m130";
  signalLabels[-30]="vbf_mass_m130";
  signalLabels[-32]="wzh_mass_m130";
  signalLabels[-31]="tth_mass_m130";
  signalLabels[-37]="ggh_mass_m125";
  signalLabels[-38]="vbf_mass_m125";
  signalLabels[-40]="wzh_mass_m125";
  signalLabels[-39]="tth_mass_m125";
  signalLabels[-25]="ggh_mass_m120";
  signalLabels[-26]="vbf_mass_m120";
  signalLabels[-28]="wzh_mass_m120";
  signalLabels[-27]="tth_mass_m120";
  signalLabels[-21]="ggh_mass_m115";
  signalLabels[-22]="vbf_mass_m115";
  signalLabels[-24]="wzh_mass_m115";
  signalLabels[-23]="tth_mass_m115";
  signalLabels[-17]="ggh_mass_m110";
  signalLabels[-18]="vbf_mass_m110";
  signalLabels[-19]="wzh_mass_m110";
  signalLabels[-20]="tth_mass_m110";
  signalLabels[-13]="ggh_mass_m105";
  signalLabels[-14]="vbf_mass_m105";
  signalLabels[-16]="wzh_mass_m105";
  signalLabels[-15]="tth_mass_m105";
  signalLabels[-69]="ggh_mass_m100";
  signalLabels[-70]="vbf_mass_m100";
  signalLabels[-72]="wzh_mass_m100";
  signalLabels[-71]="tth_mass_m100";
}
