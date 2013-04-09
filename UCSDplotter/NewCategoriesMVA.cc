#include "NewCategoriesMVA.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
void NewCategoriesMVA::setBranchAddresses(TTree* tree, std::map<std::string, std::string>& values) {
  
  tree->SetBranchAddress("itype"        , &itype);       
  tree->SetBranchAddress("full_weight"  , &full_weight); 
  tree->SetBranchAddress("mass"	        , &mass);	       
  tree->SetBranchAddress("dipho_mva"	  , &bdt);	       
  tree->SetBranchAddress("dipho_mva_cat", &bdt_cat);    

  defineLabels();
}

//----------------------------------------------------------------------

void NewCategoriesMVA::analyze(RooContainer* rooContainer) {
  
  if (bdt_cat >= 0) {
    // refill bdt_cat
    if(bdt_cat>3) {
      bdt_cat+=4;
    } else {
      if(bdt>0.95) bdt_cat=0; 
      else if(bdt>0.90) bdt_cat=1; 
      else if(bdt>0.85) bdt_cat=2; 
      else if(bdt>0.80) bdt_cat=3; 
      else if(bdt>0.65) bdt_cat=4; 
      else if(bdt>0.50) bdt_cat=5; 
      else if(bdt>0.20) bdt_cat=6; 
      else if(bdt>-.05) bdt_cat=7; 
      //bdtCategoryBoundaries=1,0.91,0.79,0.49,-0.05
    }

    //rooContainer->Fill("reco_mass", int(itype + 0.5), int(bdt_cat + 0.5), mass, full_weight);
    
    FillRooContainer(rooContainer);
    
  }
}

//----------------------------------------------------------------------

extern "C" GenericAnalysis *makeInstance() {
  return new NewCategoriesMVA();
}

//----------------------------------------------------------------------

void NewCategoriesMVA::FillRooContainer(RooContainer* rooContainer) {

  if (itype == 0 ) {
    rooContainer->InputDataPoint("data_mass", bdt_cat, mass);
  } else if (itype > 0 ) {
    rooContainer->InputDataPoint("bkg_mass",  bdt_cat, mass, full_weight);
  } else if (itype < 0) {
    rooContainer->InputDataPoint("sig_" + GetSignalLabel(itype), bdt_cat, mass, full_weight);
  }
}

//----------------------------------------------------------------------

void NewCategoriesMVA::FillRooContainerSyst(RooContainer* rooContainer) {

  if (itype < 0) {

    std::vector<double> mass_errors, weights;
    std::vector<int> categories;
    
    mass_errors.push_back(mass);
    weights.push_back(full_weight);
    categories.push_back(bdt_cat);

    rooContainer->InputSystematicSet("sig_" + GetSignalLabel(itype), *name1, categories, mass_errors, weights);
  }
}
    
//----------------------------------------------------------------------

std::string NewCategoriesMVA::GetSignalLabel(int id) {
 
  // Basically A Map of the ID (type) to the signal's name which can be filled Now:
  // For the lazy man, can return a member of the map rather than doing it yourself
  std::map<int,std::string>::iterator it = signalLabels.find(id);
  if (it!=signalLabels.end()) {
    return it->second;
  } else { 
    std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
    return "NULL";
  }  
}

//----------------------------------------------------------------------

void NewCategoriesMVA::defineLabels() {

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
