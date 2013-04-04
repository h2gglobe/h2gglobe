#include "MyAnalysis.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
void MyAnalysis::setBranchAddresses(TTree* tree, std::map<std::string, std::string>& values) {
  
  tree->SetBranchAddress("itype",       &itype);       
  tree->SetBranchAddress("eta1",        &eta1);	       
  tree->SetBranchAddress("eta2",        &eta2);	       
  tree->SetBranchAddress("et1",         &et1);
  tree->SetBranchAddress("et2",         &et2);
  tree->SetBranchAddress("muid1",       &muid1);
  tree->SetBranchAddress("muid2",       &muid2);
  tree->SetBranchAddress("cat",         &cat);
  tree->SetBranchAddress("mass",        &mass);
  tree->SetBranchAddress("full_weight", &full_weight);

  defineLabels();
}

//----------------------------------------------------------------------

void MyAnalysis::analyze(RooContainer* rooContainer) {
  
  if (muid1 == 3 && muid2==3)
    FillRooContainer(rooContainer); 
}

//----------------------------------------------------------------------

extern "C" GenericAnalysis *makeInstance() {
  return new MyAnalysis();
}

//----------------------------------------------------------------------

void MyAnalysis::FillRooContainer(RooContainer* rooContainer) {

  if (itype == 0 ) {
    rooContainer->InputDataPoint("data_mass", cat, mass);
  } else if (itype > 0 ) {
    rooContainer->InputDataPoint("bkg_mass",  cat, mass, full_weight);
  } else if (itype < 0) {
    rooContainer->InputDataPoint("sig_" + GetSignalLabel(itype), cat, mass, full_weight);
  }
}

//----------------------------------------------------------------------

void MyAnalysis::FillRooContainerSyst(RooContainer* rooContainer) {

  if (itype < 0) {

    std::vector<double> mass_errors, weights;
    std::vector<int> categories;
    
    mass_errors.push_back(mass);
    weights.push_back(full_weight);
    categories.push_back(cat);

    rooContainer->InputSystematicSet("sig_" + GetSignalLabel(itype), *name1, categories, mass_errors, weights);
  }
}
    
//----------------------------------------------------------------------

std::string MyAnalysis::GetSignalLabel(int id) {
 
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

void MyAnalysis::defineLabels() {

  // mu channel labels
  signalLabels[-113]="ggh_mass_m125";
  signalLabels[-213]="vbf_mass_m125";
  signalLabels[-513]="rsg_mass_m125";

  signalLabels[-114]="ggh_mass_m120";
  signalLabels[-214]="vbf_mass_m120";
  signalLabels[-514]="rsg_mass_m120";

  signalLabels[-115]="ggh_mass_m130";
  signalLabels[-215]="vbf_mass_m130";
  signalLabels[-515]="rsg_mass_m130";

  // e channel labels
  signalLabels[-123]="ggh_mass_m125";
  signalLabels[-223]="vbf_mass_m125";
  signalLabels[-523]="rsg_mass_m125";

  signalLabels[-124]="ggh_mass_m120";
  signalLabels[-224]="vbf_mass_m120";
  signalLabels[-524]="rsg_mass_m120";

  signalLabels[-125]="ggh_mass_m130";
  signalLabels[-225]="vbf_mass_m130";
  signalLabels[-525]="rsg_mass_m130";


}
