#include "../interface/JetResponseChange.h"

JetResponseChange::JetResponseChange(std::string filename) { 
  TFile *file = new TFile(filename.c_str());
  std::cout << "Opening Jet Response file : " << filename.c_str() << std::endl;
  hResponse       = (TH1F*) file->FindObjectAny("Eta");
  hResponse ->SetDirectory(0);
}


JetResponseChange::~JetResponseChange(){
  delete hResponse;
}


float JetResponseChange::getResponse(float eta, float lumi){
  float slope =0;
  if (hResponse)
    slope = hResponse->GetBinContent(hResponse->GetXaxis()->FindBin(fabs(eta)));
  return (1-slope*lumi);
}
