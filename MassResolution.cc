#include "MassResolution.h"

//----------------------------------------------------------//
// Project:   MassResolution
// Author:    Matt Kenzie (matthew.william.kenzie@cern.ch)
// Modified:  25/08/2011
// Admins:    Matth Kenzie (matthew.william.kenzie@cern.ch)
//---------------------------------------------------------//

/*
See MassResolution.h for instructions
*/

MassResolution::MassResolution(){}

void MassResolution::Setup(LoopAll &l, PhotonReducedInfo *leadInfo, PhotonReducedInfo *subleadInfo,int vtx_index,EnergySmearer::energySmearingParameters eSmearPars, int nR9Categories, int nEtaCategories, double beamspotSigma_in, bool usethisvtx) {
  
  beamspotSigma= beamspotSigma_in;

  leadPhoton= leadInfo;	
  subleadPhoton= subleadInfo;	
  vertex = (TVector3*)l.vtx_std_xyz->At(vtx_index);
  vtx_dxdydz = (TVector3*)l.vtx_std_dxdydz->At(vtx_index);

  //lead_sc_pos    = leadPhoton->caloPosition();
  //sublead_sc_pos = subleadPhoton->caloPosition();

  lead_Eres    = leadPhoton->corrEnergyErr();
  sublead_Eres = subleadPhoton->corrEnergyErr();

  lead_r9      = leadPhoton->r9();
  sublead_r9   = subleadPhoton->r9();

  lead_iDet    = leadPhoton->iDet();
  sublead_iDet = subleadPhoton->iDet();

  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

  higgsMass=(lead_p4+sublead_p4).M();
  _eSmearPars=eSmearPars;

/*  std::cout << "Inside the mmassreso calc -- E1 - " <<leadPhoton->energy() <<std::endl;
  std::cout << "Inside the mmassreso calc -- E2 - " <<subleadPhoton->energy() <<std::endl;
  std::cout << "Inside the mmassreso calc -- r91 - " <<lead_r9 <<std::endl;
  std::cout << "Inside the mmassreso calc -- r92 - " <<sublead_r9 <<std::endl;
  std::cout << "Inside the mmassreso calc -- mass - " <<higgsMass <<std::endl;
  std::cout << "Inside the mmassreso calc -- calopos Eta1 - " <<leadPhoton->caloPosition().Eta() <<std::endl;
  std::cout << "Inside the mmassreso calc -- calopos Eta2 - " <<subleadPhoton->caloPosition().Eta() <<std::endl;
*/


}


void MassResolution::Setup(LoopAll &l, PhotonReducedInfo *leadInfo, PhotonReducedInfo *subleadInfo,int diphoton_index,EnergySmearer::energySmearingParameters eSmearPars, int nR9Categories, int nEtaCategories, double beamspotSigma_in) 
{
  Setup(l, leadInfo, subleadInfo, l.dipho_vtxind[diphoton_index], eSmearPars, nR9Categories, nEtaCategories, beamspotSigma_in, true);

}
// return the mass resolution given correct vertex
double MassResolution::massResolutionCorrVtx(){

  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

  double lead_E = lead_p4.E();
  double sublead_E = sublead_p4.E();
  double alpha = lead_p4.Angle(sublead_p4.Vect());
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}

double MassResolution::massResolutionCorrVtxNoSmear(){

  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

  double lead_E = lead_p4.E();
  double sublead_E = sublead_p4.E();
  double alpha = lead_p4.Angle(sublead_p4.Vect());
  double lead_sig = leadPhotonResolutionNoSmear();
  double sublead_sig = subleadPhotonResolutionNoSmear();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}
// return the mass resolution wrong vertex
double MassResolution::massResolutionWrongVtx(){
  
//  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
//  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

//  double lead_E = lead_p4.E();
//  double sublead_E = sublead_p4.E();
//  double alpha = lead_p4.Angle(sublead_p4.Vect());
//  double lead_sig = leadPhotonResolution();
//  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = higgsMass*0.5*angleResolutionWrongVtx();
  
   double sigmaM = massResolutionEonly();
//  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)));
  return TMath::Sqrt((sigmaM*sigmaM)+(alpha_sig*alpha_sig));

}

double MassResolution::massResolutionWrongVtxNoSmear(){
  double alpha_sig = higgsMass*0.5*angleResolutionWrongVtx();
  
   double sigmaM = massResolutionEonlyNoSmear();
//  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)));
  return TMath::Sqrt((sigmaM*sigmaM)+(alpha_sig*alpha_sig));
}

// return energy contribution to mass resolution only
double MassResolution::massResolutionEonly() {

  double lead_E = leadPhoton->energy();
  double sublead_E = subleadPhoton->energy();
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();

  return 0.5*higgsMass*TMath::Sqrt((lead_sig*lead_sig)/(lead_E*lead_E)+(sublead_sig*sublead_sig)/(sublead_E*sublead_E));
}

double MassResolution::massResolutionEonlyNoSmear(){
  
  double lead_E = leadPhoton->energy();
  double sublead_E = subleadPhoton->energy();
  double lead_sig = leadPhotonResolutionNoSmear();
  double sublead_sig = subleadPhotonResolutionNoSmear();

  return 0.5*higgsMass*TMath::Sqrt((lead_sig*lead_sig)/(lead_E*lead_E)+(sublead_sig*sublead_sig)/(sublead_E*sublead_E));

}

double MassResolution::massResolutionAonly(){
	double aRes = angleResolution();
	return 0.5*higgsMass*aRes;
}


// return angle resolution given the vertex choice is correct
double MassResolution::angleResolutionCorrVtx() {
  return propagateDz(dzResolutionCorrVtx());
}

// return angle resolution given the vertex choice is wrong
double MassResolution::angleResolutionWrongVtx() {
  return propagateDz(dzResolutionWrongVtx());
}

// return angle resolution given a convolution of correct/wrong vertex as func of higgsPt
double MassResolution::angleResolution() {
  return propagateDz(dzResolution());
}
// return lead photon resolution without smearing
double MassResolution::leadPhotonResolutionNoSmear() {
  return lead_Eres;
}
// return sublead photon resolution without smearing
double MassResolution::subleadPhotonResolutionNoSmear() {
  return sublead_Eres;
}
// return lead photon resolution 
double MassResolution::leadPhotonResolution() {
  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  bool sphericalLeadPhoton_=leadPhoton->isSphericalPhoton();
//  std::cout << " MassResolution -- Lead special ? " << sphericalLeadPhoton_ <<std::endl;
  return getPhotonResolution(lead_p4.E(),lead_Eres, *leadPhoton);
}
// return sublead photon resolution
double MassResolution::subleadPhotonResolution() {
  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  bool sphericalSubleadPhoton_=subleadPhoton->isSphericalPhoton();
//  std::cout << " MassResolution -- SubLead special ? " << sphericalSubleadPhoton_ <<std::endl;
  return getPhotonResolution(sublead_p4.E(),sublead_Eres,*subleadPhoton);
}

// Actually compute resolution given a photon
double MassResolution::getPhotonResolution(double photonEnergy, double photonResolution, const PhotonReducedInfo &info) {
  
  // Get the photon-category sigma
  std::string myCategory = EnergySmearer::photonCategory(_eSmearPars, info);

  double categoryResolution = EnergySmearer::getSmearingSigma(_eSmearPars, myCategory, photonEnergy, 0.)*photonEnergy;
  return TMath::Sqrt(categoryResolution*categoryResolution + photonResolution*photonResolution);

}

//return dz resolution given correct vertex (used 10mm)
double MassResolution::dzResolutionCorrVtx() {
  return 0.1;
}
//return dz resolution given wrong vertex (using sqrt(2)*5.8cm)
double MassResolution::dzResolutionWrongVtx() {
  return TMath::Sqrt(2.)*beamspotSigma;
}
  
//return dz resolution from dz wrong and dz right (stored in TGraph as func of higgsPt)
double MassResolution::dzResolution() { 
  return dz;
}

// propagate error on z to error on angle
double MassResolution::propagateDz(double dz){

//  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
//  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

//  double alpha = //lead_p4.Angle(sublead_p4.Vect());
//  if (alpha!= sublead_p4.Angle(lead_p4.Vect())) std::cout << "Error: Angle between photons not consistent" << std::endl;

  TVector3 LeadPosition = (leadPhoton->caloPosition()) - *vertex;
  TVector3 SubLeadPosition = (subleadPhoton->caloPosition()) - *vertex;

/*
  double x1 = leadPhoton->caloPosition().X();
  double y1 = leadPhoton->caloPosition().Y();
  double z1 = leadPhoton->caloPosition().Z();

  double x2 = subleadPhoton->caloPosition().X();
  double y2 = subleadPhoton->caloPosition().Y();
  double z2 = subleadPhoton->caloPosition().Z();
*/
  
  double x1 = leadPhoton->caloPosition().X()-vertex->X();
  double y1 = leadPhoton->caloPosition().Y()-vertex->Y();
  double z1 = leadPhoton->caloPosition().Z()-vertex->Z();
 
  double x2 = subleadPhoton->caloPosition().X()-vertex->X();
  double y2 = subleadPhoton->caloPosition().Y()-vertex->Y();
  double z2 = subleadPhoton->caloPosition().Z()-vertex->Z();

 
  double r1 = TMath::Sqrt(x1*x1+y1*y1+z1*z1);
  double r2 = TMath::Sqrt(x2*x2+y2*y2+z2*z2);

  double cos_term = TMath::Cos(LeadPosition.Phi()-SubLeadPosition.Phi());
  double sech1 = SecH(LeadPosition.Eta());
  double sech2 = SecH(SubLeadPosition.Eta());
  double tanh1 = TanH(LeadPosition.Eta());
  double tanh2 = TanH(SubLeadPosition.Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos_term);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos_term);
  double denominator = 1. - tanh1*tanh2 - sech1*sech2*cos_term;

  double ResTerm = (-1.*dz/denominator)*(numerator1/r1 + numerator2/r2);

  //double angleResolution = ResTerm*(1.-TMath::Cos(alpha))/TMath::Sin(alpha);
  double angleResolution = ResTerm;

  return angleResolution;

}

// utility functions
double MassResolution::SecH(double x){
  return 1.0/TMath::CosH(x);
}

double MassResolution::TanH(double x){
  return TMath::TanH(x);
}

