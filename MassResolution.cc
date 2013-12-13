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

  lead_Eres    = leadPhoton->corrEnergyErr() / leadPhoton->corrEnergy();
  sublead_Eres = subleadPhoton->corrEnergyErr() / subleadPhoton->corrEnergy();

  lead_r9      = leadPhoton->r9();
  sublead_r9   = subleadPhoton->r9();

  lead_iDet    = leadPhoton->iDet();
  sublead_iDet = subleadPhoton->iDet();

  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

  higgsMass=(lead_p4+sublead_p4).M();
  _eSmearPars=eSmearPars;

}


void MassResolution::Setup(LoopAll &l, PhotonReducedInfo *leadInfo, PhotonReducedInfo *subleadInfo,int diphoton_index,EnergySmearer::energySmearingParameters eSmearPars, int nR9Categories, int nEtaCategories, double beamspotSigma_in) 
{
  Setup(l, leadInfo, subleadInfo, l.dipho_vtxind[diphoton_index], eSmearPars, nR9Categories, nEtaCategories, beamspotSigma_in, true);

}
// return the mass resolution given correct vertex
double MassResolution::relMassResolutionCorrVtx(){

  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

  double alpha = lead_p4.Angle(sublead_p4.Vect());
  double lead_sig = leadRelPhotonResolution();
  double sublead_sig = subleadRelPhotonResolution();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*TMath::Sqrt((lead_sig*lead_sig)+(sublead_sig*sublead_sig)
			 +((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}

double MassResolution::relMassResolutionCorrVtxNoSmear(){

  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

  double alpha = lead_p4.Angle(sublead_p4.Vect());
  double lead_sig = leadRelPhotonResolutionNoSmear();
  double sublead_sig = subleadRelPhotonResolutionNoSmear();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*TMath::Sqrt((lead_sig*lead_sig)+(sublead_sig*sublead_sig)
			 +((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}
// return the mass resolution wrong vertex
double MassResolution::relMassResolutionWrongVtx(){
  
  double alpha_sig = 0.5*angleResolutionWrongVtx();
  
  double sigmaM = relMassResolutionEonly();
  return TMath::Sqrt((sigmaM*sigmaM)+(alpha_sig*alpha_sig));
}

double MassResolution::relMassResolutionWrongVtxNoSmear(){
  double alpha_sig = 0.5*angleResolutionWrongVtx();
  
  double sigmaM = relMassResolutionEonlyNoSmear();
  return TMath::Sqrt((sigmaM*sigmaM)+(alpha_sig*alpha_sig));
}

// return energy contribution to mass resolution only
double MassResolution::relMassResolutionEonly() {

  double lead_sig = leadRelPhotonResolution();
  double sublead_sig = subleadRelPhotonResolution();

  return 0.5*TMath::Sqrt((lead_sig*lead_sig)+(sublead_sig*sublead_sig));
}

double MassResolution::relMassResolutionEonlyNoSmear(){
  
  double lead_sig = leadRelPhotonResolutionNoSmear();
  double sublead_sig = subleadRelPhotonResolutionNoSmear();

  return 0.5*TMath::Sqrt((lead_sig*lead_sig)+(sublead_sig*sublead_sig));

}

double MassResolution::relMassResolutionAonly(){
  double aRes = angleResolution();
  return 0.5*aRes;
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
double MassResolution::leadRelPhotonResolutionNoSmear() {
  return lead_Eres;
}
// return sublead photon resolution without smearing
double MassResolution::subleadRelPhotonResolutionNoSmear() {
  return sublead_Eres;
}
// return lead photon resolution 
double MassResolution::leadRelPhotonResolution() {
  bool sphericalLeadPhoton_=leadPhoton->isSphericalPhoton();
  return getRelPhotonResolution(lead_Eres, *leadPhoton);
}
// return sublead photon resolution
double MassResolution::subleadRelPhotonResolution() {
  bool sphericalSubleadPhoton_=subleadPhoton->isSphericalPhoton();
  return getRelPhotonResolution(sublead_Eres,*subleadPhoton);
}

// Actually compute resolution given a photon
double MassResolution::getRelPhotonResolution(double photonResolution, const PhotonReducedInfo &info) {
  
  // Get the photon-category sigma
  std::string myCategory = EnergySmearer::photonCategory(_eSmearPars, info);

  double categoryResolution = EnergySmearer::getSmearingSigma(_eSmearPars, myCategory, info.corrEnergy(), 
							      info.caloPosition().Eta(), 0.);
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

