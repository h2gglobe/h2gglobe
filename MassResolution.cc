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

MassResolution::MassResolution(){
//  dz_file = TFile::Open("../PhotonAnalysis/data/dz_vs_hpt.root");
  //dz_plot = (TGraph*)dz_file->Get("dz_vs_hpt");
  //dz_file->Close();
}

MassResolution::MassResolution(std::string ecorrmethod){
  energyCorrectionMethod = ecorrmethod;
 // dz_file = TFile::Open(fileName.c_str());
 // dz_plot = (TGraph*)dz_file->Get("dz_vs_hpt");
 // dz_file->Close();
}

void MassResolution::Setup(LoopAll &l, TLorentzVector *in_lead_p4, TLorentzVector *in_sublead_p4, int lead_index, int sublead_index, int diphoton_index, double higgsPt, double higgsM, EnergySmearer::energySmearingParameters eSmearPars, int nR9Categories, int nEtaCategories){

  lead_p4 = in_lead_p4;
  sublead_p4 = in_sublead_p4;
  lead_sc_pos = (TVector3*)(l.sc_xyz->At(l.pho_scind[lead_index]));
  sublead_sc_pos = (TVector3*)(l.sc_xyz->At(l.pho_scind[sublead_index]));
  vertex = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_index]);
  vtx_dxdydz = (TVector3*)l.vtx_std_dxdydz->At(l.dipho_vtxind[diphoton_index]);
  if (energyCorrectionMethod=="DaunceyAndKenzie"){
    lead_Eres = l.pho_residCorrResn[lead_index];
    sublead_Eres = l.pho_residCorrResn[sublead_index];
  } else if (energyCorrectionMethod=="Bendavid"){
    lead_Eres = l.pho_regr_energyerr[lead_index];
    sublead_Eres = l.pho_regr_energyerr[sublead_index];
  } else if (energyCorrectionMethod=="BendavidOTF"){
    lead_Eres = l.pho_regr_energyerr_otf[lead_index];
    sublead_Eres = l.pho_regr_energyerr_otf[sublead_index];
  }

//	lead_r9 = l.pho_r9[lead_index];
// sublead_r9 = l.pho_r9[sublead_index];
	lead_phoCat = l.PhotonCategory(lead_index,nR9Categories,nEtaCategories);
	sublead_phoCat = l.PhotonCategory(sublead_index,nR9Categories,nEtaCategories);
  _eSmearPars = eSmearPars;

  lead_iDet = (bool)l.pho_isEB[lead_index];
  sublead_iDet =(bool) l.pho_isEB[sublead_index];

//  dz = dz_plot->Eval(higgsPt);
  higgsMass = higgsM;
}
  
// return the mass resolution
/*
double MassResolution::massResolution(){
  
  double lead_E = lead_p4->E();
  double sublead_E = sublead_p4->E();
  double alpha = lead_p4->Angle(sublead_p4->Vect());
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}
*/
// return the mass resolution given correct vertex
double MassResolution::massResolutionCorrVtx(){
  
  double lead_E = lead_p4->E();
  double sublead_E = sublead_p4->E();
  double alpha = lead_p4->Angle(sublead_p4->Vect());
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}
// return the mass resolution wrong vertex
double MassResolution::massResolutionWrongVtx(){
  
  double lead_E = lead_p4->E();
  double sublead_E = sublead_p4->E();
  double alpha = lead_p4->Angle(sublead_p4->Vect());
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = angleResolutionWrongVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}

// return energy contribution to mass resolution only
double MassResolution::massResolutionEonly() {

  double lead_E = lead_p4->E();
  double sublead_E = sublead_p4->E();
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();

  return 0.5*higgsMass*TMath::Sqrt((lead_sig*lead_sig)/(lead_E*lead_E)+(sublead_sig*sublead_sig)/(sublead_E*sublead_E));
}

// return angular contribution to mass resolution only
double MassResolution::massResolutionAonly() {

  double alpha = lead_p4->Angle(sublead_p4->Vect());
  double alpha_sig = angleResolutionCorrVtx();

  return 0.5*higgsMass*(alpha_sig*TMath::Sin(alpha)/(1.-TMath::Cos(alpha)));
}

// return energy contribution with no smearing to mass resolution only
double MassResolution::massResolutionEonlyNoSmear() {

  double lead_E = lead_p4->E();
  double sublead_E = sublead_p4->E();
  double lead_sig = leadPhotonResolutionNoSmear();
  double sublead_sig = subleadPhotonResolutionNoSmear();

  return 0.5*higgsMass*TMath::Sqrt((lead_sig*lead_sig)/(lead_E*lead_E)+(sublead_sig*sublead_sig)/(sublead_E*sublead_E));
}

// return the mass resolution given correct vertex with no smearing
double MassResolution::massResolutionCorrVtxNoSmear(){
  
  double lead_E = lead_p4->E();
  double sublead_E = sublead_p4->E();
  double alpha = lead_p4->Angle(sublead_p4->Vect());
  double lead_sig = leadPhotonResolutionNoSmear();
  double sublead_sig = subleadPhotonResolutionNoSmear();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}
// return the mass resolution wrong vertex with smearing
double MassResolution::massResolutionWrongVtxNoSmear(){
  
  double lead_E = lead_p4->E();
  double sublead_E = sublead_p4->E();
  double alpha = lead_p4->Angle(sublead_p4->Vect());
  double lead_sig = leadPhotonResolutionNoSmear();
  double sublead_sig = subleadPhotonResolutionNoSmear();
  double alpha_sig = angleResolutionWrongVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

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
  return getPhotonResolution(lead_p4->E(),lead_Eres,lead_r9, lead_phoCat,lead_sc_pos->Eta(),lead_iDet);
}
// return sublead photon resolution
double MassResolution::subleadPhotonResolution() {
  return getPhotonResolution(sublead_p4->E(),sublead_Eres,sublead_r9,sublead_phoCat,sublead_sc_pos->Eta(),sublead_iDet);
}
// Actually compute resolution given a photon
double MassResolution::getPhotonResolution(double photonEnergy, double photonResolution, double r9, int phoCat, double scEta, bool iDet) {


  // Get the photon-category sigma
  std::string myCategory="";
  if (_eSmearPars.categoryType=="Automagic") 
    {
	    EnergySmearer::energySmearingParameters::phoCatVectorConstIt vit = find(_eSmearPars.photon_categories.begin(), _eSmearPars.photon_categories.end(), 
										    std::make_pair( fabs((float)scEta), (float)r9 ) );
	    if( vit ==  _eSmearPars.photon_categories.end() ) {
		    std::cerr << "Could not find energy scale correction for this photon " << (float)scEta << " " <<  (float)r9 << std::endl;
		    assert( 0 );
	    }
	    myCategory = vit->name;
    } 
  else if (_eSmearPars.categoryType=="2CatR9_EBEE")
    {
      if (iDet)
	myCategory+="EB";
      else
	myCategory+="EE";
      
      if (r9>=0.94)
	myCategory+="HighR9";
      else
	myCategory+="LowR9";
    }
  else if (_eSmearPars.categoryType=="2CatR9_EBEBm4EE")
    {
      if (iDet && fabs(scEta)      < 1.)
	myCategory+="EB";
      else if (iDet && fabs(scEta) > 1.)
	myCategory+="EBm4";
      else
	myCategory+="EE";
      
      if (r9>=0.94)
	myCategory+="HighR9";
      else
	myCategory+="LowR9";
    }
  else if (_eSmearPars.categoryType=="EBEE")
    {
      if (iDet)
	myCategory+="EB";
      else
	myCategory+="EE";
    }
  else
    {
      std::cout << "Unknown categorization. No category name is returned" << std::endl;
    }




	double categoryResolution = _eSmearPars.smearing_sigma[myCategory]*photonEnergy;	
	return TMath::Sqrt(categoryResolution*categoryResolution + photonResolution*photonResolution);
}

//return dz resolution given correct vertex (used 10mm)
double MassResolution::dzResolutionCorrVtx() {
  return 0.1;
}
//return dz resolution given wrong vertex (using sqrt(2)*5.8cm)
double MassResolution::dzResolutionWrongVtx() {
  return TMath::Sqrt(2.)*5.8;
}
  
//return dz resolution from dz wrong and dz right (stored in TGraph as func of higgsPt)
double MassResolution::dzResolution() { 
  return dz;
}

// propagate error on z to error on angle
double MassResolution::propagateDz(double dz){

  double alpha = lead_p4->Angle(sublead_p4->Vect());
  if (alpha!= sublead_p4->Angle(lead_p4->Vect())) std::cout << "Error: Angle between photons not consistent" << std::endl;
  
  double x1 = lead_sc_pos->X()-vertex->X();
  double y1 = lead_sc_pos->Y()-vertex->Y();
  double z1 = lead_sc_pos->Z()-vertex->Z();
 
  double x2 = sublead_sc_pos->X()-vertex->X();
  double y2 = sublead_sc_pos->Y()-vertex->Y();
  double z2 = sublead_sc_pos->Z()-vertex->Z();
 
  double r1 = TMath::Sqrt(x1*x1+y1*y1+z1*z1);
  double r2 = TMath::Sqrt(x2*x2+y2*y2+z2*z2);

  double cos_term = TMath::Cos(lead_p4->Phi()-sublead_p4->Phi());
  double sech1 = SecH(lead_p4->Eta());
  double sech2 = SecH(sublead_p4->Eta());
  double tanh1 = TanH(lead_p4->Eta());
  double tanh2 = TanH(sublead_p4->Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos_term);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos_term);
  double denominator = 1. - tanh1*tanh2-sech1*sech2*cos_term;

  double ResTerm = (-1.*dz/denominator)*(numerator1/r1 + numerator2/r2);

  double angleResolution = ResTerm*(1.-TMath::Cos(alpha))/TMath::Sin(alpha);

  return angleResolution;

}

// utility functions
double MassResolution::SecH(double x){
  return 1.0/TMath::CosH(x);
}

double MassResolution::TanH(double x){
  return TMath::TanH(x);
}

void MassResolution::printInfo(){
  std::cout << "MassResolution Object" << std::endl;
  std::cout << "Lead p4: (" << lead_p4->Px() << "," << lead_p4->Py() << "," << lead_p4->Pz() << "," << lead_p4->E() << std::endl;
  std::cout << "Sublead p4: (" << sublead_p4->Px() << "," << sublead_p4->Py() << "," << sublead_p4->Pz() << "," << sublead_p4->E() << std::endl;
  std::cout << "LEAD SC: (" << lead_sc_pos->X() << "," << lead_sc_pos->Y() << "," << lead_sc_pos->Z() << ")" << std::endl; 
  std::cout << "SUBL SC: (" << sublead_sc_pos->X() << "," << sublead_sc_pos->Y() << "," << sublead_sc_pos->Z() << ")" << std::endl; 
  std::cout << "VTX: (" << vertex->X() << "," << vertex->Y() << "," << vertex->Z() << ")" << std::endl;
  std::cout << "L R9: " << lead_r9 << std::endl;
  std::cout << "SL R9: " << sublead_r9 << std::endl;
  std::cout << "L cat: " << lead_phoCat << std::endl;
  std::cout << "SL cat: " << sublead_phoCat << std::endl;
}

void MassResolution::dumpInfo(std::ostream &o){
  o << "MassResolution Object" << std::endl;
  o << "Lead p4: (" << lead_p4->Px() << "," << lead_p4->Py() << "," << lead_p4->Pz() << "," << lead_p4->E() << std::endl;
  o << "Sublead p4: (" << sublead_p4->Px() << "," << sublead_p4->Py() << "," << sublead_p4->Pz() << "," << sublead_p4->E() << std::endl;
  o << "LEAD SC: (" << lead_sc_pos->X() << "," << lead_sc_pos->Y() << "," << lead_sc_pos->Z() << ")" << std::endl; 
  o << "SUBL SC: (" << sublead_sc_pos->X() << "," << sublead_sc_pos->Y() << "," << sublead_sc_pos->Z() << ")" << std::endl; 
  o << "VTX: (" << vertex->X() << "," << vertex->Y() << "," << vertex->Z() << ")" << std::endl;
  o << "L R9: " << lead_r9 << std::endl;
  o << "SL R9: " << sublead_r9 << std::endl;
  o << "L cat: " << lead_phoCat << std::endl;
  o << "SL cat: " << sublead_phoCat << std::endl;
}
