#include "Sorters.h"


void LoopAll::SetAllMVA() {
  tmvaReaderID_UCSD = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_UCSD->AddVariable("sieie",      &tmva_id_ucsd_sieie);
  tmvaReaderID_UCSD->AddVariable("goodpf_iso", &tmva_id_ucsd_goodpf_iso);
  tmvaReaderID_UCSD->AddVariable("badpf_iso",  &tmva_id_ucsd_badpf_iso);
  tmvaReaderID_UCSD->AddVariable("drtotk",     &tmva_id_ucsd_drtotk);
  tmvaReaderID_UCSD->AddVariable("hoe",        &tmva_id_ucsd_hoe);
  tmvaReaderID_UCSD->AddVariable("tkisopf",    &tmva_id_ucsd_tkisopf);
  tmvaReaderID_UCSD->AddVariable("r9",         &tmva_id_ucsd_r9);
  tmvaReaderID_UCSD->AddVariable("ptom",       &tmva_id_ucsd_ptom);
  tmvaReaderID_UCSD->AddVariable("eta",        &tmva_id_ucsd_eta);
  tmvaReaderID_UCSD->AddSpectator("isLeading", &tmva_id_ucsd_isLeading);
//  tmvaReaderID_UCSD->BookMVA("Gradient", "ID_UCSD.weights.xml");
  
  tmvaReader_dipho_UCSD = new TMVA::Reader("!Color:Silent"); 
  tmvaReader_dipho_UCSD->AddVariable("subleadptomass",   &tmva_dipho_UCSD_subleadptomass);
  tmvaReader_dipho_UCSD->AddVariable("diphoptom",        &tmva_dipho_UCSD_diphoptom);
  tmvaReader_dipho_UCSD->AddVariable("sumptom",          &tmva_dipho_UCSD_sumptom);
  tmvaReader_dipho_UCSD->AddVariable("subleadmva",       &tmva_dipho_UCSD_subleadmva);
  tmvaReader_dipho_UCSD->AddVariable("leadmva",          &tmva_dipho_UCSD_leadmva);
  tmvaReader_dipho_UCSD->AddVariable("leadeta",          &tmva_dipho_UCSD_leadeta);
  tmvaReader_dipho_UCSD->AddVariable("subleadeta",       &tmva_dipho_UCSD_subleadeta);
  tmvaReader_dipho_UCSD->AddVariable("leadr9",           &tmva_dipho_UCSD_leadr9);
  tmvaReader_dipho_UCSD->AddVariable("subleadr9",        &tmva_dipho_UCSD_subleadr9);
  tmvaReader_dipho_UCSD->AddVariable("dmom",             &tmva_dipho_UCSD_dmom);
  tmvaReader_dipho_UCSD->AddSpectator("diphocat2r92eta", &tmva_dipho_UCSD_diphocat2r92eta);
//  tmvaReader_dipho_UCSD->BookMVA("Gradient", "diphoton_UCSD.weights.xml");

  tmvaReaderID_MIT_Barrel = new TMVA::Reader("!Color:Silent"); 
  tmvaReaderID_MIT_Barrel->AddVariable("HoE",         &tmva_id_mit_hoe);
  tmvaReaderID_MIT_Barrel->AddVariable("covIEtaIEta", &tmva_id_mit_sieie);
  tmvaReaderID_MIT_Barrel->AddVariable("tIso1abs",    &tmva_id_mit_tiso1);
  tmvaReaderID_MIT_Barrel->AddVariable("tIso3abs",    &tmva_id_mit_tiso3);
  tmvaReaderID_MIT_Barrel->AddVariable("tIso2abs",    &tmva_id_mit_tiso2);
  tmvaReaderID_MIT_Barrel->AddVariable("R9",          &tmva_id_mit_r9);
  tmvaReaderID_MIT_Barrel->AddVariable("absIsoEcal",  &tmva_id_mit_ecal);
  tmvaReaderID_MIT_Barrel->AddVariable("absIsoHcal",  &tmva_id_mit_hcal);
//  tmvaReaderID_MIT_Barrel->AddVariable("RelE5x5",     &tmva_id_mit_e5x5);
//  tmvaReaderID_MIT_Barrel->AddVariable("EtaWidth",    &tmva_id_mit_etawidth);
//  tmvaReaderID_MIT_Barrel->AddVariable("PhiWidth",    &tmva_id_mit_phiwidth);
//  tmvaReaderID_MIT_Barrel->AddVariable("CoviEtaiPhi", &tmva_id_mit_sieip);
//  tmvaReaderID_MIT_Barrel->AddVariable("CoviPhiiPhi", &tmva_id_mit_sipip);
  tmvaReaderID_MIT_Barrel->AddVariable("NVertexes",   &tmva_id_mit_nvtx);
  tmvaReaderID_MIT_Barrel->AddVariable("ScEta",	      &tmva_id_mit_sceta);
  tmvaReaderID_MIT_Barrel->AddVariable("EtaWidth",    &tmva_id_mit_etawidth);
  tmvaReaderID_MIT_Barrel->AddVariable("PhiWidth",    &tmva_id_mit_phiwidth);
//  tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost", "TMVAClassificationPhotonID_Barrel_PassPreSel_Variable_6_BDTnCuts2000_BDT.weights.xml");

  tmvaReaderID_MIT_Endcap = new TMVA::Reader("!Color:Silent"); 
  tmvaReaderID_MIT_Endcap->AddVariable("HoE",         &tmva_id_mit_hoe);
  tmvaReaderID_MIT_Endcap->AddVariable("covIEtaIEta", &tmva_id_mit_sieie);
  tmvaReaderID_MIT_Endcap->AddVariable("tIso1abs",    &tmva_id_mit_tiso1);
  tmvaReaderID_MIT_Endcap->AddVariable("tIso3abs",    &tmva_id_mit_tiso3);
  tmvaReaderID_MIT_Endcap->AddVariable("tIso2abs",    &tmva_id_mit_tiso2);
  tmvaReaderID_MIT_Endcap->AddVariable("R9",          &tmva_id_mit_r9);
  tmvaReaderID_MIT_Endcap->AddVariable("absIsoEcal",  &tmva_id_mit_ecal);
  tmvaReaderID_MIT_Endcap->AddVariable("absIsoHcal",  &tmva_id_mit_hcal);
//  tmvaReaderID_MIT_Endcap->AddVariable("RelE5x5",     &tmva_id_mit_e5x5);
//  tmvaReaderID_MIT_Endcap->AddVariable("EtaWidth",    &tmva_id_mit_etawidth);
//  tmvaReaderID_MIT_Endcap->AddVariable("PhiWidth",    &tmva_id_mit_phiwidth);
//  tmvaReaderID_MIT_Endcap->AddVariable("CoviEtaiPhi", &tmva_id_mit_sieip);
//  tmvaReaderID_MIT_Endcap->AddVariable("CoviPhiiPhi", &tmva_id_mit_sipip);
  tmvaReaderID_MIT_Endcap->AddVariable("NVertexes",   &tmva_id_mit_nvtx);
  tmvaReaderID_MIT_Endcap->AddVariable("ScEta",	      &tmva_id_mit_sceta);
  tmvaReaderID_MIT_Endcap->AddVariable("EtaWidth",    &tmva_id_mit_etawidth);
  tmvaReaderID_MIT_Endcap->AddVariable("PhiWidth",    &tmva_id_mit_phiwidth);
/*
  tmvaReaderID_MIT_Endcap = new TMVA::Reader("!Color:Silent"); 
  tmvaReaderID_MIT_Endcap->AddVariable("HoE",                &tmva_id_mit_hoe);
  tmvaReaderID_MIT_Endcap->AddVariable("covIEtaIEta",        &tmva_id_mit_sieie);
  tmvaReaderID_MIT_Endcap->AddVariable("tIso1abs",           &tmva_id_mit_tiso1);
  tmvaReaderID_MIT_Endcap->AddVariable("tIso3abs",           &tmva_id_mit_tiso3);
  tmvaReaderID_MIT_Endcap->AddVariable("tIso2abs",           &tmva_id_mit_tiso2);
  tmvaReaderID_MIT_Endcap->AddVariable("R9",                 &tmva_id_mit_r9);
  tmvaReaderID_MIT_Endcap->AddVariable("absIsoEcal",         &tmva_id_mit_ecal);
  tmvaReaderID_MIT_Endcap->AddVariable("absIsoHcal",         &tmva_id_mit_hcal);
  tmvaReaderID_MIT_Endcap->AddVariable("RelE5x5",            &tmva_id_mit_e5x5);
  tmvaReaderID_MIT_Endcap->AddVariable("EtaWidth",           &tmva_id_mit_etawidth);
  tmvaReaderID_MIT_Endcap->AddVariable("PhiWidth",           &tmva_id_mit_phiwidth);
  tmvaReaderID_MIT_Endcap->AddVariable("CoviEtaiPhi",        &tmva_id_mit_sieip);
  tmvaReaderID_MIT_Endcap->AddVariable("CoviPhiiPhi",        &tmva_id_mit_sipip);
  tmvaReaderID_MIT_Endcap->AddVariable("NVertexes",          &tmva_id_mit_nvtx);
  tmvaReaderID_MIT_Endcap->AddVariable("RelPreshowerEnergy", &tmva_id_mit_preshower);
*/
//  tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost", "TMVAClassificationPhotonID_Endcap_PassPreSel_Variable_6_BDTnCuts2000_BDT.weigh  ts.xml");

  tmvaReader_dipho_MIT = new TMVA::Reader("!Color:Silent"); 
  tmvaReader_dipho_MIT->AddVariable("masserrsmeared/mass",         &tmva_dipho_MIT_dmom);
  tmvaReader_dipho_MIT->AddVariable("masserrsmearedwrongvtx/mass", &tmva_dipho_MIT_dmom_wrong_vtx);
  tmvaReader_dipho_MIT->AddVariable("vtxprob",                     &tmva_dipho_MIT_vtxprob);
  tmvaReader_dipho_MIT->AddVariable("ph1.pt/mass",                 &tmva_dipho_MIT_ptom1);
  tmvaReader_dipho_MIT->AddVariable("ph2.pt/mass",                 &tmva_dipho_MIT_ptom2);
  tmvaReader_dipho_MIT->AddVariable("ph1.eta",                     &tmva_dipho_MIT_eta1);
  tmvaReader_dipho_MIT->AddVariable("ph2.eta",                     &tmva_dipho_MIT_eta2);
  tmvaReader_dipho_MIT->AddVariable("TMath::Cos(ph1.phi-ph2.phi)", &tmva_dipho_MIT_dphi);
  tmvaReader_dipho_MIT->AddVariable("ph1.idmva",                   &tmva_dipho_MIT_ph1mva);
  tmvaReader_dipho_MIT->AddVariable("ph2.idmva",                   &tmva_dipho_MIT_ph2mva);
//  tmvaReader_dipho_MIT->BookMVA("Gradient", "HggBambu_SM_Dec9_BDTG.weights.xml");
}

Float_t LoopAll::photonIDMVA(Int_t iPhoton, Int_t vtx, TLorentzVector &p4, const char* type)  {
  
  Float_t mva = 999.;
 
  float photonEt = p4.Et();
  if (type == "UCSD") {
    Int_t cat = PhotonCategory(iPhoton);

    Float_t isomax=-99;   
    Int_t badind=0;
    for(int iv=0; iv<vtx_std_n; iv++) {
      if((*pho_pfiso_mycharged04)[iPhoton][iv]>isomax) {
	badind=iv; 
	isomax=(*pho_pfiso_mycharged04)[iPhoton][iv]; 
      }
	    }
    
    float rhofacpf[6]    = {0.075, 0.082, 0.143, 0.050, 0.091, 0.106};
    float rhofacbadpf[6] = {0.141, 0.149, 0.208, 0.135, 0.162, 0.165};
    float rhofac         = rhofacpf[cat];
    float rhofacbad      = rhofacbadpf[cat];
    
    Float_t tmva_id_ucsd_pt = photonEt;
    tmva_id_ucsd_badpf_iso = ((*pho_pfiso_mycharged04)[iPhoton][badind]+pho_pfiso_myphoton04[iPhoton]-rho*rhofacbad)*50/tmva_id_ucsd_pt;
    tmva_id_ucsd_goodpf_iso = ((*pho_pfiso_mycharged03)[iPhoton][vtx]+pho_pfiso_myphoton03[iPhoton]-rho*rhofac)*50/tmva_id_ucsd_pt;
    tmva_id_ucsd_tkisopf = (*pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_id_ucsd_pt;
    
    tmva_id_ucsd_sieie = pho_sieie[iPhoton];
    tmva_id_ucsd_drtotk = pho_drtotk_25_99[iPhoton];
    tmva_id_ucsd_hoe = pho_hoe[iPhoton];
    tmva_id_ucsd_r9 = pho_r9[iPhoton];
    tmva_id_ucsd_eta = fabs(p4.Eta());
    tmva_id_ucsd_isLeading = -1.; // not used just a spectator in the original definition
    
    mva = tmvaReaderID_UCSD->EvaluateMVA("Gradient");
  } else {
    tmva_id_mit_hoe = pho_hoe[iPhoton];
    tmva_id_mit_sieie = pho_sieie[iPhoton];
   
    float rhofacbad=0.52, rhofac=0.17;
    
    Float_t raw = sc_raw[pho_scind[iPhoton]];
//    TLorentzVector p4 = get_pho_p4(iPhoton, vtx);
//    float pho_tkiso_goodvtx = SumTrackPtInCone(&p4, vtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
//    float pho_tkiso_badvtx = (WorstSumTrackPtInCone(iPhoton, 0, 0.40, 0.02, 0.0, 1.0, 0.1)).second;

    float pho_tkiso_goodvtx = (*pho_tkiso_recvtx_030_002_0000_10_01)[iPhoton][vtx];
    float pho_tkiso_badvtx = pho_tkiso_badvtx_040_002_0000_10_01[iPhoton];
    tmva_id_mit_tiso1    = (pho_tkiso_goodvtx + pho_ecalsumetconedr03[iPhoton] + pho_hcalsumetconedr04[iPhoton] - rho*rhofac);
    tmva_id_mit_tiso2    = (pho_tkiso_badvtx + pho_ecalsumetconedr04[iPhoton] + pho_hcalsumetconedr04[iPhoton] - rho*rhofacbad);
    tmva_id_mit_tiso3    = pho_tkiso_goodvtx;
    tmva_id_mit_r9       = pho_r9[iPhoton];
    tmva_id_mit_ecal     = pho_ecalsumetconedr03[iPhoton]-rho*rhofac;
    tmva_id_mit_hcal     = pho_hcalsumetconedr04[iPhoton]-rho*rhofac;
    tmva_id_mit_e5x5     = pho_e5x5[iPhoton]/raw;
    tmva_id_mit_etawidth = sc_seta[pho_scind[iPhoton]];
    tmva_id_mit_phiwidth = sc_sphi[pho_scind[iPhoton]];
    tmva_id_mit_sieip    = pho_sieip[iPhoton];
    tmva_id_mit_sipip    = TMath::Sqrt(pho_sipip[iPhoton]);
    tmva_id_mit_nvtx      = vtx_std_n;
    tmva_id_mit_preshower = sc_pre[pho_scind[iPhoton]]/raw;
    tmva_id_mit_sceta	  = ((TVector3*)sc_xyz->At(pho_scind[iPhoton]))->Eta();

    // Print all of the variables
/*
    cout << "PHOTON VARIABLES" <<endl;
    cout << tmva_id_mit_tiso1     << endl;
    cout << tmva_id_mit_tiso2     << endl;
    cout << pho_tkiso_goodvtx << ", " << pho_tkiso_badvtx << endl;
    cout << tmva_id_mit_tiso3     << endl;
    cout << tmva_id_mit_r9        << endl;
    cout << tmva_id_mit_ecal      << endl;
    cout << tmva_id_mit_hcal      << endl;
    cout << tmva_id_mit_e5x5      << endl;
    cout << tmva_id_mit_etawidth  << endl;
    cout << tmva_id_mit_phiwidth  << endl;
    cout << tmva_id_mit_sieip     << endl;
    cout << tmva_id_mit_sipip     << endl;
    cout << tmva_id_mit_nvtx      << endl;
    cout << tmva_id_mit_preshower << endl;
    cout << "--------------------" <<endl;
*/

  if (pho_isEB[iPhoton]) 
    mva = tmvaReaderID_MIT_Barrel->EvaluateMVA("AdaBoost");
  else
    mva = tmvaReaderID_MIT_Endcap->EvaluateMVA("AdaBoost");
  }
  return mva;
}


Float_t LoopAll::diphotonMVA(Int_t leadingPho, Int_t subleadingPho, Int_t vtx, float vtxProb, TLorentzVector &leadP4, TLorentzVector &subleadP4, float sigmaMrv, float sigmaMwv, float sigmaMeonly, const char* type, float photonID_1,float photonID_2) {

  // Ok need to re-write the diphoton-mva part since the systematics won't work unless we can change the Et of the photons
  // all we have to do is to pass in the ->Et of the two photons also rather than take them from the four-vector branches
  
  Float_t mva = 99.;
  TLorentzVector Higgs = leadP4+subleadP4;
  float leadPt    = leadP4.Pt();
  float subleadPt = subleadP4.Pt();
  float mass 	  = Higgs.M();
  float diphopt   = Higgs.Pt();

  if (type == "UCSD") {
    tmva_dipho_UCSD_leadr9 = pho_r9[leadingPho];
    tmva_dipho_UCSD_subleadr9 = pho_r9[subleadingPho];
    tmva_dipho_UCSD_leadeta = fabs(leadP4.Eta());
    tmva_dipho_UCSD_subleadeta = fabs(subleadP4.Eta());
  //  tmva_dipho_UCSD_leadptomass = leadPt/mass;  
    tmva_dipho_UCSD_subleadptomass = subleadPt/mass;  
    tmva_dipho_UCSD_diphoptom = diphopt/mass;
    tmva_dipho_UCSD_sumptom = (leadPt+subleadPt)/mass;
    tmva_dipho_UCSD_subleadmva = photonIDMVA(subleadingPho, vtx,leadP4, "UCSD");
    tmva_dipho_UCSD_leadmva = photonIDMVA(leadingPho, vtx,subleadP4, "UCSD");
   // tmva_dipho_UCSD_dmom = sigmaMrv/mass;
    tmva_dipho_UCSD_dmom = sigmaMeonly/mass;
  
    mva = tmvaReader_dipho_UCSD->EvaluateMVA("Gradient");
  } else {
    tmva_dipho_MIT_dmom = sigmaMrv/mass;
    tmva_dipho_MIT_dmom_wrong_vtx = sigmaMwv/mass;
    tmva_dipho_MIT_vtxprob = vtxProb;
    tmva_dipho_MIT_ptom1 = leadPt/mass;
    tmva_dipho_MIT_ptom2 = subleadPt/mass;
//    tmva_dipho_MIT_eta1 = fabs(leadP4.Eta());
//    tmva_dipho_MIT_eta2 =  fabs(subleadP4.Eta());
    tmva_dipho_MIT_eta1 = leadP4.Eta();
    tmva_dipho_MIT_eta2 =  subleadP4.Eta();
//    tmva_dipho_MIT_dphi = TMath::Cos(((TLorentzVector*)pho_p4->At(leadingPho))->Phi() - ((TLorentzVector*)pho_p4->At(subleadingPho))->Phi());
    tmva_dipho_MIT_dphi = TMath::Cos(leadP4.Phi() - subleadP4.Phi());
      
    if (photonID_1 < -1. && photonID_2 < -1.){
    tmva_dipho_MIT_ph1mva = photonIDMVA(leadingPho,vtx, leadP4, "MIT");
    tmva_dipho_MIT_ph2mva = photonIDMVA(subleadingPho,vtx, subleadP4, "MIT");

    } else {
    tmva_dipho_MIT_ph1mva = photonID_1;//photonIDMVA(leadingPho,vtx, leadP4, "MIT");
    tmva_dipho_MIT_ph2mva = photonID_2;//photonIDMVA(subleadingPho,vtx, subleadP4, "MIT");

    }
    mva = tmvaReader_dipho_MIT->EvaluateMVA("Gradient");
  }
  
  return mva;
}

float LoopAll::getDmOverDz(Int_t pho1, Int_t pho2, Float_t* smeared) {

  TVector3* vtx_plus = new TVector3(0, 0, 0.5);
  TVector3* vtx_minus = new TVector3(0, 0, -0.5);

  TLorentzVector lead_plus = get_pho_p4(pho1, vtx_plus, smeared);
  TLorentzVector sublead_plus = get_pho_p4(pho2, vtx_minus, smeared);
  TLorentzVector diphoton_plus = lead_plus + sublead_plus;

  TLorentzVector lead_minus = get_pho_p4(pho1, vtx_minus, smeared);
  TLorentzVector sublead_minus = get_pho_p4(pho2, vtx_minus, smeared);
  TLorentzVector diphoton_minus = lead_minus + sublead_minus;
  
  float result = diphoton_plus.M() - diphoton_minus.M();

  return result;
}

Float_t LoopAll::deltaMassVtx(Int_t pho1, Int_t pho2, Float_t dz) {

  TVector3* pos1 = (TVector3*)sc_xyz->At(pho_scind[pho1]);
  TVector3* pos2 = (TVector3*)sc_xyz->At(pho_scind[pho2]);

  Float_t r1 = pos1->Mag();
  Float_t r2 = pos2->Mag();
  
  Float_t sech1 = sin(pos1->Theta());
  Float_t tanh1 = cos(pos1->Theta());
  Float_t sech2 = sin(pos2->Theta());
  Float_t tanh2 = cos(pos2->Theta());
  Float_t cos12 = cos(pos1->Phi() - pos2->Phi());
  
  Float_t rad1 = sech1*(sech1*tanh2-tanh1*sech2*cos12)/(1-tanh1*tanh2-sech1*sech2*cos12);
  Float_t rad2 = sech2*(sech2*tanh1-tanh2*sech1*cos12)/(1-tanh2*tanh1-sech2*sech1*cos12);

  return dz * 0.5*fabs(rad1/r1+rad2/r2);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::GlobeCtIsol(int mode, TLorentzVector* p4, float ptCut, float drCutMin, float drCutMax, Int_t & nIsol, Float_t & ptIsol, Float_t & angle1, Float_t & angle2, Float_t & angle3) {
  nIsol=0;
  ptIsol=0.;
  angle1=10.;
  angle2=10.;
  angle3=10.;

  //must put a track selection

  for (int i=0; i<ct_n; i++) {
    TLorentzVector * tempp4= (TLorentzVector *) ct_p4->At(i);
    if(tempp4->Et()<ptCut) continue;
    double dr=p4->DeltaR(*tempp4);
    if(dr<drCutMin) continue;
    if(dr<angle1) {
      angle3=angle2;
      angle2=angle1;
      angle1=dr;
    }
    else if (dr<angle2) {
      angle3=angle2;
      angle2=dr;
    }
    else if (dr<angle3) {
      angle3=dr;
    }
    if(dr>drCutMax) continue;
    nIsol++;
    ptIsol+=tempp4->Et();
  }
}



// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::GlobeMatchIsl(TLorentzVector* p4, Float_t & deltaR) {
  deltaR=10.;
  int imatch=-1;

  for (int i=0; i<sc_islbar_n; i++) {
    TLorentzVector * tempp4= (TLorentzVector *) sc_islbar_p4->At(i);
    double dr=p4->DeltaR(*tempp4);
    if(dr<deltaR) {
      deltaR=dr;
      imatch=i;
    }
  }
  if(imatch==-1) {
    cout<<"ERROR GlobeMatchIsl found no match!!!! "<<endl;
  }
  else if(deltaR>0.3) {
    cout<<"Strange, GlobeMatchIsl deltaR="<<deltaR<<" etapho "<<p4->Eta()<<endl;
  }
  return imatch;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
#include "eIDCuts.h"
std::pair<bool, bool> LoopAll::ElectronId(int index, eIDLevel type) { 

  std::pair<bool, bool> isoIDResult(true, true);

  //TLorentzVector* p4 = (TLorentzVector*)sc_p4->At(el_std_scind[index]);
  TLorentzVector* p4 = (TLorentzVector*)el_std_sc->At(index);
  float eta = fabs(p4->Eta());
  int eb = 0, bin = 0;
  float see = 0;

  if (p4->Et() < 20.) 
    bin = 2;
  else if (p4->Et() > 30.)
    bin = 0;
  else
    bin =1;

  //#ifndef CMSSW3
  //if (eta < 1.479) {
  //see = sc_sieie[el_std_scind[index]];
  //eb = 0;
  //} else {
  //eb = 1; 
  //see = bc_sieie[sc_bcseedind[el_std_scind[index]]];
  //}
  //#else
  if (eta < 1.479) {
    see = el_std_sieiesc[index];
    eb = 0;
  } else {
    eb = 1; 
    see = el_std_sieie[index];
  }
  //#endif

  float eseedopincor = el_std_eseedopin[index] + el_std_fbrem[index];

  if(el_std_fbrem[index]<0) 
    eseedopincor = el_std_eseedopin[index]; 

  //#ifndef CMSSW3
  //float sip = sipCalculator(index);
  //#else
  float sip = fabs(el_std_ip_gsf[index]);
  //#endif

  int cat = ElectronClassification(index);
  
  float corr_tk_iso   = el_std_tkiso03[  index];
  float corr_ecal_iso = el_std_ecaliso04[index];
  float corr_hcal_iso = el_std_hcaliso04[index];  

  corr_tk_iso   = corr_tk_iso  *pow(40/p4->Et(), 2); 
  corr_ecal_iso = corr_ecal_iso*pow(40/p4->Et(), 2);
  corr_hcal_iso = corr_hcal_iso*pow(40/p4->Et(), 2);
  
  if ((corr_tk_iso > cutisotk[bin][type][cat]) ||
      (corr_ecal_iso > cutisoecal[bin][type][cat]) ||
      (corr_hcal_iso > cutisohcal[bin][type][cat]))

    isoIDResult.first = false;
  
  if (el_std_fbrem[index] < -2) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (el_std_hoe[index]  > cuthoe[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (see > cutsee[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (fabs(el_std_dphiin[index]) > cutdphi[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;  
  }                            

  if (fabs(el_std_detain[index]) > cutdeta[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (eseedopincor < cuteopin[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }                      

  if (sip > cutip[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (el_std_hp_expin[index]  > cutmishits[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  return isoIDResult;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::ElectronClassification(int index) {
  TLorentzVector* p4 = (TLorentzVector*) el_std_sc->At(index);

  int cat = -1;
  float eta = fabs(p4->Eta());

  if (eta < 1.479) {       // BARREL
    if(el_std_fbrem[index]<0.12)
      cat=1;
    else if (el_std_eopin[index] < 1.2 && el_std_eopin[index] > 0.9) 
      cat=0;
    else 
      cat=2;
  } else {                // ENDCAP
    if(el_std_fbrem[index]<0.2)
      cat=4;
    else if (el_std_eopin[index] < 1.22 && el_std_eopin[index] > 0.82) 
      cat=3;
    else 
      cat=5;
  }

  return cat;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::sipCalculator(int index) {

  Float_t ip = 0;

  if (el_std_tkind[index] != -1) {
    TLorentzVector* tk = (TLorentzVector*)tk_p4->At(el_std_tkind[index]);
    TVector3* my_tk_pos = (TVector3*)tk_vtx_pos->At(el_std_tkind[index]);

    // FIXME to handle the case of multiple vertices

    if (vtx_std_n != 0) {
      TVector3* my_vtx_pos = (TVector3*)vtx_std_xyz->At(0);

      // this is d0 "corrected" for the vertex...
      ip = fabs((-(my_tk_pos->X()-my_vtx_pos->X())*tk->Y()+(my_tk_pos->Y()-my_vtx_pos->Y()) * tk->X())/tk->Pt());
    } else {
      ip = fabs((-(my_tk_pos->X())*tk->Y()+(my_tk_pos->Y()) * tk->X())/tk->Pt());
    }

  }

  return ip;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::eIDInfo(Int_t index, Int_t& iso_result, Int_t& id_result, Int_t eIDMaxLevel) {

  iso_result = 0;
  id_result = 0;

  // FIXME add GetEntry functions
  for(Int_t i=0; i<eIDMaxLevel; ++i) {
    std::pair<bool, bool> result = ElectronId(index, LoopAll::eIDLevel(i));

    if (result.first) 
      iso_result = i;

    if (result.second)
      id_result = i;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
// Vertex Analysis
// ---------------------------------------------------------------------------------------------------------------------------------------------
class GlobeVertexInfo : public VertexInfoAdapter
{
public:
	GlobeVertexInfo(LoopAll &);
	
	virtual int nvtx() const    { return lo_.vtx_std_n; };
	virtual int ntracks() const { return lo_.tk_n; };
	
	virtual bool hasVtxTracks() const { return true; }
	virtual const unsigned short * vtxTracks(int ii) const { return &(*lo_.vtx_std_tkind)[ii][0]; };
	virtual int vtxNTracks(int ii) const { return lo_.vtx_std_ntks[ii]; };
	virtual const float * vtxTkWeights(int ii) const { return &(*lo_.vtx_std_tkweight)[ii][0]; };

	virtual float tkpx(int ii) const { return ((TLorentzVector*)lo_.tk_p4->At(ii))->Px(); };
	virtual float tkpy(int ii) const { return ((TLorentzVector*)lo_.tk_p4->At(ii))->Py(); };
	virtual float tkpz(int ii) const { return ((TLorentzVector*)lo_.tk_p4->At(ii))->Pz(); };
	
	virtual float tkPtErr(int ii) const { return lo_.tk_pterr[ii]; };
	virtual int   tkVtxId(int ii) const { return -1; };

	virtual float tkWeight(int ii, int jj) const { return (*lo_.vtx_std_tkweight)[jj][ii]; };
	
	virtual float vtxx(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->X(); };
	virtual float vtxy(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Y(); };
	virtual float vtxz(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Z(); };

	virtual float tkd0(int ii, int jj) const { return 0.; }; // FIXME
	virtual float tkd0Err(int ii, int jj) const { return 1.; };  // FIXME

	virtual float tkdz(int ii, int jj) const { return 0.; };  // FIXME
	virtual float tkdzErr(int ii, int jj) const { return 1.; };  // FIXME

	virtual bool tkIsHighPurity(int ii) const { return ( lo_.tk_quality[ii] & (1<<2) ) >> 2; };

	virtual ~GlobeVertexInfo();
	
private:
	LoopAll & lo_;
};


// ---------------------------------------------------------------------------------------------------------------------------------------------
GlobeVertexInfo::GlobeVertexInfo(LoopAll & lo) : lo_(lo) {};

// ---------------------------------------------------------------------------------------------------------------------------------------------
GlobeVertexInfo::~GlobeVertexInfo() {};


// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::vertexAnalysis(HggVertexAnalyzer & vtxAna,  PhotonInfo pho1, PhotonInfo pho2)
{
        GlobeVertexInfo vinfo(*this); 
	//	PhotonInfo
	  //	  pho1(p1,*((TVector3*)pho_calopos->At(p1)),((TLorentzVector*)pho_p4->At(p1))->Energy()),
	  //  pho2(p2,*((TVector3*)pho_calopos->At(p2)),((TLorentzVector*)pho_p4->At(p2))->Energy());
	vtxAna.analyze(vinfo,pho1,pho2);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
PhotonInfo LoopAll::fillPhotonInfos(int p1, bool useAllConvs, float * energy) 
{
	if(LDEBUG)  cout << "  LoopAll::fillPhotonInfos with index " << p1 <<  endl;
	
	int iConv1 = useAllConvs ? matchPhotonToConversion(p1) : -1;
	if(LDEBUG)  cout << " I have matched and the index is " << iConv1 << endl;
	
	if ( iConv1 >= 0) {
		// conversions infos
		return PhotonInfo(p1,
				  // *((TVector3*)pho_calopos->At(p1)),
				  *((TVector3*)sc_xyz->At(pho_scind[p1])),
				  *((TVector3*) bs_xyz->At(0)),
				  *((TVector3*) conv_vtx->At(iConv1)),
				  *((TVector3*) conv_refitted_momentum->At(iConv1)),
				  energy == 0 ? ((TLorentzVector*)pho_p4->At(p1))->Energy() : energy[p1],
				  pho_isEB[p1],
				  conv_ntracks[iConv1],
				  conv_validvtx[iConv1],
				  conv_chi2_probability[iConv1],
				  conv_eoverp[iConv1]
			);
	} 
	//// else {
	//// 	return PhotonInfo(p1,*((TVector3*)pho_calopos->At(p1)),((TLorentzVector*)pho_p4->At(p1))->Energy());
	//// }
	
	return PhotonInfo(p1, 
			  // *((TVector3*)pho_calopos->At(p1)),                                                                                                                
			  *((TVector3*)sc_xyz->At(pho_scind[p1])),
			  *((TVector3*) bs_xyz->At(0)),                                                                                                                            
			  *((TVector3*) pho_conv_vtx->At(p1)),
			  *((TVector3*) pho_conv_refitted_momentum->At(p1)),
			  energy == 0 ? ((TLorentzVector*)pho_p4->At(p1))->Energy() : energy[p1],
			  pho_isEB[p1],                                                                                                                                     
			  pho_conv_ntracks[p1],                                                                                                                             
			  pho_conv_validvtx[p1],                                                                                                                            
			  pho_conv_chi2_probability[p1] ,                                                                                                                   
			  pho_conv_eoverp[p1]                                                                                                                               
		);
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> LoopAll::vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, 
					  PhotonInfo & pho1, PhotonInfo & pho2, std::vector<std::string> & vtxVarNames, 
					  bool useMva, TMVA::Reader * tmvaReader, std::string tmvaMethod)
{
	int p1 = pho1.id(), p2 = pho2.id();
	// assert( p1 == vtxAna.pho1() && p2 == vtxAna.pho2() );
	vtxAna.setPairID(p1,p2);
	if( useMva ) { return vtxAna.rank(*tmvaReader,tmvaMethod); }
	
	// preselect vertices : all vertices
        std::vector<int> preselAll;
        for(int i=0; i<vtx_std_n ; i++) {
          preselAll.push_back(i); 
        }

        float zconv = 0; 
        float dzconv = 0;
        std::vector<int> preselConv;

        if ( (pho1.isAConversion() || pho2.isAConversion() ) )  {
	  
          if (pho1.isAConversion()  && !pho2.isAConversion() ){
            zconv  = vtxAnaFromConv.vtxZ(pho1);
            dzconv = vtxAnaFromConv.vtxdZ(pho1);
          }
	  
          if (pho2.isAConversion() && !pho1.isAConversion()){
            zconv  = vtxAnaFromConv.vtxZ(pho2);
            dzconv = vtxAnaFromConv.vtxdZ(pho2);
          }
	  
          if ( pho1.isAConversion() && pho2.isAConversion()){
            float z1  = vtxAnaFromConv.vtxZ(pho1);
            float dz1 = vtxAnaFromConv.vtxdZ(pho1);
            
            float z2  = vtxAnaFromConv.vtxZ(pho2);
            float dz2 = vtxAnaFromConv.vtxdZ(pho2);
            
            zconv  = (z1/dz1/dz1 + z2/dz2/dz2)/(1./dz1/dz1 + 1./dz2/dz2 );  // weighted average
            dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
          }
	  
	  // preselect vertices : only vertices in a window zconv +/- dzconv
	  for(int i=0; i < vtx_std_n; i++) {
	    TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(i);
	    if ( fabs(zconv - vtxpos->Z() ) < dzconv ) 
              preselConv.push_back(i); 
          }
	  
        }
	
	// ---- METHOD 1 	
	// preselection 
	//	if ( preselConv.size()==0 )
        //  vtxAna.preselection(preselAll);
        //else
        //  vtxAna.preselection(preselConv);
	
	//std::vector<int> rankprod = vtxAna.rankprod(vtxVarNames);
	
	// ---- METHOD 1 

	
	// ---- NEW METHOD 2 (suggested by MarcoP) : first use ranking, then conversions info, e.g. on the N vtxs with best rank
	// preselection  : all vtxs
	std::vector<int> rankprodAll = useMva ? vtxAna.rank(*tmvaReader,tmvaMethod) : vtxAna.rankprod(vtxVarNames);
	int iClosestConv = -1;
	float dminconv = 9999999;
	
	TLorentzVector dipho = get_pho_p4( p1, 0 ) + get_pho_p4( p2, 0 ) ;
	
	int nbest ;
	if (  dipho.Pt() < 30 ) nbest = 5;
	else nbest = 3; 
	if (rankprodAll.size() < nbest ) nbest = rankprodAll.size();

	for (int ii = 0; ii < nbest; ii++ ){
	   TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(rankprodAll[ii]);
	   if ( fabs( vtxpos->Z()-zconv ) < dzconv && fabs(vtxpos->Z() - zconv ) < dminconv){
	     iClosestConv = rankprodAll[ii];
	     dminconv = fabs(vtxpos->Z()-zconv );
	   }
	}
	std::vector<int> rankprod;
	rankprod.clear();
	if (iClosestConv!=-1 ) rankprod.push_back(iClosestConv);

	//for (int kk = 0; kk  < nbest; kk++ ){
	for (int kk = 0; kk  < rankprodAll.size(); kk++ ){
	  if ( iClosestConv == rankprodAll[kk] ) continue;
	  else rankprod.push_back(rankprodAll[kk]);
	}
	// ---- METHOD 2


 
	return rankprod;
}

vector<double> LoopAll::generate_flat10_weights(TH1D* data_npu_estimated){
    // see
    // SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py;
    // copy and paste from there:
    double npu_probs[25] =
   {0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,
    0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,
    0.0698146584,0.0630151648,0.0526654164,0.0402754482,0.0292988928,
    0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322, 
    0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05};
    vector<double> result(25);
    double s = 0.0;

    for (int npu=0; npu<25; ++npu){
  
        s += npu_probs[npu];
    } 
    for (int npu=0; npu<25; ++npu){
  
        npu_probs[npu] *= (1./s);
    } 

    data_npu_estimated->Scale(1./data_npu_estimated->Integral());
    for (int npu=0; npu<25; ++npu){
        double npu_estimated = 
            data_npu_estimated->GetBinContent(
            data_npu_estimated->GetXaxis()->FindBin(npu));     
        result[npu] = npu_estimated / npu_probs[npu];
        //s += npu_estimated;
    }
    return result;
}


double LoopAll::phiNorm(float &phi) {

  const float pi = 3.1415927;
  const float twopi = 2.0*pi;

  if(phi >  pi) {phi = phi - twopi;}
  if(phi < -pi) {phi = phi + twopi;}

  return phi;
}


double LoopAll::etaTransformation(  float EtaParticle , float Zvertex)  {

  //---Definitions
  const float pi = 3.1415927;

  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479; 
   
  //---ETA correction

  float Theta = 0.0  ; 
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+pi ;
  double ETA = - log(tan(0.5*Theta));
         
  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle); 
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+pi ;
      ETA = - log(tan(0.5*Theta));		      
    } 
  //---Return the result
  return ETA;
  //---end
}


int  LoopAll::matchPhotonToConversion( int lpho) {

  int result=-99;
  double conv_eta=-999.;
  double conv_phi=-999.;
  
 // float sc_eta  = ((TVector3 *) pho_calopos->At(lpho))->Eta();
  float sc_eta  = ((TVector3 *) sc_xyz->At(pho_scind[lpho]))->Eta();
  float phi  = ((TVector3 *) sc_xyz->At(pho_scind[lpho]))->Phi();
//  float  phi  = ((TVector3 *) pho_calopos->At(lpho))->Phi();
  double sc_phi = phiNorm(phi);
  
  TLorentzVector * p4 = (TLorentzVector *) pho_p4->At(lpho);
  float et = p4->Energy() / cosh(sc_eta);
  //  cout << " photon index " << lpho << " eta " <<sc_eta<< " phi " << sc_phi << " et " << et << endl; 
  
  float detaMin=999.;
  float dphiMin=999.;   
  float dRMin = 999.;

  float mconv_pt=-999999;
  int iMatch=-1;     
  float conv_pt = -9999;

  if(LDEBUG)  cout << "   LoopAll::matchPhotonToConversion conv_n " << conv_n << endl; 
  for(int iconv=0; iconv<conv_n; iconv++) {
    TVector3 refittedPairMomentum= *((TVector3*) conv_refitted_momentum->At(iconv));
    conv_pt =  refittedPairMomentum.Pt();
    if (conv_pt < 1 ) continue;    
    if ( !conv_validvtx[iconv] || conv_ntracks[iconv]!=2 || conv_chi2_probability[iconv]<0.000001) continue; // Changed back based on meeting on 21.03.2012

    phi  = ((TVector3 *) conv_refitted_momentum->At(iconv))->Phi();
    conv_phi  = phiNorm(phi);
    float eta  = ((TVector3 *) conv_refitted_momentum->At(iconv))->Eta();
    conv_eta = etaTransformation(eta, conv_zofprimvtxfromtrks[iconv] );

    //    cout << " conversion index " << iconv << " eta " <<conv_eta<<  " norm phi " << conv_phi << " PT " << conv_pt << endl; 

    //  float dPhi =conv_phi - sc_phi;       
    /// double delta_phi = conv_phi - sc_phi;       
    double delta_phi = acos( cos(conv_phi - sc_phi) );       
    double delta_eta = conv_eta - sc_eta;

    //// assert( delta_phi >= 0. && delta_phi <= TMath::Pi() ); 
 
    //cout << " delta_eta " << delta_eta << " delta_phi " << delta_phi << endl;
    delta_phi*=delta_phi;
    delta_eta*=delta_eta;
    float dR = sqrt( delta_phi + delta_eta ); 
    
    //if ( fabs(delta_eta) < detaMin && fabs(delta_phi) < dphiMin ) {
     if ( dR < dRMin ) {
      detaMin=  fabs(delta_eta);
      dphiMin=  fabs(delta_phi);
      dRMin=dR;
      iMatch=iconv;
      mconv_pt = conv_pt;
    }
    
  }
  
  //  cout << " minimized conversion index " << iMatch << " eta " <<conv_eta<< " phi " << conv_phi <<endl; 

  //if ( detaMin < 0.1 && dphiMin < 0.1 ) {
  if ( dRMin< 0.1 ) {
    if(LDEBUG)    cout << " matched conversion index " << iMatch << " eta " <<conv_eta<< " phi " << conv_phi << " pt " << mconv_pt << endl; 	
    result = iMatch;
  } else {
    result = -1;
  }
  
  return result;
  

}


// ---------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector LoopAll::get_pho_p4(int ipho, int ivtx, float * energy)
{
	/// /// PhotonInfo p(ipho, *((TVector3*)sc_xyz->At(pho_scind[ipho])),
	/// PhotonInfo p(ipho, *((TVector3*)pho_calopos->At(ipho)),
	/// 	     energy != 0 ? energy[ipho] : ((TLorentzVector*)pho_p4->At(ipho))->Energy() );
	/// TVector3 * vtx = (TVector3*) vtx_std_xyz->At(ivtx);
	/// return p.p4( vtx->X(), vtx->Y(), vtx->Z() );
	TVector3 * vtx = (TVector3*) vtx_std_xyz->At(ivtx);
	return get_pho_p4(ipho,vtx,energy);
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector LoopAll::get_pho_p4(int ipho, TVector3 * vtx, float * energy)
{
	/// PhotonInfo p(ipho, *((TVector3*)sc_xyz->At(pho_scind[ipho])),
	PhotonInfo p(ipho, *((TVector3*)sc_xyz->At(pho_scind[ipho])),
		     energy != 0 ? energy[ipho] : ((TLorentzVector*)pho_p4->At(ipho))->Energy() );
	return p.p4( vtx->X(), vtx->Y(), vtx->Z() );
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::set_pho_p4(int ipho, int ivtx, float *pho_energy_array)
{
	*((TLorentzVector*)pho_p4->At(ipho)) = get_pho_p4(ipho,ivtx,pho_energy_array);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::FillCICInputs()
{
	pho_tkiso_recvtx_030_002_0000_10_01->clear(); pho_tkiso_recvtx_030_002_0000_10_01->resize(pho_n,std::vector<float>(vtx_std_n,0.));
	
	for(int ipho=0;ipho<pho_n;++ipho){
		// TLorentzVector * phop4 = (TLorentzVector*)pho_p4->At(ipho);
		std::pair<Int_t, Float_t> worse_iso = WorstSumTrackPtInCone(ipho, 0,0, 0.40, 0.02, 0.0, 1.0, 0.1); 
		pho_tkiso_badvtx_040_002_0000_10_01[ipho] = worse_iso.second;
		pho_tkiso_badvtx_id[ipho] = worse_iso.first;
		pho_drtotk_25_99[ipho] = DeltaRToTrack(ipho, vtx_std_sel, 2.5, 99.);
		
		for(int ivtx=0;ivtx<vtx_std_n;++ivtx) {
			TLorentzVector p4 = get_pho_p4( ipho, ivtx );
			(*pho_tkiso_recvtx_030_002_0000_10_01)[ipho][ivtx] = SumTrackPtInCone(&p4, ivtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::FillCIC()
{
	// 6 categories
	pho_cic6passcuts_lead->clear(); pho_cic6passcuts_lead->resize( pho_n, std::vector<std::vector<UInt_t> >(vtx_std_n, std::vector<UInt_t> (phoNCUTLEVELS,0) ) ); 
	pho_cic6cutlevel_lead->clear(); pho_cic6cutlevel_lead->resize( pho_n, std::vector<Short_t>(vtx_std_n,0) );
	pho_cic6passcuts_sublead->clear(); pho_cic6passcuts_sublead->resize( pho_n, std::vector<std::vector<UInt_t> >(vtx_std_n,std::vector<UInt_t>(phoNCUTLEVELS,0)) ); 
	pho_cic6cutlevel_sublead->clear(); pho_cic6cutlevel_sublead->resize( pho_n, std::vector<Short_t>(vtx_std_n,0) );
	std::vector<std::vector<bool> > cic6_passcut_lead, cic6_passcut_sublead;
	// 4 categories
	pho_cic4passcuts_lead->clear(); pho_cic4passcuts_lead->resize( pho_n, std::vector<std::vector<UInt_t> >(vtx_std_n, std::vector<UInt_t>(phoNCUTLEVELS,0) ) ); 
	pho_cic4cutlevel_lead->clear(); pho_cic4cutlevel_lead->resize( pho_n, std::vector<Short_t>(vtx_std_n,0) );
	pho_cic4passcuts_sublead->clear(); pho_cic4passcuts_sublead->resize( pho_n, std::vector<std::vector<UInt_t> >(vtx_std_n, std::vector<UInt_t>(phoNCUTLEVELS,0) ) ); 
	pho_cic4cutlevel_sublead->clear(); pho_cic4cutlevel_sublead->resize( pho_n, std::vector<Short_t>(vtx_std_n,0) );
	std::vector<std::vector<bool> > cic4_passcut_lead, cic4_passcut_sublead;
	
	for(int ipho=0;ipho<pho_n;++ipho){
		for(int ivtx=0;ivtx<vtx_std_n;++ivtx){
		// 6 categories
			int cic6_level_lead = PhotonCiCSelectionLevel(ipho, ivtx, cic6_passcut_lead, 6, 0);
			int cic6_level_sublead = PhotonCiCSelectionLevel(ipho, ivtx, cic6_passcut_sublead, 6, 1);
			(*pho_cic6cutlevel_lead)[ipho][ivtx] = cic6_level_lead;
			(*pho_cic6cutlevel_sublead)[ipho][ivtx] = cic6_level_sublead;
			// 4 categories
			int cic4_level_lead = PhotonCiCSelectionLevel(ipho, ivtx, cic4_passcut_lead, 4, 0);
			int cic4_level_sublead = PhotonCiCSelectionLevel(ipho, ivtx, cic4_passcut_sublead, 4, 1);
			(*pho_cic4cutlevel_lead)[ipho][ivtx] = cic4_level_lead;
			(*pho_cic4cutlevel_sublead)[ipho][ivtx] = cic4_level_sublead;
		
			for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {
				// 6 categories
				UInt_t cic6_leadw=0, cic6_subleadw=0;
				for(size_t icut=0; icut<cic6_passcut_lead.size(); ++icut) {
					cic6_leadw |= ( (!cic6_passcut_lead[iCUTLEVEL][icut] & 0x1) << icut);
				}
				(*pho_cic6passcuts_lead)[ipho][ivtx][iCUTLEVEL] = cic6_leadw;
				for(size_t icut=0; icut<cic6_passcut_sublead.size(); ++icut) {
					cic6_subleadw |= ( (!cic6_passcut_sublead[iCUTLEVEL][icut] & 0x1) << icut);
				}
				(*pho_cic6passcuts_sublead)[ipho][ivtx][iCUTLEVEL] = cic6_subleadw;
				// 4 categories
				UInt_t cic4_leadw=0, cic4_subleadw=0;
				for(size_t icut=0; icut<cic4_passcut_lead.size(); ++icut) {
					cic4_leadw |= ( (!cic4_passcut_lead[iCUTLEVEL][icut] & 0x1) << icut);
				}
				(*pho_cic4passcuts_lead)[ipho][ivtx][iCUTLEVEL] = cic4_leadw;
				for(size_t icut=0; icut<cic4_passcut_sublead.size(); ++icut) {
					cic4_subleadw |= ( (!cic4_passcut_sublead[iCUTLEVEL][icut] & 0x1) << icut);
				}
				(*pho_cic4passcuts_sublead)[ipho][ivtx][iCUTLEVEL] = cic4_subleadw;
			}
		}
	}
}


// CiC SELECTION CODE BEGIN - SSIMON
// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cic6_allcuts_lead, float * cic6_allcuts_sublead, float * cic4_allcuts_lead, float * cic4_allcuts_sublead) {

  //thresholds are in this order below
  // isosumoet[]
  // isosumoetbad[]
  // trkisooetom[]
  // sieie[]
  // hovere[]
  // r9[]
  // drtotk_25_99[]
  // pixel[]

// 6 categories
// phoNOCUTS      - thresholds so all photons pass
// phoLOOSE       - sob value=0.0002         - iteration 8 - eff=0.947448  fake=0.0783937
// phoMEDIUM      - sob value=0.0004         - iteration 6 - eff=0.928017  fake=0.0572683
// phoTIGHT       - sob value=0.0008         - iteration 6 - eff=0.895238  fake=0.0392572
// phoSUPERTIGHT  - sob value=0.0016         - iteration 6 - eff=0.849812  fake=0.0256949
// phoHYPERTIGHT1 - sob value=0.0032         - iteration 6 - eff=0.784283  fake=0.016346
// phoHYPERTIGHT2 - sob value=0.00625        - iteration 6 - eff=0.699128  fake=0.00991561
// phoHYPERTIGHT3 - sob value=0.0125         - iteration 6 - eff=0.573171  fake=0.00520159
// phoHYPERTIGHT4 - sob value=0.025          - iteration 6 - eff=0.41176   fake=0.00217666


// 4 categories
// phoNOCUTS      - thresholds so all photons pass
// phoLOOSE       - sob value=0.0002         - iteration 8 - eff=0.939229  fake=0.0815158
// phoMEDIUM      - sob value=0.0004         - iteration 6 - eff=0.91754   fake=0.0581047
// phoTIGHT       - sob value=0.0008         - iteration 6 - eff=0.886869  fake=0.041063
// phoSUPERTIGHT  - sob value=0.0016         - iteration 6 - eff=0.844314  fake=0.0286033
// phoHYPERTIGHT1 - sob value=0.0032         - iteration 6 - eff=0.774552  fake=0.0191603
// phoHYPERTIGHT2 - sob value=0.00625        - iteration 6 - eff=0.67859   fake=0.0121262
// phoHYPERTIGHT3 - sob value=0.0125         - iteration 6 - eff=0.521328  fake=0.00626966
// phoHYPERTIGHT4 - sob value=0.025          - iteration 6 - eff=0.381192  fake=0.00357649



  const unsigned int ncuts = 8;
  const unsigned int ncat_cic6 = 6;
  const unsigned int ncat_cic4 = 4;
  switch(cutlevel) {

    std::cout << "FFFSSSSSSSSSSSSSSSSSS!!!!!!!!!!!!" << cutlevel << std::endl;

    case(phoNOCUTS) : {
                        float cic6_allcuts_temp_lead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                        float cic6_allcuts_temp_sublead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                        float cic4_allcuts_temp_lead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5 };
                        float cic4_allcuts_temp_sublead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5 };
                        for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                          cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                          for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                      } break;
    case(phoLOOSE) : {
                       float cic6_allcuts_temp_lead[] = { 
                         14.1278,     11.7187,     9.78826,     10.9814,     9.21945,     8.89621,
                         72.5178,     59.1506,     85.1822,     93.8969,     74.2109,     14.4058,
                         7.89015,     5.61652,     4.45536,     5.87563,     4.24725,     2.96206,
                         0.0114196,   0.0109898,   0.0100549,    0.029265,   0.0290002,   0.0279397,
                         0.0907646,   0.0791189,   0.0835245,    0.102617,   0.0596196,    0.098899,
                         0.94,    0.899976,    0.262285,    0.94,    0.90,    0.276953,
                         12.0314,     98.0038,  0.00968623,  0.00636153,  0.00476398,  0.00610842,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                       float cic6_allcuts_temp_sublead[] = { 
                         14.1278,     11.7187,     9.78826,     10.9814,     9.21945,     8.89621,
                         72.5178,     59.1506,     85.1822,     93.8969,     74.2109,     14.4058,
                         7.89015,     5.61652,     4.45536,     5.87563,     4.24725,     2.96206,
                         0.0114196,   0.0109898,   0.0100549,    0.029265,   0.0290002,   0.0279397,
                         0.0907646,   0.0791189,   0.0835245,    0.102617,   0.0596196,    0.098899,
                         0.94,    0.899976,    0.262285,    0.94,    0.90,    0.276953,
                         12.0314,     98.0038,  0.00968623,  0.00636153,  0.00476398,  0.00610842,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                       float cic4_allcuts_temp_lead[] = { 
			       8.2,       4.1,       5.4,       2.6,
			       67,        69,        85,       7.2,
			       7.5,       4.5,       5.2,       2.5,
			       0.0112,    0.0102,     0.029,     0.028,
			       0.09,     0.089,     0.101,     0.073,
			       0.94,      0.31,      0.92,      0.29,
			       0.26,     0.029,    0.0062,    0.0055,
			       1.5,         1.5,         1.5,         1.5};
                       float cic4_allcuts_temp_sublead[] = { 
				       8.2,       4.1,       5.4,       2.6,
				       67,        69,        85,       7.2,
				       7.5,       4.5,       5.2,       2.5,
				       0.0112,    0.0102,     0.029,     0.028,
				       0.09,     0.089,     0.101,     0.073,
				       0.94,      0.31,      0.92,      0.29,
				       0.26,     0.029,    0.0062,    0.0055,
				       1.5,         1.5,         1.5,         1.5};
		       for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                         cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                         for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                           cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                     } break;
    case(phoMEDIUM) : {
                        float cic6_allcuts_temp_lead[] = {  
                          12.5084,      10.156,     9.23141,     10.0482,     8.34498,     8.73704,
                          70.9011,     50.0742,     21.9926,     24.2436,     18.7884,     12.6882,
                          6.58797,     4.68564,     4.38815,     5.67876,     2.41162,     2.19991,
                          0.0110266,   0.0106749,  0.00983011,   0.0287021,   0.0286817,   0.0272739,
                          0.0891215,   0.0763711,   0.0798623,   0.0911974,   0.0511163,   0.0627764,
                          0.94,    0.90,    0.274434,    0.94,    0.90,    0.276953,
                          96.5654,     98.9721,   0.0119942,   0.0111399,  0.00855448,    0.012159,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                        float cic6_allcuts_temp_sublead[] = {  
                          12.5084,      10.156,     9.23141,     10.0482,     8.34498,     8.73704,
                          70.9011,     50.0742,     21.9926,     24.2436,     18.7884,     12.6882,
                          6.58797,     4.68564,     4.38815,     5.67876,     2.41162,     2.19991,
                          0.0110266,   0.0106749,  0.00983011,   0.0287021,   0.0286817,   0.0272739,
                          0.0891215,   0.0763711,   0.0798623,   0.0911974,   0.0511163,   0.0627764,
                          0.94,    0.90,    0.274434,    0.94,    0.90,    0.276953,
                          96.5654,     98.9721,   0.0119942,   0.0111399,  0.00855448,    0.012159,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                        float cic4_allcuts_temp_lead[] = {  
				6.4,       3.2,       3.4,       2.2,
				64,      10.8,        13,       3.5,
				6.4,       3.4,       3.8,       2.1,
				0.0109,      0.01,     0.029,     0.028,
				0.089,     0.079,      0.09,     0.061,
				0.94,      0.32,      0.94,      0.29,
				0.98,     0.029,    0.0109,    0.0111,
				1.5,         1.5,         1.5,         1.5 };
                        float cic4_allcuts_temp_sublead[] = {  
				6.4,       3.2,       3.4,       2.2,
				64,      10.8,        13,       3.5,
				6.4,       3.4,       3.8,       2.1,
				0.0109,      0.01,     0.029,     0.028,
				0.089,     0.079,      0.09,     0.061,
				0.94,      0.32,      0.94,      0.29,
				0.98,     0.029,    0.0109,    0.0111,
				1.5,         1.5,         1.5,         1.5 };
                        for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                          cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                          for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                      } break;
    case(phoTIGHT) : {
                       float cic6_allcuts_temp_lead[] = { 
                         11.1845,     9.28445,     8.98759,     9.19055,     7.94171,     8.16991,
                         70.7835,     16.7873,     13.7361,     15.6259,     13.2407,     10.3932,
                         5.76122,     3.97439,     2.89137,     4.62749,     2.34848,      1.9302,
                         0.010781,   0.0104673,  0.00965497,   0.0284936,    0.028082,   0.0270328,
                         0.0844869,   0.0703749,    0.060775,   0.0881813,   0.0502974,    0.060877,
                         0.94,    0.90,       0.321,    0.94,    0.90,    0.279956,
                         98.9318,     98.9992,   0.0146256,   0.0207672,     34.1809,   0.0261029,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                       float cic6_allcuts_temp_sublead[] = { 
                         11.1845,     9.28445,     8.98759,     9.19055,     7.94171,     8.16991,
                         70.7835,     16.7873,     13.7361,     15.6259,     13.2407,     10.3932,
                         5.76122,     3.97439,     2.89137,     4.62749,     2.34848,      1.9302,
                         0.010781,   0.0104673,  0.00965497,   0.0284936,    0.028082,   0.0270328,
                         0.0844869,   0.0703749,    0.060775,   0.0881813,   0.0502974,    0.060877,
                         0.94,    0.90,       0.321,    0.94,    0.90,    0.279956,
                         98.9318,     98.9992,   0.0146256,   0.0207672,     34.1809,   0.0261029,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                       float cic4_allcuts_temp_lead[] = { 
			       4.7,       2.8,       2.5,      1.46,
			       62,       5.2,       7.3,       2.5,
			       4.7,       2.9,       3.8,      1.63,
			       0.0107,    0.0099,     0.028,     0.027,
			       0.087,     0.065,     0.087,      0.05,
			       0.94,      0.34,      0.94,      0.29,
			       1,     0.029,     0.021,     0.028,
			       1.5,         1.5,         1.5,         1.5 };
                       float cic4_allcuts_temp_sublead[] = { 
			       4.7,       2.8,       2.5,      1.46,
			       62,       5.2,       7.3,       2.5,
			       4.7,       2.9,       3.8,      1.63,
			       0.0107,    0.0099,     0.028,     0.027,
			       0.087,     0.065,     0.087,      0.05,
			       0.94,      0.34,      0.94,      0.29,
			       1,     0.029,     0.021,     0.028,
			       1.5,         1.5,         1.5,         1.5 };
                       for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                         cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                         for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                           cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                     } break;
    case(phoSUPERTIGHT) : {
                            float cic6_allcuts_temp_lead[] = { 
                              10.0171,     8.81037,     8.74909,     8.47393,     7.94171,     7.47883,
                              54.9366,     14.3545,     11.5208,      12.939,     10.2496,      9.7095,
                              4.11252,     3.35092,     2.49296,     2.05592,     1.67021,     1.66678,
                              0.0106315,   0.0101656,  0.00950936,   0.0283215,   0.0276216,   0.0263378,
                              0.0823828,   0.0598641,   0.0494497,   0.0706222,   0.0502974,    0.060877,
                              0.94,    0.90,       0.321,    0.94,    0.90,    0.282153,
                              98.9981,          99,   0.0216484,     96.2292,     97.1855,     96.2294,
                              1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                            float cic6_allcuts_temp_sublead[] = { 
                              10.0171,     8.81037,     8.74909,     8.47393,     7.94171,     7.47883,
                              54.9366,     14.3545,     11.5208,      12.939,     10.2496,      9.7095,
                              4.11252,     3.35092,     2.49296,     2.05592,     1.67021,     1.66678,
                              0.0106315,   0.0101656,  0.00950936,   0.0283215,   0.0276216,   0.0263378,
                              0.0823828,   0.0598641,   0.0494497,   0.0706222,   0.0502974,    0.060877,
                              0.94,    0.90,       0.321,    0.94,    0.90,    0.282153,
                              98.9981,          99,   0.0216484,     96.2292,     97.1855,     96.2294,
                              1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                            float cic4_allcuts_temp_lead[] = { 
                              3.8,     2.2,     1.77,    1.29,
                              11.7,    3.4,     3.9,     1.84,
                              3.5,     2.2,     2.3,     1.45,
                              0.0106,  0.0097,  0.028,   0.027,
                              0.082,   0.062,   0.065,   0.048,
                              0.94,    0.36,    0.94,    0.32,
                              1.,      0.062,   0.97,    0.97,
                              1.5,     1.5,     1.5,     1.5 };
                            float cic4_allcuts_temp_sublead[] = { 
                              3.8,     2.2,     1.77,    1.29,
                              11.7,    3.4,     3.9,     1.84,
                              3.5,     2.2,     2.3,     1.45,
                              0.0106,  0.0097,  0.028,   0.027,
                              0.082,   0.062,   0.065,   0.048,
                              0.94,    0.36,    0.94,    0.32,
                              1.,      0.062,   0.97,    0.97,
                              1.5,     1.5,     1.5,     1.5 };
                            for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                              cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                              for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                          } break;
    case(phoHYPERTIGHT1) : {
                             float cic6_allcuts_temp_lead[] = { 
                               9.14323,     8.13617,     7.43416,     7.97795,     5.88227,     6.60691,
                               16.4126,     10.7813,     10.1764,     11.3829,     8.63128,     8.75289,
                               3.49873,     2.93013,     2.00419,     1.60673,     1.36163,     1.36132,
                               0.0105033,  0.00999387,  0.00946607,   0.0282088,   0.0273334,   0.0256399,
                               0.0782034,   0.0598641,   0.0273668,   0.0553324,   0.0502974,   0.0465477,
                               0.94,    0.90,    0.347653,    0.94,    0.90,    0.301546,
                               98.9999,          99,     1.92089,     98.9224,     98.9492,     98.9224,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               9.14323,     8.13617,     7.43416,     7.97795,     5.88227,     6.60691,
                               16.4126,     10.7813,     10.1764,     11.3829,     8.63128,     8.75289,
                               3.49873,     2.93013,     2.00419,     1.60673,     1.36163,     1.36132,
                               0.0105033,  0.00999387,  0.00946607,   0.0282088,   0.0273334,   0.0256399,
                               0.0782034,   0.0598641,   0.0273668,   0.0553324,   0.0502974,   0.0465477,
                               0.94,    0.90,    0.347653,    0.94,    0.90,    0.301546,
                               98.9999,          99,     1.92089,     98.9224,     98.9492,     98.9224,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_lead[] = { 
				     3.2,      1.76,      1.39,      1.18,
				     6.1,       2.7,       2.8,      0.66,
				     3.4,      1.86,      1.67,      1.44,
				     0.0104,    0.0094,     0.028,     0.025,
				     0.076,      0.03,     0.047,     0.046,
				     0.94,      0.41,      0.94,      0.34,
				     1,      0.97,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = { 
				     3.2,      1.76,      1.39,      1.18,
				     6.1,       2.7,       2.8,      0.66,
				     3.4,      1.86,      1.67,      1.44,
				     0.0104,    0.0094,     0.028,     0.025,
				     0.076,      0.03,     0.047,     0.046,
				     0.94,      0.41,      0.94,      0.34,
				     1,      0.97,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
    case(phoHYPERTIGHT2) : {
                             float cic6_allcuts_temp_lead[] = { 
                               8.57184,     6.64014,     6.82022,     7.13109,     5.88011,      6.2565,
                               13.4065,     10.4316,     9.18551,     9.30193,     7.51729,     7.30382,
                               2.73319,     2.93013,     1.55723,     1.54876,     1.05254,     1.36132,
                               0.0103615,  0.00978982,  0.00940152,   0.0279141,   0.0260354,   0.0241246,
                               0.0572816,   0.0232443,   0.0173437,   0.0553324,   0.0365276,   0.0465477,
                               0.94,    0.90,    0.367082,    0.94,    0.90,    0.579434,
                               99,          99,     96.2824,     98.9978,     98.9986,     98.9978,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               8.57184,     6.64014,     6.82022,     7.13109,     5.88011,      6.2565,
                               13.4065,     10.4316,     9.18551,     9.30193,     7.51729,     7.30382,
                               2.73319,     2.93013,     1.55723,     1.54876,     1.05254,     1.36132,
                               0.0103615,  0.00978982,  0.00940152,   0.0279141,   0.0260354,   0.0241246,
                               0.0572816,   0.0232443,   0.0173437,   0.0553324,   0.0365276,   0.0465477,
                               0.94,    0.90,    0.367082,    0.94,    0.90,    0.579434,
                               99,          99,     96.2824,     98.9978,     98.9986,     98.9978,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_lead[] = { 
				     2.6,      1.31,      1.33,      0.82,
				     5.1,      1.62,      1.38,  -0.224864,
				     2.9,       1.6,      1.55,      1.44,
				     0.0101,    0.0093,     0.027,     0.023,
				     0.048,    0.0189,     0.032,    0.0085,
				     0.94,      0.47,      0.94,      0.52,
				     1,         1,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = { 
				     2.6,      1.31,      1.33,      0.82,
				     5.1,      1.62,      1.38,  -0.224864,
				     2.9,       1.6,      1.55,      1.44,
				     0.0101,    0.0093,     0.027,     0.023,
				     0.048,    0.0189,     0.032,    0.0085,
				     0.94,      0.47,      0.94,      0.52,
				     1,         1,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
    case(phoHYPERTIGHT3) : {
                             float cic6_allcuts_temp_lead[] = { 
                               7.97897,     6.64014,     6.60332,     5.14765,     5.02192,     5.72775,
                               11.3476,     8.93788,     8.36279,     7.88566,     5.83093,     6.66771,
                               2.348,     2.59173,     1.55158,     1.54876,     0.98618,     1.06927,
                               0.0100676,  0.00971589,  0.00932669,   0.0279141,    0.025781,   0.0229432,
                               0.0372854,   0.0215628,   0.0132992,   0.0412051,   0.0322458,   0.0465477,
                               0.94,    0.90,    0.375623,    0.94,    0.90,    0.579434,
                               99,          99,     98.9239,     98.9999,     98.9997,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               7.97897,     6.64014,     6.60332,     5.14765,     5.02192,     5.72775,
                               11.3476,     8.93788,     8.36279,     7.88566,     5.83093,     6.66771,
                               2.348,     2.59173,     1.55158,     1.54876,     0.98618,     1.06927,
                               0.0100676,  0.00971589,  0.00932669,   0.0279141,    0.025781,   0.0229432,
                               0.0372854,   0.0215628,   0.0132992,   0.0412051,   0.0322458,   0.0465477,
                               0.94,    0.90,    0.375623,    0.94,    0.90,    0.579434,
                               99,          99,     98.9239,     98.9999,     98.9997,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_lead[] = { 
				     1.85,      0.96,      1.21,  -0.028513,
				     3.7,      0.97,      1.38,  -0.880416,
				     1.93,       1.4,      1.48,     0.056,
				     0.0099,    0.0092,     0.027,     0.023,
				     0.042,    0.0173,     0.023,    0.0085,
				     0.94,      0.69,      0.97,      0.52,
				     1,         1,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = { 
				     1.85,      0.96,      1.21,  -0.028513,
				     3.7,      0.97,      1.38,  -0.880416,
				     1.93,       1.4,      1.48,     0.056,
				     0.0099,    0.0092,     0.027,     0.023,
				     0.042,    0.0173,     0.023,    0.0085,
				     0.94,      0.69,      0.97,      0.52,
				     1,         1,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
    case(phoHYPERTIGHT4) : {
                             float cic6_allcuts_temp_lead[] = { 
                               6.53539,     6.07874,     5.51521,     4.78731,     5.00511,     4.90969,
                               9.30747,      8.0574,     7.70153,     7.43339,     5.44326,     6.66771,
                               1.96543,     1.78829,    0.819072,     1.54876,     0.98618,    0.255192,
                               0.0100676,  0.00919753,  0.00911379,   0.0278098,   0.0249354,   0.0221531,
                               0.03099,   0.0153957,   0.0132992,   0.0214415,   0.0322458,   0.0138186,
                               0.94,    0.90,    0.397401,    0.94,    0.90,     0.68715,
                               99,          99,     98.9979,          99,          99,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               6.53539,     6.07874,     5.51521,     4.78731,     5.00511,     4.90969,
                               9.30747,      8.0574,     7.70153,     7.43339,     5.44326,     6.66771,
                               1.96543,     1.78829,    0.819072,     1.54876,     0.98618,    0.255192,
                               0.0100676,  0.00919753,  0.00911379,   0.0278098,   0.0249354,   0.0221531,
                               0.03099,   0.0153957,   0.0132992,   0.0214415,   0.0322458,   0.0138186,
                               0.94,    0.90,    0.397401,    0.94,    0.90,     0.68715,
                               99,          99,     98.9979,          99,          99,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_lead[] = { 
				     1.31,       0.3,      1.15,  -0.028513,
				     1.72,      0.69,      1.14,  -0.880416,
				     1.42,      0.76,      1.48,     0.056,
				     0.0098,     0.009,     0.026,     0.023,
				     0.037,   0.00049,    0.0198,   0.00024,
				     0.94,      0.69,      0.97,      0.73,
				     1,         1,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = { 
				     1.31,       0.3,      1.15,  -0.028513,
				     1.72,      0.69,      1.14,  -0.880416,
				     1.42,      0.76,      1.48,     0.056,
				     0.0098,     0.009,     0.026,     0.023,
				     0.037,   0.00049,    0.0198,   0.00024,
				     0.94,      0.69,      0.97,      0.73,
				     1,         1,         1,         1,
				     1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
    default:std::cout << "UNKNOWN phoCiCIDLevel: " << cutlevel << std::endl;

  }
}



// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::DiphotonCiCSelection( phoCiCIDLevel LEADCUTLEVEL, phoCiCIDLevel SUBLEADCUTLEVEL, 
				   Float_t leadPtMin, Float_t subleadPtMin, int ncategories, bool applyPtoverM, 
				   float *pho_energy_array, bool split) {

  //rho=0;// CAUTION SETTING RHO TO 0 FOR 2010 DATA FILES (RHO ISN'T IN THESE FILES)
  int selected_lead_index = -1;
  int selected_sublead_index = -1;
  float selected_lead_pt = -1;
  float selected_sublead_pt = -1;
  
  std::vector<int> passing_dipho;
  std::vector<float> passing_sumpt;
  for(int idipho = 0; idipho < dipho_n; ++idipho ) {
	  int ivtx = dipho_vtxind[idipho];
	  int lead = dipho_leadind[idipho];
	  int sublead = dipho_subleadind[idipho];
	  
	  if( lead == sublead ) { continue; }

	  TLorentzVector lead_p4 = get_pho_p4(lead,ivtx,pho_energy_array); 
	  TLorentzVector sublead_p4 = get_pho_p4(sublead,ivtx,pho_energy_array); 
	  
	  float leadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[lead]))->Eta());
	  float subleadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[sublead]))->Eta());
	  float m_gamgam = (lead_p4+sublead_p4).M();

	  if( leadEta > 2.5 || subleadEta > 2.5 || 
	      ( leadEta > 1.4442 && leadEta < 1.566 ) ||
	      ( subleadEta > 1.4442 && subleadEta < 1.566 ) ) { continue; }

	  float leadpt = lead_p4.Pt() > sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt();
      	  float subleadpt = lead_p4.Pt() < sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt(); 	  

	  if( applyPtoverM ) {
	    if ( leadpt/m_gamgam < leadPtMin/120. || subleadpt/m_gamgam < subleadPtMin/120. ) { continue; }
	  } else {
	    if ( leadpt < leadPtMin || subleadpt < subleadPtMin ) { continue; }
	  }
	  
	  if (subleadpt > leadpt){ // Swap them
		int tmp = lead;
		lead = sublead;
		sublead =tmp;
		dipho_leadind[idipho] = lead;
		dipho_subleadind[idipho] = sublead;
	  }

	  std::vector<std::vector<bool> > ph_passcut;
	  if( PhotonCiCSelectionLevel(lead, ivtx, ph_passcut, ncategories, 0, pho_energy_array ) < LEADCUTLEVEL ) { continue; }
	  if( PhotonCiCSelectionLevel(sublead, ivtx, ph_passcut, ncategories, 1, pho_energy_array ) < SUBLEADCUTLEVEL ) { continue; }
	  
	  passing_dipho.push_back(idipho);
	  passing_sumpt.push_back(lead_p4.Et()+sublead_p4.Et());
  }
  
  if( passing_dipho.empty() ) { return -1; }

  std::sort(passing_dipho.begin(),passing_dipho.end(),
	    SimpleSorter<float,std::greater<double> >(&passing_sumpt[0]));

  return passing_dipho[0];

  ///// std::vector<std::vector<bool> > ph_passcut;
  ///// for(int ipho=0;ipho!=pho_n;++ipho) {
  /////   TLorentzVector * iphop4 = (TLorentzVector*)pho_p4_array->At(ipho);
  /////   /// float scEta = fabs(((TVector3 *)pho_calopos->At(ipho))->Eta());
  /////   float scEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[ipho]))->Eta());
  /////   if(iphop4->Et() < leadPtMin || scEta > 2.5 || ( scEta > 1.4442 && scEta < 1.566 ) )continue;
  ///// 
  /////   if(PhotonCiCSelectionLevel(ipho, ivtx, ph_passcut, ncategories, 0) < LEADCUTLEVEL)continue;
  ///// 
  /////   for(int iipho=0;iipho!=pho_n;++iipho) {
  /////     if(iipho == ipho)continue;
  /////     TLorentzVector * iiphop4 = (TLorentzVector*)pho_p4_array->At(iipho);
  /////     /// float iiscEta = fabs(((TVector3 *)pho_calopos->At(iipho))->Eta());
  /////     float iiscEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[iipho]))->Eta());
  /////     if(iiphop4->Et() < subleadPtMin || iiscEta > 2.5 || ( iiscEta > 1.4442 && iiscEta < 1.566 ) )continue;
  /////     if(iiphop4->Et() > iphop4->Et())continue;
  /////     float m_gamgam = (*iphop4+*iiphop4).M();
  /////     float L_ptom = iphop4->Et()/m_gamgam;
  /////     float S_ptom = iiphop4->Et()/m_gamgam;
  /////     if(applyPtoverM && (L_ptom < 0.33 || S_ptom<0.25))continue;
  ///// 
  /////     if(PhotonCiCSelectionLevel(iipho, ivtx, ph_passcut, ncategories, 1) < SUBLEADCUTLEVEL )continue;
  /////     // if here, diphoton passed all cuts.
  /////     //std::cout << "FOUND DIPHOTON" << std::endl;
  /////     if( (iphop4->Et()>selected_lead_pt) || (ipho==selected_lead_index&&iiphop4->Et()>selected_sublead_pt) ) {
  ///// 	      selected_lead_pt = iphop4->Et();
  ///// 	      selected_sublead_pt = iiphop4->Et();
  ///// 	      selected_lead_index = ipho;
  ///// 	      selected_sublead_index = iipho;
  /////     }
  /////     
  /////   }// end photon loop (iipho), aka sublead
  ///// }// end photon loop (ipho), aka lead
  ////  
  ////  std::pair<int,int> dipho_inds(selected_lead_index,selected_sublead_index);
  ////  return dipho_inds;

}

int LoopAll::DiphotonMITPreSelection(Float_t leadPtMin, Float_t subleadPtMin,bool applyPtoverM,float *pho_energy_array) {

  //rho=0;// CAUTION SETTING RHO TO 0 FOR 2010 DATA FILES (RHO ISN'T IN THESE FILES)
  int selected_lead_index = -1;
  int selected_sublead_index = -1;
  float selected_lead_pt = -1;
  float selected_sublead_pt = -1;
  
  std::vector<int> passing_dipho;
  std::vector<float> passing_sumpt;
  for(int idipho = 0; idipho < dipho_n; ++idipho ) {
	  int ivtx = dipho_vtxind[idipho];
	  int lead = dipho_leadind[idipho];
	  int sublead = dipho_subleadind[idipho];
	  
	  if( lead == sublead ) { continue; }

	  TLorentzVector lead_p4 = get_pho_p4(lead,ivtx,pho_energy_array); 
	  TLorentzVector sublead_p4 = get_pho_p4(sublead,ivtx,pho_energy_array); 
	  
	  float leadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[lead]))->Eta());
	  float subleadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[sublead]))->Eta());
	  float m_gamgam = (lead_p4+sublead_p4).M();

	  if( leadEta > 2.5 || subleadEta > 2.5 || 
	      ( leadEta > 1.4442 && leadEta < 1.566 ) ||
	      ( subleadEta > 1.4442 && subleadEta < 1.566 ) ) { continue; }



	  float leadpt = lead_p4.Pt() > sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt();
      	  float subleadpt = lead_p4.Pt() < sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt(); 	  

   if (run==173439 && lumis==155 && event==236661419) {
 std::cout <<	"Inside diphotonSelection now" <<std::endl;
 std::cout <<	"mass - " <<m_gamgam<<std::endl;
 std::cout <<	"lead ET - " <<leadpt<<std::endl;
 std::cout <<	"sublead ET - " <<subleadpt<<std::endl;
}
	  if( applyPtoverM ) {
	    if ( leadpt/m_gamgam < leadPtMin/120. || subleadpt/m_gamgam < subleadPtMin/120. ) { continue; }
	  } else {
	    if ( leadpt < leadPtMin || subleadpt < subleadPtMin ) { continue; }
	  }
	  
	  if (subleadpt > leadpt){ // Swap them
		int tmp = lead;
		lead = sublead;
		sublead =tmp;
		dipho_leadind[idipho] = lead;
		dipho_subleadind[idipho] = sublead;
	  }

	  std::vector<std::vector<bool> > ph_passcut;
          if (!( PhotonMITPreSelection(lead, ivtx, pho_energy_array ) && PhotonMITPreSelection(sublead, ivtx,  pho_energy_array ))) continue; 
	  
	  passing_dipho.push_back(idipho);
	  passing_sumpt.push_back(lead_p4.Et()+sublead_p4.Et());
  }
  
  if( passing_dipho.empty() ) { return -1; }

  std::sort(passing_dipho.begin(),passing_dipho.end(),
	    SimpleSorter<float,std::greater<double> >(&passing_sumpt[0]));
  int selected_dipho_ind     = passing_dipho[0];
  int selected_dipho_lead    = dipho_leadind[selected_dipho_ind];
  int selected_dipho_sublead = dipho_subleadind[selected_dipho_ind];
  int selected_dipho_vtx     = dipho_vtxind[selected_dipho_ind];
  TLorentzVector selected_lead_p4 = get_pho_p4(selected_dipho_lead,selected_dipho_vtx,pho_energy_array); 
  TLorentzVector selected_sublead_p4 = get_pho_p4(selected_dipho_sublead,selected_dipho_vtx,pho_energy_array);
   if (run==173439 && lumis==155 && event==236661419) {
 std::cout <<	"Inside diphotonSelection now" <<std::endl;
std::cout << photonIDMVA(selected_dipho_lead,selected_dipho_vtx,selected_lead_p4,"MIT")<<std::endl;
std::cout << photonIDMVA(selected_dipho_sublead,selected_dipho_vtx,selected_sublead_p4,"MIT")<<std::endl;

}
 
  if ( photonIDMVA(selected_dipho_lead,selected_dipho_vtx,selected_lead_p4,"MIT") <= -0.3
    || photonIDMVA(selected_dipho_sublead,selected_dipho_vtx,selected_sublead_p4,"MIT")	<= -0.3
     ) {return -1;}
   if (run==173439 && lumis==155 && event==236661419) {
 std::cout <<	"Diphoton selected index -- "<< selected_dipho_ind <<std::endl;
}
  
  return selected_dipho_ind;

}
// Define newfunction to calculate MIT (Pre-)Selection                                                      
bool LoopAll::PhotonMITPreSelection( int photon_index, int vertex_index, float *pho_energy_array ) {


   int r9_category = (int) (pho_r9[photon_index] <= 0.9);                                                      
   int photon_category = r9_category + 2*PhotonEtaCategory(photon_index,2);                                 
   int photon_cic_category = PhotonCategory(photon_index,2,2);
   
   float mitCuts_hoe[4]                 = {0.082,0.075,0.075,0.075};                                        
   float mitCuts_sieie[4]               = {0.014,0.014,0.034,0.034};                                        
   float mitCuts_ecaliso[4]             = {50,4,50,4};                                                      
   float mitCuts_hcaliso[4]             = {50,4,50,4};                                                      
   float mitCuts_trkiso[4]              = {50,4,50,4};                                                      
   float mitCuts_hcalecal[4]            = {3,3,3,3};                                                        
   float mitCuts_abstrkiso[4]           = {2.8,2.8,2.8,2.8};                                                
   float mitCuts_trkiso_hollow03[4]     = {4,4,4,4};                                                       
   float mitCuts_drtotk_25_99[4]	= {0.26,0.029,0.0062,0.0055};

   TLorentzVector phop4 = get_pho_p4( photon_index, vertex_index, pho_energy_array  );                      
   TLorentzVector phop4_badvtx = get_pho_p4( photon_index, pho_tkiso_badvtx_id[photon_index], pho_energy_array  );

   float rhofac=0.17;
   float val_hoe        = pho_hoe[photon_index];
   float val_sieie      = pho_sieie[photon_index];                                                          
   float val_ecaliso = pho_ecalsumetconedr03[photon_index] - 0.012*phop4.Et();                              
   float val_hcaliso = pho_hcalsumetconedr03[photon_index] - 0.005*phop4.Et(); 
                             
   float val_trkiso  = pho_trksumpthollowconedr03[photon_index] - 0.002*phop4.Et();                          
   float val_hcalecal   = (pho_ecalsumetconedr03[photon_index]+pho_hcalsumetconedr03[photon_index]-rho*rhofac);                                             
   float val_abstrkiso  = (*pho_tkiso_recvtx_030_002_0000_10_01)[photon_index][vertex_index];                
   float val_trkiso_hollow03 = pho_trksumpthollowconedr03[photon_index];                                    
//   float val_drtotk_25_99 = pho_drtotk_25_99[photon_index];
   int   val_pho_isconv = pho_isconv[photon_index];

   //if (run==173439 && lumis==155 && event==236661419) {

//	std::cout << "I will lose this photon, but MIT will not" << std::endl;
	
//	std::cout << "rho " << rho <<std::endl;
//	std::cout << "pho_n " << pho_n <<std::endl;
//	std::cout << "pho_index " << photon_index <<std::endl;
//	std::cout << "pho_et " << phop4.Et() <<std::endl;
//	std::cout << "hoe " << val_hoe <<std::endl;
//	std::cout << "sieie " << val_sieie <<std::endl;
//	std::cout << "ecaliso " << val_ecaliso <<std::endl;
//	std::cout << "hcaliso " << val_hcaliso <<std::endl;
//	std::cout << "hcalecal " << val_hcalecal <<std::endl;
//	std::cout << "abstrkiso " << val_abstrkiso <<std::endl;
//	std::cout << "trkiso " << val_trkiso_hollow03 <<std::endl;
//	std::cout << "isconv " << val_pho_isconv <<std::endl;

//	std::cout << "r9 " << pho_r9[photon_index] <<std::endl;
//	std::cout << "r9 categoty " << r9_category <<std::endl;
//	std::cout << "eta scxyz" << fabs(((TVector3*)sc_xyz->At(pho_scind[photon_index]))->Eta()) <<std::endl;
//	std::cout << "category" << photon_category <<std::endl;

//   }

   if (val_hoe             >= mitCuts_hoe[photon_category]         ) return false;                                           
   if (val_sieie           >= mitCuts_sieie[photon_category]       ) return false;
   if (val_ecaliso         >= mitCuts_ecaliso[photon_category]     ) return false;
   if (val_hcaliso         >= mitCuts_hcaliso[photon_category]     ) return false;                                           
   if (val_trkiso          >= mitCuts_trkiso[photon_category]      ) return false;
   if (val_hcalecal        >= mitCuts_hcalecal[photon_category]    ) return false;
   if (val_abstrkiso       >= mitCuts_abstrkiso[photon_category]   ) return false;                   
//   if (val_drtotk_25_99    <  mitCuts_drtotk_25_99[photon_category]   ) return false; // Electron Rejection based on CiC for now
   if (!val_pho_isconv						   ) return false; // Electron Rejection based Conversion Safe Veto
   if (val_trkiso_hollow03 >= mitCuts_trkiso_hollow03[photon_category]) return false;                                        

   return true;

}

// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::PhotonCiCSelectionLevel( int photon_index, int vertex_index, std::vector<std::vector<bool> > & ph_passcut, int ncategories, 
				      int doSublead, float *pho_energy_array ) {

  if( ! runCiC ) {
    switch(ncategories) {
    case (4):
      return doSublead ? (*pho_cic4cutlevel_sublead)[photon_index][vertex_index] : (*pho_cic4cutlevel_lead)[photon_index][vertex_index] ;
    case (6):
      return doSublead ? (*pho_cic6cutlevel_sublead)[photon_index][vertex_index] : (*pho_cic6cutlevel_lead)[photon_index][vertex_index] ;
    }
  }
  int cutlevelpassed = -1;

  int n_r9_categories = -1;
  int n_eta_categories = -1;
  if(ncategories==6) {
    n_r9_categories = 3;
    n_eta_categories = 2;
  } else if(ncategories==4) {
    n_r9_categories = 2;
    n_eta_categories = 2;
  } else {
    std::cout << "UNKNOWN ncategories must be 4 or 6, not " << ncategories << std::endl;
  }
  int photon_category = PhotonCategory(photon_index,n_r9_categories,n_eta_categories);

  TLorentzVector phop4 = get_pho_p4( photon_index, vertex_index, pho_energy_array  );
  TLorentzVector phop4_badvtx = get_pho_p4( photon_index, pho_tkiso_badvtx_id[photon_index], pho_energy_array  );

  float val_tkiso = (*pho_tkiso_recvtx_030_002_0000_10_01)[photon_index][vertex_index];
  float val_ecaliso = pho_ecalsumetconedr03[photon_index];
  float val_hcaliso = pho_hcalsumetconedr04[photon_index];
  float val_ecalisobad = pho_ecalsumetconedr04[photon_index];
  float val_hcalisobad = pho_hcalsumetconedr04[photon_index];
  float val_tkisobad = pho_tkiso_badvtx_040_002_0000_10_01[photon_index];
  float val_sieie = pho_sieie[photon_index];
  float val_hoe = pho_hoe[photon_index];
  float val_r9 = pho_r9[photon_index];
  float val_drtotk_25_99 = pho_drtotk_25_99[photon_index];
  float val_pixel = (float)pho_haspixseed[photon_index];

  float isosumconst = 5.;
  float isosumconstbad = 7.;
  if(ncategories==4) {
    isosumconst = 0.;
    isosumconstbad = 0.;
  }
  
  /// float rhofacbad=0.40, rhofac=0.05;
  float rhofacbad=0.52, rhofac=0.17;
  float val_isosumoet=(val_tkiso+val_ecaliso+val_hcaliso+isosumconst-rho*rhofac)*50./phop4.Et();
  float val_isosumoetbad=(val_tkisobad+val_ecalisobad+val_hcalisobad+isosumconstbad-rho*rhofacbad)*50./phop4_badvtx.Et();
  float val_trkisooet=(val_tkiso)*50./phop4.Et();

  ph_passcut.clear();
  ph_passcut.resize(phoNCUTLEVELS,std::vector<bool>(8,true) );
  if(!doSublead) {
    for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {
      switch(ncategories) {
        case(6) :
          ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic6_cut_lead_isosumoet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic6_cut_lead_isosumoetbad[iCUTLEVEL][photon_category]  );
          ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic6_cut_lead_trkisooet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic6_cut_lead_sieie[iCUTLEVEL][photon_category]         );
          ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic6_cut_lead_hovere[iCUTLEVEL][photon_category]        );
          ph_passcut[iCUTLEVEL][5] = (val_r9             >=     cic6_cut_lead_r9[iCUTLEVEL][photon_category]            );// gt cut
          ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic6_cut_lead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
          ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic6_cut_lead_pixel[iCUTLEVEL][photon_category]         );
          break;
        case(4) :
          ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4_cut_lead_isosumoet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4_cut_lead_isosumoetbad[iCUTLEVEL][photon_category]  );
          ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4_cut_lead_trkisooet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4_cut_lead_sieie[iCUTLEVEL][photon_category]         );
          ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4_cut_lead_hovere[iCUTLEVEL][photon_category]        );
          ph_passcut[iCUTLEVEL][5] = (val_r9             >=     cic4_cut_lead_r9[iCUTLEVEL][photon_category]            );// gt cut
          ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic4_cut_lead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
          ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic4_cut_lead_pixel[iCUTLEVEL][photon_category]         );
          break;
      }
      bool ph_passcut_all = true;
      for(int icut=0;icut!=8;++icut) {
	ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
      }
      if(ph_passcut_all) {
	if( cutlevelpassed != iCUTLEVEL - 1 ) {
	  std::cerr << "photon " << photon_index << " (category " << photon_category << ") in run/event " << run << "/" << event << " passed CiC cut level " 
		    << iCUTLEVEL << " but not "  << iCUTLEVEL - 1 << ". Did you load your cut values correctly? "<< std::endl;
	  /// assert( 0 );
	}
	cutlevelpassed=iCUTLEVEL;
      }
    }
  } else if(doSublead) {
    for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {
      switch(ncategories) {
        case(6) :
          ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic6_cut_sublead_isosumoet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic6_cut_sublead_isosumoetbad[iCUTLEVEL][photon_category]  );
          ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic6_cut_sublead_trkisooet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic6_cut_sublead_sieie[iCUTLEVEL][photon_category]         );
          ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic6_cut_sublead_hovere[iCUTLEVEL][photon_category]        );
          ph_passcut[iCUTLEVEL][5] = (val_r9             >=     cic6_cut_sublead_r9[iCUTLEVEL][photon_category]            );// gt cut
          ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic6_cut_sublead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
          ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic6_cut_sublead_pixel[iCUTLEVEL][photon_category]         );
          break;
        case(4) :
          ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4_cut_sublead_isosumoet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4_cut_sublead_isosumoetbad[iCUTLEVEL][photon_category]  );
          ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4_cut_sublead_trkisooet[iCUTLEVEL][photon_category]     );
          ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4_cut_sublead_sieie[iCUTLEVEL][photon_category]         );
          ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4_cut_sublead_hovere[iCUTLEVEL][photon_category]        );
          ph_passcut[iCUTLEVEL][5] = (val_r9             >=     cic4_cut_sublead_r9[iCUTLEVEL][photon_category]            );// gt cut
          ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic4_cut_sublead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
          ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic4_cut_sublead_pixel[iCUTLEVEL][photon_category]         );
          break;
      }
      bool ph_passcut_all = true;
      for(int icut=0;icut!=8;++icut) {
	ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
      }
      if(ph_passcut_all) {
	if( cutlevelpassed != iCUTLEVEL - 1 ) {
	  std::cerr << "photon " << photon_index << " (category " << photon_category << ") in run/event " << run << "/" << event << " passed CiC cut level " 
		    << iCUTLEVEL << " but not "  << iCUTLEVEL - 1 << ". Did you load your cut values correctly? " << std::endl;
	  //// assert( 0 );
	}
	cutlevelpassed=iCUTLEVEL;
      }
    }
  }

  return cutlevelpassed;

}


// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::DeltaRToTrack(Int_t photonind, Int_t vtxind, Float_t PtMin, Float_t dzmax, Float_t dxymax, int maxlosthits){
  if(LDEBUG)std::cout << "DeltaRToTrack BEGIN" << std::endl;
  int elind = -1;
  float eldr = 99.;
  for(int iel=0;iel!=el_std_n;++iel) {
    if(el_std_hp_expin[iel]>maxlosthits // || ((TLorentzVector*)el_std_p4->At(iel))->Pt()<PtMin 
       ) continue;
    if(el_std_scind[iel] == pho_scind[photonind]) {
      elind = iel;
      break;
    }
  }
  if(elind >= 0)
    eldr = sqrt(el_std_detain[elind]*el_std_detain[elind] + el_std_dphiin[elind]*el_std_dphiin[elind]);
  return eldr;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::IsoEcalHitsSumEtNumCrystal( TVector3 *calopos, Float_t innerConeDR, Float_t outerConeDR, Float_t stripEtaHalfWidth, Float_t stripHalfLength) {
  if(LDEBUG)std::cout << "IsoEcalHitsSumEtNumCrystal begin" << std::endl;
  Float_t ecalhitsSumEt=0.;
  for(int i=0; i!= ecalhit_n; ++i) {
    //if(ecalhit_type[i]==2)continue;// don't use preshower (still want this?? - crashes without it)
    TLorentzVector * rechitp4 = (TLorentzVector *) ecalhit_p4->At(i);
    Float_t etaclus = calopos->Eta();
    Float_t phiclus = calopos->Phi();
    Float_t eta = rechitp4->Eta();
    Float_t phi = rechitp4->Phi();
    Float_t etaDiff = eta - etaclus;
    Float_t phiDiff = phi - phiclus;
    if(phiDiff>TMath::Pi())phiDiff = TMath::TwoPi() - phiDiff;
    Float_t deltaR = calopos->DeltaR(rechitp4->Vect());
    Float_t deltaEta = fabs(rechitp4->Eta() - calopos->Eta());
    if(sqrt(etaDiff*etaDiff + phiDiff*phiDiff)>outerConeDR)continue;
    if(fabs(etaclus) < 1.479) {
      if(fabs(etaDiff) < 0.0174*stripEtaHalfWidth)continue;
      if(sqrt(etaDiff*etaDiff + phiDiff*phiDiff) < 0.0174*innerConeDR)continue;
    } else {
      if(fabs(etaDiff) < 0.00864*fabs(sinh(eta))*stripEtaHalfWidth)continue;
      if(sqrt(etaDiff*etaDiff + phiDiff*phiDiff) < 0.00864*fabs(sinh(eta))*innerConeDR)continue;
    }
    if(fabs(rechitp4->E())>0.08 && fabs(rechitp4->Et()) > 0.)ecalhitsSumEt+=rechitp4->Et();
  }
  if(LDEBUG)std::cout << "IsoEcalHitsSumEtNumCrystal end" << std::endl;
  return ecalhitsSumEt;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
std::pair<Int_t, Float_t> LoopAll::WorstSumTrackPtInCone(Int_t ipho, Int_t returnvtxind, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {

  Int_t worstvtxind = -1;
  Float_t maxisosum = -100;
  for(int ivtx=0;ivtx!=vtx_std_n;++ivtx) {
    TLorentzVector photon_p4 = get_pho_p4( ipho, ivtx );
    Float_t thisvtxisosum = SumTrackPtInCone(&photon_p4, ivtx, PtMin, OuterConeRadius, InnerConeRadius, EtaStripHalfWidth, dzmax, dxymax);
    if(thisvtxisosum > maxisosum) {
      maxisosum = thisvtxisosum;
      worstvtxind = ivtx;
    }
  }

  /// if(returnvtxind == 1) {
  ///   return 0.5+(float)worstvtxind;
  /// } else {
  ///   return maxisosum;
  /// }
  return std::make_pair(worstvtxind,maxisosum); 
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::SumTrackPtInCone(TLorentzVector *photon_p4, Int_t vtxind, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {
  // TRACKER Isolation
  if(vtxind<0)return -99;
  TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(vtxind);
  float SumTrackPt=0;
  for(unsigned int itk=0; itk!=tk_n; itk++) {
    TLorentzVector * tkp4= (TLorentzVector *) tk_p4->At(itk);
    if(tkp4->Pt() < PtMin)continue;
    TVector3 * tkpos= (TVector3 *) tk_vtx_pos->At(itk);
    /// double deltaz = fabs(vtxpos->Z() - tkpos->Z()); 
    double deltaz = fabs( (tkpos->Z()-vtxpos->Z()) - ( (tkpos->X()-vtxpos->X())*tkp4->Px() + (tkpos->Y()-vtxpos->Y())*tkp4->Py() )/tkp4->Pt() * tkp4->Pz()/tkp4->Pt() );
    if(deltaz > dzmax)continue;
    double dxy = ( -(tkpos->X() - vtxpos->X())*tkp4->Py() + (tkpos->Y() - vtxpos->Y())*tkp4->Px()) / tkp4->Pt();
    if(fabs(dxy) > dxymax)continue;
    double tk_eta = tkp4->Eta();
    double tk_phi = tkp4->Phi();
    double deta = fabs(photon_p4->Eta() - tk_eta);
    double dphi = fabs(photon_p4->Phi() - tk_phi);
    if(dphi > TMath::Pi())dphi = TMath::TwoPi() - dphi;
    double deltaR = sqrt(deta*deta + dphi*dphi);
    if(deltaR < OuterConeRadius && deltaR >= InnerConeRadius && deta >= EtaStripHalfWidth)SumTrackPt+=tkp4->Pt();
  }
  return SumTrackPt;
}


bool LoopAll::CheckSphericalPhoton(int phoid){

  TVector3 *phoCalo = (TVector3*)sc_xyz->At(pho_scind[phoid]);
  if (pho_r9[phoid]<0.94 || fabs(phoCalo->Eta())>1.) return false;

  TLorentzVector *bcpos   = (TLorentzVector*)bc_p4->At(sc_bcseedind[pho_scind[phoid]]);
  //TVector3 bcxyz = bcpos->Vect();

  double minDR=999.;
  int closestHit=-1;
  for (int i=0;i<ecalhit_n;i++){
	TLorentzVector *xtalpos = (TLorentzVector*)ecalhit_p4->At(i);
  	//TVector3 xtalxyz = xtalpos->Vect();
	//double dR = (xtalxyz-bcxyz).Mag();
	double dR = xtalpos->DeltaR(*bcpos);
	if(dR<minDR){
	 	closestHit = i;
		minDR = dR;
	}
  }
  if (closestHit<0) std::cout << "Fishy !!!!!!" <<std::endl;
 
  int detid = ecalhit_detid[closestHit];
  //int detid = ecalhit_detid[bc_seed[sc_bcseedind[pho_scind[phoid]]]];
  int ieta  = (detid>>9)&0x7F; 
  int iphi  = detid&0x1FF; 

 //int ieta=IEta(bcpos->Eta());
 //int iphi=IPhi(bcpos->Phi());
//int ieta = IEta( ((TLorentzVector*)ecalhit_p4->At(closestHit))->Eta()); 
//int iphi = IPhi( ((TLorentzVector*)ecalhit_p4->At(closestHit))->Phi()); 

  if ((iphi %20)<=5 || (iphi%20)>=16){
   return false;
  }

  int ietaTT=(std::abs(ieta)-1)/5+1;
  if
    (
     (ietaTT>= 2&&     ietaTT<    5 ) ||
     (ietaTT>= 7&&     ietaTT<    9 ) ||
     (ietaTT>= 11&&    ietaTT<    13) ||
     (ietaTT>= 15&&    ietaTT<    17)
    ){	
      return true; 
   } 


  return false;

}
// CiC SELECTION CODE END - SSIMON

// Functions moved from Tools.h
// ---------------------------------------------------------------------------------------------------------------------------------------------
double LoopAll::DeltaPhi(double phi1, double phi2) {
  double deltaphi;
  if(phi1<0) phi1+=TWOPI;
  if(phi2<0) phi2+=TWOPI;
  deltaphi=fabs(phi1-phi2);
  if(deltaphi>TWOPI) deltaphi-=TWOPI;
  if(deltaphi>PI) deltaphi=TWOPI-deltaphi;
  return deltaphi;
}
//
// Generate dictionary entries for branches from GeneralFunctions_h 
//
// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::DefineUserBranches() 
{
  runCiC = true;
#ifndef __CINT__
  BRANCH_DICT(gh_gen2reco1);
  BRANCH_DICT(gh_gen2reco2);
  BRANCH_DICT(gh_vbfq1_pdgid);
  BRANCH_DICT(gh_vbfq2_pdgid);
  BRANCH_DICT(gh_vh_pdgid);
  BRANCH_DICT(gh_vh1_pdgid);
  BRANCH_DICT(gh_vh2_pdgid);
  BRANCH_DICT(gh_higgs_p4);
  BRANCH_DICT(gh_pho1_p4);
  BRANCH_DICT(gh_pho2_p4);
  BRANCH_DICT(gh_vbfq1_p4);
  BRANCH_DICT(gh_vbfq2_p4);
  BRANCH_DICT(gh_vh1_p4);
  BRANCH_DICT(gh_vh2_p4);

	BRANCH_DICT(gv_n  );
	BRANCH_DICT(gv_pos);
	BRANCH_DICT(pu_n);
	BRANCH_DICT(pu_zpos);
	BRANCH_DICT(pu_sumpt_lowpt);
	BRANCH_DICT(pu_sumpt_highpt);
	BRANCH_DICT(pu_ntrks_lowpt);
	BRANCH_DICT(pu_ntrks_highpt);
	
	BRANCH_DICT(rho);
	
	BRANCH_DICT(vtx_std_sel);
	BRANCH_DICT(vtx_std_ranked_list);
	BRANCH_DICT(vtx_std_evt_mva);

	BRANCH_DICT(pho_tkiso_recvtx_030_002_0000_10_01);
	BRANCH_DICT(pho_tkiso_badvtx_040_002_0000_10_01);
	BRANCH_DICT(pho_tkiso_badvtx_id);
	BRANCH_DICT(pho_drtotk_25_99);

	BRANCH_DICT(pho_cic6cutlevel_lead);
	BRANCH_DICT(pho_cic6passcuts_lead);
	BRANCH_DICT(pho_cic6cutlevel_sublead);
	BRANCH_DICT(pho_cic6passcuts_sublead);

	BRANCH_DICT(pho_cic4cutlevel_lead);
	BRANCH_DICT(pho_cic4passcuts_lead);
	BRANCH_DICT(pho_cic4cutlevel_sublead);
	BRANCH_DICT(pho_cic4passcuts_sublead);

	BRANCH_DICT(pho_cutlevel_lead);
	BRANCH_DICT(pho_passcuts_lead);
	BRANCH_DICT(pho_cutlevel_sublead);
	BRANCH_DICT(pho_passcuts_sublead);

	BRANCH_DICT(pho_matchingConv);
	
	BRANCH_DICT(dipho_n);
	BRANCH_DICT(dipho_leadind);
	BRANCH_DICT(dipho_subleadind);
	BRANCH_DICT(dipho_vtxind);
	BRANCH_DICT(dipho_sumpt);
	BRANCH_DICT(pho_genmatched);
	BRANCH_DICT(pho_regr_energy_otf);
	BRANCH_DICT(pho_regr_energyerr_otf);
	//// BRANCH_DICT(dipho_leadet);
	//// BRANCH_DICT(dipho_subleadet);
	//// BRANCH_DICT(dipho_leadeta);
	//// BRANCH_DICT(dipho_subleadeta);
	//// BRANCH_DICT(dipho_leadci6cindex);
	//// BRANCH_DICT(dipho_subleadci6cindex);
	//// BRANCH_DICT(dipho_leadci4cindex);
	//// BRANCH_DICT(dipho_subleadci4cindex);
	//// BRANCH_DICT(dipho_mass);
	//// BRANCH_DICT(dipho_pt);
	//// BRANCH_DICT(dipho_eta);
	//// BRANCH_DICT(dipho_phi);
	//// BRANCH_DICT(dipho_cts);

	BRANCH_DICT(pho_matchingConv);

#endif
}

int  LoopAll::RescaleJetEnergy() {
  for (int i = 0; i<jet_algoPF1_n; i++) {
    TLorentzVector * thisjet = (TLorentzVector *) jet_algoPF1_p4->At(i);
    *thisjet*=jet_algoPF1_erescale[i];
  }
  return 1;
}

std::pair<int, int> LoopAll::Select2HighestPtJets(TLorentzVector& leadpho, TLorentzVector& subleadpho, float jtLMinPt, float jtTMinPt){

  std::pair<int, int> myJets(-1,-1);
  std::pair<int, int> fail(-1,-1);

  std::pair<float, float> myJetspt(-1.,-1.);

  float dr2pho = 0.5;
  float dr2jet = 0.5;

  TLorentzVector* j1p4;
  TLorentzVector* j2p4;
  float j1pt=-1;
  float j2pt=-1;

  // select highest pt jets

  for(int j1_i=0; j1_i<jet_algoPF1_n; j1_i++){
    j1p4 = (TLorentzVector*) jet_algoPF1_p4->At(j1_i);
    if(fabs(j1p4->Eta()) > 4.7) continue;
    if(j1p4->DeltaR(leadpho) < dr2pho) continue;
    if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
    j1pt=j1p4->Pt();


    if(j1pt>myJetspt.first) {
      myJets.second=myJets.first;
      myJetspt.second=myJetspt.first;
      myJetspt.first=j1pt;
      myJets.first=j1_i;
    }
    else if(j1pt>myJetspt.second) {
      myJetspt.second=j1pt;
      myJets.second=j1_i;
    }
  }

  if(myJets.second!=-1&&myJetspt.first>jtTMinPt&&myJetspt.second>jtTMinPt) {
    if(myJets!=myJets) {
      j1p4 = (TLorentzVector*) jet_algoPF1_p4->At(myJets.first);
      j2p4 = (TLorentzVector*) jet_algoPF1_p4->At(myJets.second);
      float dr=j2p4->DeltaR(*j1p4);
    }
  }

  return myJets;
}

int LoopAll::MuonSelection(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind){
  int mymu = -1;

  if(run == 170249 && (lumis>= 37 && lumis<= 191 )) return mymu;
  if(run == 170249 && (lumis>= 191 && lumis<= 507 )) return mymu;
  if(run == 170255 && (lumis>= 1 && lumis<= 387 )) return mymu;
  if(run == 170286 && (lumis>= 77 && lumis<= 260 )) return mymu;
  if(run == 170292 && (lumis>= 1 && lumis<= 258 )) return mymu;
  if(run == 170298 && (lumis>= 1 && lumis<= 178 )) return mymu;
  if(run == 170354 && (lumis>= 1 && lumis<= 308)) return mymu;
  if(run == 170397 && (lumis>= 1 && lumis<= 345)) return mymu;
  if(run == 170406 && (lumis>= 1 && lumis<= 171)) return mymu;
  if(run == 170452 && (lumis>= 72 && lumis<= 110)) return mymu;
  if(run == 170527 && (lumis>= 49 && lumis<= 92)) return mymu;
  
  TLorentzVector* thismu;
  float thiseta = -100;
  float thispt = -100;
  float thisiso =1000;

  int passingMu = 0;

  for( int indmu=0; indmu<mu_glo_n; indmu++){
    thismu = (TLorentzVector*) mu_glo_p4->At(indmu);
    thiseta = fabs(thismu->Eta());
    if(thiseta>2.4) continue;
    thispt = thismu->Pt();
    if(thispt<20) continue;
    if(mu_glo_type[indmu]<1100) continue;  // global and tracker
    if(mu_glo_chi2[indmu]/mu_glo_dof[indmu]>=10) continue;
    if(mu_glo_pixelhits[indmu]<=0) continue;
    if(mu_glo_validhits[indmu]<=10) continue;
    if(mu_glo_validChmbhits[indmu]<=0) continue;
    if(mu_glo_nmatches[indmu]<=1) continue;
    thisiso=mu_glo_ecaliso03[indmu]+mu_glo_hcaliso03[indmu]+mu_glo_tkiso03[indmu] - rho*3.1415926*0.09;
    if(thisiso/thispt>=0.1) continue;
    if(std::min(pho1.DeltaR(*thismu),pho2.DeltaR(*thismu)) < 1) continue;
    
    // need to calculate d0, dz wrt chosen vtx
    if(fabs(mu_glo_D0Vtx[indmu][vtxind]) > 0.02) continue;
    if(fabs(mu_glo_DZVtx[indmu][vtxind]) > 0.1)  continue;


    passingMu++;

    std::cout << setprecision(4) << "Run = " << run << "  LS = " << lumis << "  Event = " << event << "  SelVtx = " << vtxind << " muEta = " << thiseta << "  muPhi = " << thismu->Phi() <<  "  muEt = " << thismu->Et() << endl;

    mymu = indmu;
  }

  if(passingMu>1) std::cout<<"There are "<<passingMu<<" passing muons!!"<<std::endl;


  /////////////////fabs(muD0Vtx[i][vtx]) < 0.02 &&       
  /////////////////fabs(muDzVtx[i][vtx]) < 0.1 
  /////////////////(here D0 and DZ are wrt vertex selected by the mva vertexing)

  return mymu;
}


int LoopAll::ElectronSelection(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind){
  int myel = -1;

  if(run == 170249 && (lumis>= 37 && lumis<= 191 ))return myel;
  if(run == 170249 && (lumis>= 191 && lumis<= 507 )) return myel;
  if(run == 170255 && (lumis>= 1 && lumis<= 387 )) return myel;
  if(run == 170286 && (lumis>= 77 && lumis<= 260 )) return myel;
  if(run == 170292 && (lumis>= 1 && lumis<= 258 )) return myel;
  if(run == 170298 && (lumis>= 1 && lumis<= 178 )) return myel;
  if(run == 170354 && (lumis>= 1 && lumis<= 308)) return myel;
  if(run == 170397 && (lumis>= 1 && lumis<= 345)) return myel;
  if(run == 170406 && (lumis>= 1 && lumis<= 171)) return myel;
  if(run == 170452 && (lumis>= 72 && lumis<= 110)) return myel;
  if(run == 170527 && (lumis>= 49 && lumis<= 92)) return myel;
  
  TLorentzVector* thisel;
  TLorentzVector* thissc;
  float thiseta = -100;
  float thispt = -100;
  float thisiso =1000;

  int passingEl = 0;

  for( int indel=0; indel<el_std_n; indel++){
    if(el_std_hp_expin[indel]!=0) continue;

    thisel = (TLorentzVector*) el_std_p4->At(indel);
    thissc = (TLorentzVector*) el_std_sc->At(indel);
    thiseta = fabs(thissc->Eta());
    if(thiseta>2.5 || (thiseta>1.442 && thiseta<1.566)) continue;
    thispt = thisel->Pt();
    if(thispt<20) continue;
    if(thiseta<1.442) {   // EB cuts
      if(el_std_sieie[indel]>=0.01) continue; 
      if(fabs(el_std_dphiin[indel])>=0.039) continue;
      if(fabs(el_std_detain[indel])>=0.005) continue;
      thisiso = el_std_tkiso03[indel] + std::max(0.,(double)el_std_ecaliso03[indel]-1.)
              + el_std_hcaliso03[indel] + (el_std_hoe[indel]*thissc->Energy()*sin(thissc->Theta())) - rho*3.1415926*0.09;
      if(thisiso/thispt>=0.053) continue; 
    } else {  // EE cuts
      if(el_std_sieie[indel]>=0.03) continue; 
      if(fabs(el_std_dphiin[indel])>=0.028) continue;
      if(fabs(el_std_detain[indel])>=0.007) continue;
      thisiso = el_std_tkiso03[indel] + el_std_ecaliso03[indel]
              + el_std_hcaliso03[indel] + (el_std_hoe[indel]*thissc->Energy()*sin(thissc->Theta())) 
              - rho*3.1415926*0.09;
      if(thisiso/thispt>=0.042) continue; 
    }

    // conversion rejection
    if(fabs(el_std_dcot[indel])<0.02 && fabs(el_std_dist[indel])<0.02) continue;

    if(std::min( pho1.DeltaR(*thisel), pho2.DeltaR(*thisel))<=1.) continue;

    TLorentzVector elpho1 = *thisel + pho1;
    if( fabs(elpho1.M() - 91.19) <= 5) continue;

    TLorentzVector elpho2 = *thisel + pho2;
    if( fabs(elpho2.M() - 91.19) <= 5) continue;

    // need to calculate d0, dz wrt chosen vtx
    if(fabs(el_std_D0Vtx[indel][vtxind]) > 0.02) continue;
    if(fabs(el_std_DZVtx[indel][vtxind]) > 0.1)  continue;

    passingEl++;

    std::cout << setprecision(4) << "Run = " << run << "  LS = " << lumis << "  Event = " << event << "  SelVtx = " << vtxind << " elEta = " << thiseta << "  elPhi = " << thisel->Phi() <<  "  elEt = " << thisel->Et() << endl;

    myel = indel;
  }

  if(passingEl>1) std::cout<<"There are "<<passingEl<<" passing electrons!!"<<std::endl;


  /////////////////fabs(elD0Vtx[i][vtx]) < 0.02 &&       
  /////////////////fabs(elDzVtx[i][vtx]) < 0.1 
  /////////////////(here D0 and DZ are wrt vertex selected by the mva vertexing)

  return myel;
}


//--- RECO-MC PHOTONS MATCHING --------------------------------------------------------------------------------------------------------
bool LoopAll::FindMCHiggsPhotons(int& higgsind, int& mc1, int& mc2, int& i1, int& i2 )
{
  bool is_mcmatched = false;
  
  int higgsstat2 = -1;
  int higgsstat3 = -1;

  int pho1stat3 = -1;
  int pho2stat3 = -1;
  
  for ( int i = 0; i< gp_n ; i++ ){
    int pid        = gp_pdgid[i];
    int status     = gp_status[i];
    
    if (pid == 25) {
      if(status==2) higgsstat2=i;
      if(status==3) higgsstat3=higgsind=i;
      continue;
    }
   
    
    if (pid == 22 && status==3 && gp_mother[i]==higgsstat3) {
      if(pho1stat3==-1) pho1stat3=i;
      else if(pho2stat3==-1) pho2stat3=i;
      continue;
    }

    if (pid != 22 || status!=1) continue;
   
    if ( (gp_mother[i]==pho1stat3 || gp_mother[i]==pho2stat3) && mc1 < 0 ) { mc1  = i; }
    else if ( (gp_mother[i]==pho1stat3 || gp_mother[i]==pho2stat3) && mc2 < 0 ) { mc2  = i; break; }
  
  }
  
  if(higgsstat2 == -1) return is_mcmatched;

  TLorentzVector* higgs2 = (TLorentzVector*)gp_p4->At(higgsstat2);
  TLorentzVector* higgs3 = (TLorentzVector*)gp_p4->At(higgsstat3);
  

  if (mc1 < 0 || mc2 < 0 ) return is_mcmatched;
  
  TLorentzVector *mcV1 = (TLorentzVector*)gp_p4->At(mc1);  
  TLorentzVector *mcV2 = (TLorentzVector*)gp_p4->At(mc2);  
 
  TLorentzVector higgs_gg = *mcV1 + *mcV2;


  int index1 = -100, index2 = -100;
  double dr1min = 0.3;
  double dr2min = 0.3;
  
  for (int i = 0; i < pho_n; i++){
  
    TVector3 *scpos = (TVector3*) pho_calopos -> At(i);
    
    double dr1 = scpos->DeltaR(mcV1->Vect());
    if (dr1 < dr1min) { dr1min = dr1; index1 = i; }
    
    double dr2 = scpos->DeltaR(mcV2->Vect());
    if (dr2 < dr2min) { dr2min = dr2; index2 = i; }
  }
  
  if(index1!=-100 && index2!=-100){
    TLorentzVector *photon1 = (TLorentzVector*)pho_p4->At(index1);
    TLorentzVector *photon2 = (TLorentzVector*)pho_p4->At(index2);


    // order photons by pt
    if (photon1->Pt() > photon2->Pt()) {
      i1 = index1;
      i2 = index2;
    } else {
      i1 = index2;
      i2 = index1;
    }
    
    is_mcmatched   = (dr1min < 0.15) && (dr2min< 0.15);
  } else if(index1!=-100 || index2!=-100) {
    i1 = std::max(index1,index2);
  }
    

  return is_mcmatched;
}


//--- GEN-VBF MATCHING --------------------------------------------------------------------------------------------------------
bool LoopAll::FindMCVBF(int higgsind, int& vbfq1, int& vbfq2 )
{
  bool matchfound = false;
  int higgsmom = gp_mother[higgsind];

  for ( int i=higgsind+1; i< gp_n ; i++ ){
    if(gp_status[i]!=3) continue;
    //std::cout<<"vbfq1 vbfq2 i gp_mother[i] gp_pdgid[i] gp_status[i]  "<<vbfq1<<"   "<<vbfq2<<"   "<<i<<"   "<<gp_mother[i]<<"   "<<gp_pdgid[i]<<"   "<<gp_status[i]<<std::endl;
    if(higgsmom==gp_mother[i]&&gp_pdgid[i]!=21&&abs(gp_pdgid[i])!=23&&abs(gp_pdgid[i])!=24){
      if(vbfq1==-100) vbfq1=i;
      else if(vbfq2==-100) {
        vbfq2=i;
        matchfound=true; 
        break;
      }
    }
  }
  

  return matchfound;
}


//--- GEN-VH MATCHING --------------------------------------------------------------------------------------------------------
bool LoopAll::FindMCVH(int higgsind, int& vh, int& vh1, int& vh2 )
{
  bool matchfound = false;
  int higgsmom = gp_mother[higgsind];

  for ( int i=higgsind-1; i< gp_n ; i++ ){
    if(gp_status[i]!=3) continue;
    if(higgsmom==gp_mother[i] && i!=higgsind && vh==-100) {
      if(abs(gp_pdgid[i])==23 || abs(gp_pdgid[i])==24) vh=i;
      continue;
    }

    //std::cout<<"vh vh1 vh2 i gp_mother[i] gp_pdgid[i] gp_status[i]  "<<vh<<"   "<<vh1<<"   "<<vh2<<"   "<<i<<"   "<<gp_mother[i]<<"   "<<gp_pdgid[i]<<"   "<<gp_status[i]<<std::endl;

    if(gp_mother[i]==vh && vh1==-100) vh1=i;
    else if(gp_mother[i]==vh && vh2==-100) { 
      vh2=i; 
      matchfound=true; 
      break; 
    }
  }
  
  return matchfound;
}


#ifdef NewFeatures
#include "Marco/plotInteractive_cc.h"
#endif
