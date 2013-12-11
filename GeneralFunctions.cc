#include "LoopAll.h"
#include "Sorters.h"
#include "TRandom3.h"
#define GFDEBUG 0

float LoopAll::pfTkIsoWithVertex(int phoindex, int vtxInd, float dRmax, float dRvetoBarrel, float dRvetoEndcap, 
                                 float ptMin, float dzMax, float dxyMax, int pfToUse) {
  
    float dRveto;
    if (pho_isEB[phoindex])
        dRveto = dRvetoBarrel;
    else
        dRveto = dRvetoEndcap;
  
    TLorentzVector photonDirectionWrtVtx = get_pho_p4(phoindex, vtxInd, 0);
  
    float sum = 0;
    // Loop over the PFCandidates
    for(unsigned i=0; i<pfcand_n; i++) {
    
        //require that PFCandidate is a charged hadron
        if (pfcand_pdgid[i] == pfToUse) {
      
            TLorentzVector* pfc = (TLorentzVector*)pfcand_p4->At(i);
      
            if (pfc->Pt() < ptMin)
                continue;
    
            TVector3* vtx = (TVector3*)vtx_std_xyz->At(vtxInd);
            TVector3* pfCandVtx = (TVector3*)pfcand_posvtx->At(i);

            float dz = fabs(pfCandVtx->Z() - vtx->Z());
      
            if (dz > dzMax) 
                continue;

            double dxy = (-(pfCandVtx->X() - vtx->X())*pfc->Py() + (pfCandVtx->Y() - vtx->Y())*pfc->Px()) / pfc->Pt();
            if(fabs(dxy) > dxyMax) 
                continue;
      
            float dR = photonDirectionWrtVtx.DeltaR(*pfc);
            if(dR > dRmax || dR < dRveto) 
                continue;
      
            sum += pfc->Pt();
        }
    }
  
    return sum;
}

float LoopAll::pfEcalIso(int phoindex, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStripBarrel, 
                         float etaStripEndcap, float thrBarrel, float thrEndcaps, int pfToUse) {
  
    float dRVeto, etaStrip, thr;
    if (pho_isEB[phoindex]) {
        dRVeto = dRVetoBarrel;
        etaStrip = etaStripBarrel;
        thr = thrBarrel;
    } else {
        dRVeto = dRVetoEndcap;
        etaStrip = etaStripEndcap;
        thr = thrEndcaps;
    }

    float sum = 0;
    for(unsigned i=0; i<pfcand_n; i++) {
    
        if (pfcand_pdgid[i] == pfToUse) {
      
            // FIXME questo non so come implementarlo...
            //if(pfc.superClusterRef().isNonnull() && localPho->superCluster().isNonnull()) {
            //if (pfc.superClusterRef() == localPho->superCluster()) 
            //  continue;
            //}
      
            TVector3* pfvtx = (TVector3*)pfcand_posvtx->At(i);
            TVector3* phoEcalPos = (TVector3*)sc_xyz->At(pho_scind[phoindex]);

            TVector3 photonDirectionWrtVtx = TVector3(phoEcalPos->X() - pfvtx->X(),
                                                      phoEcalPos->Y() - pfvtx->Y(),
                                                      phoEcalPos->Z() - pfvtx->Z());

            TLorentzVector* pfc = (TLorentzVector*)pfcand_p4->At(i);

            if( pfc->Pt() < thr ) 
                continue;

            float dEta = fabs(photonDirectionWrtVtx.Eta() - pfc->Eta());
            float dR = photonDirectionWrtVtx.DeltaR(pfc->Vect());
      
            if (dEta < etaStrip)
                continue;
      
            if(dR > dRmax || dR < dRVeto)
                continue;
      
            sum += pfc->Pt();
        }
    }
  
    return sum;
}



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
  tmvaReaderID_MIT_Barrel->AddVariable("ScEta",        &tmva_id_mit_sceta);
  tmvaReaderID_MIT_Barrel->AddVariable("EtaWidth",    &tmva_id_mit_etawidth);
  tmvaReaderID_MIT_Barrel->AddVariable("PhiWidth",    &tmva_id_mit_phiwidth);
  
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
  tmvaReaderID_MIT_Endcap->AddVariable("ScEta",        &tmva_id_mit_sceta);
  tmvaReaderID_MIT_Endcap->AddVariable("EtaWidth",    &tmva_id_mit_etawidth);
  tmvaReaderID_MIT_Endcap->AddVariable("PhiWidth",    &tmva_id_mit_phiwidth);
  
  tmva_dipho_MIT_buf.resize(10,0.);
  tmva_dipho_MIT_dmom = &tmva_dipho_MIT_buf[0];
  tmva_dipho_MIT_dmom_wrong_vtx= &tmva_dipho_MIT_buf[1];
  tmva_dipho_MIT_vtxprob= &tmva_dipho_MIT_buf[2];
  tmva_dipho_MIT_ptom1= &tmva_dipho_MIT_buf[3];
  tmva_dipho_MIT_ptom2= &tmva_dipho_MIT_buf[4];
  tmva_dipho_MIT_eta1= &tmva_dipho_MIT_buf[5];
  tmva_dipho_MIT_eta2= &tmva_dipho_MIT_buf[6];
  tmva_dipho_MIT_dphi= &tmva_dipho_MIT_buf[7];
  tmva_dipho_MIT_ph1mva= &tmva_dipho_MIT_buf[8];
  tmva_dipho_MIT_ph2mva= &tmva_dipho_MIT_buf[9];
  if( funcReader_dipho_MIT != 0 ) {
	  tmvaReader_dipho_MIT = 0;
	  //// masserr,masserrwrong,vtxprob,pt1,pt2,eta1,eta2,dphi,idmva1,idmva2
	  funcReader_dipho_MIT->bookVariable("masserr",         tmva_dipho_MIT_dmom);
	  funcReader_dipho_MIT->bookVariable("masserrwrongvtx", tmva_dipho_MIT_dmom_wrong_vtx);
	  funcReader_dipho_MIT->bookVariable("vtxprob",         tmva_dipho_MIT_vtxprob);
	  funcReader_dipho_MIT->bookVariable("pt1",             tmva_dipho_MIT_ptom1);
	  funcReader_dipho_MIT->bookVariable("pt2",             tmva_dipho_MIT_ptom2);
	  funcReader_dipho_MIT->bookVariable("eta1",            tmva_dipho_MIT_eta1);
	  funcReader_dipho_MIT->bookVariable("eta2",            tmva_dipho_MIT_eta2);
	  funcReader_dipho_MIT->bookVariable("dphi",            tmva_dipho_MIT_dphi);
	  funcReader_dipho_MIT->bookVariable("idmva1",          tmva_dipho_MIT_ph1mva);
	  funcReader_dipho_MIT->bookVariable("idmva2",          tmva_dipho_MIT_ph2mva);
  } else { 
	  tmvaReader_dipho_MIT = new TMVA::Reader("!Color:Silent"); 
	  tmvaReader_dipho_MIT->AddVariable("masserrsmeared/mass",         tmva_dipho_MIT_dmom);
	  tmvaReader_dipho_MIT->AddVariable("masserrsmearedwrongvtx/mass", tmva_dipho_MIT_dmom_wrong_vtx);
	  tmvaReader_dipho_MIT->AddVariable("vtxprob",                     tmva_dipho_MIT_vtxprob);
	  tmvaReader_dipho_MIT->AddVariable("ph1.pt/mass",                 tmva_dipho_MIT_ptom1);
	  tmvaReader_dipho_MIT->AddVariable("ph2.pt/mass",                 tmva_dipho_MIT_ptom2);
	  tmvaReader_dipho_MIT->AddVariable("ph1.eta",                     tmva_dipho_MIT_eta1);
	  tmvaReader_dipho_MIT->AddVariable("ph2.eta",                     tmva_dipho_MIT_eta2);
	  tmvaReader_dipho_MIT->AddVariable("TMath::Cos(ph1.phi-ph2.phi)", tmva_dipho_MIT_dphi);
	  tmvaReader_dipho_MIT->AddVariable("ph1.idmva",                   tmva_dipho_MIT_ph1mva);
	  tmvaReader_dipho_MIT->AddVariable("ph2.idmva",                   tmva_dipho_MIT_ph2mva);
  }

  tmvaReaderID_Single_Barrel = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_Single_Barrel->AddVariable("ph.r9",   &tmva_photonid_r9 );
  tmvaReaderID_Single_Barrel->AddVariable("ph.sigietaieta",   &tmva_photonid_sieie );
  tmvaReaderID_Single_Barrel->AddVariable("ph.scetawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_Single_Barrel->AddVariable("ph.scphiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_CoviEtaiPhi",   &tmva_photonid_sieip );
  tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_GammaIso",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_ChargedIso_selvtx",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_ChargedIso_worstvtx",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_Single_Barrel->AddVariable("ph.sceta",   &tmva_photonid_sceta );
  tmvaReaderID_Single_Barrel->AddVariable("rho",   &tmva_photonid_eventrho );
  
  tmvaReaderID_Single_Endcap = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_Single_Endcap->AddVariable("ph.r9",   &tmva_photonid_r9 );
  tmvaReaderID_Single_Endcap->AddVariable("ph.sigietaieta",   &tmva_photonid_sieie );
  tmvaReaderID_Single_Endcap->AddVariable("ph.scetawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_Single_Endcap->AddVariable("ph.scphiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_CoviEtaiPhi",   &tmva_photonid_sieip );
  tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_GammaIso",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_ChargedIso_selvtx",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_ChargedIso_worstvtx",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_Single_Endcap->AddVariable("ph.sceta",   &tmva_photonid_sceta );
  tmvaReaderID_Single_Endcap->AddVariable("rho",   &tmva_photonid_eventrho );
  tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &tmva_photonid_ESEffSigmaRR );
  
  tmvaReaderID_2013_Barrel = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_2013_Barrel->AddVariable("ph.scrawe",   &tmva_photonid_scrawe );
  tmvaReaderID_2013_Barrel->AddVariable("ph.r9",   &tmva_photonid_r9 );
  tmvaReaderID_2013_Barrel->AddVariable("ph.sigietaieta",   &tmva_photonid_sieie );
  tmvaReaderID_2013_Barrel->AddVariable("ph.scetawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_2013_Barrel->AddVariable("ph.scphiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_2013_Barrel->AddVariable("ph.idmva_CoviEtaiPhi",   &tmva_photonid_sieip );
  tmvaReaderID_2013_Barrel->AddVariable("ph.idmva_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_2013_Barrel->AddVariable("ph.idmva_GammaIso",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_2013_Barrel->AddVariable("ph.idmva_ChargedIso_selvtx",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_2013_Barrel->AddVariable("ph.idmva_ChargedIso_worstvtx",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_2013_Barrel->AddVariable("ph.sceta",   &tmva_photonid_sceta );
  tmvaReaderID_2013_Barrel->AddVariable("rho",   &tmva_photonid_eventrho );
  
  tmvaReaderID_2013_Endcap = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_2013_Endcap->AddVariable("ph.scrawe",   &tmva_photonid_scrawe );
  tmvaReaderID_2013_Endcap->AddVariable("ph.r9",   &tmva_photonid_r9 );
  tmvaReaderID_2013_Endcap->AddVariable("ph.sigietaieta",   &tmva_photonid_sieie );
  tmvaReaderID_2013_Endcap->AddVariable("ph.scetawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_2013_Endcap->AddVariable("ph.scphiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_2013_Endcap->AddVariable("ph.idmva_CoviEtaiPhi",   &tmva_photonid_sieip );
  tmvaReaderID_2013_Endcap->AddVariable("ph.idmva_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_2013_Endcap->AddVariable("ph.idmva_GammaIso",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_2013_Endcap->AddVariable("ph.idmva_ChargedIso_selvtx",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_2013_Endcap->AddVariable("ph.idmva_ChargedIso_worstvtx",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_2013_Endcap->AddVariable("ph.sceta",   &tmva_photonid_sceta );
  tmvaReaderID_2013_Endcap->AddVariable("rho",   &tmva_photonid_eventrho );
  tmvaReaderID_2013_Endcap->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &tmva_photonid_ESEffSigmaRR ); 

  tmvaReaderID_2013_7TeV_MIT_Barrel = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.scrawe",   &tmva_photonid_scrawe );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.r9",   &tmva_photonid_r9 );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.sigietaieta",   &tmva_photonid_sieie );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.scetawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.scphiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.idmva_CoviEtaiPhi",   &tmva_photonid_sieip );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.idmva_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.idmva_GammaIso",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.idmva_ChargedIso_selvtx",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.idmva_ChargedIso_worstvtx",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("ph.sceta",   &tmva_photonid_sceta );
  tmvaReaderID_2013_7TeV_MIT_Barrel->AddVariable("rho",   &tmva_photonid_eventrho );
  
  tmvaReaderID_2013_7TeV_MIT_Endcap = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.scrawe",   &tmva_photonid_scrawe );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.r9",   &tmva_photonid_r9 );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.sigietaieta",   &tmva_photonid_sieie );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.scetawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.scphiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.idmva_CoviEtaiPhi",   &tmva_photonid_sieip );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.idmva_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.idmva_GammaIso",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.idmva_ChargedIso_selvtx",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.idmva_ChargedIso_worstvtx",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.sceta",   &tmva_photonid_sceta );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("rho",   &tmva_photonid_eventrho );
  tmvaReaderID_2013_7TeV_MIT_Endcap->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &tmva_photonid_ESEffSigmaRR );
}

Float_t LoopAll::photonIDMVA2013(Int_t iPhoton, Int_t vtx, TLorentzVector &p4, const char* type)  {

    Float_t mva = 999.;

    double pfchargedisobad03=0.;
    for(int ivtx=0; ivtx<vtx_std_n; ivtx++) {
        pfchargedisobad03=(*pho_pfiso_mycharged03)[iPhoton][ivtx]>pfchargedisobad03?(*pho_pfiso_mycharged03)[iPhoton][ivtx]:pfchargedisobad03;
    }

    tmva_photonid_pfchargedisogood03 = (*pho_pfiso_mycharged03)[iPhoton][vtx];
    tmva_photonid_pfchargedisobad03  = pfchargedisobad03;
    tmva_photonid_pfphotoniso03      = pho_pfiso_myphoton03[iPhoton];
    tmva_photonid_pfneutraliso03     = pho_pfiso_myneutral03[iPhoton]; 
  
    tmva_photonid_sieie        = pho_sieie[iPhoton];
    tmva_photonid_sieip        = pho_sieip[iPhoton];
    tmva_photonid_etawidth     = pho_etawidth[iPhoton];
    int scind=pho_scind[iPhoton];
    tmva_photonid_scrawe       = sc_raw[scind];
    tmva_photonid_phiwidth     = sc_sphi[scind]; //pho_etawidth[iPhoton]*pho_brem[iPhoton]; //sc_sphi[pho_scind[ipho]]
    tmva_photonid_r9           = pho_r9[iPhoton];
    tmva_photonid_lambdaratio  = pho_lambdaratio[iPhoton];
  
    //  tmva_photonid_s4ratio  = pho_e2x2[iPhoton]/pho_e5x5[iPhoton];
  tmva_photonid_s4ratio  = pho_s4ratio[iPhoton];
  tmva_photonid_eventrho = rho_algo1;
  tmva_photonid_sceta    = ((TVector3*)sc_xyz->At(pho_scind[iPhoton]))->Eta(); 
  tmva_photonid_ESEffSigmaRR = pho_ESEffSigmaRR[iPhoton];

  if (pho_isEB[iPhoton]) {
    mva = tmvaReaderID_2013_Barrel->EvaluateMVA("AdaBoost");
  } else {
    mva = tmvaReaderID_2013_Endcap->EvaluateMVA("AdaBoost");
  }

    return mva;
}

Float_t LoopAll::photonIDMVA2013_7TeV(Int_t iPhoton, Int_t vtx, TLorentzVector &p4, const char* type)  {

    Float_t mva = 999.;

    double pfchargedisobad03=0.;
    for(int ivtx=0; ivtx<vtx_std_n; ivtx++) {
        pfchargedisobad03=(*pho_pfiso_mycharged03)[iPhoton][ivtx]>pfchargedisobad03?(*pho_pfiso_mycharged03)[iPhoton][ivtx]:pfchargedisobad03;
    }

    tmva_photonid_pfchargedisogood03 = (*pho_pfiso_mycharged03)[iPhoton][vtx];
    tmva_photonid_pfchargedisobad03  = pfchargedisobad03;
    tmva_photonid_pfphotoniso03      = pho_pfiso_myphoton03[iPhoton];
    tmva_photonid_pfneutraliso03     = pho_pfiso_myneutral03[iPhoton]; 
  
    tmva_photonid_sieie        = pho_sieie[iPhoton];
    tmva_photonid_sieip        = pho_sieip[iPhoton];
    tmva_photonid_etawidth     = pho_etawidth[iPhoton];
    int scind=pho_scind[iPhoton];
    tmva_photonid_scrawe       = sc_raw[scind];
    tmva_photonid_phiwidth     = sc_sphi[scind]; //pho_etawidth[iPhoton]*pho_brem[iPhoton]; //sc_sphi[pho_scind[ipho]]
    tmva_photonid_r9           = pho_r9[iPhoton];
    tmva_photonid_lambdaratio  = pho_lambdaratio[iPhoton];
  
    tmva_photonid_s4ratio  = pho_s4ratio[iPhoton];
    tmva_photonid_eventrho = rho_algo1;
    tmva_photonid_sceta    = ((TVector3*)sc_xyz->At(pho_scind[iPhoton]))->Eta(); 
    tmva_photonid_ESEffSigmaRR = pho_ESEffSigmaRR[iPhoton];
    
    //std::cout << tmva_photonid_pfchargedisogood03 << " "  << tmva_photonid_pfchargedisobad03 << " " 
    //	      << tmva_photonid_pfphotoniso03 << " "<< tmva_photonid_sieie << " " << tmva_photonid_sieip << " " 
    //	      << tmva_photonid_etawidth << " " << tmva_photonid_scrawe << " " << tmva_photonid_phiwidth << " " 
    //	      << tmva_photonid_lambdaratio << " " <<  tmva_photonid_s4ratio  << " " << tmva_photonid_eventrho << " "
    //	      << tmva_photonid_ESEffSigmaRR << std::endl;
      

    if (pho_isEB[iPhoton]) {
      mva = tmvaReaderID_2013_7TeV_MIT_Barrel->EvaluateMVA("AdaBoost");
    } else {
      mva = tmvaReaderID_2013_7TeV_MIT_Endcap->EvaluateMVA("AdaBoost");
    }
    
    return mva;
}

Float_t LoopAll::photonIDMVA2012(Int_t iPhoton, Int_t vtx, TLorentzVector &p4, const char* type)  {

    Float_t mva = 999.;

    double pfchargedisobad03=0.;
    for(int ivtx=0; ivtx<vtx_std_n; ivtx++) {
        pfchargedisobad03=(*pho_pfiso_mycharged03)[iPhoton][ivtx]>pfchargedisobad03?(*pho_pfiso_mycharged03)[iPhoton][ivtx]:pfchargedisobad03;
    }

    tmva_photonid_pfchargedisogood03 = (*pho_pfiso_mycharged03)[iPhoton][vtx];
    tmva_photonid_pfchargedisobad03  = pfchargedisobad03;
    tmva_photonid_pfphotoniso03      = pho_pfiso_myphoton03[iPhoton];
    tmva_photonid_pfneutraliso03     = pho_pfiso_myneutral03[iPhoton]; 
  
    tmva_photonid_sieie        = pho_sieie[iPhoton];
    tmva_photonid_sieip        = pho_sieip[iPhoton];
    tmva_photonid_etawidth     = pho_etawidth[iPhoton];
    int scind=pho_scind[iPhoton];
    tmva_photonid_phiwidth     = sc_sphi[scind]; //pho_etawidth[iPhoton]*pho_brem[iPhoton]; //sc_sphi[pho_scind[ipho]]
    tmva_photonid_r9           = pho_r9[iPhoton];
    tmva_photonid_lambdaratio  = pho_lambdaratio[iPhoton];
  
//  tmva_photonid_s4ratio  = pho_e2x2[iPhoton]/pho_e5x5[iPhoton];
  tmva_photonid_s4ratio  = pho_s4ratio[iPhoton];
  tmva_photonid_eventrho = rho_algo1;
  tmva_photonid_sceta    = ((TVector3*)sc_xyz->At(pho_scind[iPhoton]))->Eta(); 
  tmva_photonid_ESEffSigmaRR = pho_ESEffSigmaRR[iPhoton];

  if (pho_isEB[iPhoton]) {
    mva = tmvaReaderID_Single_Barrel->EvaluateMVA("AdaBoost");
  }
  else {
    mva = tmvaReaderID_Single_Endcap->EvaluateMVA("AdaBoost");
  }

    return mva;
}

Float_t LoopAll::photonIDMVA(Int_t iPhoton, Int_t vtx, TLorentzVector &p4, const char* type)  {
  if( pho_idmva_cached && pho_idmva[iPhoton][vtx] > -2 ) {
    return pho_idmva[iPhoton][vtx];
  }
  TString t(type);
  if( t == "MIT" ) {
    return photonIDMVA2013(iPhoton,vtx,p4,"MIT");
  } else if( t == "Moriond2013" ) {
    return photonIDMVA2012(iPhoton,vtx,p4,"MIT");
  } else if( t == "Old7TeV" ) {
    return photonIDMVA2011(iPhoton,vtx,p4,"MIT");
  } else if( t == "2013_7TeV") {
    return photonIDMVA2013_7TeV(iPhoton, vtx, p4, "MIT");
  } else {
    std::cerr << "Uknown BDT type " << t << std::endl;
    assert(0);
  }
}

Float_t LoopAll::photonIDMVA2011(Int_t iPhoton, Int_t vtx, TLorentzVector &p4, const char* type)  {
  
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
        tmva_id_mit_sceta    = ((TVector3*)sc_xyz->At(pho_scind[iPhoton]))->Eta();

        if (pho_isEB[iPhoton]) 
            mva = tmvaReaderID_MIT_Barrel->EvaluateMVA("AdaBoost");
        else
            mva = tmvaReaderID_MIT_Endcap->EvaluateMVA("AdaBoost");
    }
    return mva;
}

Float_t LoopAll::diphotonMVA(Int_t diphoton_id, Int_t leadingPho, Int_t subleadingPho, Int_t vtx, float vtxProb, TLorentzVector &leadP4, TLorentzVector &subleadP4, float sigmaMrv, float sigmaMwv, float sigmaMeonly, const char* idType, const char* bdtType, float photonID_1,float photonID_2) {

    // Ok need to re-write the diphoton-mva part since the systematics won't work unless we can change the Et of the photons
    // all we have to do is to pass in the ->Et of the two photons also rather than take them from the four-vector branches
  
    Float_t mva = 99.;
    TLorentzVector Higgs = leadP4+subleadP4;
    float leadPt    = leadP4.Pt();
    float subleadPt = subleadP4.Pt();
    float mass     = Higgs.M();
    float diphopt   = Higgs.Pt();

    if (idType == "UCSD") {
        tmva_dipho_UCSD_leadr9 = pho_r9[leadingPho];
        tmva_dipho_UCSD_subleadr9 = pho_r9[subleadingPho];
        tmva_dipho_UCSD_leadeta = fabs(leadP4.Eta());
        tmva_dipho_UCSD_subleadeta = fabs(subleadP4.Eta());
        //  tmva_dipho_UCSD_leadptomass = leadPt/mass;  
        tmva_dipho_UCSD_subleadptomass = subleadPt/mass;  
        tmva_dipho_UCSD_diphoptom = diphopt/mass;
        tmva_dipho_UCSD_sumptom = (leadPt+subleadPt)/mass;
        tmva_dipho_UCSD_subleadmva = photonIDMVA2011(subleadingPho, vtx,leadP4, "UCSD");
        tmva_dipho_UCSD_leadmva = photonIDMVA2011(leadingPho, vtx,subleadP4, "UCSD");
        // tmva_dipho_UCSD_dmom = sigmaMrv/mass;
        tmva_dipho_UCSD_dmom = sigmaMeonly/mass;
  
        mva = tmvaReader_dipho_UCSD->EvaluateMVA("Gradient");
    } else {
        *tmva_dipho_MIT_dmom = sigmaMrv/mass;
        *tmva_dipho_MIT_dmom_wrong_vtx = sigmaMwv/mass;
        *tmva_dipho_MIT_vtxprob = vtxProb;
        *tmva_dipho_MIT_ptom1 = leadPt/mass;
        *tmva_dipho_MIT_ptom2 = subleadPt/mass;
        *tmva_dipho_MIT_eta1 = leadP4.Eta();
        *tmva_dipho_MIT_eta2 =  subleadP4.Eta();
        *tmva_dipho_MIT_dphi = TMath::Cos(leadP4.Phi() - subleadP4.Phi());
      
        if (photonID_1 < -1. && photonID_2 < -1.) {
          *tmva_dipho_MIT_ph1mva = photonIDMVA(leadingPho,vtx, leadP4, bdtType);
          *tmva_dipho_MIT_ph2mva = photonIDMVA(subleadingPho,vtx, subleadP4, bdtType);
        } else {
          *tmva_dipho_MIT_ph1mva = photonID_1;
          *tmva_dipho_MIT_ph2mva = photonID_2;
        }
	
	//std::cout << *tmva_dipho_MIT_dmom << " " << *tmva_dipho_MIT_dmom_wrong_vtx << " " << *tmva_dipho_MIT_vtxprob
	//	  << *tmva_dipho_MIT_ptom1 << " " << *tmva_dipho_MIT_ptom2 << " " << *tmva_dipho_MIT_eta1 << " " << *tmva_dipho_MIT_dphi 
	//	  << *tmva_dipho_MIT_ph1mva << " " << *tmva_dipho_MIT_ph2mva << std::endl;

        tmva_dipho_MIT_cache[diphoton_id] = tmva_dipho_MIT_buf;
        mva = ( funcReader_dipho_MIT != 0 ? funcReader_dipho_MIT->eval() : tmvaReader_dipho_MIT->EvaluateMVA("Gradient") );
	//std::cout << mva << std::endl;
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
    virtual int ntracks() const { return lo_.tk_p4->GetEntries(); } // return lo_.tk_n; };
  
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
    //  PhotonInfo
    //    pho1(p1,*((TVector3*)pho_calopos->At(p1)),((TLorentzVector*)pho_p4->At(p1))->Energy()),
    //  pho2(p2,*((TVector3*)pho_calopos->At(p2)),((TLorentzVector*)pho_p4->At(p2))->Energy());
    vtxAna.analyze(vinfo,pho1,pho2);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
PhotonInfo LoopAll::fillPhotonInfos(int p1, int useAllConvs, float * energy) 
{
  
    int iConv1 = useAllConvs>0 ? matchPhotonToConversion(p1,useAllConvs) : -1;
    
    if ( iConv1 >= 0) {
        // conversions infos
        return PhotonInfo(p1,
                          // *((TVector3*)pho_calopos->At(p1)),
                          *((TVector3*)sc_xyz->At(pho_scind[p1])),
                          *((TVector3*) bs_xyz->At(0)),
                          *((TVector3*) conv_vtx->At(iConv1)),
                          conv_ntracks[iConv1] == 1 ? *((TVector3*) conv_singleleg_momentum->At(iConv1)) : *((TVector3*) conv_refitted_momentum->At(iConv1)),
                          energy == 0 ? ((TLorentzVector*)pho_p4->At(p1))->Energy() : energy[p1],
                          pho_isEB[p1],
                          conv_ntracks[iConv1],
                          conv_validvtx[iConv1],
                          conv_chi2_probability[iConv1],
                          conv_eoverp[iConv1]
            );
    } 
    //// else {
    ////   return PhotonInfo(p1,*((TVector3*)pho_calopos->At(p1)),((TLorentzVector*)pho_p4->At(p1))->Energy());
    //// }
  
    return PhotonInfo(p1, 
                      // *((TVector3*)pho_calopos->At(p1)),                                                                                                                
                      *((TVector3*)sc_xyz->At(pho_scind[p1])),
                      *((TVector3*) bs_xyz->At(0)),
                      *((TVector3*) pho_conv_vtx->At(p1)),
                      *((TVector3*) pho_conv_refitted_momentum->At(p1)),
                      energy == 0 ? ((TLorentzVector*)pho_p4->At(p1))->Energy() : energy[p1],
                      pho_isEB[p1],
                      0,
                      0,
                      0,
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

        // at least one of the photons is identified as a conversion
    
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
    
    } // end if at least one photon is a conversion
  
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

//----------------------------------------------------------------------

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

int  LoopAll::matchPhotonToConversion( int lpho, int useAllConvs) {

    int iMatch=-1;
    TVector3 Photonxyz = *((TVector3*) sc_xyz->At(pho_scind[lpho]));

    float detaMin=999.;
    float dphiMin=999.;   
    float dRMin = 999.;

    for(int iconv=0; iconv<conv_n; iconv++) {

        TVector3 refittedPairMomentum= conv_ntracks[iconv]==1 ? *((TVector3*) conv_singleleg_momentum->At(iconv)) : *((TVector3*) conv_refitted_momentum->At(iconv));

        //Conversion Selection
        if ( refittedPairMomentum.Pt() < 10 ) continue;
        if ( useAllConvs==1 && conv_ntracks[iconv]!=1 ) continue;
        if ( useAllConvs==2 && conv_ntracks[iconv]!=2 ) continue;
        if ( useAllConvs==3 && conv_ntracks[iconv]!=1 && conv_ntracks[iconv]!=2 ) continue;
        if ( conv_ntracks[iconv]==2 && (!conv_validvtx[iconv] || conv_chi2_probability[iconv]<0.000001) ) continue; // Changed back based on meeting on 21.03.2012

        //New matching technique from meeting on 06.08.12
        TVector3 ConversionVertex = *((TVector3*) conv_vtx->At(iconv));
        TVector3 NewPhotonxyz = Photonxyz-ConversionVertex;
        double dR = NewPhotonxyz.DeltaR(refittedPairMomentum);
        double delta_eta = NewPhotonxyz.Eta()-refittedPairMomentum.Eta();
        double delta_phi = NewPhotonxyz.DeltaPhi(refittedPairMomentum);

        if ( dR < dRMin ) {
            detaMin=fabs(delta_eta);
            dphiMin=fabs(delta_phi);
            dRMin=dR;
            iMatch=iconv;
        }
    
    }

    if ( dRMin< 0.1 ) return iMatch;
    else return -1;
    
}

//-Applying MET correction-------------
//correctMETinRED
// pfjet resolutions. taken from AN-2010-371
// --------------------------------------------------------------------------
double LoopAll::ErrEt( double Et, double Eta) {
  
    double InvPerr2;
  
    double N, S, C, m;
    if(fabs(Eta) < 0.5 ) {
        N = 3.96859;
        S = 0.18348;
        C = 0.;
        m = 0.62627;
    } else if( fabs(Eta) < 1. ) {
        N = 3.55226;
        S = 0.24026;
        C = 0.;
        m = 0.52571;
    } else if( fabs(Eta) < 1.5 ) {
        N = 4.54826;
        S = 0.22652;
        C = 0.;
        m = 0.58963;
    } else if( fabs(Eta) < 2. ) {
        N = 4.62622;
        S = 0.23664;
        C = 0.;
        m = 0.48738;
    } else if( fabs(Eta) < 3. ) {
        N = 2.53324;
        S = 0.34306;
        C = 0.;
        m = 0.28662;
    } else if( fabs(Eta) < 5. ) {
        N = 2.95397;
        S = 0.11619;
        C = 0.;
        m = 0.96086;
    }
  
    // this is the absolute resolution (squared), not sigma(pt)/pt
    // so have to multiply by pt^2, thats why m+1 instead of m-1
    InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;
  
  
    return sqrt(InvPerr2)/Et;

}


TLorentzVector LoopAll::shiftMet(TLorentzVector *uncormet, bool isMC) {

    TLorentzVector correctedMet;
    
    // correction for METx, METy bias
    double px(0), py(0), e(0);
    
////2011
//if(!isMC){
//px = uncormet->Pt()*cos(uncormet->Phi())-0.00563109*met_sumet_pfmet+0.959742;
//py = uncormet->Pt()*sin(uncormet->Phi())+0.00586162*met_sumet_pfmet-0.540137;
//}else{
//px = uncormet->Pt()*cos(uncormet->Phi())-0.00069992*met_sumet_pfmet+0.430059;
//py = uncormet->Pt()*sin(uncormet->Phi())+0.00262869*met_sumet_pfmet+0.210784;
//}
    
    //2012
    // data
    if(!isMC){
        px = uncormet->Pt()*cos(uncormet->Phi())-0.006239*met_sumet_pfmet+0.662;
        py = uncormet->Pt()*sin(uncormet->Phi())+0.004613*met_sumet_pfmet-0.673;
    // MC
    }else{
        px = uncormet->Pt()*cos(uncormet->Phi())+0.00135*met_sumet_pfmet-0.021;
        py = uncormet->Pt()*sin(uncormet->Phi())+0.00371*met_sumet_pfmet-0.826;
    }
    e = sqrt(px*px+py*py);
    correctedMet.SetPxPyPzE(px,py,0,e);
    
    return correctedMet;
}


//met at analysis step
TLorentzVector LoopAll::correctMet_Simple( TLorentzVector & pho_lead, TLorentzVector & pho_sublead,TLorentzVector *uncormet, bool smearing, bool scale) {

    TRandom3 jSmearRan(event);
    TLorentzVector jetSumSmeared;
    jetSumSmeared.SetXYZT(0.,0.,0.,0);
    TLorentzVector jetSumUnsmeared;
    jetSumUnsmeared.SetXYZT(0.,0.,0.,0);
    
    TLorentzVector p4_jet;

    //associating reco - gen met
    for(int i=0; i<jet_algoPF1_n; i++){
      // TLorentzVector * p4_jet = (TLorentzVector *) jet_algoPF1_p4->At(i);
	p4_jet = (TLorentzVector&) *( jet_algoPF1_p4->At(i) );
        if( version >= 13 ) {
	    //// p4_jet = (TLorentzVector *) p4_jet->Clone();
            //// *p4_jet = (*p4_jet) * (1/jet_algoPF1_erescale[i]);
	    p4_jet *= 1. / jet_algoPF1_erescale[i];
        }
        bool isJet_LeadPho = false;
        bool isJet_SubLeadPho = false;
        
        double dR_jet_PhoLead = p4_jet.DeltaR(pho_lead);
        if( dR_jet_PhoLead<0.5 ) isJet_LeadPho = true;
        
        double dR_jet_PhoSubLead = p4_jet.DeltaR(pho_sublead);
        if( dR_jet_PhoSubLead<0.5 ) isJet_SubLeadPho = true;
        
        if( isJet_LeadPho || isJet_SubLeadPho ) continue;
        
        double ptJet_pfakt5 = p4_jet.Pt();
        double eJet_pfakt5 = p4_jet.Energy();
        double etaJet_pfakt5 = p4_jet.Eta();
        double phiJet_pfakt5 = p4_jet.Phi();
        double ptCorrJet_pfakt5 = ptJet_pfakt5*jet_algoPF1_erescale[i];
        
        //smearing via association with genjets
        double DRmin(999.);    
        double expres = ErrEt(ptCorrJet_pfakt5,etaJet_pfakt5);
        
        if (jet_algoPF1_genMatched[i] && (ptCorrJet_pfakt5-jet_algoPF1_genPt[i])/ptCorrJet_pfakt5 < 5. * expres) {
            DRmin=jet_algoPF1_genDr[i];
        }
        
        if (jet_algoPF1_genMatched[i]) {
            if(DRmin > 0.1 + 0.3 * exp(-0.05*(jet_algoPF1_genPt[i]-10)))  { jet_algoPF1_genMatched[i]=false; }
        }
        
        //smearing for non-associated jets, using expected resolutions
        float smear = -999.;
        if (fabs(etaJet_pfakt5)<=1.1)                            smear = 1.06177;
        if (fabs(etaJet_pfakt5)<=1.7 && fabs(etaJet_pfakt5)>1.1) smear = 1.08352;
        if (fabs(etaJet_pfakt5)<=2.3 && fabs(etaJet_pfakt5)>1.7) smear = 1.02911;
        if (fabs(etaJet_pfakt5)>2.3)                             smear = 1.15288;
        
        double shift(0);
        
        if(jet_algoPF1_genMatched[i]) {    
            shift = (smear-1) * (ptCorrJet_pfakt5 - jet_algoPF1_genPt[i])/ptCorrJet_pfakt5; }
        else {
            double expres = ErrEt(ptJet_pfakt5, etaJet_pfakt5);
            double relsmear = expres * sqrt(smear*smear-1);
            jSmearRan.SetSeed(event+(Int_t)(etaJet_pfakt5*1000));
            shift = jSmearRan.Gaus(0.,relsmear);
        }
        
        float ptSmeared  = ptJet_pfakt5;
        float eneSmeared = eJet_pfakt5;
        
        if(smearing && shift>-1 && shift < 2) {
            ptSmeared  *= 1 + shift;
            eneSmeared *= 1 + shift;
        }
        
        //JEC scaling to correct for residual jet corrections
        if(scale) {
            double factor(1);
            if(TMath::Abs(etaJet_pfakt5)<1.5) factor = 1.015;
            else if(TMath::Abs(etaJet_pfakt5)<3) factor = 1.04;
            else factor = 1.15;
            ptSmeared  *= factor;
            eneSmeared *= factor;
        }
        
        TLorentzVector thisJetSmeared;
        thisJetSmeared.SetPtEtaPhiE(ptSmeared,etaJet_pfakt5,phiJet_pfakt5,eneSmeared);
        
        TLorentzVector thisJetUnsmeared;
        thisJetUnsmeared.SetPtEtaPhiE(ptJet_pfakt5,etaJet_pfakt5,phiJet_pfakt5,eJet_pfakt5);
        
        if (ptJet_pfakt5>10 && TMath::Abs(etaJet_pfakt5)<4.7) {
            jetSumSmeared   += thisJetSmeared;
            jetSumUnsmeared += thisJetUnsmeared;
        }
    
    }
    
    TLorentzVector correctedMet;
    correctedMet = (*uncormet) + jetSumUnsmeared - jetSumSmeared;
    return correctedMet;
}

TLorentzVector LoopAll::METCorrection2012B(TLorentzVector lead_p4, TLorentzVector sublead_p4, bool moriond2013MetCorrection){
  
  // corrected met
  static TLorentzVector finalCorrMET;
  static int lastEvent=-1, lastLumi=-1, lastRun=-1;
  static bool lastIsMC;
  static TLorentzVector last_lead, last_sublead;
  
  bool isMC = itype[current]!=0;

  if( event == lastEvent && lumis == lastLumi && run == lastRun && lastIsMC == isMC &&
      last_lead == lead_p4 && last_sublead == sublead_p4 ) {
	  return finalCorrMET;
  }

  // uncorrected PF met
  TLorentzVector unpfMET;
  unpfMET.SetPxPyPzE (met_pfmet*cos(met_phi_pfmet),met_pfmet*sin(met_phi_pfmet),0,
		      sqrt(met_pfmet*cos(met_phi_pfmet) * met_pfmet*cos(met_phi_pfmet) 
			   + met_pfmet*sin(met_phi_pfmet) * met_pfmet*sin(met_phi_pfmet))); 
  
  if (isMC) {
    // smear raw met
    TLorentzVector smearMET_corr = correctMet_Simple( lead_p4, sublead_p4 , &unpfMET, true, false);
    // then shift smeared met
    float px  = smearMET_corr.Pt()*cos(smearMET_corr.Phi())-0.000685*met_sumet_pfmet+0.084403;
    float py  = smearMET_corr.Pt()*sin(smearMET_corr.Phi())+0.003950*met_sumet_pfmet-0.415907;
    if(moriond2013MetCorrection){
      px  = smearMET_corr.Pt()*cos(smearMET_corr.Phi())+0.00135*met_sumet_pfmet-0.021;
      py  = smearMET_corr.Pt()*sin(smearMET_corr.Phi())+0.00371*met_sumet_pfmet-0.826;
    }
    float ene = sqrt(px*px+py*py);
    finalCorrMET.SetPxPyPzE(px,py,0,ene);
  } else {
    // shifted met for data
    float px  = unpfMET.Pt()*cos(unpfMET.Phi())-0.005117*met_sumet_pfmet+0.830150;
    float py  = unpfMET.Pt()*sin(unpfMET.Phi())+0.002813*met_sumet_pfmet-0.384595;
    if(moriond2013MetCorrection){
      float px  = unpfMET.Pt()*cos(unpfMET.Phi())-0.006239*met_sumet_pfmet+0.662;
      float py  = unpfMET.Pt()*sin(unpfMET.Phi())+0.004613*met_sumet_pfmet-0.673;
    }
    float ene = sqrt(px*px+py*py);
    TLorentzVector shiftedMET;
    shiftedMET.SetPxPyPzE(px,py,0,ene);
    // scale shifted met
    TLorentzVector shiftscaleMET_corr = correctMet_Simple( lead_p4, sublead_p4 , &shiftedMET, false , true);
    finalCorrMET = shiftscaleMET_corr;
  }

  lastEvent = event;
  lastLumi  = lumis;
  lastRun   = run;
  lastIsMC  = isMC;
  last_lead = lead_p4;
  last_sublead = sublead_p4;
  
  return finalCorrMET;
}

bool LoopAll::METAnalysis2012B(TLorentzVector lead_p4, TLorentzVector sublead_p4, bool useUncorr, bool doMETCleaning, bool moriond2013MetCorrection){

  bool tag=false;

  TLorentzVector myMet;
  if (useUncorr) {
    myMet.SetPxPyPzE (met_pfmet*cos(met_phi_pfmet),met_pfmet*sin(met_phi_pfmet),0, sqrt(met_pfmet*cos(met_phi_pfmet) * met_pfmet*cos(met_phi_pfmet) + met_pfmet*sin(met_phi_pfmet) * met_pfmet*sin(met_phi_pfmet))); 
  } else {
    myMet = METCorrection2012B(lead_p4, sublead_p4, moriond2013MetCorrection);
  }		   

  // TLorentzVector TwoPhoton_p4 = lead_p4 + sublead_p4;
  // TVector3 lead_p3      = lead_p4.Vect();
  // TVector3 TwoPhoton_p3 = TwoPhoton_p4.Vect();

  // additional angular cuts
  // TVector3 met_p3;
  // met_p3.SetPtEtaPhi(met_pfmet,0,met_phi_pfmet);
  // float dPhiMetGG   = fabs(180./3.1415927 * TwoPhoton_p3.DeltaPhi(met_p3));
  // float dPhiMetLead = fabs(180./3.1415927 * lead_p3.DeltaPhi(met_p3));
  // if ( dPhiMetGG<40)             return tag;
  // if ( dPhiMetLead<40)           return tag;
  // if( l.correctedpfMET > 70 )    tag = true;

  // EBEB only
  float leadEta    = lead_p4.Eta();
  float subleadEta = sublead_p4.Eta();
  if (fabs(leadEta)>1.5)    return tag;
  if (fabs(subleadEta)>1.5) return tag;
 
  if(doMETCleaning){
    bool passcleaning = METCleaning2012B(lead_p4, sublead_p4, myMet);
    if(!passcleaning) return tag;
  }

  if ( myMet.Pt() > 70 ) tag = true;

  return tag;
}

bool LoopAll::METCleaning2012B(TLorentzVector& lead_p4, TLorentzVector& sublead_p4, TLorentzVector& myMet){

    bool pass=false;
    if(GFDEBUG) std::cout<<"cleaning met"<<std::endl;
    TLorentzVector diphoton_p4 = lead_p4 + sublead_p4;
    float dPhiMetGG   = fabs(diphoton_p4.DeltaPhi(myMet));
    if( dPhiMetGG <= 2.1 ) return pass;

    float jetptmin=50;
    float jethighpt=0;
    TLorentzVector* hiptjet;
    for(int ijet=0; ijet<jet_algoPF1_n; ijet++){
        TLorentzVector* thisjet=(TLorentzVector*) jet_algoPF1_p4->At(ijet);
        if( thisjet->DeltaR(lead_p4)<0.5) continue;
        if( thisjet->DeltaR(sublead_p4)<0.5) continue;
        if(thisjet->Pt() > jethighpt){
            jethighpt= thisjet->Pt();
            hiptjet=thisjet;
        }
    }

    if(jethighpt>jetptmin){
        if(fabs(hiptjet->DeltaPhi(myMet))>=2.7){
          return pass;
        }
    }

    pass=true;
    return pass;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector LoopAll::get_pho_p4(int ipho, int ivtx, const float * energy) const
{
    /// /// PhotonInfo p(ipho, *((TVector3*)sc_xyz->At(pho_scind[ipho])),
    /// PhotonInfo p(ipho, *((TVector3*)pho_calopos->At(ipho)),
    ///        energy != 0 ? energy[ipho] : ((TLorentzVector*)pho_p4->At(ipho))->Energy() );
    /// TVector3 * vtx = (TVector3*) vtx_std_xyz->At(ivtx);
    /// return p.p4( vtx->X(), vtx->Y(), vtx->Z() );
    TVector3 * vtx = (TVector3*) vtx_std_xyz->At(ivtx);
    return get_pho_p4(ipho,vtx,energy);
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector LoopAll::get_pho_p4(int ipho, TVector3 * vtx, const float * energy) const
{
    /// PhotonInfo p(ipho, *((TVector3*)sc_xyz->At(pho_scind[ipho])),
    if(GFDEBUG) std::cout<<"General Functions::get p4 -- ipho energy p4energy"<<ipho<<" "<<energy[ipho]<<" "<<((TLorentzVector*)pho_p4->At(ipho))->Energy()<<std::endl;
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
void LoopAll::FillCICPFInputs()
{
    for(int ipho=0; ipho<pho_n; ++ipho) {
        float neu01 = pfEcalIso(ipho, 0.1, 0., 0., 0., 0., 0., 0., 5);
        float neu02 = pfEcalIso(ipho, 0.2, 0., 0., 0., 0., 0., 0., 5);
        float neu03 = pfEcalIso(ipho, 0.3, 0., 0., 0., 0., 0., 0., 5);
        float neu04 = pfEcalIso(ipho, 0.4, 0., 0., 0., 0., 0., 0., 5); 
        float neu05 = pfEcalIso(ipho, 0.5, 0., 0., 0., 0., 0., 0., 5); 
        float neu06 = pfEcalIso(ipho, 0.6, 0., 0., 0., 0., 0., 0., 5); 
        if( GFDEBUG ) {
            if( ( pho_pfiso_myneutral03[ipho] != neu03 || 
                  pho_pfiso_myneutral04[ipho] != neu04   )
                ) { std::cout << "Fishy... "; }
            std::cout << "neu03: " << pho_pfiso_myneutral03[ipho] << " " << neu03 
                      << " neu04: " << pho_pfiso_myneutral04[ipho] << " " << neu04 
                      << " " << ((TLorentzVector*)pho_p4->At(ipho))->Pt() 
                      << " " << ((TLorentzVector*)pho_p4->At(ipho))->Eta() 
                      << std::endl;
        }
        pho_pfiso_myneutral01[ipho] = neu01;
        pho_pfiso_myneutral02[ipho] = neu02;
        pho_pfiso_myneutral03[ipho] = neu03;
        pho_pfiso_myneutral04[ipho] = neu04;
        pho_pfiso_myneutral05[ipho] = neu05;
        pho_pfiso_myneutral06[ipho] = neu06;


        float pho01 = pfEcalIso(ipho, 0.1, 0., 0.070, 0.015, 0., 0., 0.);
        float pho02 = pfEcalIso(ipho, 0.2, 0., 0.070, 0.015, 0., 0., 0.);
        float pho03 = pfEcalIso(ipho, 0.3, 0., 0.070, 0.015, 0., 0., 0.);
        float pho04 = pfEcalIso(ipho, 0.4, 0., 0.070, 0.015, 0., 0., 0.); 
        float pho05 = pfEcalIso(ipho, 0.5, 0., 0.070, 0.015, 0., 0., 0.); 
        float pho06 = pfEcalIso(ipho, 0.6, 0., 0.070, 0.015, 0., 0., 0.); 
        ///// float pho03 = pfEcalIso(ipho, 0.3, 0.045, 0.070, 0.015, 0.015, 0.08, 0.1);
        ///// float pho04 = pfEcalIso(ipho, 0.4, 0.045, 0.070, 0.015, 0.015, 0.08, 0.1); 
        if( GFDEBUG ) {
            if( ( pho_pfiso_myphoton03[ipho] != pho03 || 
                  pho_pfiso_myphoton04[ipho] != pho04   )
                ) { std::cout << "Fishy... "; }
            std::cout << "pho03: " << pho_pfiso_myphoton03[ipho] << " " << pho03 
                      << " pho04: " << pho_pfiso_myphoton04[ipho] << " " << pho04 
                      << " " << ((TLorentzVector*)pho_p4->At(ipho))->Pt() 
                      << " " << ((TLorentzVector*)pho_p4->At(ipho))->Eta() 
                      << std::endl;
        }
        pho_pfiso_myphoton01[ipho] = pho01;
        pho_pfiso_myphoton02[ipho] = pho02;
        pho_pfiso_myphoton03[ipho] = pho03;
        pho_pfiso_myphoton04[ipho] = pho04;
        pho_pfiso_myphoton05[ipho] = pho05;
        pho_pfiso_myphoton06[ipho] = pho06;

        int badvtx = 0;
        float badiso = 0.;
        for(int ivtx=0; ivtx<vtx_std_n; ++ivtx) {
            float ch01 = pfTkIsoWithVertex(ipho,ivtx,0.1,0.02,0.02,0.0,0.2,0.1);
            float ch02 = pfTkIsoWithVertex(ipho,ivtx,0.2,0.02,0.02,0.0,0.2,0.1);
            float ch03 = pfTkIsoWithVertex(ipho,ivtx,0.3,0.02,0.02,0.0,0.2,0.1);
            float ch04 = pfTkIsoWithVertex(ipho,ivtx,0.4,0.02,0.02,0.0,0.2,0.1);
            float ch05 = pfTkIsoWithVertex(ipho,ivtx,0.5,0.02,0.02,0.0,0.2,0.1);
            float ch06 = pfTkIsoWithVertex(ipho,ivtx,0.6,0.02,0.02,0.0,0.2,0.1);
            ///// float ch03 = pfTkIsoWithVertex(ipho,ivtx,0.3,0.02,0.02,1.0,0.2,0.1);
            ///// float ch04 = pfTkIsoWithVertex(ipho,ivtx,0.4,0.02,0.02,1.0,0.2,0.1);
            if( GFDEBUG ) {
                if( ( pho_pfiso_mycharged03->at(ipho).at(ivtx) != ch03 ||
                      pho_pfiso_mycharged04->at(ipho).at(ivtx) != ch04   )
                    )  { std::cout << "Fishy... "; }
                std::cout << "ch03: " << pho_pfiso_mycharged03->at(ipho).at(ivtx) << " " << ch03 
                          << " ch04: " << pho_pfiso_mycharged04->at(ipho).at(ivtx) << " " << ch04 
                          << " " << ((TLorentzVector*)pho_p4->At(ipho))->Pt() 
                          << " " << ((TLorentzVector*)pho_p4->At(ipho))->Eta() 
                          << std::endl;
            }
            pho_pfiso_mycharged01->at(ipho).at(ivtx) = ch01;
            pho_pfiso_mycharged02->at(ipho).at(ivtx) = ch02;
            pho_pfiso_mycharged03->at(ipho).at(ivtx) = ch03;
            pho_pfiso_mycharged04->at(ipho).at(ivtx) = ch04;
            pho_pfiso_mycharged05->at(ipho).at(ivtx) = ch05;
            pho_pfiso_mycharged06->at(ipho).at(ivtx) = ch06;
            if( ch04 > badiso ) {
                badiso = ch04;
                badvtx = ivtx;
            }
        }
        pho_tkiso_badvtx_id[ipho] = badvtx;
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::FillCICInputs()
{
    pho_tkiso_recvtx_030_002_0000_10_01->clear(); pho_tkiso_recvtx_030_002_0000_10_01->resize(pho_n,std::vector<float>(vtx_std_n,0.));
  
    pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01->clear(); pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01->resize(pho_n,std::vector<float>(vtx_std_n,0.));
  
    for(int ipho=0;ipho<pho_n;++ipho){
        // TLorentzVector * phop4 = (TLorentzVector*)pho_p4->At(ipho);
        std::pair<Int_t, Float_t> worse_iso = WorstSumTrackPtInCone(ipho, 0,0, 0.40, 0.02, 0.0, 1.0, 0.1); 
        pho_tkiso_badvtx_040_002_0000_10_01[ipho] = worse_iso.second;
        pho_tkiso_badvtx_id[ipho] = worse_iso.first;
        pho_drtotk_25_99[ipho] = DeltaRToTrack(ipho, vtx_std_sel, 2.5, 99.);
    
        std::pair<Int_t, Float_t> worse_ZeeVal_iso = WorstSumTrackPtInCone(ipho, 0,0, 0.40, 0.02, 0.0, 1.0, 0.1, true); 
        pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01[ipho] = worse_iso.second;
        pho_ZeeVal_tkiso_badvtx_id[ipho] = worse_iso.first;
    
        for(int ivtx=0;ivtx<vtx_std_n;++ivtx) {
            TLorentzVector p4 = get_pho_p4( ipho, ivtx );
            (*pho_tkiso_recvtx_030_002_0000_10_01)[ipho][ivtx] = SumTrackPtInCone(&p4, ivtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
            (*pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01)[ipho][ivtx] = SumTrackPtInCone(&p4, ivtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1,true,ipho);
      
            // Need to fill the CiCpf here
            float largestIso= -1;
            if((*pho_pfiso_mycharged04)[ipho][ivtx] > largestIso){
                pho_pfiso_charged_badvtx_04[ipho]=(*pho_pfiso_mycharged04)[ipho][ivtx];
                pho_pfiso_charged_badvtx_id[ipho]=ivtx;
            }

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
  
    // PF - 4 categories
    pho_cic4pfpasscuts_lead->clear(); pho_cic4pfpasscuts_lead->resize( pho_n, std::vector<std::vector<UInt_t> >(vtx_std_n, std::vector<UInt_t>(phoNCUTLEVELS,0) ) ); 
    pho_cic4pfcutlevel_lead->clear(); pho_cic4pfcutlevel_lead->resize( pho_n, std::vector<Short_t>(vtx_std_n,0) );
    pho_cic4pfpasscuts_sublead->clear(); pho_cic4pfpasscuts_sublead->resize( pho_n, std::vector<std::vector<UInt_t> >(vtx_std_n, std::vector<UInt_t>(phoNCUTLEVELS,0) ) ); 
    pho_cic4pfcutlevel_sublead->clear(); pho_cic4pfcutlevel_sublead->resize( pho_n, std::vector<Short_t>(vtx_std_n,0) );
    std::vector<std::vector<bool> > cic4pf_passcut_lead, cic4pf_passcut_sublead;
  
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
    
            // PF - 4 categories
            int cic4pf_level_lead = PhotonCiCPFSelectionLevel(ipho, ivtx, cic4pf_passcut_lead, 4, 0);
            if(GFDEBUG) std::cout<<"cic4pf_level_lead  "<<cic4pf_level_lead<<std::endl;
            int cic4pf_level_sublead = PhotonCiCPFSelectionLevel(ipho, ivtx, cic4pf_passcut_sublead, 4, 1);
            if(GFDEBUG) std::cout<<"cic4pf_level_sublead  "<<cic4pf_level_sublead<<std::endl;
            (*pho_cic4pfcutlevel_lead)[ipho][ivtx] = cic4pf_level_lead;
            (*pho_cic4pfcutlevel_sublead)[ipho][ivtx] = cic4pf_level_sublead;
    
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
        
                // PF - 4 categories
                UInt_t cic4pf_leadw=0, cic4pf_subleadw=0;
                for(size_t icut=0; icut<cic4pf_passcut_lead.size(); ++icut) {
                    cic4pf_leadw |= ( (!cic4pf_passcut_lead[iCUTLEVEL][icut] & 0x1) << icut);
                }
                (*pho_cic4pfpasscuts_lead)[ipho][ivtx][iCUTLEVEL] = cic4pf_leadw;
                for(size_t icut=0; icut<cic4pf_passcut_sublead.size(); ++icut) {
                    cic4pf_subleadw |= ( (!cic4pf_passcut_sublead[iCUTLEVEL][icut] & 0x1) << icut);
                }
                (*pho_cic4pfpasscuts_sublead)[ipho][ivtx][iCUTLEVEL] = cic4pf_subleadw;
            }
        }
    }
}


// CiC SELECTION CODE BEGIN - SSIMON
// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cic6_allcuts_lead, float * cic6_allcuts_sublead, 
                                        float * cic4_allcuts_lead, float * cic4_allcuts_sublead,
                                        float * cic4pf_allcuts_lead, float * cic4pf_allcuts_sublead) {

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

        float cic4pf_allcuts_temp_lead[] = { 
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            1.5,         1.5,         1.5,         1.5 };
        float cic4pf_allcuts_temp_sublead[] = { 
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            1.5,         1.5,         1.5,         1.5 };
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
          
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            8.9,       6.3,       9.8,       6.8,
            43,      19.4,        24,       7.9,
            6.2,       4.3,         5,       4.3,
            0.0117,    0.0105,     0.031,     0.031,
            0.137,      0.14,     0.145,     0.143,
            0.94,      0.25,      0.93,      0.24,
            1,    0.0136,    0.0138,    0.0122,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            8.9,       6.3,       9.8,       6.8,
            43,      19.4,        24,       7.9,
            6.2,       4.3,         5,       4.3,
            0.0117,    0.0105,     0.031,     0.031,
            0.137,      0.14,     0.145,     0.143,
            0.94,      0.25,      0.93,      0.24,
            1,    0.0136,    0.0138,    0.0122,
            1.5,         1.5,         1.5,         1.5};
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            8.4,       5.4,       6.1,       6.2,
            43,       8.6,       8.9,       6.1,
            5.5,       3.8,       3.8,       3.4,
            0.0116,    0.0104,     0.031,     0.029,
            0.137,     0.103,     0.145,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.021,    0.0138,    0.0162, 
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            8.4,       5.4,       6.1,       6.2,
            43,       8.6,       8.9,       6.1,
            5.5,       3.8,       3.8,       3.4,
            0.0116,    0.0104,     0.031,     0.029,
            0.137,     0.103,     0.145,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.021,    0.0138,    0.0162, 
            1.5,         1.5,         1.5,         1.5};

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
           
        float cic4pf_allcuts_temp_lead[] = {
            6.6,       4.9,       5.6,       4.2,
            11.5,       7.2,       8.4,       5.5,
            4.4,       2.7,       3.3,       2.5,
            0.0111,    0.0103,     0.028,     0.028,
            0.124,     0.103,     0.142,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.032,     0.024,    0.0173,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            6.6,       4.9,       5.6,       4.2,
            11.5,       7.2,       8.4,       5.5,
            4.4,       2.7,       3.3,       2.5,
            0.0111,    0.0103,     0.028,     0.028,
            0.124,     0.103,     0.142,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.032,     0.024,    0.0173,
            1.5,         1.5,         1.5,         1.5};

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
           
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

            // category from PhotonCategory(..)
            // cat 0      cat 1    cat 2    cat 3 
            // high R9    low R9   high R9  low R9
            // barrel     barrel   endcap   endcap

            3.8,     2.2,     1.77,    1.29,   // isosumoet       sum of isolation cone energies divided by Et (?)   (upper cut, maximum value)
            11.7,    3.4,     3.9,     1.84,   // isosumoetbad    same as isosumoet but for 'bad' (worst ?) vertex   (upper cut, maximum value)
            3.5,     2.2,     2.3,     1.45,   // trkisooetom     tracker isolation cone only (divided by Et)        (upper cut, maximum value)
            0.0106,  0.0097,  0.028,   0.027,  // sieie           sigma(ieta,ieta)                                   (upper cut, maximum value)
            0.082,   0.062,   0.065,   0.048,  // hovere          H/E                                                (upper cut, maximum value)
            0.94,    0.36,    0.94,    0.32,   // r9              R9                                               (lower cut, minimum value) 
            1.,      0.062,   0.97,    0.97,   // drtotk_25_99    Delta R to track ?                               (lower cut, minimum value)
            1.5,     1.5,     1.5,     1.5 };  // pixel           has pixel seed ?                                   (upper cut, maximum value)
        float cic4_allcuts_temp_sublead[] = { 
            3.8,     2.2,     1.77,    1.29,
            11.7,    3.4,     3.9,     1.84,
            3.5,     2.2,     2.3,     1.45,
            0.0106,  0.0097,  0.028,   0.027,
            0.082,   0.062,   0.065,   0.048,
            0.94,    0.36,    0.94,    0.32,
            1.,      0.062,   0.97,    0.97,
            1.5,     1.5,     1.5,     1.5 };
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            6,       4.7,       5.6,       3.6,
            10,       6.5,       5.6,       4.4,
            3.8,       2.5,       3.1,       2.2,
            0.0108,    0.0102,     0.028,     0.028,
            0.124,     0.092,     0.142,     0.063,
            0.94,      0.28,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf_allcuts_temp_sublead[] = { 
            6,       4.7,       5.6,       3.6,
            10,       6.5,       5.6,       4.4,
            3.8,       2.5,       3.1,       2.2,
            0.0108,    0.0102,     0.028,     0.028,
            0.124,     0.092,     0.142,     0.063,
            0.94,      0.28,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf8tev_allcuts_temp_lead[] = {     
            6.3,       5.6,       5.8,       5.1,
            18.9,         8,        10,       6.2,
            4.5,       2.8,         4,      1.62,
            0.0125,    0.0103,     0.029,     0.028,
            0.141,     0.138,      0.12,     0.091,
            0.94,      0.33,      0.94,      0.37,
            1,     0.051,     0.054,     0.064,
            1.5,         1.5,         1.5,         1.5};
    
        float cic4pf8tev_allcuts_temp_sublead[] = {  
            6.3,       5.6,       5.8,       5.1,
            18.9,         8,        10,       6.2,
            4.5,       2.8,         4,      1.62,
            0.0125,    0.0103,     0.029,     0.028,
            0.141,     0.138,      0.12,     0.091,
            0.94,      0.33,      0.94,      0.37,
            1,     0.051,     0.054,     0.064,
            1.5,         1.5,         1.5,         1.5};
    
        float cic4pfichep_allcuts_temp_lead[] = {     
            7.7,       4.1,       1.8,       2.1,
            8,       6.3,       8.6,         4,
            5.7,       3.4,       2.3,       2.4,
            0.0191,    0.0101,     0.033,     0.025,
            0.23,      0.38,     0.168,     0.055,
            0.94,      0.35,      0.95,      0.41,    
            1,      0.31,      0.85,      0.99,
            1.5,         1.5,         1.5,         1.5};
    
        float cic4pfichep_allcuts_temp_sublead[] = {  
            7.7,       4.1,       1.8,       2.1,
            8,       6.3,       8.6,         4,
            5.7,       3.4,       2.3,       2.4,
            0.0191,    0.0101,     0.033,     0.025,
            0.23,      0.38,     0.168,     0.055,
            0.94,      0.35,      0.95,      0.41,    
            1,      0.31,      0.85,      0.99,
            1.5,         1.5,         1.5,         1.5};

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i]; 
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf8tev_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf8tev_allcuts_temp_sublead[i]; 
            }
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            5.6,       4.3,       5.5,       2.8,
            9.5,         6,       5.1,       4.3,
            3.3,       2.3,       2.2,       1.2,
            0.0107,    0.0101,     0.028,     0.028,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.38,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf_allcuts_temp_sublead[] = { 
            5.6,       4.3,       5.5,       2.8,
            9.5,         6,       5.1,       4.3,
            3.3,       2.3,       2.2,       1.2,
            0.0107,    0.0101,     0.028,     0.028,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.38,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf8tev_allcuts_temp_lead[] = { 
            5.9,       4.8,       5.3,       3.9,
            11.2,         8,        10,       5.8,
            4.2,       2.8,       2.8,      1.41,
            0.0125,    0.0101,     0.028,     0.028,
            0.141,     0.138,      0.12,     0.058,
            0.94,      0.33,      0.94,      0.39,
            1,     0.095,      0.77,       0.1,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf8tev_allcuts_temp_sublead[] = {   
            5.9,       4.8,       5.3,       3.9,
            11.2,         8,        10,       5.8,
            4.2,       2.8,       2.8,      1.41,
            0.0125,    0.0101,     0.028,     0.028,
            0.141,     0.138,      0.12,     0.058,
            0.94,      0.33,      0.94,      0.39,
            1,     0.095,      0.77,       0.1,
            1.5,         1.5,         1.5,         1.5};


        float cic4pfichep_allcuts_temp_lead[] = {     
            4.3,       2.5,      1.59,      1.53,
            7.1,       4.9,       4.3,      1.05,
            3.4,       2.4,       2.3,      1.02,
            0.019,     0.01,     0.027,     0.023,
            0.23,      0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,         0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};
        float cic4pfichep_allcuts_temp_sublead[] = {  
            4.3,       2.5,      1.59,      1.53,
            7.1,       4.9,       4.3,      1.05,
            3.4,       2.4,       2.3,      1.02,
            0.019,      0.01,     0.027,     0.023,
            0.23,      0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i]; 
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf8tev_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf8tev_allcuts_temp_sublead[i]; 
            }
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            5.1,       3.8,       3.3,       2.6,
            7.5,         5,       4.7,         4,
            3,       2.2,      1.35,      0.38,
            0.0107,      0.01,     0.028,     0.027,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,      0.99,      0.06,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            5.1,       3.8,       3.3,       2.6,
            7.5,         5,       4.7,         4,
            3,       2.2,      1.35,      0.38,
            0.0107,      0.01,     0.028,     0.027,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,      0.99,      0.06,
            1.5,         1.5,         1.5,         1.5};

        float cic4pfichep_allcuts_temp_lead[] = {     
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};    
        float cic4pfichep_allcuts_temp_sublead[] = {  
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};    


        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i]; 
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
            }
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            4.5,       2.7,       2.6,       2.4,
            5.9,       4.8,       4.6,       3.7,
            2.5,      1.48,      0.87,      0.38,
            0.0106,      0.01,     0.028,     0.027,
            0.099,     0.092,     0.104,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.78,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            4.5,       2.7,       2.6,       2.4,
            5.9,       4.8,       4.6,       3.7,
            2.5,      1.48,      0.87,      0.38,
            0.0106,      0.01,     0.028,     0.027,
            0.099,     0.092,     0.104,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.78,
            1.5,         1.5,         1.5,         1.5};

        float cic4pfichep_allcuts_temp_lead[] = {     
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};    
        float cic4pfichep_allcuts_temp_sublead[] = {  
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};    


        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i]; 
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
            }
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2.7,       2.4,       2.4,       2.3,
            5.9,       3.7,       4.6,       2.7,
            1.28,      1.15,      0.65,      0.38,
            .0106,    0.0099,     0.028,     0.027,
            0.095,     0.092,     0.065,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.7,       2.4,       2.4,       2.3,
            5.9,       3.7,       4.6,       2.7,
            1.28,      1.15,      0.65,      0.38,
            .0106,     0.0099,     0.028,     0.027,
            0.095,     0.092,     0.065,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;
        float cic4pfichep_allcuts_temp_lead[] = {
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};
        float cic4pfichep_allcuts_temp_sublead[] = {
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};


        for(int i=0;i!=ncuts*ncat_cic4;++i) {
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i];
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            }
        }
    } break;
    
    case(phoHYPERTIGHT5) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }


        float cic4pf_allcuts_temp_lead[] = { 
            2.5,       2.2,       2.3,       2.2,
            5.9,         3,       3.2,       2.7,
            1.05,      1.15,    0.0065,    0.0038,
            0.0106,    0.0099,     0.027,     0.027,
            0.086,     0.055,     0.065,    0.0032,
            0.94,      0.39,      0.95,      0.78,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.5,       2.2,       2.3,       2.2,
            5.9,         3,       3.2,       2.7,
            1.05,      1.15,    0.0065,    0.0038,
            0.0106,    0.0099,     0.027,     0.027,
            0.086,     0.055,     0.065,    0.0032,
            0.94,      0.39,      0.95,      0.78,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;


        float cic4pfichep_allcuts_temp_lead[] = {
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};
        float cic4pfichep_allcuts_temp_sublead[] = {
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};


        for(int i=0;i!=ncuts*ncat_cic4;++i) {
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i];
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            }
        }

    } break;
    
    case(phoHYPERTIGHT6) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2.3,         2,       2.2,      1.88,
            5.5,       2.2,       2.1,      1.74,
            0.96,      0.99,   6.6e-05,   3.9e-05,
            0.0106,    0.0098,     0.026,     0.027,
            0.032,     0.035,      0.03,   3.2e-05,
            0.94,      0.39,      0.96,      0.86,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.3,         2,       2.2,      1.88,
            5.5,       2.2,       2.1,      1.74,
            0.96,      0.99,   6.6e-05,   3.9e-05,
            0.0106,    0.0098,     0.026,     0.027,
            0.032,     0.035,      0.03,   3.2e-05,
            0.94,      0.39,      0.96,      0.86,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;
        float cic4pfichep_allcuts_temp_lead[] = {
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};
        float cic4pfichep_allcuts_temp_sublead[] = {
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};


        for(int i=0;i!=ncuts*ncat_cic4;++i) {
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i];
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            }
        }

    } break;
    
    case(phoHYPERTIGHT7) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2.2,      1.64,         2,      1.88,
            4.2,       2.2,       2.1,      1.39,
            0.96,      0.01,     1e-06,         0,
            0.0101,     0.009,     0.026,     0.027,
            0.017,   0.00156,     0.023,   3.2e-05,
            0.94,      0.68,      0.96,      0.86,
            1,         1,         1,      0.99,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.2,      1.64,         2,      1.88,
            4.2,       2.2,       2.1,      1.39,
            0.96,      0.01,     1e-06,         0,
            0.0101,     0.009,     0.026,     0.027,
            0.017,   0.00156,     0.023,   3.2e-05,
            0.94,      0.68,      0.96,      0.86,
            1,         1,         1,      0.99,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pfichep_allcuts_temp_lead[] = {
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};
        float cic4pfichep_allcuts_temp_sublead[] = {
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};


        for(int i=0;i!=ncuts*ncat_cic4;++i) {
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i];
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            }
        }

    } break;
    
    case(phoHYPERTIGHT8) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2,      1.62,      1.87,      1.88,
            4.2,      1.57,      1.38,      1.39,
            0.0097,  0.000101,         0,         0,
            0.0098,    0.0088,     0.025,     0.027,
            0.0164,   0.00156,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2,      1.62,      1.87,      1.88,
            4.2,      1.57,      1.38,      1.39,
            0.0097,  0.000101,         0,         0,
            0.0098,    0.0088,     0.025,     0.027,
            0.0164,   0.00156,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pfichep_allcuts_temp_lead[] = {
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};
        float cic4pfichep_allcuts_temp_sublead[] = {
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};


        for(int i=0;i!=ncuts*ncat_cic4;++i) {
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i];
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            }
        }

    } break;
    
    case(phoHYPERTIGHT9) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2,      1.47,      1.87,      1.88,
            3.1,      1.37,      1.38,      1.39,
            8e-05,     1e-06,         0,         0,
            0.0098,    0.0087,     0.024,     0.027,
            0.000165,   0.00093,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2,      1.47,      1.87,      1.88,
            3.1,      1.37,      1.38,      1.39,
            8e-05,     1e-06,         0,         0,
            0.0098,    0.0087,     0.024,     0.027,
            0.000165,   0.00093,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pfichep_allcuts_temp_lead[] = {
            4.3,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,    0.027,     0.023,
            0.23,     0.06,     0.041,     0.033,
            0.94,      0.36,      0.95,      0.45,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};
        float cic4pfichep_allcuts_temp_sublead[] = {
            4.7,       2.3,       1.4,     0.069,
            6.4,       2.9,       4.8,      0.25,
            3.3,       2.3,      1.73,      1.02,
            0.0188,      0.01,     0.028,     0.023,
            0.23,     0.169,     0.041,     0.041,
            0.94,      0.46,      0.97,      0.59,
            1,      0.32,         1,         1,
            1.5,         1.5,         1.5,         1.5};


        for(int i=0;i!=ncuts*ncat_cic4;++i) {
            if (cicVersion == "7TeV") {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            } else if (cicVersion == "ichep") {
                cic4pf_allcuts_lead[i]    = cic4pfichep_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pfichep_allcuts_temp_sublead[i];
            } else {
                cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
                cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i];
            }
        }

    } break;
    
    
    
    
    default:std::cout << "UNKNOWN phoCiCIDLevel: " << cutlevel << std::endl;

    }
}



// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::DiphotonCiCSelection( phoCiCIDLevel LEADCUTLEVEL, phoCiCIDLevel SUBLEADCUTLEVEL, 
                   Float_t leadPtMin, Float_t subleadPtMin, int ncategories, bool applyPtoverM, 
                   float *pho_energy_array, bool split, int fixedvtx, std::vector<bool> veto_indices, 
                   std::vector<int> cutsbycat) {

  //rho=0;// CAUTION SETTING RHO TO 0 FOR 2010 DATA FILES (RHO ISN'T IN THESE FILES)
  int g = -1;
  int selected_sublead_index = -1;
  float selected_lead_pt = -1;
  float selected_sublead_pt = -1;

  if( ! cutsbycat.empty() ) {
      assert( cutsbycat.size() == 4 );
      /// std::cout << "cutsbycat " << cutsbycat.size() <<std::endl;
  }
  
  std::vector<int> passing_dipho;
  std::vector<float> passing_sumpt;
  for(int idipho = 0; idipho < dipho_n; ++idipho ) {
    if( idipho >= MAX_DIPHOTONS-1 ) { 
      std::cout << "Warning diphoton index exceeds array capacity. Throwing even away " << idipho << " " << MAX_DIPHOTONS <<  dipho_n << " " << run << " " << lumis << " " << event << " " << std::endl;
      if( itype[current] == 0 ) { assert( 0 ); }
      return -1;
    }
    int ivtx = (fixedvtx==-1) ? dipho_vtxind[idipho] : fixedvtx;
    int lead = dipho_leadind[idipho];
    int sublead = dipho_subleadind[idipho];
    
    if( lead == sublead ) { continue; }

    if(veto_indices.size()!=0) {
        if(veto_indices[lead]) continue;
        if(veto_indices[sublead]) continue;
    }

    if( usePFCiC 
	&& ( ! PhotonMITPreSelection(lead, ivtx, pho_energy_array) || 
	     ! PhotonMITPreSelection(sublead, ivtx,  pho_energy_array) ) ) { continue; }
    
    TLorentzVector lead_p4 = get_pho_p4(lead,ivtx,pho_energy_array); 
    TLorentzVector sublead_p4 = get_pho_p4(sublead,ivtx,pho_energy_array); 
    
    if (sublead_p4.Pt() > lead_p4.Pt()){ // Swap them but also swap the indeces
            int tmp = lead;
            lead = sublead;
            sublead =tmp;
            dipho_leadind[idipho] = lead;
            dipho_subleadind[idipho] = sublead;
    }
    
    float leadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[lead]))->Eta());
    float subleadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[sublead]))->Eta());
    float m_gamgam = (lead_p4+sublead_p4).M();
    
    if( leadEta > 2.5 || subleadEta > 2.5 || 
	( leadEta > 1.4442 && leadEta < 1.566 ) ||
	( subleadEta > 1.4442 && subleadEta < 1.566 ) ) { continue; }
    
    float leadpt = lead_p4.Pt() > sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt();
    float subleadpt = lead_p4.Pt() < sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt();       
    // Exclusive modes cut smoothly on lead pt/M but on straight pt on sublead to save sig eff and avoid HLT turn-on  
    if(createCS_){//for control sample only cut at 33 and 25                                                                                                                    
      if(!(leadpt>=33 && subleadpt>=25)) continue;
    }else{

    if(split){   
            if ( leadpt/m_gamgam < leadPtMin/120. || subleadpt< subleadPtMin ) { continue; }  
    }else{
            if( applyPtoverM ) {
		    if ( leadpt/m_gamgam < leadPtMin/120. || subleadpt/m_gamgam < subleadPtMin/120. ||
			 leadpt < 100./3. || subleadpt < 100./4.) { continue; } 
            } else {
		    if ( leadpt < leadPtMin || subleadpt < subleadPtMin ) { continue; }
            }
    }
    }

    std::vector<std::vector<bool> > ph_passcut;
    if( ! cutsbycat.empty() ) {
            int leadCat = PhotonCategory(lead,2,2);
            int subleadCat = PhotonCategory(sublead,2,2);
            LEADCUTLEVEL    = (phoCiCIDLevel)cutsbycat[leadCat];
            SUBLEADCUTLEVEL = (phoCiCIDLevel)cutsbycat[subleadCat];
    }

    if(createCS_){
      if( (PhotonCiCSelectionLevel(lead, ivtx, ph_passcut, ncategories, 0, pho_energy_array ) > LEADCUTLEVEL && PhotonCiCSelectionLevel(sublead, ivtx, ph_passcut, ncategories, 1, pho_energy_array ) > SUBLEADCUTLEVEL) ||(PhotonCiCSelectionLevel(lead, ivtx, ph_passcut, ncategories, 0, pho_energy_array ) < LEADCUTLEVEL && PhotonCiCSelectionLevel(sublead, ivtx, ph_passcut, ncategories, 1, pho_energy_array ) < SUBLEADCUTLEVEL) ) { continue; }
  
    }else{

    if( PhotonCiCSelectionLevel(lead, ivtx, ph_passcut, ncategories, 0, pho_energy_array ) < LEADCUTLEVEL ) { continue; }
    if( PhotonCiCSelectionLevel(sublead, ivtx, ph_passcut, ncategories, 1, pho_energy_array ) < SUBLEADCUTLEVEL ) { continue; }

    }
    
    passing_dipho.push_back(idipho);
    passing_sumpt.push_back(leadpt+subleadpt);
  }
  
  if( passing_dipho.empty() ) { return -1; }
  
  //// std::vector<int> passing_dipho_places;
  //// for (int counting=0;counting<passing_dipho.size();counting++){
  //// 	  passing_dipho_places.push_back(counting); // This is very weird, but needed for later  
  //// }
  //// std::sort(passing_dipho_places.begin(),passing_dipho_places.end(),
  //// 	    SimpleSorter<float,std::greater<float> >(&passing_sumpt[0]));
  //// return passing_dipho[passing_dipho_places[0]];
  std::sort(passing_dipho.begin(),passing_dipho.end(),
	    SimpleSorter<float,std::greater<double> >(&passing_sumpt[0]));
  
  return passing_dipho[0];
}


int LoopAll::DiphotonMITPreSelection(const char * type, Float_t leadPtMin, Float_t subleadPtMin, Float_t phoidMvaCut, bool applyPtoverM, float *pho_energy_array, bool vetodipho, bool kinonly, float dipho_BDT_Cut,int fixedvtx, bool split, std::vector<bool> veto_indices) {

    //rho=0;// CAUTION SETTING RHO TO 0 FOR 2010 DATA FILES (RHO ISN'T IN THESE FILES)
    int selected_lead_index = -1;
    int selected_sublead_index = -1;
    float selected_lead_pt = -1;
    float selected_sublead_pt = -1;
    
    std::vector<int> passing_dipho;
    std::vector<float> passing_sumpt;
    for(int idipho = 0; idipho < dipho_n; ++idipho ) {
        if( idipho >= MAX_DIPHOTONS-1 ) { 
            std::cout << "Warning diphoton index exceeds array capacity. Throwing event away " << idipho << " " << MAX_DIPHOTONS <<  dipho_n << " " << run << " " << lumis << " " << event << " " << std::endl;
            if( itype[current] == 0 ) { assert( 0 ); }
            return -1;
        }
        
        if(vetodipho && dipho_sel[idipho]!=true) continue;
        if(dipho_BDT[idipho]<dipho_BDT_Cut) continue;

        float sumpt = DiphotonMITPreSelectionPerDipho(type, idipho, leadPtMin, subleadPtMin, phoidMvaCut, applyPtoverM, pho_energy_array, fixedvtx, split, kinonly, veto_indices);

        if(sumpt!=-99){
            passing_dipho.push_back(idipho);
            passing_sumpt.push_back(sumpt); // need to use reordered pt!
        }
    }
  
    if( passing_dipho.empty() ) { return -1; }

    std::vector<int> passing_dipho_places;
    for (int counting=0;counting<passing_dipho.size();counting++){
        passing_dipho_places.push_back(counting); // This is very weird, but needed for later  
    }
    std::sort(passing_dipho_places.begin(),passing_dipho_places.end(),
              SimpleSorter<float,std::greater<float> >(&passing_sumpt[0]));

    int selected_dipho_ind     = passing_dipho[passing_dipho_places[0]];
    
    // MOVED INTO DiphotonMITPreSelectionPerDipho
    //int selected_dipho_lead    = dipho_leadind[selected_dipho_ind];
    //int selected_dipho_sublead = dipho_subleadind[selected_dipho_ind];
    //int selected_dipho_vtx     = dipho_vtxind[selected_dipho_ind];
    //TLorentzVector selected_lead_p4 = get_pho_p4(selected_dipho_lead,selected_dipho_vtx,pho_energy_array); 
    //TLorentzVector selected_sublead_p4 = get_pho_p4(selected_dipho_sublead,selected_dipho_vtx,pho_energy_array);
    //
    //if (!kinonly) {
    //    if( version >= 13 ) {
    //        if ( photonIDMVANew(selected_dipho_lead,selected_dipho_vtx,selected_lead_p4,"MIT") <= phoidMvaCut
    //            || photonIDMVANew(selected_dipho_sublead,selected_dipho_vtx,selected_sublead_p4,"MIT")  <= phoidMvaCut
	  //            ) {return -1;}
    //    } else {
	  //        if ( photonIDMVA(selected_dipho_lead,selected_dipho_vtx,selected_lead_p4,"MIT") <= phoidMvaCut
	  //            || photonIDMVA(selected_dipho_sublead,selected_dipho_vtx,selected_sublead_p4,"MIT") <= phoidMvaCut
	  //            ) {return -1;}
    //    }
    //}
    
    return selected_dipho_ind;
}

float LoopAll::DiphotonMITPreSelectionPerDipho(const char * type, int idipho, Float_t leadPtMin, Float_t subleadPtMin, Float_t phoidMvaCut, bool applyPtoverM, float *pho_energy_array, int fixedvtx, bool split, bool kinonly, std::vector<bool> veto_indices) {
    
    int ivtx = (fixedvtx==-1) ? dipho_vtxind[idipho] : fixedvtx;
    int lead = dipho_leadind[idipho];
    int sublead = dipho_subleadind[idipho];
    
    if( lead == sublead ) { return -99; }
    
    if(veto_indices.size()!=0) {
        if(veto_indices[lead]) return -99;
        if(veto_indices[sublead]) return -99;
    }
    
    TLorentzVector lead_p4 = get_pho_p4(lead,ivtx,pho_energy_array); 
    TLorentzVector sublead_p4 = get_pho_p4(sublead,ivtx,pho_energy_array);
    
    if (sublead_p4.Pt() > lead_p4.Pt()){ // Swap them but also swap the indeces
        int tmp = lead;
        lead = sublead;
        sublead =tmp;
        dipho_leadind[idipho] = lead;
        dipho_subleadind[idipho] = sublead;
    }
    
    float leadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[lead]))->Eta());
    float subleadEta = fabs(((TVector3 *)sc_xyz->At(pho_scind[sublead]))->Eta());
    float m_gamgam = (lead_p4+sublead_p4).M();
    
    if( leadEta > 2.5 || subleadEta > 2.5 || 
        ( leadEta > 1.4442 && leadEta < 1.566 ) ||
        ( subleadEta > 1.4442 && subleadEta < 1.566 ) ) { return -99; }
    
    float leadpt = lead_p4.Pt() > sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt();
    float subleadpt = lead_p4.Pt() < sublead_p4.Pt() ? lead_p4.Pt() : sublead_p4.Pt();     
    
    if( split ) {
        if ( leadpt/m_gamgam < leadPtMin/120. || subleadpt < 25. ||
             leadpt < 100./3. ) { return -99; } 
    } else if( applyPtoverM ) {
        if ( leadpt/m_gamgam < leadPtMin/120. || subleadpt/m_gamgam < subleadPtMin/120. ||
             leadpt < 100./3. || subleadpt < 100./4.) { return -99; } 
    } else {
        if ( leadpt < leadPtMin || subleadpt < subleadPtMin ) { return -99; }
    }
    
    if ( runZeeValidation && m_gamgam<70. ) { return -99; }
    
    std::vector<std::vector<bool> > ph_passcut;
    if (!kinonly) {
        if( version >= 13 ) {
            if (!( PhotonMITPreSelection(lead, ivtx, pho_energy_array ) && PhotonMITPreSelection(sublead, ivtx,  pho_energy_array ))) return -99; 
        } else {
            if (!( PhotonMITPreSelection2011(lead, ivtx, pho_energy_array ) && PhotonMITPreSelection2011(sublead, ivtx,  pho_energy_array ))) return -99; 
        }
    }

    if (!kinonly) {
	    /// std::cout << "Diphoton preselection " << dipho_leadind[idipho] << " " << dipho_subleadind[idipho] << " " 
	    /// 	      << photonIDMVA(lead,ivtx,lead_p4,type) << " " << photonIDMVA(sublead,ivtx,sublead_p4,type) << std::endl;
	if ( photonIDMVA(lead,ivtx,lead_p4,type) <= phoidMvaCut
	     || photonIDMVA(sublead,ivtx,sublead_p4,type)  <= phoidMvaCut
	    ) {return -99;}
    }
    
    return (leadpt+subleadpt);
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
    //float mitCuts_hcalecal[4]            = {3,3,3,3};                                                        
    //float mitCuts_abstrkiso[4]           = {2.8,2.8,2.8,2.8};                                                
    //float mitCuts_trkiso_hollow03[4]     = {4,4,4,4};                                                       
    //float mitCuts_drtotk_25_99[4]  = {0.26,0.029,0.0062,0.0055};
    float mitCuts_pfiso[4]               = {4,4,4,4}; // WARN if depends on category should change below

    TLorentzVector phop4 = get_pho_p4( photon_index, vertex_index, pho_energy_array);                      
    //TLorentzVector phop4_badvtx = get_pho_p4( photon_index, pho_tkiso_badvtx_id[photon_index], pho_energy_array  );

    //float rhofac=0.17;
    float val_hoe        = pho_hoe[photon_index];
    float val_sieie      = pho_sieie[photon_index];                                                          
    float val_ecaliso = pho_ecalsumetconedr03[photon_index] - 0.012*phop4.Et();                              
    float val_hcaliso = pho_hcalsumetconedr03[photon_index] - 0.005*phop4.Et(); 
    float val_trkiso  = pho_trksumpthollowconedr03[photon_index] - 0.002*phop4.Et();                          \

    //float val_hcalecal   = (pho_ecalsumetconedr03[photon_index]+pho_hcalsumetconedr03[photon_index]-rho_algo1*rhofac);                                             
    //float val_abstrkiso  = (*pho_tkiso_recvtx_030_002_0000_10_01)[photon_index][vertex_index];                
    //float val_trkiso_hollow03 = pho_trksumpthollowconedr03[photon_index];                                    
    //float val_drtotk_25_99 = pho_drtotk_25_99[photon_index];
    int val_pho_isconv = pho_isconv[photon_index];
    float val_pfiso02 = (*pho_pfiso_mycharged02)[photon_index][vertex_index];

    /*
      if (run==170397 && lumis==279 &&event==304405242){
     
      std::cout << "rho " << rho <<std::endl;
      std::cout << "pho_n " << pho_n <<std::endl;
      std::cout << "pho_index " << photon_index <<std::endl;
      std::cout << "pho_et " << phop4.Et() <<std::endl;
      std::cout << "hoe " << val_hoe <<std::endl;
      std::cout << "sieie " << val_sieie <<std::endl;
      std::cout << "ecaliso " << val_ecaliso <<std::endl;
      std::cout << "hcaliso " << val_hcaliso <<std::endl;
      std::cout << "hcalecal " << val_hcalecal <<std::endl;
      std::cout << "abstrkiso " << val_abstrkiso <<std::endl;
      std::cout << "trkiso " << val_trkiso_hollow03 <<std::endl;
      std::cout << "isconv " << val_pho_isconv <<std::endl;
     
      std::cout << "r9 " << pho_r9[photon_index] <<std::endl;
      std::cout << "r9 categoty " << r9_category <<std::endl;
      std::cout << "eta scxyz" << fabs(((TVector3*)sc_xyz->At(pho_scind[photon_index]))->Eta()) <<std::endl;
      std::cout << "category" << photon_category <<std::endl;
      }
    */
  
    // can't apply cuts in categories at reduction level because of the shape rescaling
    if( itype[current] == 0 || typerun != kReduce ) {
        if (val_hoe             >= mitCuts_hoe[photon_category]         ) return false;                                           
        if (val_sieie           >= mitCuts_sieie[photon_category]       ) return false;
        if (val_hcaliso         >= mitCuts_hcaliso[photon_category]     ) return false;                                           
        if (val_trkiso          >= mitCuts_trkiso[photon_category]      ) return false;
        //if (val_hcalecal        >= mitCuts_hcalecal[photon_category]    ) return false;
        //if (val_abstrkiso       >= mitCuts_abstrkiso[photon_category]   ) return false;                   
        //if (val_trkiso_hollow03 >= mitCuts_trkiso_hollow03[photon_category]) return false;                                        
    }
    if( typerun != kReduce ) {
        // if (val_drtotk_25_99    <  mitCuts_drtotk_25_99[photon_category]   ) return false; // Electron Rejection based on CiC for now
        if ((!val_pho_isconv && !runZeeValidation) || (runZeeValidation && val_pho_isconv) ) return false; // Electron Rejection based Conversion Safe Veto
    }
   
    // this does not depend on R9
    if( typerun != kReduce ) {
        if (val_pfiso02 >= mitCuts_pfiso[photon_category]) return false;            
	if (applyEcalIsoPresel && (val_ecaliso         >= mitCuts_ecaliso[photon_category])     ) return false;

    }
   
    return true;
}

// Define newfunction to calculate MIT (Pre-)Selection                                                      
bool LoopAll::PhotonMITPreSelection2011( int photon_index, int vertex_index, float *pho_energy_array ) 
{
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
   int   val_pho_isconv = pho_isconv[photon_index];


   if (val_hoe             >= mitCuts_hoe[photon_category]         ) return false;                                           
   if (val_sieie           >= mitCuts_sieie[photon_category]       ) return false;
   if (val_ecaliso         >= mitCuts_ecaliso[photon_category]     ) return false;
   if (val_hcaliso         >= mitCuts_hcaliso[photon_category]     ) return false;                                           
   if (val_trkiso          >= mitCuts_trkiso[photon_category]      ) return false;
   if (val_hcalecal        >= mitCuts_hcalecal[photon_category]    ) return false;
   if (val_abstrkiso       >= mitCuts_abstrkiso[photon_category]   ) return false;                   
   if ((!val_pho_isconv && !runZeeValidation) || (runZeeValidation && val_pho_isconv) ) return false; // Electron Rejection based Conversion Safe Veto
   if (val_trkiso_hollow03 >= mitCuts_trkiso_hollow03[photon_category]) return false;                                        

   return true;

}

// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::PhotonCiCPFSelectionLevel( int photon_index, int vertex_index, std::vector<std::vector<bool> > & ph_passcut, int ncategories, 
                                        int doSublead, float *pho_energy_array ) {
    assert( version >= 13 );
    if( ! runCiC ) {
        switch(ncategories) {
        case (4):
            return doSublead ? (*pho_cic4pfcutlevel_sublead)[photon_index][vertex_index] : (*pho_cic4pfcutlevel_lead)[photon_index][vertex_index] ;
        default:
            std::cout << "CiCPF is available only in four categories." << std::endl;
            return -1;
        }
    }
  
    int cutlevelpassed = -1;

    int n_r9_categories = -1;
    int n_eta_categories = -1;
    if(ncategories == 4) {
        n_r9_categories = 2;
        n_eta_categories = 2;
    } else {
        std::cout << "UNKNOWN ncategories must be 4, not " << ncategories << std::endl;
    } 
    int photon_category = PhotonCategory(photon_index,n_r9_categories,n_eta_categories);

    TLorentzVector phop4 = get_pho_p4( photon_index, vertex_index, pho_energy_array  );
    TLorentzVector phop4_badvtx = get_pho_p4( photon_index, pho_tkiso_badvtx_id[photon_index], pho_energy_array  );

    // quantities which are used for the cuts

    float val_tkisobad = -99;
    for(int iv=0; iv < vtx_std_n; iv++) {
        if((*pho_pfiso_mycharged04)[photon_index][iv] > val_tkisobad) {
            val_tkisobad = (*pho_pfiso_mycharged04)[photon_index][iv];
        }
    }

    float val_tkiso        = (*pho_pfiso_mycharged03)[photon_index][vertex_index];
    float val_ecaliso      = pho_pfiso_myphoton03[photon_index];
    float val_ecalisobad   = pho_pfiso_myphoton04[photon_index];
    float val_sieie        = pho_sieie[photon_index];
    float val_hoe          = pho_hoe[photon_index];
    float val_r9           = pho_r9_cic[photon_index];
    float val_conv         = pho_isconv[photon_index];

    float rhofacbad=0.23, rhofac=0.09;

    float val_isosumoet    = (val_tkiso    + val_ecaliso    + pfisoOffset - rho_algo1 * rhofac )   * 50. / phop4.Et();
    float val_isosumoetbad = (val_tkisobad + val_ecalisobad + pfisoOffset - rho_algo1 * rhofacbad) * 50. / phop4_badvtx.Et();

    // tracker isolation cone energy divided by Et
    float val_trkisooet    = (val_tkiso) * 50. / phop4.Et();

    ph_passcut.clear();
    ph_passcut.resize(phoNCUTLEVELS,std::vector<bool>(8,true) );
    if(!doSublead) {
        for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {
            switch(ncategories) {
            case(4) :
                if(GFDEBUG) cout<<"iCUTLEVEL phoNCUTLEVELS "<<iCUTLEVEL<<" "<<phoNCUTLEVELS<<endl;
                ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4pf_cut_lead_isosumoet[iCUTLEVEL][photon_category]     );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][0] <<"   "<< val_isosumoet        <<"   "<<   cic4pf_cut_lead_isosumoet[iCUTLEVEL][photon_category]  <<std::endl;
                ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4pf_cut_lead_isosumoetbad[iCUTLEVEL][photon_category]  );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][1] <<"   "<< val_isosumoetbad     <<"   "<<   cic4pf_cut_lead_isosumoetbad[iCUTLEVEL][photon_category]<<std::endl;
                ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4pf_cut_lead_trkisooet[iCUTLEVEL][photon_category]     );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][2] <<"   "<< val_trkisooet        <<"   "<<   cic4pf_cut_lead_trkisooet[iCUTLEVEL][photon_category]  <<std::endl;
                ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4pf_cut_lead_sieie[iCUTLEVEL][photon_category]         );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][3] <<"   "<< val_sieie            <<"   "<<   cic4pf_cut_lead_sieie[iCUTLEVEL][photon_category]      <<std::endl;
                ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4pf_cut_lead_hovere[iCUTLEVEL][photon_category]        );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][4] <<"   "<< val_hoe              <<"   "<<   cic4pf_cut_lead_hovere[iCUTLEVEL][photon_category]     <<std::endl;
                ph_passcut[iCUTLEVEL][5] = (val_r9               >=   cic4pf_cut_lead_r9[iCUTLEVEL][photon_category]            );// gt cut
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][5] <<"   "<< val_r9               <<"   "<<   cic4pf_cut_lead_r9[iCUTLEVEL][photon_category]         <<std::endl;
                if(GFDEBUG) cout<<"val_conv:  "<<val_conv<<endl;
                if( runZeeValidation ) { 
                    ph_passcut[iCUTLEVEL][6] = !val_conv;
                } else {
                    ph_passcut[iCUTLEVEL][6] = val_conv;
                }
                break;
            }
            bool ph_passcut_all = true;
            for(int icut=0;icut!=8;++icut) {
                ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
            }
            if(ph_passcut_all) {
                if( cutlevelpassed != iCUTLEVEL - 1 && ( cicVersion != "ichep" || cutlevelpassed >4 ) ) {
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
            case(4) :
                if(GFDEBUG) cout<<"iCUTLEVEL phoNCUTLEVELS "<<iCUTLEVEL<<" "<<phoNCUTLEVELS<<endl;
                ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4pf_cut_sublead_isosumoet[iCUTLEVEL][photon_category]     );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][0] <<"   "<< val_isosumoet        <<"   "<<   cic4pf_cut_sublead_isosumoet[iCUTLEVEL][photon_category]  <<std::endl;
                ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4pf_cut_sublead_isosumoetbad[iCUTLEVEL][photon_category]  );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][1] <<"   "<< val_isosumoetbad     <<"   "<<   cic4pf_cut_sublead_isosumoetbad[iCUTLEVEL][photon_category]<<std::endl;
                ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4pf_cut_sublead_trkisooet[iCUTLEVEL][photon_category]     );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][2] <<"   "<< val_trkisooet        <<"   "<<   cic4pf_cut_sublead_trkisooet[iCUTLEVEL][photon_category]  <<std::endl;
                ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4pf_cut_sublead_sieie[iCUTLEVEL][photon_category]         );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][3] <<"   "<< val_sieie            <<"   "<<   cic4pf_cut_sublead_sieie[iCUTLEVEL][photon_category]      <<std::endl;
                ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4pf_cut_sublead_hovere[iCUTLEVEL][photon_category]        );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][4] <<"   "<< val_hoe              <<"   "<<   cic4pf_cut_sublead_hovere[iCUTLEVEL][photon_category]     <<std::endl;
                ph_passcut[iCUTLEVEL][5] = (val_r9               >=   cic4pf_cut_sublead_r9[iCUTLEVEL][photon_category]            );
                if(GFDEBUG) cout<<"pass, val, cut:  "<<ph_passcut[iCUTLEVEL][5] <<"   "<< val_r9               <<"   "<<   cic4pf_cut_sublead_r9[iCUTLEVEL][photon_category]         <<std::endl;
                if(GFDEBUG) cout<<"val_conv:  "<<val_conv<<endl;
                if( runZeeValidation ) { 
                    ph_passcut[iCUTLEVEL][6] = !val_conv;
                } else {
                    ph_passcut[iCUTLEVEL][6] = val_conv;
                }
                break;
            }
            bool ph_passcut_all = true;
            for(int icut=0;icut!=8;++icut) {
                ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
            }
            if(ph_passcut_all) {
                if( cutlevelpassed != iCUTLEVEL - 1  && ( cicVersion != "ichep" || cutlevelpassed >4 ) ) {
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
int LoopAll::PhotonCiCSelectionLevel( int photon_index, int vertex_index, std::vector<std::vector<bool> > & ph_passcut, int ncategories, 
                                      int doSublead, float *pho_energy_array ) {
  
    if( usePFCiC && ncategories==4 ) {
        return PhotonCiCPFSelectionLevel( photon_index, vertex_index, ph_passcut, ncategories, doSublead, pho_energy_array );
    }
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

    // quantities which are used for the cuts
    float val_tkiso        = (*pho_tkiso_recvtx_030_002_0000_10_01)[photon_index][vertex_index];
    float val_ecaliso      = pho_ecalsumetconedr03[photon_index];
    float val_hcaliso      = pho_hcalsumetconedr04[photon_index];
    float val_ecalisobad   = pho_ecalsumetconedr04[photon_index];
    float val_hcalisobad   = pho_hcalsumetconedr04[photon_index];
    float val_tkisobad     = pho_tkiso_badvtx_040_002_0000_10_01[photon_index];
    float val_sieie        = pho_sieie[photon_index];
    float val_hoe          = pho_hoe[photon_index];
    float val_r9           = pho_r9[photon_index];
    float val_drtotk_25_99 = pho_drtotk_25_99[photon_index];
    float val_pixel        = (float)pho_haspixseed[photon_index];
    int val_pho_isconv     = pho_isconv[photon_index];
  
    float isosumconst = 5.;
    float isosumconstbad = 7.;
    if(ncategories==4) {
        isosumconst = 0.;
        isosumconstbad = 0.;
    }
  
    /// float rhofacbad=0.40, rhofac=0.05;
    float rhofacbad=0.52, rhofac=0.17;

    // isolation cone energies divided by Et
    //                        tracker iso    ecal iso         hcal iso         offset ?!        mean energy
    //                                                                                          per unit
    //                                                                                          area in the event
    if( version >= 13 ) rho = rho_algo1;
    float val_isosumoet    = (val_tkiso    + val_ecaliso    + val_hcaliso    + isosumconst    - rho * rhofac )   * 50. / phop4.Et();
    float val_isosumoetbad = (val_tkisobad + val_ecalisobad + val_hcalisobad + isosumconstbad - rho * rhofacbad) * 50. / phop4_badvtx.Et();

    // tracker isolation cone energy divided by Et
    float val_trkisooet    = (val_tkiso) * 50. / phop4.Et();

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
                if( runZeeValidation ) { 
                    ph_passcut[iCUTLEVEL][6] = ! val_pho_isconv;
                } else {
                    ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic6_cut_lead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
                    ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic6_cut_lead_pixel[iCUTLEVEL][photon_category]         );
                }
                break;
            case(4) :
                ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4_cut_lead_isosumoet[iCUTLEVEL][photon_category]     );
                ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4_cut_lead_isosumoetbad[iCUTLEVEL][photon_category]  );
                ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4_cut_lead_trkisooet[iCUTLEVEL][photon_category]     );
                ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4_cut_lead_sieie[iCUTLEVEL][photon_category]         );
                ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4_cut_lead_hovere[iCUTLEVEL][photon_category]        );
                ph_passcut[iCUTLEVEL][5] = (val_r9             >=     cic4_cut_lead_r9[iCUTLEVEL][photon_category]            );// gt cut
                if( runZeeValidation ) { 
                    ph_passcut[iCUTLEVEL][6] = ! val_pho_isconv;
                } else {
                    ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic4_cut_lead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
                    ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic4_cut_lead_pixel[iCUTLEVEL][photon_category]         );
                }
                break;
            }
            bool ph_passcut_all = true;
            for(int icut=0;icut!=8;++icut) {
                ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
            }
            if(ph_passcut_all) {
                if( cutlevelpassed != iCUTLEVEL - 1 && ! runZeeValidation) {
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
                if( runZeeValidation ) { 
                    ph_passcut[iCUTLEVEL][6] = ! val_pho_isconv;
                } else {
                    ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic6_cut_sublead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
                    ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic6_cut_sublead_pixel[iCUTLEVEL][photon_category]         );
                }
                break;
            case(4) :
                ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4_cut_sublead_isosumoet[iCUTLEVEL][photon_category]     );
                ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4_cut_sublead_isosumoetbad[iCUTLEVEL][photon_category]  );
                ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4_cut_sublead_trkisooet[iCUTLEVEL][photon_category]     );
                ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4_cut_sublead_sieie[iCUTLEVEL][photon_category]         );
                ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4_cut_sublead_hovere[iCUTLEVEL][photon_category]        );
                ph_passcut[iCUTLEVEL][5] = (val_r9             >=     cic4_cut_sublead_r9[iCUTLEVEL][photon_category]            );// gt cut
                if( runZeeValidation ) { 
                    ph_passcut[iCUTLEVEL][6] = ! val_pho_isconv;
                } else {
                    ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >=     cic4_cut_sublead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
                    ph_passcut[iCUTLEVEL][7] = (val_pixel            <=   cic4_cut_sublead_pixel[iCUTLEVEL][photon_category]         );
                }
                break;
            }
            bool ph_passcut_all = true;
            for(int icut=0;icut!=8;++icut) {
                ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
            }
            if(ph_passcut_all) {
                if( cutlevelpassed != iCUTLEVEL - 1 && ! runZeeValidation ) {
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
    return ecalhitsSumEt;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
std::pair<Int_t, Float_t> LoopAll::WorstSumTrackPtInCone(Int_t ipho, Int_t returnvtxind, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax, bool Zee_validation) {

    Int_t worstvtxind = -1;
    Float_t maxisosum = -100;
    for(int ivtx=0;ivtx!=vtx_std_n;++ivtx) {
        TLorentzVector photon_p4 = get_pho_p4( ipho, ivtx );
        Float_t thisvtxisosum = SumTrackPtInCone(&photon_p4, ivtx, PtMin, OuterConeRadius, InnerConeRadius, EtaStripHalfWidth, dzmax, dxymax, Zee_validation, ipho);
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
Float_t LoopAll::SumTrackPtInCone(TLorentzVector *photon_p4, Int_t vtxind, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax, bool Zee_validation, Int_t pho_ind) {
    // TRACKER Isolation
    if(vtxind<0)return -99;
    TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(vtxind);
    float SumTrackPt=0;
  
    int electronMatch =-1;

    if(Zee_validation && pho_ind!=-1) {
        for(int iel=0; iel<el_std_n; iel++){
            if(el_std_scind[iel]==pho_scind[pho_ind]) {
                electronMatch=iel;
                break;
            }
        }
    }
  
    for(unsigned int itk=0; itk!=tk_n; itk++) {
        // remove electron track for Zee validation when there is a match
        if(Zee_validation && !pho_isconv[pho_ind] && electronMatch!=-1) {
            if(itk==el_std_tkind[electronMatch]) continue;
        }  

        TLorentzVector * tkp4= (TLorentzVector *) tk_p4->At(itk);
        if(tkp4->Pt() < PtMin)continue;
	
	if( itk >= tk_vtx_pos->GetEntries() ) {
		std::cout << "WARNING: track position at vertex not available " << vtxind << " " << itk << " " << tkp4->Pt() << std::endl;
	} else {
		TVector3 * tkpos= (TVector3 *) tk_vtx_pos->At(itk);
		/// double deltaz = fabs(vtxpos->Z() - tkpos->Z()); 
		double deltaz = fabs( (tkpos->Z()-vtxpos->Z()) - ( (tkpos->X()-vtxpos->X())*tkp4->Px() + (tkpos->Y()-vtxpos->Y())*tkp4->Py() )/tkp4->Pt() * tkp4->Pz()/tkp4->Pt() );
		if(deltaz > dzmax)continue;
		
		double dxy = ( -(tkpos->X() - vtxpos->X())*tkp4->Py() + (tkpos->Y() - vtxpos->Y())*tkp4->Px()) / tkp4->Pt();
		if(fabs(dxy) > dxymax)continue;
	}
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


void LoopAll::getIetaIPhi(int phoid, int & ieta, int & iphi ) const
{
    TLorentzVector *bcpos   = (TLorentzVector*)bc_p4->At(sc_bcseedind[pho_scind[phoid]]);
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
    //// int detid = ecalhit_detid[bc_seed[sc_bcseedind[pho_scind[phoid]]]];
    ieta  = (detid>>9)&0x7F; 
    iphi  = detid&0x1FF; 
}

bool LoopAll::CheckSphericalPhoton(int ieta, int iphi) const 
{
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

bool LoopAll::CheckSphericalPhoton(int phoid) const 
{

    TVector3 *phoCalo = (TVector3*)sc_xyz->At(pho_scind[phoid]);
    if (pho_r9[phoid]<0.94 || fabs(phoCalo->Eta())>1.) return false;

    int ieta, iphi;
    getIetaIPhi(phoid,ieta,iphi);
  
    return CheckSphericalPhoton(ieta,iphi);

}
// CiC SELECTION CODE END - SSIMON

// Functions moved from Tools.h
// ---------------------------------------------------------------------------------------------------------------------------------------------
double LoopAll::DeltaPhi(double phi1, double phi2) {
    double deltaphi;
    if(phi1<0) phi1+=2*TMath::Pi();
    if(phi2<0) phi2+=2*TMath::Pi();
    deltaphi=fabs(phi1-phi2);
    if(deltaphi>2*TMath::Pi()) deltaphi-=2*TMath::Pi();
    if(deltaphi>TMath::Pi()) deltaphi=2*TMath::Pi()-deltaphi;
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
//  BRANCH_DICT(METcorrected); //met at analysis step
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
  
    BRANCH_DICT(vtx_std_sel);
    BRANCH_DICT(vtx_std_ranked_list);
    BRANCH_DICT(vtx_std_evt_mva);
  
    BRANCH_DICT(pho_mitmva);
 
    BRANCH_DICT(pho_tkiso_recvtx_030_002_0000_10_01);
    BRANCH_DICT(pho_tkiso_badvtx_040_002_0000_10_01);
    BRANCH_DICT(pho_pfiso_charged_badvtx_04);
    BRANCH_DICT(pho_pfiso_charged_badvtx_id);
    BRANCH_DICT(pho_tkiso_badvtx_id);
    BRANCH_DICT(pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01);
    BRANCH_DICT(pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01);
    BRANCH_DICT(pho_ZeeVal_tkiso_badvtx_id);
    BRANCH_DICT(pho_drtotk_25_99);

    BRANCH_DICT(pho_cic6cutlevel_lead);
    BRANCH_DICT(pho_cic6passcuts_lead);
    BRANCH_DICT(pho_cic6cutlevel_sublead);
    BRANCH_DICT(pho_cic6passcuts_sublead);

    BRANCH_DICT(pho_cic4cutlevel_lead);
    BRANCH_DICT(pho_cic4passcuts_lead);
    BRANCH_DICT(pho_cic4cutlevel_sublead);
    BRANCH_DICT(pho_cic4passcuts_sublead);

    BRANCH_DICT(pho_cic4pfcutlevel_lead);
    BRANCH_DICT(pho_cic4pfpasscuts_lead);
    BRANCH_DICT(pho_cic4pfcutlevel_sublead);
    BRANCH_DICT(pho_cic4pfpasscuts_sublead);

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


  //correctMETinRED
  BRANCH_DICT(shiftMET_pt);
  BRANCH_DICT(shiftMET_phi);
  BRANCH_DICT(smearMET_pt);
  BRANCH_DICT(smearMET_phi);
  BRANCH_DICT(shiftscaleMET_pt);
  BRANCH_DICT(shiftscaleMET_phi);
  BRANCH_DICT(shiftsmearMET_pt);
  BRANCH_DICT(shiftsmearMET_phi);
  BRANCH_DICT(correctedpfMET);
  BRANCH_DICT(correctedpfMET_phi);

  BRANCH_DICT(shiftMET_e);
  BRANCH_DICT(shiftMET_eta);
  BRANCH_DICT(shiftscaleMET_e);
  BRANCH_DICT(shiftscaleMET_eta);
  
  BRANCH_DICT(mc_et);
  BRANCH_DICT(mc_eta);
  BRANCH_DICT(mc_phi);
  BRANCH_DICT(fsr_et);
  BRANCH_DICT(fsr_eta);
  BRANCH_DICT(fsr_phi);
  BRANCH_DICT(higgs);

#endif
}

int  LoopAll::RescaleJetEnergy(bool force) {
    static int nprints = 1;
    if( version >= 13 && ! force ) {
        if( nprints-- > 0 ) {
            std::cout << "Data format version is " << version << " jets are already corrected " << std::endl;
        }
        return 1;
    }
    for (int i = 0; i<jet_algoPF1_n; i++) {
        TLorentzVector * thisjet = (TLorentzVector *) jet_algoPF1_p4->At(i);
        *thisjet*=jet_algoPF1_erescale[i];
    }
    return 1;
}
   


void LoopAll::PhotonsToVeto(TLorentzVector* veto_p4, float drtoveto, std::vector<bool>& vetos, bool drtotkveto,float drgsftoveto){
    vetos.clear();

    for(int ipho=0; ipho<pho_n; ipho++){
        TLorentzVector* photon = (TLorentzVector*) pho_p4->At(ipho);
        if(GFDEBUG) std::cout<<"ipho "<<ipho<<std::endl;
        if(GFDEBUG) std::cout<<"eta "<<photon->Eta()<<std::endl;
        if(GFDEBUG) std::cout<<"dr "<<photon->DeltaR(*veto_p4)<<std::endl;
        if(GFDEBUG) std::cout<<"drtotk "<<pho_drtotk_25_99[ipho]<<std::endl;
        if(photon->DeltaR(*veto_p4)<drtoveto) {
            vetos.push_back(true);
        } else if (drtotkveto  &&  pho_drtotk_25_99[ipho]<drgsftoveto) {
            vetos.push_back(true);
        } else {
            vetos.push_back(false);
	    float drLepPho=photon->DeltaR(*veto_p4);
	    float drGsf=pho_drtotk_25_99[ipho];
        }
    }

}

std::pair<int, int> LoopAll::Select2HighestPtJets(TLorentzVector& leadpho, TLorentzVector& subleadpho, Bool_t * jetid_flags)
{
    std::pair<int, int> myJets(-1,-1);
    std::pair<int, int> fail(-1,-1);

    std::pair<float, float> myJetspt(-1.,-1.);
    if(GFDEBUG) std::cout<<"Select2HighestPtJets--start"<<std::endl;
    float dr2pho = 0.5;
    float dr2jet = 0.5;

    TLorentzVector* j1p4;
    TLorentzVector* j2p4;
    float j1pt=-1;
    float j2pt=-1;

    // select highest pt jets

    for(int j1_i=0; j1_i<jet_algoPF1_n; j1_i++){
        j1p4 = (TLorentzVector*) jet_algoPF1_p4->At(j1_i);
        if(GFDEBUG) std::cout<<"jet "<<j1_i<<std::endl;
        if(GFDEBUG) std::cout<<"passing pu id?"<<std::endl;
        if(jetid_flags != 0 && !jetid_flags[j1_i]) continue; 
        if(GFDEBUG) std::cout<<"within eta 4.7?"<<std::endl;
        if(fabs(j1p4->Eta()) > 4.7) continue;
        if(GFDEBUG) std::cout<<"close to leadpho?"<<std::endl;
        if(j1p4->DeltaR(leadpho) < dr2pho) continue;
        if(GFDEBUG) std::cout<<"close to subleadpho?"<<std::endl;
        if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
        j1pt=j1p4->Pt();
        if(GFDEBUG) std::cout<<"passing all single jet requirements with pt"<<j1pt<<std::endl;


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
    passingMu++;mymu = indmu;
  }
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

        //std::cout << setprecision(4) << "Run = " << run << "  LS = " << lumis << "  Event = " << event << "  SelVtx = " << vtxind << " elEta = " << thiseta << "  elPhi = " << thisel->Phi() <<  "  elEt = " << thisel->Et() << endl;

        myel = indel;
    }

    //if(passingEl>1) std::cout<<"There are "<<passingEl<<" passing electrons!!"<<std::endl;


    /////////////////fabs(elD0Vtx[i][vtx]) < 0.02 &&       
    /////////////////fabs(elDzVtx[i][vtx]) < 0.1 
    /////////////////(here D0 and DZ are wrt vertex selected by the mva vertexing)

    return myel;
}

int LoopAll::MuonSelection2012(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind){
    int mymu = -1;  
    TLorentzVector* thismu;
    float thiseta = -100;
    float thispt = -100;
    float thisiso =1000;
    int passingMu = 0;

    for( int indmu=0; indmu<mu_glo_n; indmu++){

        thismu = (TLorentzVector*) mu_glo_p4->At(indmu);
        thiseta = fabs(thismu->Eta());
        thispt = thismu->Pt();

        if(thiseta>2.4) continue;
        if(thispt<20) continue;

        if(!MuonTightID2012(indmu,vtxind)) continue;
        if(!MuonIsolation2012(indmu, thispt)) continue;
        if(!MuonPhotonCuts2012(pho1, pho2, thismu)) continue;
    
        passingMu++;
        mymu = indmu;
        std::cout << setprecision(4) << "MUON EVENT -> Run = " << run << "  Lumis = " << lumis << "  Event = " << event << "  SelVtx = " << vtxind << " muEta = " << thiseta << "  muPhi = " << thismu->Phi() <<  "  muPt = " << thismu->Pt() <<   "  muCHIso = " << mu_glo_chhadiso04[indmu] <<endl;

    }
    return mymu;
}

bool LoopAll::MuonLooseID2012(int indmu){
    if(mu_glo_type[indmu]%10000==0) return false;  // PF Muon
    if(mu_glo_type[indmu]%1000==0 && mu_glo_type[indmu]%100==0) return false;  // not global, not tracker
    return true;
}

bool LoopAll::MuonTightID2012(int indmu, int vtxind){
    
    if(mu_glo_type[indmu]<11000) return false;  // global and PF Muon
    if(mu_glo_chi2[indmu]/mu_glo_dof[indmu]>=10) return false;
    if(mu_glo_validChmbhits[indmu]<=0) return false;
    if(mu_glo_nmatches[indmu]<=1) return false;

    // need to calculate d0, dz wrt chosen vtx
    if(vtxind!=-1){
        if(fabs(mu_glo_D0Vtx[indmu][vtxind]) > 0.2) return false;
        if(fabs(mu_glo_DZVtx[indmu][vtxind]) > 0.5)  return false;
    }
    if(mu_glo_pixelhits[indmu]<=0) return false;  
    if(mu_tkLayers[indmu]<=5) return false;

    return true;
}

bool LoopAll::MuonIsolation2012(int indmu, float mupt, bool isTight){

    float thisiso=((mu_glo_nehadiso04[indmu]+mu_glo_photiso04[indmu])>mu_dbCorr[indmu]) ?
      mu_glo_chhadiso04[indmu]+mu_glo_nehadiso04[indmu]+mu_glo_photiso04[indmu]-mu_dbCorr[indmu] : mu_glo_chhadiso04[indmu];    
    if ((thisiso/mupt)>0.2) return false;
    if(isTight){
      if ((thisiso/mupt)>0.12) return false;
    }

    return true;
}

bool LoopAll::MuonPhotonCuts2012(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector* thismu){
    if(std::min(pho1.DeltaR(*thismu),pho2.DeltaR(*thismu)) < 1) return false;   

    return true;
}


int LoopAll::MuonSelection2012B(float muptcut){

    int mymu = -1;  
    TLorentzVector* thismu;
    float thiseta = -100;
    float thispt = -100;
    float thisiso =1000;
    int passingMu = 0;
    float bestpt = -2.0;

    if(GFDEBUG) std::cout<<"mu_glo_n "<<mu_glo_n<<std::endl;
    for( int indmu=0; indmu<mu_glo_n; indmu++){

        thismu = (TLorentzVector*) mu_glo_p4->At(indmu);
        thiseta = fabs(thismu->Eta());
        thispt = thismu->Pt();

        if(thiseta>2.4) continue;
        if(thispt<muptcut) continue;

        if(!MuonTightID2012(indmu)) continue;
        if(!MuonIsolation2012(indmu, thispt)) continue;
    
	if(bestpt<thispt) {
	  bestpt=thispt;
	  mymu = indmu;
	}
        if(GFDEBUG) std::cout<<"new mymu "<<mymu<<std::endl;
    }
    if(GFDEBUG) std::cout<<"final mymu "<<mymu<<std::endl;
        
    return mymu;
}



bool LoopAll::MuonPhotonCuts2012B(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector* thismu,float deltaRcut){
  //  cout<<"mu:"<<pho1.DeltaR(*thismu)<<" "<<pho2.DeltaR(*thismu)<<endl;
    if(pho1.DeltaR(*thismu)<deltaRcut) return false;
    if(pho2.DeltaR(*thismu)<deltaRcut) return false;   

    return true;
}


int LoopAll::FindMuonVertex(int mu_ind){
    int vtx_ind=-1;
    float vtx_dz=10000;
    for(int ivtx=0; ivtx<vtx_std_n; ivtx++){
      if(vtx_dz>fabs(mu_glo_DZVtx[mu_ind][ivtx])) {
	vtx_dz=fabs(mu_glo_DZVtx[mu_ind][ivtx]);
            vtx_ind=ivtx;
        }
    }
    return vtx_ind;
}


int LoopAll::ElectronSelection2012(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind, bool phodepend){

  int myel = -1;
  
  //Loose CUT-BASED ELECTRON ID
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification

  TLorentzVector* thisel;
  TLorentzVector* thissc;
  float thiseta = -100;
  float thispt = -100;
  float thisiso =1000;
  int passingEl = 0;
  
  for( int indel=0; indel<el_std_n; indel++){
    
    if(!ElectronLooseEGammaID(indel,vtxind)) continue;
    
    thisel = (TLorentzVector*) el_std_p4->At(indel);
    thissc = (TLorentzVector*) el_std_sc->At(indel);
    thispt = thisel->Pt();
    thiseta = fabs(thissc->Eta());
    if(thispt<20) continue;
    if(thiseta>2.5 || (thiseta>1.442 && thiseta<1.566)) continue;
    
    if(phodepend){
        if(!ElectronPhotonCuts(pho1, pho2, *thisel)) continue;
    }

    passingEl++;
    myel = indel;

    //std::cout << setprecision(4) << "ELECTRON EVENT -> Run = " << run << "  Lumis = " << lumis << "  Event = " << event << "  SelVtx = " << vtxind << " elEta = " << thiseta << "  elPhi = " << thisel->Phi() <<  "  elPt = " << thispt <<   "  elIso = " << thisiso/thispt << " pfiso_charged = " <<el_std_pfiso_charged[indel]<<"  pfiso_neutral = "<<el_std_pfiso_neutral[indel]<<"  pfiso_photon = "<<el_std_pfiso_photon[indel] << "  rho = "<<rho<<endl;

  }
  return myel;
}

bool LoopAll::ElectronLooseEGammaID(int electron_index, int vtxind){

    bool pass=false;

    if(electron_index<0 || electron_index>=el_std_n){
        std::cout<<"LoopAll::ElectronLooseEGammaID:  electron_index "<<electron_index<<" is out of bounds."<<std::endl;
        std::cout<<"el_std_n "<<el_std_n<<std::endl;
        return pass;
    }

    TLorentzVector* thisel = (TLorentzVector*) el_std_p4->At(electron_index);
    TLorentzVector* thissc = (TLorentzVector*) el_std_sc->At(electron_index);
    
    float thiseta = fabs(thissc->Eta());
    float thispt = thisel->Pt();

    double Aeff=0.;
    if(thiseta<1.0)                   Aeff=0.10;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.12;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.085;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.11;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.12;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.12;
    if(thiseta>=2.4)                  Aeff=0.13;

    
    //EE-EB common cuts
    float overE_overP=fabs((1/el_std_pin[electron_index])-(1/(el_std_pin[electron_index]*el_std_eopin[electron_index])));
    if(vtxind!=-1){
        if(fabs(el_std_D0Vtx[electron_index][vtxind]) > 0.02) return pass;
        if(fabs(el_std_DZVtx[electron_index][vtxind]) > 0.2)  return pass;
    }
    if(overE_overP>0.05)return pass;    
    if(el_std_hp_expin[electron_index]>1) return pass;
    if(el_std_conv[electron_index]==0) return pass;
    float thisiso=el_std_pfiso_charged[electron_index]+std::max(el_std_pfiso_neutral[electron_index]+el_std_pfiso_photon[electron_index]-rho*Aeff,0.);
    if (thisiso/thispt >0.15) return pass;  
    
    if(thiseta<1.442) {   // EB cuts
      if(fabs(el_std_detain[electron_index])>=0.007) return pass;
      if(fabs(el_std_dphiin[electron_index])>=0.15) return pass;
      if(el_std_sieie[electron_index]>=0.01) return pass;
      if(el_std_hoe[electron_index]>=0.12) return pass;
    } else {  // EE cuts
      if(fabs(el_std_detain[electron_index])>=0.009) return pass;
      if(fabs(el_std_dphiin[electron_index])>=0.10) return pass;
      if(el_std_sieie[electron_index]>=0.03) return pass; 
      if(el_std_hoe[electron_index]>=0.10) return pass;
      if(thispt<10){
        if (thisiso/thispt >0.1) return pass;  
      } 
    }

    pass=true;
    return pass;

}

bool LoopAll::ElectronTightEGammaID(int electron_index, int vtxind){

    bool pass=false;

    if(electron_index<0 || electron_index>=el_std_n){
        std::cout<<"LoopAll::ElectronTightEGammaID:  electron_index "<<electron_index<<" is out of bounds."<<std::endl;
        std::cout<<"el_std_n "<<el_std_n<<std::endl;
        return pass;
    }

    TLorentzVector* thisel = (TLorentzVector*) el_std_p4->At(electron_index);
    TLorentzVector* thissc = (TLorentzVector*) el_std_sc->At(electron_index);
    
    float thiseta = fabs(thissc->Eta());
    float thispt = thisel->Pt();

    double Aeff=0.;
    if(thiseta<1.0)                   Aeff=0.10;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.12;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.085;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.11;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.12;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.12;
    if(thiseta>=2.4)                  Aeff=0.13;

    
    //EE-EB common cuts
    float overE_overP=fabs((1/el_std_pin[electron_index])-(1/(el_std_pin[electron_index]*el_std_eopin[electron_index])));
    if(vtxind!=-1){
        if(fabs(el_std_D0Vtx[electron_index][vtxind]) > 0.02) return pass;
        if(fabs(el_std_DZVtx[electron_index][vtxind]) > 0.1)  return pass;
    }
    if(overE_overP>0.05)return pass;    
    if(el_std_hp_expin[electron_index]>0) return pass;
    if(el_std_conv[electron_index]==0) return pass;
    float thisiso=el_std_pfiso_charged[electron_index]+std::max(el_std_pfiso_neutral[electron_index]+el_std_pfiso_photon[electron_index]-rho*Aeff,0.);
    if (thisiso/thispt >0.1) return pass;  
    
    if(thiseta<1.442) {   // EB cuts
      if(fabs(el_std_detain[electron_index])>=0.004) return pass;
      if(fabs(el_std_dphiin[electron_index])>=0.03) return pass;
      if(el_std_sieie[electron_index]>=0.01) return pass;
      if(el_std_hoe[electron_index]>=0.12) return pass;
    } else {  // EE cuts
      if(fabs(el_std_detain[electron_index])>=0.005) return pass;
      if(fabs(el_std_dphiin[electron_index])>=0.02) return pass;
      if(el_std_sieie[electron_index]>=0.03) return pass; 
      if(el_std_hoe[electron_index]>=0.10) return pass;
      if(thispt<10){
        if (thisiso/thispt >0.07) return pass;  
      } 
    }

    pass=true;
    return pass;
}

bool LoopAll::ElectronPhotonCuts(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector& ele){
    bool pass=false;
    if(std::min( pho1.DeltaR(ele), pho2.DeltaR(ele))<=1) return pass;
    TLorentzVector elpho1=ele + pho1;
    TLorentzVector elpho2=ele + pho2;
    if( fabs(elpho1.M() - 91.19) <= 5) return pass;
    if( fabs(elpho2.M() - 91.19) <= 5) return pass;
    
    pass=true;
    return pass;
}


int LoopAll::ElectronSelectionMVA2012(float elptcut){
    
    int el_ind=-1;
    float bestmvaval=-2;

    for(int iel=0; iel<el_std_n; iel++){
        int vtx_ind = FindElectronVertex(iel);
	if(ElectronMVACuts(iel,vtx_ind)){
	  if(GFDEBUG) std::cout<<"passing mva "<<std::endl;
            TLorentzVector* thiselp4 = (TLorentzVector*) el_std_p4->At(iel);
            if(GFDEBUG) std::cout<<"passing eta "<<thiselp4->Eta()<<std::endl;
            if(elptcut<thiselp4->Pt()){
                if(GFDEBUG) std::cout<<"passing pt "<<std::endl;
                if(bestmvaval<el_std_mva_nontrig[iel]) {
                    bestmvaval=el_std_mva_nontrig[iel];
                    el_ind=iel;
                }
            }
        }
    }
   
    if(GFDEBUG) std::cout<<"final el_ind "<<el_ind<<std::endl;
    return el_ind;
}

bool LoopAll::ElectronMVACuts(int el_ind, int vtx_ind){
    bool pass=false;

    if(GFDEBUG) std::cout<<"Is el in bounds?  el el_std_n "<<el_ind<<" "<<el_std_n<<std::endl;
    if(el_ind<0 || el_ind>=el_std_n) return pass;

    if(GFDEBUG) std::cout<<"Passes el mva?  mva "<<el_std_mva_nontrig[el_ind]<<std::endl;
    if(el_std_mva_nontrig[el_ind]<0.9) return pass;

    if(GFDEBUG) std::cout<<"Passes el iso/pt?   "<<el_std_mva_nontrig[el_ind]<<std::endl;
    TLorentzVector* thisel = (TLorentzVector*) el_std_p4->At(el_ind);
    TLorentzVector* thissc = (TLorentzVector*) el_std_sc->At(el_ind);
    float thiseta = fabs(thissc->Eta());
    float thispt = thisel->Pt();

    if(thiseta>2.5 || (thiseta>1.442 && thiseta<1.566)) return pass;

    double Aeff=0.;
    /*
    if(thiseta<1.0)                   Aeff=0.10;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.12;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.085;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.11;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.12;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.12;
    if(thiseta>=2.4)                  Aeff=0.13;
    */
    if(thiseta<1.0)                   Aeff=0.135;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.168;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.068;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.116;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.162;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.241;
    if(thiseta>=2.4)                  Aeff=0.23;
    float thisiso=el_std_pfiso_charged[el_ind]+std::max(el_std_pfiso_neutral[el_ind]+el_std_pfiso_photon[el_ind]-rho*Aeff,0.);
    
    if(GFDEBUG) std::cout<<"Passes el iso/pt?  iso pt "<<thisiso<<" "<<thispt<<std::endl;
    if (thisiso/thispt >0.15) return pass;  

    if(vtx_ind!=-1){
        if(GFDEBUG) std::cout<<"Passes d0 and dZ cuts?  d0 dZ "<<el_std_D0Vtx[el_ind][vtx_ind]<<" "<<el_std_DZVtx[el_ind][vtx_ind]<<std::endl;
        if(fabs(el_std_D0Vtx[el_ind][vtx_ind]) > 0.02) return pass;
        if(fabs(el_std_DZVtx[el_ind][vtx_ind]) > 0.2)  return pass;
    }

    if(el_std_hp_expin[el_ind]>1) return pass;
    if(el_std_conv[el_ind]==0)    return pass;

    pass=true;
    return pass;
}



bool LoopAll::ElectronPhotonCuts2012B(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector& ele, bool includeVHlepPlusMet, float deltaRcut){
    bool pass=false;
    if(GFDEBUG)     std::cout<<"dreg1 dreg2 "<<pho1.DeltaR(ele)<<" "<<pho2.DeltaR(ele)<<std::endl;
    if( pho1.DeltaR(ele) <= deltaRcut) return pass;
    if( pho2.DeltaR(ele) <= deltaRcut) return pass;
    TLorentzVector elpho1=ele + pho1;
    TLorentzVector elpho2=ele + pho2;
    if(GFDEBUG) std::cout<<"dMeg1 dMeg2 "<<fabs(elpho1.M() - 91.19)<<" "<<fabs(elpho2.M() - 91.19)<<std::endl;
    if(!includeVHlepPlusMet){
    if( fabs(elpho1.M() - 91.19) <= 10) return pass;
    if( fabs(elpho2.M() - 91.19) <= 10) return pass;
    }
    pass=true;
    return pass;
}

int LoopAll::FindElectronVertex(int el_ind){
    int vtx_ind=-1;
    float vtx_dz=10000;
    for(int ivtx=0; ivtx<vtx_std_n; ivtx++){
      if(vtx_dz>fabs(el_std_DZVtx[el_ind][ivtx])) {
	vtx_dz=fabs(el_std_DZVtx[el_ind][ivtx]);
            vtx_ind=ivtx;
        }
    }
    return vtx_ind;
}



//--- RECO-MC JET MATCHING --------------------------------------------------------------------------------------------------------
void LoopAll::doJetMatching(TClonesArray & reco, TClonesArray & gen, 
                            Bool_t  * match_flag, Bool_t * match_vbf_flag,  Bool_t * match_bjet_flag, Bool_t * match_cjet_flag, Bool_t * match_ljet_flag,
                            Float_t * match_pt,   Float_t * match_dr,      Float_t maxDr )
{
    int ngen  = gen.GetEntries();
    int nreco = reco.GetEntries();

    std::vector<TLorentzVector*> bs;
    std::vector<TLorentzVector*> cs;
    std::vector<TLorentzVector*> ls;
    for(int ipart=0; ipart<gp_n; ++ipart) {
	    if( abs(gp_pdgid[ipart]) == 5 && gp_status[ipart] == 3 ) {
		    TLorentzVector * gp4 = (TLorentzVector*)gp_p4->At(ipart);
		    if( gp4->Pt() > 0. ) {
			    bs.push_back(gp4);
		    }
	    }else if( abs(gp_pdgid[ipart]) == 4 && gp_status[ipart] == 3 ) {
	      TLorentzVector * gp4 = (TLorentzVector*)gp_p4->At(ipart);
	      if( gp4->Pt() > 0. ) {
		cs.push_back(gp4);
	      }
	    }else if( abs(gp_pdgid[ipart]) <= 3 && gp_status[ipart] == 3 ) {
	      TLorentzVector * gp4 = (TLorentzVector*)gp_p4->At(ipart);
	      if( gp4->Pt() > 0. ) {
		ls.push_back(gp4);
	      }
	    }
    } 
    for(int ipart=0; ipart<gp_n; ++ipart) {
	    if( abs(gp_pdgid[ipart]) == 5 && (gp_status[ipart] == 2 || gp_status[ipart] == 1) ) {
		    TLorentzVector * gp4 = (TLorentzVector*)gp_p4->At(ipart);
		    if( gp4->Pt() < 15 ) { continue; }
		    bool duplicate = false;
		    for( size_t ib=0; ib<bs.size(); ++ib ) {
			    if( gp4->DeltaR( *(bs[ib]) ) < 0.3 ) { 
				    duplicate = true;
				    break; 
			    }
		    }
		    if( ! duplicate ) {
			    bs.push_back(gp4);
		    }
	    }else if( abs(gp_pdgid[ipart]) == 4 && gp_status[ipart] == 3 ) {
	      TLorentzVector * gp4 = (TLorentzVector*)gp_p4->At(ipart);
	      if( gp4->Pt() > 0. ) {
		cs.push_back(gp4);
	      }
	    }else if( abs(gp_pdgid[ipart]) <= 3 && gp_status[ipart] == 3 ) {
	      TLorentzVector * gp4 = (TLorentzVector*)gp_p4->At(ipart);
	      if( gp4->Pt() > 0. ) {
		ls.push_back(gp4);
	      }
	    }
    }
    


    for(int ir=0; ir<nreco; ++ir) {
        match_flag[ir] = false;
        match_vbf_flag[ir] = false;
        match_bjet_flag[ir] = false;
        match_cjet_flag[ir] = false;
        match_ljet_flag[ir] = false;
        match_pt[ir] = 0.;
        match_dr[ir] = 999.;
        TLorentzVector & recop4 = *(TLorentzVector*)reco.At(ir);
        for(int ig=0; ig<ngen; ++ig) {
            /// std::cerr << "ir "  << ir << " ig " << ig << std::endl; 
            TLorentzVector & genp4 = *(TLorentzVector*)gen.At(ig);
            Float_t dR = recop4.DeltaR(genp4);
            if( dR < maxDr && dR < match_dr[ir] ) {
                match_flag[ir] = true;
                match_pt[ir]   = genp4.Pt();
                match_dr[ir]   = dR;
                if( ! match_vbf_flag[ir] ) {
                    match_vbf_flag[ir] = ( gh_vbfq1_pdgid != 0 && recop4.DeltaR( *(TLorentzVector*)gh_vbfq1_p4->At(0) ) < maxDr ||
                                           gh_vbfq2_pdgid != 0 && recop4.DeltaR( *(TLorentzVector*)gh_vbfq2_p4->At(0) ) < maxDr );
                }
            }
        }
	for( size_t ib=0; ib<bs.size(); ++ib ) {
		if( recop4.DeltaR( *(bs[ib]) ) < 0.3 ) { 
			match_bjet_flag[ir] = true;
		}
	}
	for( size_t ic=0; ic<cs.size(); ++ic ) {
		if( recop4.DeltaR( *(cs[ic]) ) < 0.3 ) { 
			match_cjet_flag[ir] = true;
		}
	}
	for( size_t il=0; il<ls.size(); ++il ) {
		if( recop4.DeltaR( *(ls[il]) ) < 0.3 ) { 
			match_ljet_flag[ir] = true;
		}
	}
    }
}


//--- RECO-MC PHOTONS MATCHING --------------------------------------------------------------------------------------------------------
bool LoopAll::FindMCLeptons(int index, int& mc1, int& mc2, int& pho, int type) {

  TLorentzVector* lep = 0;
  if (abs(type) == 11) 
    lep = (TLorentzVector*) el_std_p4->At(index);
  else
    lep = (TLorentzVector*) mu_glo_p4->At(index);

  float drmin = 0.2;
  for (int i=0; i<gp_n; i++) {
    if (abs(gp_pdgid[i]) == type && gp_status[i] == 2) {
      TLorentzVector* p4 = (TLorentzVector*) gp_p4->At(i);
      float dr = lep->DeltaR(*p4);
      if (dr < drmin) {
	drmin = dr;
	mc1 = i;
      }
    }
  }
  
  if (mc1 != -1) {
    for (int i=0; i<gp_n;  i++) {
      if (gp_mother[i] == mc1 && gp_status[i] == 1) {
	if (gp_pdgid[i] == gp_pdgid[mc1]) {
	  mc2 = i;
	}
	
	if (gp_pdgid[i] == 22) {
	  pho = i;
	}
      }
    }
  }
  
  return true;
  
  // 23 2 
  // 11 3
  // -11 3  
  // 11 2
  // -11 2
  // 22 1
}

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

void LoopAll::FillMuonGsfTracks() {
  for( int indmu=0; indmu<mu_glo_n; indmu++){
    bool hasgsftrack = false;
    for (int indel=0; indel<el_std_n; indel++) if (mu_glo_tkind[indmu]==el_std_tkind[indel]) hasgsftrack=true;
    mu_glo_hasgsftrack[indmu]=int(hasgsftrack);
  }
}

void LoopAll::VHNewLeptonCategorization(bool & VHlep1event, bool & VHlep2event, int diphotonVHlep_id, int vertex, bool VHelevent_prov, bool VHmuevent_prov, int el_ind, int mu_ind, float* smeared_pho_energy, float METcut, bool moriond2013MetCorrection){
  TLorentzVector lead_p4 = get_pho_p4( dipho_leadind[diphotonVHlep_id], vertex, &smeared_pho_energy[0]);
  TLorentzVector sublead_p4 = get_pho_p4( dipho_subleadind[diphotonVHlep_id], vertex, &smeared_pho_energy[0]); 
  TLorentzVector MET;
  MET = METCorrection2012B(lead_p4, sublead_p4, moriond2013MetCorrection);
  if(MET.Pt()>=METcut) VHlep1event=true;
  else{
    if(VHelevent_prov){
      TLorentzVector* myelectron   = (TLorentzVector*) el_std_p4->At(el_ind);
      TLorentzVector elpho1=*myelectron + lead_p4;
      TLorentzVector elpho2=*myelectron + sublead_p4;
      if(fabs(elpho1.M() - 91.19)>10 && fabs(elpho2.M() - 91.19)>10) VHlep2event=true;
    }		    
    if(VHmuevent_prov) VHlep2event=true;
  }
}

void LoopAll::VHTwoMuonsEvents(bool & VHlep1event, bool & VHlep2event, int & diphotonVHlep_id, int & muVtx, float* smeared_pho_energy, float leadEtVHlepCut, float subleadEtVHlepCut, bool applyPtoverM, bool mvaselection, float diphobdt_output_Cut_VHLep, float phoidMvaCut, bool vetodipho, bool kinonly, const char * type){
  int mu_ind_1 = MuonSelection2012B(10);
  if(mu_ind_1!=-1) {
    TLorentzVector* mymu_1 = (TLorentzVector*) mu_glo_p4->At(mu_ind_1);
    int muVtx_1 = FindMuonVertex(mu_ind_1);
    std::vector<bool> veto_indices; veto_indices.clear();
    PhotonsToVeto(mymu_1, 0.5, veto_indices, false);
    int diphotonVHlep_id_1 =  -1;
    if(mvaselection) {
      diphotonVHlep_id_1 = DiphotonMITPreSelection(type,leadEtVHlepCut,subleadEtVHlepCut,phoidMvaCut,
						   applyPtoverM, &smeared_pho_energy[0], vetodipho, kinonly, diphobdt_output_Cut_VHLep, -1, false, veto_indices);
    } else {
      diphotonVHlep_id_1 = DiphotonCiCSelection( phoSUPERTIGHT, phoSUPERTIGHT, leadEtVHlepCut,subleadEtVHlepCut, 4,
						 applyPtoverM, &smeared_pho_energy[0], true, -1, veto_indices);
    }
    if(diphotonVHlep_id_1!=-1){
      TLorentzVector lead_p4_1    = get_pho_p4( dipho_leadind[diphotonVHlep_id_1],    dipho_vtxind[diphotonVHlep_id_1], &smeared_pho_energy[0]);
      TLorentzVector sublead_p4_1 = get_pho_p4( dipho_subleadind[diphotonVHlep_id_1], dipho_vtxind[diphotonVHlep_id_1], &smeared_pho_energy[0]);
      TLorentzVector* thismu;
      int mu_ind_2 = -1; float bestpt = -2.0;
      for( int indmu=0; indmu<mu_glo_n; indmu++){
	thismu = (TLorentzVector*) mu_glo_p4->At(indmu);
	if(indmu==mu_ind_1) continue;
	if(fabs(thismu->Eta())>2.4) continue;
	if((thismu->Pt())<10) continue;
	if(!MuonTightID2012(indmu)) continue;
	if(!MuonIsolation2012(indmu, (thismu->Pt()))) continue;
	if(bestpt<(thismu->Pt())) {
	  bestpt=thismu->Pt();
	  mu_ind_2 = indmu;
	}
      }
      if(mu_ind_2!=-1){
	TLorentzVector* mymu_2 = (TLorentzVector*) mu_glo_p4->At(mu_ind_2);
	if(mymu_2->DeltaR(lead_p4_1)>0.5 && mymu_2->DeltaR(sublead_p4_1)>0.5 && (*mymu_1+*mymu_2).M()<110 && (*mymu_1+*mymu_2).M()>70){
	  VHlep1event=true;
	  VHlep2event=false;
	  diphotonVHlep_id = diphotonVHlep_id_1;
	  muVtx = muVtx_1;
	}
      }
    }
  }
}

void LoopAll::VHTwoElectronsEvents(bool & VHlep1event, bool & VHlep2event, int & diphotonVHlep_id, int & elVtx, float* smeared_pho_energy, float leadEtVHlepCut, float subleadEtVHlepCut, bool applyPtoverM, bool mvaselection, float diphobdt_output_Cut_VHLep, float phoidMvaCut, bool vetodipho, bool kinonly, const char * type){
  int el_ind_1=ElectronSelectionMVA2012(10);
  if(el_ind_1!=-1) {
    TLorentzVector* myel_1 = (TLorentzVector*) el_std_p4->At(el_ind_1);
    TLorentzVector* mysc_1 = (TLorentzVector*) el_std_sc->At(el_ind_1);
    int elVtx_1 = FindElectronVertex(el_ind_1);
    std::vector<bool> veto_indices; veto_indices.clear();
    PhotonsToVeto(mysc_1, 0.5, veto_indices, true);
    int diphotonVHlep_id_1 =  -1;
    if(mvaselection) {
      diphotonVHlep_id_1 = DiphotonMITPreSelection(type,leadEtVHlepCut,subleadEtVHlepCut,phoidMvaCut,
						   applyPtoverM, &smeared_pho_energy[0], vetodipho, kinonly, diphobdt_output_Cut_VHLep, -1, false, veto_indices);
    } else {
      diphotonVHlep_id_1 = DiphotonCiCSelection( phoSUPERTIGHT, phoSUPERTIGHT, leadEtVHlepCut,subleadEtVHlepCut, 4,
						 applyPtoverM, &smeared_pho_energy[0], true, -1, veto_indices);
    }
    if(diphotonVHlep_id_1!=-1 && (ElectronMVACuts(el_ind_1, elVtx_1))==true){
      TLorentzVector lead_p4_1 = get_pho_p4( dipho_leadind[diphotonVHlep_id_1], dipho_vtxind[diphotonVHlep_id_1], &smeared_pho_energy[0]);
      TLorentzVector sublead_p4_1 = get_pho_p4( dipho_subleadind[diphotonVHlep_id_1], dipho_vtxind[diphotonVHlep_id_1], &smeared_pho_energy[0]);
      int el_ind_2=-1; float bestmvaval=-2;
      for(int iel=0; iel<el_std_n; iel++){
	TLorentzVector* thiselp4 = (TLorentzVector*) el_std_p4->At(iel);
	if(iel==el_ind_1) continue;
	if(el_std_mva_nontrig[iel]<0.9) continue;
	if(thiselp4->Eta()>2.5 || (thiselp4->Eta()>1.442 && thiselp4->Eta()<1.566)) continue;
	if(fabs(el_std_D0Vtx[iel][elVtx_1]) > 0.02) continue;
	if(fabs(el_std_DZVtx[iel][elVtx_1]) > 0.2)  continue;
	if(el_std_hp_expin[iel]>1) continue;
	if(el_std_conv[iel]==0)    continue;
	if(ElectronMVACuts(iel) && thiselp4->Pt()>10 && bestmvaval<el_std_mva_nontrig[iel]){
	  bestmvaval=el_std_mva_nontrig[iel];
	  el_ind_2=iel;
	}
      }
      if(el_ind_2!=-1){
	TLorentzVector* myel_2 = (TLorentzVector*) el_std_p4->At(el_ind_2);
	TLorentzVector* mysc_2 = (TLorentzVector*) el_std_sc->At(el_ind_2);
	if(myel_2->DeltaR(lead_p4_1)>0.5 && myel_2->DeltaR(sublead_p4_1)>0.5 && (*myel_1+*myel_2).M()<110 && (*myel_1+*myel_2).M()>70){
	  VHlep1event=true;
	  VHlep2event=false;
	  diphotonVHlep_id = diphotonVHlep_id_1;
	  elVtx = elVtx_1;
	}
      }
    }
  }
}



#ifdef NewFeatures
#include "Marco/plotInteractive_cc.h"
#endif
