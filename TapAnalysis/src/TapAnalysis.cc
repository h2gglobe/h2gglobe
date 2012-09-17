
#include "TapAnalysis/interface/TapAnalysis.h"
#include "TapAnalysis/interface/PhotonIDCuts.h"

#include "TRegexp.h"
#include <iostream>

#define TapAnalysisDEBUG 0

using namespace std;

void TapAnalysis::FillFlatTree(LoopAll& l, Int_t type, Int_t ipho1, Int_t ipho2, Int_t iele1, Int_t iele2,
			       Float_t mass, Float_t weight, Float_t id) {

  std::cout << l.run << std::endl;
  l.FillTree("run", l.run);
  l.FillTree("lumis", l.lumis);
  l.FillTree("event", l.event);
  l.FillTree("type", type);
  l.FillTree("charge", l.el_std_charge[iele1]);
  l.FillTree("chargeTag", l.el_std_charge[iele2]);
  TLorentzVector* p4_tag = (TLorentzVector*) l.el_std_p4->At(iele1);
  TLorentzVector* p4_probe = (TLorentzVector*) l.el_std_p4->At(iele2);
  l.FillTree("eta", p4_probe->Eta());
  l.FillTree("etatag", p4_tag->Eta());
  l.FillTree("mass", mass);
  l.FillTree("weight", weight);
  l.FillTree("fbrem", l.el_std_fbrem[iele2]);  
  l.FillTree("r9tag", l.pho_r9[ipho1]);
  l.FillTree("r9", l.pho_r9[ipho2]);

  TLorentzVector* p4_pho_tag = (TLorentzVector*) l.pho_p4->At(ipho1);
  TLorentzVector* p4_pho_probe = (TLorentzVector*) l.pho_p4->At(ipho2);
  l.FillTree("ettag", p4_pho_tag->Eta()*sin(p4_tag->Theta()));
  l.FillTree("et",  p4_pho_probe->Eta()*sin(p4_probe->Theta()));
  l.FillTree("nvtx", l.vtx_std_n);
  l.FillTree("id", id);
}

// ----------------------------------------------------------------------------------------------------
TapAnalysis::TapAnalysis()  : 
  name_("TapAnalysis")
{}

// ----------------------------------------------------------------------------------------------------
TapAnalysis::~TapAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void TapAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void TapAnalysis::Init(LoopAll& l) {
  if(TapAnalysisDEBUG) 
    cout << "InitRealTapAnalysis START"<<endl;
  
  l.tmvaReaderID_Single_Barrel = new TMVA::Reader("!Color:Silent");
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.r9",   &l.tmva_photonid_r9 );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.sigietaieta",   &l.tmva_photonid_sieie );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.scetawidth",   &l.tmva_photonid_etawidth );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.scphiwidth",   &l.tmva_photonid_phiwidth );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_CoviEtaiPhi",   &l.tmva_photonid_sieip );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_s4ratio",   &l.tmva_photonid_s4ratio );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_GammaIso",   &l.tmva_photonid_pfphotoniso03 );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_ChargedIso_selvtx",   &l.tmva_photonid_pfchargedisogood03 );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.idmva_ChargedIso_worstvtx",   &l.tmva_photonid_pfchargedisobad03 );
  l.tmvaReaderID_Single_Barrel->AddVariable("ph.sceta",   &l.tmva_photonid_sceta );
  l.tmvaReaderID_Single_Barrel->AddVariable("rho",   &l.tmva_photonid_eventrho );
  l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost", "aux/2012ICHEP_PhotonID_Barrel_BDT.weights.xml");

  l.tmvaReaderID_Single_Endcap = new TMVA::Reader("!Color:Silent");
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.r9",   &l.tmva_photonid_r9 );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.sigietaieta",   &l.tmva_photonid_sieie );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.scetawidth",   &l.tmva_photonid_etawidth );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.scphiwidth",   &l.tmva_photonid_phiwidth );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_CoviEtaiPhi",   &l.tmva_photonid_sieip );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_s4ratio",   &l.tmva_photonid_s4ratio );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_GammaIso",   &l.tmva_photonid_pfphotoniso03 );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_ChargedIso_selvtx",   &l.tmva_photonid_pfchargedisogood03 );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_ChargedIso_worstvtx",   &l.tmva_photonid_pfchargedisobad03 );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.sceta",   &l.tmva_photonid_sceta );
  l.tmvaReaderID_Single_Endcap->AddVariable("rho",   &l.tmva_photonid_eventrho );
  l.tmvaReaderID_Single_Endcap->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &l.tmva_photonid_ESEffSigmaRR );
  l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost", "aux/2012ICHEP_PhotonID_Endcap_BDT.weights_PSCorr.xml");
  
  if(TapAnalysisDEBUG) 
    cout <<"InitRealTapAnalysis END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
bool TapAnalysis::Analysis(LoopAll& l, Int_t jentry) {
  if(TapAnalysisDEBUG) 
    cout <<"Analysis START"<<endl;
   if (puHist != "" && puHist != "auto" ) {
        if(TapAnalysisDEBUG) 
            cout << "Opening PU file"<<endl;
        TFile* puFile = TFile::Open( puHist );
        if (puFile) {
        TH1 * target = 0;
        
        if( puTarget != "" ) {
        TFile * puTargetFile = TFile::Open( puTarget ); 
        assert( puTargetFile != 0 );
        target = (TH1*)puTargetFile->Get("pileup");
        if( target == 0 ) { target = (TH1*)puTargetFile->Get("target_pileup"); }
        target->Scale( 1. / target->Integral() );
        }
        
            if( puMap != "" ) {
                loadPuMap(puMap, puFile, target); 
            } else {
                loadPuWeights(0, puFile, target);
            }
            puFile->Close();
        }
        else {
            cout<<"Error opening " <<puHist<<" pileup reweighting histogram, using 1.0"<<endl; 
            weights[0].resize(50);
            for (unsigned int i=0; i<weights[0].size(); i++) weights[0][i] = 1.0;
        }
        if(TapAnalysisDEBUG) 
            cout << "Opening PU file END"<<endl;
    } else if ( puHist == "auto" ) {
    TFile * puTargetFile = TFile::Open( puTarget ); 
    assert( puTargetFile != 0 );
    puTargetHist = (TH1*)puTargetFile->Get("pileup");
    if( puTargetHist == 0 ) { puTargetHist = (TH1*)puTargetFile->Get("target_pileup"); }
    puTargetHist = (TH1*)puTargetHist->Clone();
    puTargetHist->SetDirectory(0);
    puTargetHist->Scale( 1. / puTargetHist->Integral() );
    puTargetFile->Close();
    }

    //apply pileup reweighting
    unsigned int n_pu = l.pu_n;
    float weight =1.;
    if (l.itype[l.current] !=0 && puHist != "") {
        std::vector<double> & puweights = weights.find( l.itype[l.current] ) != weights.end() ? weights[ l.itype[l.current] ] : weights[0]; 
        if(n_pu<puweights.size()){
            weight *= puweights[n_pu];
        }    
        else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
            cout<<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<"), event will not be reweighted for pileup"<<endl;
        }
    }

  // move ate reduction level
  bool hlTrigger = true; //checkElectronHLT(iElectron, mp, indexfiles, itype, jentry);
  if(!hlTrigger)
    return 0;

  // MET CUT
  if (l.met_pfmet > cutPFMET) 
    return 0;

  ////////////////
  // TAG ELECTRONS
  ////////////////

  // Vector of protoTag (tags without eID selection)
  std::vector<int> selectedTags; selectedTags.clear();
  
  for(unsigned int iElectron = 0; iElectron < l.el_std_n; iElectron++) {
    int iIndexSuperCluster = l.el_std_scind[iElectron];
 
    if (iIndexSuperCluster == -1)
      continue;
 
    TLorentzVector p4_sc = *((TLorentzVector *) l.el_std_sc->At(iElectron));
 
    if(fabs(p4_sc.Eta()) > 1.4442 && fabs(p4_sc.Eta()) < 1.566)
      continue;
    if (fabs(p4_sc.Et()) < cutETTag)
      continue;
 
    Int_t phoindex = -1;
    for (int k=0; k < l.pho_n; k++) {
      if (l.pho_scind[k] == l.el_std_scind[iElectron])
	phoindex = k;
    }
 
    if (phoindex != -1) {
       if ((applyPreselection > 0) && PhotonId(l, iElectron, phoindex, std::string("Presel"), 1)) {
	 if (PhotonId(l, iElectron, phoindex, selectionTypeTag, cutSelectionTag)) {
	  selectedTags.push_back(iElectron);
	}
      }
    }
  }

  //std::cout << "tags" << selectedTags.size() << std::endl;
    
  ////////////////////////////
  // ELECTRON PROBE CANDIDATES
  ////////////////////////////
  // Vector of probe electrons without eID requirements
  std::vector<int> selectedProbes; selectedProbes.clear();
  
  // Define Probe Candidates
  for (unsigned int iEl = 0; iEl<l.el_std_n; iEl++){
    int iIndexSuperCluster = l.el_std_scind[iEl];
    
    // CHECK ELECTRON TYPE
    if(iIndexSuperCluster == -1)
      continue;

    TLorentzVector p4_sc = *((TLorentzVector *) l.el_std_sc->At(iEl));
    if( fabs(p4_sc.Eta()) > 1.4442 && fabs(p4_sc.Eta()) < 1.566)
      continue;
    
    if(p4_sc.Et() < cutETProbe)
      continue;
	  
    Int_t phoindex = -1;
    for (int k=0; k < l.pho_n; k++) {
      if (l.pho_scind[k] == l.el_std_scind[iEl])
	phoindex = k;
    }
	  
    if (phoindex != -1) {
      if ((applyPreselection > 0) && PhotonId(l, iEl, phoindex, std::string("Presel"), 1)) {
	if (PhotonId(l, iEl, phoindex, selectionTypeProbe, cutSelectionProbe))
	  selectedProbes.push_back(iEl);
      }
    }
  }

  //std::cout << "probes" << selectedProbes.size() << std::endl;

  std::vector<std::pair<int, int> > phoEffPairs; phoEffPairs.clear();
  phoEffPairs = TPPairs(l, selectedTags, selectedProbes, 1, 0);

  for (unsigned int tpPair=0; tpPair<phoEffPairs.size(); tpPair++) {
    
    TLorentzVector p4_tag = *((TLorentzVector*) l.el_std_p4->At(phoEffPairs[tpPair].first));
    TLorentzVector* p4_probe = (TLorentzVector*) l.el_std_p4->At(phoEffPairs[tpPair].second);

    // APPLY Energy corrections 
    /*
      if (itype == 0) {
      if (fabs(p4_tag.Eta()) < 1.479) {
	p4_tag.SetE(p4_tag.E() *  (double)GetCutValue("value_EnergyCorrection_EB",0));
	p4_tag.SetPx(p4_tag.Px() *  (double)GetCutValue("value_EnergyCorrection_EB",0));
	p4_tag.SetPy(p4_tag.Py() *  (double)GetCutValue("value_EnergyCorrection_EB",0));
	p4_tag.SetPz(p4_tag.Pz() *  (double)GetCutValue("value_EnergyCorrection_EB",0));
      } else {
	p4_tag.SetE(p4_tag.E() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
	p4_tag.SetPx(p4_tag.Px() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
	p4_tag.SetPy(p4_tag.Py() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
	p4_tag.SetPz(p4_tag.Pz() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
      }
      
      if (fabs(p4_probe->Eta()) < 1.479) {
	p4_probe->SetE(p4_probe->E() *  (double)GetCutValue("value_EnergyCorrection_EB",0));	
	p4_probe->SetPx(p4_probe->Px() *  (double)GetCutValue("value_EnergyCorrection_EB",0));
	p4_probe->SetPy(p4_probe->Py() *  (double)GetCutValue("value_EnergyCorrection_EB",0));
	p4_probe->SetPz(p4_probe->Pz() *  (double)GetCutValue("value_EnergyCorrection_EB",0));
      } else {
	p4_probe->SetE(p4_probe->E() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
	p4_probe->SetPx(p4_probe->Px() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
	p4_probe->SetPy(p4_probe->Py() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
	p4_probe->SetPz(p4_probe->Pz() *  (double)GetCutValue("value_EnergyCorrection_EE",0));
      }
    }
    */

    double zMass = (p4_tag + *p4_probe).Mag();
    //cout<<l.run << " zmass = "<<zMass<<endl;
	
    if ((zMass > minZMass) && (zMass < maxZMass)) {
      
      int iElTag = phoEffPairs[tpPair].first;
      int iEl = phoEffPairs[tpPair].second;
      
      //std::cout << "Pair" << std::endl;
      bool pass = false;
      Int_t phoindex_tag = -1;
      Int_t phoindex = -1;
      for (int k=0; k < l.pho_n; k++) {
	if (l.pho_scind[k] == l.el_std_scind[iEl])
	  phoindex = k;
	if (l.pho_scind[k] == l.el_std_scind[iElTag])
	  phoindex_tag = k;
      }
      
      if (phoindex != -1) {
	int type = l.itype[l.current];
	Float_t id = PhotonId(l, iEl, phoindex, selectionTypeToMeasure, cutSelectionToMeasure);
	
	FillFlatTree(l, type, phoindex_tag, phoindex, iElTag, iEl, zMass, weight, id);
      }
    }
  }    
  if(TapAnalysisDEBUG) 
    cout<<"Analysis END"<<endl;
}

std::vector<std::pair<int, int> > TapAnalysis::TPPairs(LoopAll& l, std::vector<int> tags, std::vector<int> probes, int type, int chargePairing) {
  
  std::vector<std::pair<int, int> > tpPairs;

  if (type == 0)
    chargePairing = 0;
  if (chargePairing > 0)
    chargePairing = 1;
  if (chargePairing < 0)
    chargePairing = -1;
  
  // PRODUCE TAG AND PROBE PAIRS
  /*
  for(unsigned int i=0; i<tags.size(); i++) {
    for(unsigned int j=0; j<probes.size(); j++) {
      
      
      // CHECK TAG != PROBE
      if (type == 1) {
	if(probes[j] == tags[i])
	  continue;
      } else { // FOR SC PROBES
	if (probes[j] == el_std_scind[tags[i]])
	  continue;
      }
      
      if (chargePairing == 0)
	tpPairs.push_back( std::pair<int, int>(tags[i], probes[j]));
      else if (chargePairing == (el_std_charge[tags[i]]*el_std_charge[probes[j]]))
	tpPairs.push_back( std::pair<int, int>(tags[i], probes[j]));
    }
  }
  */
  if (tags.size() > 0) {
    for(unsigned int i=0; i<tags.size()-1; i++) {
      for(unsigned int j=i+1; j<tags.size(); j++) {
	//std::cout << "COPPIE" <<  i << " " << j << std::endl;
	
      // CHECK TAG != PROBE
	if (type == 1) {
	  if(tags[j] == tags[i])
	    continue;
	} else { // FOR SC PROBES
	  if (tags[j] == l.el_std_scind[tags[i]])
	    continue;
	}
	
	if (chargePairing == 0)
	  tpPairs.push_back( std::pair<int, int>(tags[i], tags[j]));
	else if (chargePairing == (l.el_std_charge[tags[i]]*l.el_std_charge[tags[j]]))
	  tpPairs.push_back( std::pair<int, int>(tags[i], tags[j]));
      }
    }
  }

  /*
  if (chargePairing == 0) {
    for (int i=0; i<tpPairs.size(); i++) {
      std::cout << tpPairs[i].first << " " << tpPairs[i].second << std::endl; 
    }
    std::cout << "__________" << std::endl;
  }
  */
  return tpPairs;
}

// ----------------------------------------------------------------------------------------------------
void TapAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{}

// ----------------------------------------------------------------------------------------------------
void TapAnalysis::FillReductionVariables(LoopAll& l, int jentry) {

  // FIX PILEUP
  for (int ipho=0;ipho<l.pho_n;ipho++){
    l.pho_s4ratio[ipho]  = l.pho_e2x2[ipho]/l.pho_e5x5[ipho];
    float rr2=l.pho_eseffsixix[ipho]*l.pho_eseffsixix[ipho]+l.pho_eseffsiyiy[ipho]*l.pho_eseffsiyiy[ipho];
    l.pho_ESEffSigmaRR[ipho] = 0.0; 
    if(rr2>0. && rr2<999999.) {
      l.pho_ESEffSigmaRR[ipho] = sqrt(rr2);
    }
  }

  // Data driven MC corrections to cluster shape variables and energy resolution estimate
  if (l.itype[l.current] !=0)
    rescaleClusterVariables(l);

  for(int ipho=0; ipho<l.pho_n; ++ipho) {
    std::vector<float> temp;

    //float neu01 = l.pfEcalIso(ipho, 0.1, 0., 0., 0., 0., 0., 0., 5);
    //float neu02 = l.pfEcalIso(ipho, 0.2, 0., 0., 0., 0., 0., 0., 5);
    float neu03 = l.pfEcalIso(ipho, 0.3, 0., 0., 0., 0., 0., 0., 5);
    float neu04 = l.pfEcalIso(ipho, 0.4, 0., 0., 0., 0., 0., 0., 5); 
    //float neu05 = l.pfEcalIso(ipho, 0.5, 0., 0., 0., 0., 0., 0., 5); 
    //float neu06 = l.pfEcalIso(ipho, 0.6, 0., 0., 0., 0., 0., 0., 5); 
    
    //l.pho_pfiso_myneutral01[ipho] = neu01;
    //l.pho_pfiso_myneutral02[ipho] = neu02;
    l.pho_pfiso_myneutral03[ipho] = neu03;
    l.pho_pfiso_myneutral04[ipho] = neu04;
    //l.pho_pfiso_myneutral05[ipho] = neu05;
    //l.pho_pfiso_myneutral06[ipho] = neu06;
    
    //float pho01 = l.pfEcalIso(ipho, 0.1, 0., 0.070, 0.015, 0., 0., 0.);
    //float pho02 = l.pfEcalIso(ipho, 0.2, 0., 0.070, 0.015, 0., 0., 0.);
    float pho03 = l.pfEcalIso(ipho, 0.3, 0., 0.070, 0.015, 0., 0., 0.);
    float pho04 = l.pfEcalIso(ipho, 0.4, 0., 0.070, 0.015, 0., 0., 0.); 
    //float pho05 = l.pfEcalIso(ipho, 0.5, 0., 0.070, 0.015, 0., 0., 0.); 
    //float pho06 = l.pfEcalIso(ipho, 0.6, 0., 0.070, 0.015, 0., 0., 0.); 
    
    //l.pho_pfiso_myphoton01[ipho] = pho01;
    //l.pho_pfiso_myphoton02[ipho] = pho02;
    l.pho_pfiso_myphoton03[ipho] = pho03;
    l.pho_pfiso_myphoton04[ipho] = pho04;
    //l.pho_pfiso_myphoton05[ipho] = pho05;
    //l.pho_pfiso_myphoton06[ipho] = pho06;
	
    float badiso = 0;
    for(int ivtx=0; ivtx<l.vtx_std_n; ++ivtx) {
      
      //float ch01 = l.pfTkIsoWithVertex(ipho,ivtx,0.1,0.02,0.02,0.0,0.2,0.1);
      float ch02 = l.pfTkIsoWithVertex(ipho,ivtx,0.2,0.02,0.02,0.0,0.2,0.1);
      float ch03 = l.pfTkIsoWithVertex(ipho,ivtx,0.3,0.02,0.02,0.0,0.2,0.1);
      float ch04 = l.pfTkIsoWithVertex(ipho,ivtx,0.4,0.02,0.02,0.0,0.2,0.1);
      //float ch05 = l.pfTkIsoWithVertex(ipho,ivtx,0.5,0.02,0.02,0.0,0.2,0.1);
      //float ch06 = l.pfTkIsoWithVertex(ipho,ivtx,0.6,0.02,0.02,0.0,0.2,0.1);
      
      //l.pho_pfiso_mycharged01->at(ipho).at(ivtx) = ch01;
      l.pho_pfiso_mycharged02->at(ipho).at(ivtx) = ch02;
      l.pho_pfiso_mycharged03->at(ipho).at(ivtx) = ch03;
      l.pho_pfiso_mycharged04->at(ipho).at(ivtx) = ch04;
      //l.pho_pfiso_mycharged05->at(ipho).at(ivtx) = ch05;
      //l.pho_pfiso_mycharged06->at(ipho).at(ivtx) = ch06;
      
      TLorentzVector* p4 = (TLorentzVector*)l.pho_p4->At(ipho);
      temp.push_back(l.photonIDMVANew(ipho, ivtx, *p4, ""));
      
      if( ch04 > badiso ) {
	badiso = ch04;
	//badvtx = ivtx;
      }
    }
   
    //pho_tkiso_badvtx_id[ipho] = badvtx;
    l.pho_pfiso_charged_badvtx_04[ipho] = badiso;
    l.pho_mitmva->push_back(temp);
  }
}

// ----------------------------------------------------------------------------------------------------
bool TapAnalysis::SelectEventsReduction(LoopAll& l, int jentry) {
  
  std::vector<int> eleIndex; 
  
  int goodEl = 0;
  for(int i =0; i<l.el_std_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.el_std_sc->At(i);
    if (p4->Et() > 20.)
      eleIndex.push_back(i);
  }
  
  Float_t mass = 0.;
  if (eleIndex.size() > 1) {
    for (unsigned int i=0; i<eleIndex.size()-1; i++) {
      TLorentzVector* p1 = (TLorentzVector*)l.el_std_sc->At(eleIndex[i]);
      
      for (unsigned int j=i+1; j<eleIndex.size(); j++) {
	TLorentzVector* p2 = (TLorentzVector*)l.el_std_sc->At(eleIndex[j]);
	Float_t mass_temp = ((*p1)+(*p2)).M();
	if (mass < mass_temp) {
	  mass = mass_temp;
	}
      }
    }
    
    if (mass > 60.)
      return 1;
  }
    
  return 0;
  
  if(TapAnalysisDEBUG)  cout << " ****************** SelectEventsReduction " << endl;
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool TapAnalysis::SkimEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool TapAnalysis::SelectEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
void TapAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) {
  l.pho_mitmva = new std::vector<std::vector<float> >();

  if (outputTree) {
    l.Branch_pho_mitmva(outputTree);
    l.Branch_pho_pfiso_charged_badvtx_04(outputTree);
  }
}

// ----------------------------------------------------------------------------------------------------
void TapAnalysis::ResetAnalysis()
{}

bool TapAnalysis::PhotonId(LoopAll& l, Int_t eleIndex, Int_t phoIndex, std::string type, Float_t selection) {

  if (selection < -10.)
    return true;
 
  bool result = false;

  int thevertexind = ChooseVertex(l, eleIndex, true);

  if (type == "CiC") {
    result = CiCPhotonIDPF(l, 4, phoIndex, thevertexind, (int)selection);
  } else if (type == "MVA") {
    result = ((*l.pho_mitmva)[phoIndex][thevertexind] > selection);
  } else if (type == "Presel") {
    result = l.PhotonMITPreSelection(phoIndex, thevertexind, 0);
  }

  return result;
}

bool TapAnalysis::checkEventHLT(LoopAll& l, std::vector<std::string> paths) {

  bool result = false;
  std::vector<unsigned short> hltNumbers;

  for (unsigned int i=0; i<paths.size(); i++) { 
    TRegexp e(TString(paths[i].c_str()));
    for (unsigned int j=0; j<l.hlt_path_names_HLT1->size(); j++) {
      TString str1((*l.hlt_path_names_HLT1)[j].c_str());
      //std::cout << (*hlt_path_names_HLT1)[j] << std::endl;
      if (str1.Contains(e)) {	
	//std::cout << (*hlt_path_names_HLT1)[j] << std::endl;
	hltNumbers.push_back(j);
      }
    }
  }

  for (int j=0; j< hltNumbers.size(); j++) {
    for (int i=0; i<(*l.hlt1_bit).size(); i++) {
      if (hltNumbers[j] == (*l.hlt1_bit)[i]) {
	result = true;
	break;
      }
    }
  }

  return result;
}

int TapAnalysis::ChooseVertex(LoopAll& l, int iEl, bool isGood) {

  TVector3* el_pos = (TVector3*)l.el_std_posvtx->At(iEl); 

  float dist[100];
  int index[100];

  for (int i=0; i<l.vtx_std_n; i++) {
    TVector3* pos = (TVector3*)l.vtx_std_xyz->At(i); 
    dist[i] = ((*pos) - (*el_pos)).Mag();
    index[i] = i;
  }
  
  for (int i=0; i<l.vtx_std_n-1; i++) {
    for (int j=i+1; j<l.vtx_std_n; j++) {
      if (dist[i] >= dist[j]) {
	float ftemp = dist[j]; 
	dist[j] = dist[i];
	dist[i] = ftemp;

	int itemp = index[j]; 
	index[j] = index[i];
	index[i] = itemp;
      }
    }
  }
  
  return index[0];
  /*
  // TOGLI IL VERTICE dell'elettrone e prendi il best sum pT

  int secondGoodIndex = index[0];
  float maxPt = 0;
  for(int i=1; i<l.vtx_std_n; i++) {
    float temp = SumPt2(jentry, index[i]);
    if (temp > maxPt) {
      maxPt = temp;
      secondGoodIndex = index[i];
    }
  }
    
  if (isGood) 
    return index[0];
  else 
    return secondGoodIndex;

  return -1;
  */
}	

TLorentzVector TapAnalysis::get_pho_p4(LoopAll& l, Int_t ipho, int ivtx) {

  TVector3 * vtx = (TVector3*) l.vtx_std_xyz->At(ivtx);
  TVector3 direction = *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])) - *vtx;
  TVector3 p = direction.Unit() * ((TLorentzVector*)l.pho_p4->At(ipho))->Energy();
  TLorentzVector p4(p.x(),p.y(),p.z(),((TLorentzVector*)l.pho_p4->At(ipho))->Energy());
  return p4;
}

bool TapAnalysis::CiCPhotonIDPF(LoopAll& l, int nCategories, int photonindex, int chosenVtx, int IDlevel) {
  
  bool passes = true;

  float vars[10];
  //float rhofacbad=0.52, rhofac=0.17;
  // need pho_calopos and pho_r9 branches here
  int thisCat = PhotonIDCategory(l, photonindex, nCategories);  
  TLorentzVector phop4 = get_pho_p4(l, photonindex, chosenVtx);

  float rhofac = 0.09;
  float rhofacbad = 0.23;
  float isosumconst    = 2.8;
  float isosumconstbad = 4.8;

  float val_isosumoet    = ((*l.pho_pfiso_mycharged03)[photonindex][chosenVtx] + l.pho_pfiso_myphoton03[photonindex] + isosumconst - l.rho_algo1*rhofac)*50./phop4.Et();
  float val_isosumoetbad = (l.pho_pfiso_myphoton04[photonindex] + l.pho_pfiso_charged_badvtx_04[photonindex] + isosumconstbad - l.rho_algo1*rhofacbad)*50./phop4.Et();
  float val_trkisooet    = ((*l.pho_pfiso_mycharged03)[photonindex][chosenVtx])*50./phop4.Et();

  vars[0] = val_isosumoet;
  vars[1] = val_isosumoetbad;
  vars[2] = val_trkisooet;

  vars[3] = l.pho_sieie[photonindex];
  vars[4] = l.pho_hoe[photonindex];
  vars[5] = l.pho_r9[photonindex];
  vars[6] = 99.; //DeltaRToTrackHgg(jentry, photonindex, chosenVtx, 2.5, 1., 0.1,0);
  vars[7] = 0;//(float)pho_haspixseed[photonindex];
  
  for(int var=0; var<NVARS; var++){
    if (nCategories == 6) {
      if(var == 5 || var == 6){   //cuts from below
	if(vars[var] < phoIDcuts6catPF[IDlevel][var][thisCat]) {
	  passes =false;
	  break;
	}
      } else {                    //cuts from above
	if (vars[var] > phoIDcuts4catPF[IDlevel][var][thisCat]) {
	  passes =false;
	  break;
	}
      }
    }
  } 
  
  //std::cout << passes << std::endl;
  return passes;
}

int TapAnalysis::PhotonIDCategory(LoopAll& l, int photonindex, int nCategories) {

  int thisCat = -1;
  
  float eta = fabs(((TVector3*)l.sc_xyz->At(l.pho_scind[photonindex]))->Eta());
  float r9  = l.pho_r9[photonindex];

  if (nCategories == 4) {
    bool etaCat = eta>1.479;
    bool r9Cat  = r9<0.94;
    
    thisCat = 2*(etaCat) + r9Cat;
  } else if (nCategories == 6) {
    bool r9Cat = (Int_t)(r9<0.94) + (r9<0.9); // 0, 1, or 2 (high r9 --> low r9)
    bool etaCat = eta>1.479;

    thisCat = 3*(etaCat) + r9Cat;
  } else {
    std::cout << "Wrong number of categories." << std::endl;
  }
    
  return thisCat;
}
