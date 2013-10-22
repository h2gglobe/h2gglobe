#include "../interface/PhotonJetAnalysis.h"

#define PADEBUG 0

using namespace std;

PhotonJetAnalysis::PhotonJetAnalysis(){}
PhotonJetAnalysis::~PhotonJetAnalysis(){}

// ----------------------------------------------------------------------------------------------------
bool PhotonJetAnalysis::SelectEventsReduction(LoopAll& l, int jentry)
{

  if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
  // require at least two reconstructed photons to store the event
  //if( pho_acc.size() < 2 || l.get_pho_p4( pho_acc[0], 0, &corrected_pho_energy[0] ).Pt() < presel_scet1 ) { return false; }
  if( pho_acc.size() < 1 || l.get_pho_p4( pho_acc[0], 0, &corrected_pho_energy[0] ).Pt() < presel_scet1 ) { return false; }

  vtxAna_.clear();
  l.vtx_std_ranked_list->clear();
  l.dipho_vtx_std_sel->clear();
  l.vtx_std_ranked_list->clear();
  l.vtx_std_evt_mva->clear();
  l.vtx_std_sel=0;
  float maxSumPt = 0.;
  l.dipho_n = 0;
  bool oneKinSelected = false;

  //Jet Stuff
  TLorentzVector MaxPtJet(0,0,0,0);
  for (unsigned int ijet=0; ijet<(unsigned int) l.jet_algoPF1_n; ijet++) {
    TLorentzVector JetP4 = *((TLorentzVector*) l.jet_algoPF1_p4->At(ijet));
    if (fabs(JetP4.Eta())>2.4 || l.jet_algoPF1_ntk[ijet]<3 || JetP4.Pt()<30) continue;
    TLorentzVector JetTrackSumP4(0,0,0,0);
    for (unsigned int itk=0; itk<(unsigned int) l.jet_algoPF1_ntk[ijet]; itk++) {
      if (l.jet_algoPF1_tkind->at(ijet).at(itk)>=10000) continue;
      TLorentzVector trackp4 = *((TLorentzVector*) l.tk_p4->At(l.jet_algoPF1_tkind->at(ijet).at(itk)));
      if (trackp4.Pt()>1.0) JetTrackSumP4 += trackp4;
    }
    if (JetP4.Pt()>30.0 && JetTrackSumP4.Pt()>30.0 && JetP4.Pt()>MaxPtJet.Pt()) MaxPtJet=JetP4;
  }

  // fill ID variables
  if( forcedRho >= 0. ) {
    l.rho = forcedRho;
  } else if ( l.rho == 0. ) {
    l.rho = l.rho_algo1;
  }
  l.FillCICInputs();
  if(reComputeCiCPF) { l.FillCICPFInputs(); }
  l.FillCIC();
  l.FillMuonGsfTracks();

  if(l.itype[l.current]<0) {
    bool foundHiggs=FindHiggsObjects(l);
    if(PADEBUG)  cout << " foundHiggs? "<<foundHiggs<<std::endl;
  } else {
    SetNullHiggs(l);
  }
  /// Jet matching
  // pfJets ak5
  l.doJetMatching(*l.jet_algoPF1_p4,*l.genjet_algo1_p4,l.jet_algoPF1_genMatched,l.jet_algoPF1_vbfMatched,l.jet_algoPF1_bgenMatched,l.jet_algoPF1_cgenMatched,l.jet_algoPF1_lgenMatched,l.jet_algoPF1_genPt,l.jet_algoPF1_genDr);
  // pfJets ak7
  //l.doJetMatching(*l.jet_algoPF2_p4,*l.genjet_algo2_p4,l.jet_algoPF2_genMatched,l.jet_algoPF2_vbfMatched,l.jet_algoPF2_genPt,l.jet_algoPF2_genDr);
  // CHS ak5
  l.doJetMatching(*l.jet_algoPF3_p4,*l.genjet_algo1_p4,l.jet_algoPF3_genMatched,l.jet_algoPF3_vbfMatched,l.jet_algoPF3_bgenMatched,l.jet_algoPF3_cgenMatched,l.jet_algoPF3_lgenMatched,l.jet_algoPF3_genPt,l.jet_algoPF3_genDr);


  if( pho_presel.size() < 1 ) {
    // zero or one photons, can't determine a vertex based on photon pairs
    l.vtx_std_ranked_list->push_back( std::vector<int>() );
    for(int ii=0;ii<l.vtx_std_n; ++ii) { l.vtx_std_ranked_list->back().push_back(ii); }
    l.vtx_std_sel = 0;
  } else {
    // fully combinatorial vertex selection
    std::vector<std::pair<int,int> > diphotons;
    for(size_t ip=0; ip<pho_presel.size(); ++ip) {
      for(size_t jp=ip+1; jp<pho_presel.size(); ++jp) {
        diphotons.push_back( std::make_pair( pho_presel[ip], pho_presel[jp] ) );
      }
    }
    for(size_t id=0; id<l.pho_n; ++id ) {

      int ipho1 = id;
      int ipho2 = id;

      if(PADEBUG)        cout << " SelectEventsReduction going to fill photon info " << endl;
      PhotonInfo pho1=l.fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions,&corrected_pho_energy[0]);
      //PhotonInfo pho2=l.fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions);
      PhotonInfo pho2=PhotonInfo(-1,MaxPtJet.Vect(),MaxPtJet.E());
      pho2.isFake(true);
      if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;

      l.vertexAnalysis(vtxAna_, pho1, pho2 );
      std::vector<int> vtxs = l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames, mvaVertexSelection,
                                                tmvaPerVtxReader_, tmvaPerVtxMethod);

      TLorentzVector lead_p4 = l.get_pho_p4( ipho1, vtxs[0], &corrected_pho_energy[0] ).Pt();
      TLorentzVector sublead_p4 = MaxPtJet;

      //           if(sublead_p4.Pt()  > lead_p4.Pt() ) {
      //             std::swap( diphotons[id].first,  diphotons[id].second );
      //             std::swap( lead_p4,  sublead_p4 );
      //           }

      //           if( lead_p4.Pt() < presel_scet1 || sublead_p4.Pt() < presel_scet2 ||
      //               fabs(lead_p4.Eta()) > presel_maxeta || fabs(sublead_p4.Eta()) > presel_maxeta ) {
      //             vtxAna_.discardLastDipho();
      //          continue;
      //           }
      oneKinSelected = true;

      //           if( ! l.PhotonMITPreSelection(ipho1, vtxs[0], &corrected_pho_energy[0] )
      //               || ! l.PhotonMITPreSelection(ipho2, vtxs[0], &corrected_pho_energy[0] ) ) {
      //          vtxAna_.discardLastDipho();
      //          continue;
      //           }

      l.vtx_std_ranked_list->push_back(vtxs);
      if( tmvaPerEvtReader_ ) {
        float vtxEvtMva = vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back() );
        l.vtx_std_evt_mva->push_back(vtxEvtMva);
      }
      if( l.vtx_std_ranked_list->back().size() != 0 && ! useDefaultVertex ) {
        l.dipho_vtx_std_sel->push_back( (l.vtx_std_ranked_list)->back()[0] );
      } else {
        l.dipho_vtx_std_sel->push_back(0);
        std::cerr << "NO VERTEX SELECTED " << l.event << " " << l.run << " " << id << " " << id << std::endl;
      }

      l.dipho_leadind[id] = id;
      l.dipho_subleadind[id] = id;
      l.dipho_vtxind[id] = l.dipho_vtx_std_sel->back();

      l.dipho_sumpt[id] = lead_p4.Pt() + MaxPtJet.Pt();

      if( l.dipho_sumpt[id] > maxSumPt ) {
        l.vtx_std_sel = l.dipho_vtx_std_sel->back();
        maxSumPt = l.dipho_sumpt[id];
      }

      // make sure that vertex analysis indexes are in synch
      assert( id == vtxAna_.pairID(ipho1,-1) );

      //l.dipho_n++;
    }

    MetCorrections2012( l );
  }

  // Post-process jets and compute beta variables for missing vertexes if needed.
  int highestVtx = ( ! l.dipho_vtx_std_sel->empty() ?
                     *std::max_element(l.dipho_vtx_std_sel->begin(), l.dipho_vtx_std_sel->end()) + 1
                     : 1 );
  for(int ivtx = 0; ivtx<highestVtx; ++ivtx ) {
    postProcessJets(l,ivtx);
  }

  return oneKinSelected;
}

bool PhotonJetAnalysis::SkimEvents(LoopAll& l, int jentry)
{

  if( run7TeV4Xanalysis ) { l.version=12; }
  else { l.version=13; }

  l.b_pho_n->GetEntry(jentry);
  if( l.pho_n < 1 ) { return false; }

  if( skimOnDiphoN && l.typerun == l.kFill ) {
    l.b_dipho_n->GetEntry(jentry);
    if( l.dipho_n < 1 ) {
      return false;
    }
  }

  //Jet Stuff
  vector <TLorentzVector> GoodJetVector;
  for (unsigned int ijet=0; ijet<(unsigned int) l.jet_algoPF1_n; ijet++) {
    TLorentzVector JetP4 = *((TLorentzVector*) l.jet_algoPF1_p4->At(ijet));
    if (fabs(JetP4.Eta())>2.4 || l.jet_algoPF1_ntk[ijet]<3 || JetP4.Pt()<30) continue;
    TLorentzVector JetTrackSumP4(0,0,0,0);
    for (unsigned int itk=0; itk<(unsigned int) l.jet_algoPF1_ntk[ijet]; itk++) {
      if (l.jet_algoPF1_tkind->at(ijet).at(itk)>=10000) continue;
      TLorentzVector trackp4 = *((TLorentzVector*) l.tk_p4->At(l.jet_algoPF1_tkind->at(ijet).at(itk)));
      if (trackp4.Pt()>1.0) JetTrackSumP4 += *((TLorentzVector*) l.tk_p4->At(l.jet_algoPF1_tkind->at(ijet).at(itk)));
    }
    if (JetP4.Pt()>30. && JetTrackSumP4.Pt()>30.) GoodJetVector.push_back(JetP4);
  }

  if (GoodJetVector.size()<1) return false;

  // do not run trigger selection on MC
  int filetype = l.itype[l.current];
  bool skipTrigger = !doTriggerSelection || ( filetype != 0 && !l.runZeeValidation ) || triggerSelections.empty();
  if( ! skipTrigger ) {
    // get the trigger selection for this run
    l.b_run->GetEntry(jentry);
    std::vector<TriggerSelection>::iterator isel = find(triggerSelections.begin(), triggerSelections.end(), l.run );
    if(isel == triggerSelections.end() ) {
      std::cerr << "No trigger selection for run " << l.run << "defined" << std::endl;
      return true;
    }

    // get the trigger data
    if( l.version < 13 ) {
      l.b_hlt1_bit->GetEntry(jentry);
      l.b_hlt_path_names_HLT1->GetEntry(jentry);
      if( !  isel->pass( *(l.hlt_path_names_HLT1), *(l.hlt1_bit) ) ) {
        return false;
      }
    } else {
      l.b_hlt_bit->GetEntry(jentry);
      l.b_hlt_path_names_HLT->GetEntry(jentry);
      if( !  isel->pass( *(l.hlt_path_names_HLT), *(l.hlt_bit) ) ) {
        return false;
      }
    }
    //l.countersred[trigCounter_]++;
  }

  if( l.typerun == l.kReduce || l.typerun == l.kFillReduce ) {
    //// if( filetype == 2 ) { // photon+jet
    ////    l.b_process_id->GetEntry(jentry);
    ////    if( l.process_id == 18 ) {
    ////        return false;
    ////    }
    //// }

    if(selectprocess){
      if(processtoselect!=l.process_id){
        return false;
      }
    }

    if( filetype != 0 && ! (keepPP && keepPF && keepFF) ) {
      l.b_gp_n->GetEntry(jentry);
      l.b_gp_mother->GetEntry(jentry);
      l.b_gp_status->GetEntry(jentry);
      l.b_gp_pdgid->GetEntry(jentry);
      l.b_gp_p4->GetEntry(jentry);

      int np = 0;
      for(int ip=0;ip<l.gp_n;++ip) {
        if( l.gp_status[ip] != 1 || l.gp_pdgid[ip] != 22 ) {
          continue;
        }
        TLorentzVector * p4 = (TLorentzVector*) l.gp_p4->At(ip);
        if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
        int mother_id = abs( l.gp_pdgid[ l.gp_mother[ip] ] );
        if( mother_id <= 25 ) { ++np; }
        if( np >= 2 ) { break; }
      }
      if( np >= 2 && ! keepPP ) { return false; }
      if( np == 1 && ! keepPF ) { return false; }
      if( np == 0 && ! keepFF ) { return false; }
    }
  }

  return true;
}
