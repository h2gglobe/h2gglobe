#include "../interface/PhotonJetAnalysis.h"

#define PADEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
bool PhotonJetAnalysis::SelectEventsReduction(LoopAll& l, int jentry) 
{

    if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
    // require at least two reconstructed photons to store the event
    //if( pho_acc.size() < 2 || l.get_pho_p4( pho_acc[0], 0, &corrected_pho_energy[0] ).Pt() < presel_scet1 ) { return false; }
    if (pho_acc.size() < 1) {return false;}
    
    vtxAna_.clear();
    l.vtx_std_ranked_list->clear();
    l.dipho_vtx_std_sel->clear();
    l.vtx_std_ranked_list->clear();
    l.vtx_std_evt_mva->clear();
    l.vtx_std_sel=0;
    float maxSumPt = 0.;
    l.dipho_n = 0;

    //Jet Stuff
    vector <TLorentzVector> GoodJetVector;
    for (unsigned int ijet=0; ijet<(unsigned int) l.jet_algoPF1_n; ijet++) {
      TLorentzVector JetP4 = *((TLorentzVector*) l.jet_algoPF1_p4->At(ijet));
      TLorentzVector TrackSumP4(0,0,0,0);
      if (fabs(JetP4.Eta())>2.4 || l.jet_algoPF1_ntk[ijet]<3) continue;
      for (unsigned int itrack=0; itrack<(unsigned int) l.jet_algoPF1_ntk[ijet]; itrack++) {
        unsigned short trackindex = l.jet_algoPF1_tkind->at(ijet).at(itrack);
        if (trackindex>=3000) continue;
        TLorentzVector TrackP4 = *((TLorentzVector*) l.tk_p4->At(trackindex));
        if (TrackP4.Pt()>1.0) TrackSumP4+=TrackP4;
      }
      if (JetP4.Pt()>30 && TrackSumP4.Pt()>30) {
        GoodJetVector.push_back(JetP4);
      }
    }
    
    // fill ID variables
    if( forcedRho >= 0. ) {
      l.rho = forcedRho;
    } else if ( l.rho == 0. ) {
      l.rho = l.rho_algo1;
    }
    l.FillCICInputs();
    l.FillCIC();

    if(l.itype[l.current]<0) {
        bool foundHiggs=FindHiggsObjects(l);
        if(PADEBUG)  cout << " foundHiggs? "<<foundHiggs<<std::endl;
    } else {
        SetNullHiggs(l);
    }
    /// Jet matching
    // pfJets ak5
    l.doJetMatching(*l.jet_algoPF1_p4,*l.genjet_algo1_p4,l.jet_algoPF1_genMatched,l.jet_algoPF1_vbfMatched,l.jet_algoPF1_genPt,l.jet_algoPF1_genDr);
    // pfJets ak7
    l.doJetMatching(*l.jet_algoPF2_p4,*l.genjet_algo2_p4,l.jet_algoPF2_genMatched,l.jet_algoPF2_vbfMatched,l.jet_algoPF2_genPt,l.jet_algoPF2_genDr);
    // CHS ak5
    l.doJetMatching(*l.jet_algoPF3_p4,*l.genjet_algo1_p4,l.jet_algoPF3_genMatched,l.jet_algoPF3_vbfMatched,l.jet_algoPF3_genPt,l.jet_algoPF3_genDr);
    
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
            PhotonInfo pho2=PhotonInfo(-1,GoodJetVector[0].Vect(),GoodJetVector[0].E());
            pho2.isFake(true);
            if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;
            
            l.vertexAnalysis(vtxAna_, pho1, pho2 );
            // make sure that vertex analysis indexes are in synch 
            assert( (int)id == vtxAna_.pairID(ipho1,ipho2) );
            
            l.vtx_std_ranked_list->push_back( l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames, mvaVertexSelection, 
                                                                tmvaPerVtxReader_, tmvaPerVtxMethod) );
            if( tmvaPerEvtReader_ ) {
                float vtxEvtMva = vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back() );
                l.vtx_std_evt_mva->push_back(vtxEvtMva);
            }
            if( l.vtx_std_ranked_list->back().size() != 0 && ! useDefaultVertex ) {  
                l.dipho_vtx_std_sel->push_back( (l.vtx_std_ranked_list)->back()[0] );
            } else {
                l.dipho_vtx_std_sel->push_back(0);
                std::cerr << "NO VERTEX SELECTED " << l.event << " " << l.run << " " << diphotons[id].first << " " << diphotons[id].second << std::endl;
            }
            l.dipho_n = id+1;
            l.dipho_leadind[id] = id;
            l.dipho_subleadind[id] = id;
            l.dipho_vtxind[id] = l.dipho_vtx_std_sel->back();
            
            TLorentzVector lead_p4 = l.get_pho_p4( id, l.dipho_vtxind[id], &corrected_pho_energy[0] );
            TLorentzVector sublead_p4 = GoodJetVector[0];
            l.dipho_sumpt[id] = lead_p4.Pt() + sublead_p4.Pt();
            
            if( l.dipho_sumpt[id] > maxSumPt ) {
                l.vtx_std_sel = l.dipho_vtx_std_sel->back();
                maxSumPt = l.dipho_sumpt[id];
            }
        }
    }
    

    return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonJetAnalysis::SkimEvents(LoopAll& l, int jentry)
{
    l.b_pho_n->GetEntry(jentry);
    if( l.pho_n < 1 ) {
        return false;
    }

    //Jet Stuff
    vector <TLorentzVector> GoodJetVector;
    for (unsigned int ijet=0; ijet<(unsigned int) l.jet_algoPF1_n; ijet++) {
      TLorentzVector JetP4 = *((TLorentzVector*) l.jet_algoPF1_p4->At(ijet));
      TLorentzVector TrackSumP4(0,0,0,0);
      if (fabs(JetP4.Eta())>2.4 || l.jet_algoPF1_ntk[ijet]<3) continue;
      for (unsigned int itrack=0; itrack<(unsigned int) l.jet_algoPF1_ntk[ijet]; itrack++) {
        unsigned short trackindex = l.jet_algoPF1_tkind->at(ijet).at(itrack);
        if (trackindex>=3000) continue;
        TLorentzVector TrackP4 = *((TLorentzVector*) l.tk_p4->At(trackindex));
        if (TrackP4.Pt()>1.0) TrackSumP4+=TrackP4;
      }
      if (JetP4.Pt()>30 && TrackSumP4.Pt()>30) {
        GoodJetVector.push_back(JetP4);
      }
    }

    if (GoodJetVector.size()!=1) return false;
    
    return true;
}

// ----------------------------------------------------------------------------------------------------
void PhotonJetAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    vtxAna_.branches(outputTree,"vtx_std_");    

    l.pho_matchingConv = new  std::vector<int>();
    l.Branch_pho_matchingConv(outputTree);
    
    l.vtx_std_evt_mva = new std::vector<float>();
    l.vtx_std_ranked_list = new std::vector<std::vector<int> >();
    l.pho_tkiso_recvtx_030_002_0000_10_01 = new std::vector<std::vector<float> >();
    l.pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01 = new std::vector<std::vector<float> >();
    l.pho_cic6cutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic6passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic6cutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic6passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4cutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4cutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.dipho_vtx_std_sel =  new std::vector<int>();

    l.Branch_vtx_std_evt_mva(outputTree);
    l.Branch_vtx_std_ranked_list(outputTree);
    l.Branch_vtx_std_sel(outputTree);
    l.Branch_pho_tkiso_recvtx_030_002_0000_10_01(outputTree);
    l.Branch_pho_tkiso_badvtx_040_002_0000_10_01(outputTree);
    l.Branch_pho_tkiso_badvtx_id(outputTree);
    l.Branch_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01(outputTree);
    l.Branch_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01(outputTree);
    l.Branch_pho_ZeeVal_tkiso_badvtx_id(outputTree);
    l.Branch_pho_drtotk_25_99(outputTree);

    l.Branch_dipho_n(outputTree);
    l.Branch_dipho_leadind(outputTree);
    l.Branch_dipho_subleadind(outputTree);
    l.Branch_dipho_vtxind(outputTree);
    l.Branch_dipho_sumpt(outputTree);
    
    l.Branch_pho_cic6cutlevel_lead( outputTree );
    l.Branch_pho_cic6passcuts_lead( outputTree );
    l.Branch_pho_cic6cutlevel_sublead( outputTree );
    l.Branch_pho_cic6passcuts_sublead( outputTree );
    l.Branch_pho_cic4cutlevel_lead( outputTree );
    l.Branch_pho_cic4passcuts_lead( outputTree );
    l.Branch_pho_cic4cutlevel_sublead( outputTree );
    l.Branch_pho_cic4passcuts_sublead( outputTree );

    l.Branch_pho_genmatched(outputTree);
    l.Branch_pho_regr_energy_otf(outputTree);
    l.Branch_pho_regr_energyerr_otf(outputTree);
    
    l.Branch_jet_algoPF1_genMatched(outputTree);
    l.Branch_jet_algoPF1_vbfMatched(outputTree);
    l.Branch_jet_algoPF1_genPt(outputTree);
    l.Branch_jet_algoPF1_genDr(outputTree);
    
    l.Branch_jet_algoPF2_genMatched(outputTree);
    l.Branch_jet_algoPF2_vbfMatched(outputTree);
    l.Branch_jet_algoPF2_genPt(outputTree);
    l.Branch_jet_algoPF2_genDr(outputTree);
    
    l.Branch_jet_algoPF3_genMatched(outputTree);
    l.Branch_jet_algoPF3_vbfMatched(outputTree);
    l.Branch_jet_algoPF3_genPt(outputTree);
    l.Branch_jet_algoPF3_genDr(outputTree);
    
    l.gh_higgs_p4 = new TClonesArray("TLorentzVector", 1); 
    l.gh_higgs_p4->Clear();
    ((*l.gh_higgs_p4)[0]) = new TLorentzVector();

    l.gh_pho1_p4 = new TClonesArray("TLorentzVector", 1); 
    l.gh_pho1_p4->Clear();
    ((*l.gh_pho1_p4)[0]) = new TLorentzVector();

    l.gh_pho2_p4 = new TClonesArray("TLorentzVector", 1); 
    l.gh_pho2_p4->Clear();
    ((*l.gh_pho2_p4)[0]) = new TLorentzVector();

    l.gh_vbfq1_p4 = new TClonesArray("TLorentzVector", 1); 
    l.gh_vbfq1_p4->Clear();
    ((*l.gh_vbfq1_p4)[0]) = new TLorentzVector();

    l.gh_vbfq2_p4 = new TClonesArray("TLorentzVector", 1); 
    l.gh_vbfq2_p4->Clear();
    ((*l.gh_vbfq2_p4)[0]) = new TLorentzVector();

    l.gh_vh1_p4 = new TClonesArray("TLorentzVector", 1); 
    l.gh_vh1_p4->Clear();
    ((*l.gh_vh1_p4)[0]) = new TLorentzVector();

    l.gh_vh2_p4 = new TClonesArray("TLorentzVector", 1); 
    l.gh_vh2_p4->Clear();
    ((*l.gh_vh2_p4)[0]) = new TLorentzVector();
    l.Branch_gh_gen2reco1( outputTree );
    l.Branch_gh_gen2reco2( outputTree );
    l.Branch_gh_vbfq1_pdgid( outputTree );
    l.Branch_gh_vbfq2_pdgid( outputTree );
    l.Branch_gh_vh_pdgid( outputTree );
    l.Branch_gh_vh1_pdgid( outputTree );
    l.Branch_gh_vh2_pdgid( outputTree );
    l.Branch_gh_higgs_p4( outputTree );
    l.Branch_gh_pho1_p4( outputTree );
    l.Branch_gh_pho2_p4( outputTree );
    l.Branch_gh_vbfq1_p4( outputTree );
    l.Branch_gh_vbfq2_p4( outputTree );
    l.Branch_gh_vh1_p4( outputTree );
    l.Branch_gh_vh2_p4( outputTree );
}

