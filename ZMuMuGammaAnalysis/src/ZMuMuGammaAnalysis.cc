#include "../interface/ZMuMuGammaAnalysis.h"

#include "PhotonReducedInfo.h"
#include "Sorters.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "JetAnalysis/interface/JetHandler.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

#define PADEBUG 0

using namespace std;

ZMuMuGammaAnalysis::ZMuMuGammaAnalysis(){}
ZMuMuGammaAnalysis::~ZMuMuGammaAnalysis(){}

// ----------------------------------------------------------------------------------------------------
void ZMuMuGammaAnalysis::Init(LoopAll& l)
{
    PhotonAnalysis::Init(l);

    l.SetAllMVA();
    cout << "Weights file is: " << photonLevel2012IDMVA_EB.c_str() << endl;
    l.tmvaReaderID_2013_Barrel->BookMVA("AdaBoost",photonLevel2012IDMVA_EB.c_str());
    l.tmvaReaderID_2013_Endcap->BookMVA("AdaBoost",photonLevel2012IDMVA_EE.c_str());
    
}

//----------------------------------------------------------------------------------------------------
bool ZMuMuGammaAnalysis::SkimEvents(LoopAll& l, int jentry)
{
    if( run7TeV4Xanalysis ) { l.version=12; }
    else { l.version=13; }

    return true;
}

// ----------------------------------------------------------------------------------------------------
bool ZMuMuGammaAnalysis::SelectEventsReduction(LoopAll& l, int jentry)
{

    if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
    // require at least two reconstructed photons to store the event

    //if( pho_acc.size() < 2 ) { return false; }

    vtxAna_.clear();
    l.vtx_std_ranked_list->clear();
    l.dipho_vtx_std_sel->clear();
    l.vtx_std_ranked_list->clear();
    l.vtx_std_evt_mva->clear();
    l.vtx_std_sel=0;
    l.pho_mitmva->clear();
    float maxSumPt = 0.;
    l.dipho_n = 0;
    bool oneKinSelected = true;

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

    //Calculate cluster shape variables prior to shape rescaling
    for (int ipho=0;ipho<l.pho_n;ipho++){
      l.pho_s4ratio[ipho] = l.pho_e2x2[ipho]/l.bc_s25[l.sc_bcseedind[l.pho_scind[ipho]]];
      float rr2=l.pho_eseffsixix[ipho]*l.pho_eseffsixix[ipho]+l.pho_eseffsiyiy[ipho]*l.pho_eseffsiyiy[ipho];
      l.pho_ESEffSigmaRR[ipho] = 0.0;
      if(rr2>0. && rr2<999999.) {
        l.pho_ESEffSigmaRR[ipho] = sqrt(rr2);
      }
    }

    for (int ipho=0; ipho<l.pho_n; ipho++) {
      vector<float> MVAValues;
      for(int ivtx=0; ivtx<l.vtx_std_n; ++ivtx) {
        TLorentzVector pho_p4 = l.get_pho_p4(ipho, ivtx);
        MVAValues.push_back(l.photonIDMVA(ipho, ivtx, pho_p4, bdtTrainingType.c_str()));
      }
      l.pho_mitmva->push_back(MVAValues);
    }
        
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

    if( pho_presel.size() < 2 ) {
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
        l.dipho_n = 0;
        for(size_t id=0; id<diphotons.size(); ++id ) {

	    if( l.dipho_n >= MAX_DIPHOTONS-1 ) { continue; }
            int ipho1 = diphotons[id].first;
            int ipho2 = diphotons[id].second;

            if(PADEBUG)        cout << " SelectEventsReduction going to fill photon info " << endl;
            PhotonInfo pho1=l.fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions,&corrected_pho_energy[0]);
            PhotonInfo pho2=l.fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions,&corrected_pho_energy[0]);
            if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;

            l.vertexAnalysis(vtxAna_, pho1, pho2 );
            std::vector<int> vtxs = l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames, mvaVertexSelection,
                              tmvaPerVtxReader_, tmvaPerVtxMethod);

            TLorentzVector lead_p4 = l.get_pho_p4( ipho2, vtxs[0], &corrected_pho_energy[0] ).Pt();
            TLorentzVector sublead_p4 = l.get_pho_p4( ipho1, vtxs[0], &corrected_pho_energy[0] ).Pt();

            if(sublead_p4.Pt()  > lead_p4.Pt() ) {
                std::swap( diphotons[id].first,  diphotons[id].second );
                std::swap( lead_p4,  sublead_p4 );
            }

            if( lead_p4.Pt() < presel_scet1 || sublead_p4.Pt() < presel_scet2 ||
                fabs(lead_p4.Eta()) > presel_maxeta || fabs(sublead_p4.Eta()) > presel_maxeta ) {
                vtxAna_.discardLastDipho();
                continue;
            }
	    oneKinSelected = true;

            if( ! l.PhotonMITPreSelection(ipho1, vtxs[0], &corrected_pho_energy[0] )
                || ! l.PhotonMITPreSelection(ipho2, vtxs[0], &corrected_pho_energy[0] ) ) {
                vtxAna_.discardLastDipho();
                continue;
            }

            l.vtx_std_ranked_list->push_back(vtxs);
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

            l.dipho_leadind[l.dipho_n] = diphotons[id].first;
            l.dipho_subleadind[l.dipho_n] = diphotons[id].second;
            l.dipho_vtxind[l.dipho_n] = l.dipho_vtx_std_sel->back();

            l.dipho_sumpt[l.dipho_n] = lead_p4.Pt() + sublead_p4.Pt();

            if( l.dipho_sumpt[l.dipho_n] > maxSumPt ) {
                l.vtx_std_sel = l.dipho_vtx_std_sel->back();
                maxSumPt = l.dipho_sumpt[l.dipho_n];
            }

            // make sure that vertex analysis indexes are in synch
            assert( l.dipho_n == vtxAna_.pairID(ipho1,ipho2) );

            l.dipho_n++;
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

// ----------------------------------------------------------------------------------------------------
void ZMuMuGammaAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree)
{
    if( outputTree ) { vtxAna_.branches(outputTree,"vtx_std_"); }

    l.pho_matchingConv = new  std::vector<int>();
    if( outputTree ) { l.Branch_pho_matchingConv(outputTree); }

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
    l.pho_cic4pfcutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4pfpasscuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4pfcutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4pfpasscuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_mitmva = new std::vector<std::vector<float> >();
    l.dipho_vtx_std_sel =  new std::vector<int>();

    if( outputTree ) {

	l.Branch_vtx_std_evt_mva(outputTree);
	l.Branch_vtx_std_ranked_list(outputTree);
	l.Branch_vtx_std_sel(outputTree);
	l.Branch_pho_tkiso_recvtx_030_002_0000_10_01(outputTree);
	l.Branch_pho_tkiso_badvtx_040_002_0000_10_01(outputTree);
	l.Branch_pho_tkiso_badvtx_id(outputTree);
	l.Branch_pho_pfiso_charged_badvtx_04(outputTree);
	l.Branch_pho_pfiso_charged_badvtx_id(outputTree);
	l.Branch_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01(outputTree);
	l.Branch_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01(outputTree);
	l.Branch_pho_ZeeVal_tkiso_badvtx_id(outputTree);
	l.Branch_pho_mitmva(outputTree);
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
	l.Branch_pho_cic4pfcutlevel_lead( outputTree );
	l.Branch_pho_cic4pfpasscuts_lead( outputTree );
	l.Branch_pho_cic4pfcutlevel_sublead( outputTree );
	l.Branch_pho_cic4pfpasscuts_sublead( outputTree );

	l.Branch_pho_genmatched(outputTree);
	l.Branch_pho_regr_energy_otf(outputTree);
	l.Branch_pho_regr_energyerr_otf(outputTree);

	l.Branch_jet_algoPF1_genMatched(outputTree);
	l.Branch_jet_algoPF1_bgenMatched(outputTree);
	l.Branch_jet_algoPF1_cgenMatched(outputTree);
	l.Branch_jet_algoPF1_lgenMatched(outputTree);
	l.Branch_jet_algoPF1_vbfMatched(outputTree);
	l.Branch_jet_algoPF1_genPt(outputTree);
	l.Branch_jet_algoPF1_genDr(outputTree);

	//l.Branch_jet_algoPF2_genMatched(outputTree);
	//l.Branch_jet_algoPF2_vbfMatched(outputTree);
	//l.Branch_jet_algoPF2_genPt(outputTree);
	//l.Branch_jet_algoPF2_genDr(outputTree);

	l.Branch_jet_algoPF3_genMatched(outputTree);
	l.Branch_jet_algoPF3_bgenMatched(outputTree);
	l.Branch_jet_algoPF3_cgenMatched(outputTree);
	l.Branch_jet_algoPF3_lgenMatched(outputTree);
	l.Branch_jet_algoPF3_vbfMatched(outputTree);
	l.Branch_jet_algoPF3_genPt(outputTree);
	l.Branch_jet_algoPF3_genDr(outputTree);

	//correctMETinRED
	l.Branch_shiftMET_pt(outputTree);
	l.Branch_shiftMET_phi(outputTree);
	l.Branch_smearMET_pt(outputTree);
	l.Branch_smearMET_phi(outputTree);
	l.Branch_shiftsmearMET_pt(outputTree);
	l.Branch_shiftsmearMET_phi(outputTree);
	l.Branch_shiftscaleMET_pt(outputTree);
	l.Branch_shiftscaleMET_phi(outputTree);
	l.Branch_shiftMET_eta(outputTree);
	l.Branch_shiftMET_e(outputTree);
	l.Branch_shiftscaleMET_eta(outputTree);
	l.Branch_shiftscaleMET_e(outputTree);
    }

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

//    l.METcorrected = new TClonesArray("TLorentzVector", 1);     //met at analysis step
//    l.METcorrected->Clear();                    //met at analysis step
//    ((*l.METcorrected)[0]) = new TLorentzVector();        //met at analysis step


    if( outputTree ) {
    l.Branch_gh_gen2reco1( outputTree );
    l.Branch_gh_gen2reco2( outputTree );
    l.Branch_gh_vbfq1_pdgid( outputTree );
    l.Branch_gh_vbfq2_pdgid( outputTree );
    l.Branch_gh_vh_pdgid( outputTree );
    l.Branch_gh_vh1_pdgid( outputTree );
    l.Branch_gh_vh2_pdgid( outputTree );
//    l.Branch_METcorrected( outputTree );  //met at analysis step
    l.Branch_gh_higgs_p4( outputTree );
    l.Branch_gh_pho1_p4( outputTree );
    l.Branch_gh_pho2_p4( outputTree );
    l.Branch_gh_vbfq1_p4( outputTree );
    l.Branch_gh_vbfq2_p4( outputTree );
    l.Branch_gh_vh1_p4( outputTree );
    l.Branch_gh_vh2_p4( outputTree );
    }
    l.Branch_mu_glo_hasgsftrack(outputTree);
}
