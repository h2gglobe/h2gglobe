#include "../interface/VbfAnalysis.h"

#include "../interface/JetHandler.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#include "CMGTools/External/interface/PileupJetIdentifier.h"

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
VbfAnalysis::VbfAnalysis()  
{
    name_ = "VbfAnalysis";

    recomputeJetId = false;
    expoMatching   = false;
    dumpFlatTree   = false;
    
    flatTree_ = 0;
    outputFile_ = 0;
}

// ----------------------------------------------------------------------------------------------------
VbfAnalysis::~VbfAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void VbfAnalysis::Term(LoopAll& l) 
{
    /// l.outputFile->cd();
    if( dumpFlatTree ) {
	if( outputFile_ ) {
	    outputFile_->cd();
	} else {
	    l.outputFile->cd();
	}
	flatTree_->Write();
	if( outputFile_ ) {
	    outputFile_->Close();
	}
    }

}

// ----------------------------------------------------------------------------------------------------
void VbfAnalysis::Init(LoopAll& l) 
{
    MassFactorizedMvaAnalysis::Init(l);
    if( jetHandler_ == 0 ) {
    	jetHandler_ = new JetHandler(jetHandlerCfg, l);
    }

    if( dumpFlatTree ) {
	if( flatTree_ == 0 ) { 
	    outputFile_ = TFile::Open("vbfAnalysisTree_"+l.histFileName,"recreate");
	    flatTree_ = new TTree("flatTree","flatTree");
	    tree_entry = 0;
	}
	flatTree_->Branch( "entry",         &tree_entry );
	flatTree_->Branch( "pho1pt",        &tree_pho1pt );
	flatTree_->Branch( "pho2pt",        &tree_pho2pt );
	flatTree_->Branch( "diphopt",       &tree_diphopt );
	flatTree_->Branch( "diphoM",        &tree_diphoM );
	flatTree_->Branch( "diphoEta",      &tree_diphoEta );
	flatTree_->Branch( "dijetEta",      &tree_dijetEta );
	flatTree_->Branch( "jet1isMatched", &tree_jet1isMatched );
	flatTree_->Branch( "jet2isMatched", &tree_jet2isMatched );
	flatTree_->Branch( "jet1genPt",     &tree_jet1genPt );
	flatTree_->Branch( "jet2genPt",     &tree_jet2genPt );
	flatTree_->Branch( "jet1genDr",     &tree_jet1genDr );
	flatTree_->Branch( "jet2genDr",     &tree_jet2genDr );
	flatTree_->Branch( "jet1pt",        &tree_jet1pt );
	flatTree_->Branch( "jet2pt",        &tree_jet2pt );
	flatTree_->Branch( "jet1eta",       &tree_jet1eta );
	flatTree_->Branch( "jet2eta",       &tree_jet2eta );
	flatTree_->Branch( "zepp",          &tree_zepp );
	flatTree_->Branch( "mj1j2",         &tree_mj1j2 );
	flatTree_->Branch( "dphi",          &tree_dphi);
	flatTree_->Branch( "dphiJJ",        &tree_dphiJJ);
	flatTree_->Branch( "dphiJJ2",       &tree_dphiJJ2);
	flatTree_->Branch( "deltaEta3",     &tree_deltaEta3);
	flatTree_->Branch( "jet1PileupID",  &tree_jet1PileupID );
	flatTree_->Branch( "jet2PileupID",  &tree_jet2PileupID );
	flatTree_->Branch( "isSignal",      &tree_isSignal );
	flatTree_->Branch( "mctype",        &tree_mctype );
   }
}

// ----------------------------------------------------------------------------------------------------
void VbfAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    dumpFlatTree=true;
    flatTree_ = new TTree("flatTree","flatTree");
}

// ----------------------------------------------------------------------------------------------------
void VbfAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
}
   
// ----------------------------------------------------------------------------------------------------
bool VbfAnalysis::SelectEventsReduction(LoopAll&, int)
{
    return true;
}

//void switchJetIdVertex(LoopAll &l, int ivtx) ;

// ----------------------------------------------------------------------------------------------------
bool VbfAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, 
				 int & category, int & diphoton_id,
				 bool & isCorrectVertex,
				 float &kinematic_bdtout,
				 bool isSyst, 
				 float syst_shift, bool skipSelection,
				 BaseGenLevelSmearer *genSys, BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys)
{
    assert( isSyst || ! skipSelection );

    int cur_type = l.itype[l.current];
    float sampleweight = l.sampleContainer[l.current_sample_index].weight;
    /// diphoton_id = -1;
    
    std::pair<int,int> diphoton_index;
   
    // do gen-level dependent first (e.g. k-factor); only for signal
    genLevWeight=1.;
    if(cur_type!=0 ) {
	applyGenLevelSmearings(genLevWeight,gP4,l.pu_n,cur_type,genSys,syst_shift);
    }

    // event selection
    if( ! skipSelection ) {
	
	// first apply corrections and smearing on the single photons 
	smeared_pho_energy.clear(); smeared_pho_energy.resize(l.pho_n,0.); 
	smeared_pho_r9.clear();     smeared_pho_r9.resize(l.pho_n,0.); 
	smeared_pho_weight.clear(); smeared_pho_weight.resize(l.pho_n,1.);
	applySinglePhotonSmearings(smeared_pho_energy, smeared_pho_r9, smeared_pho_weight, cur_type, l, energyCorrected, energyCorrectedError,
				   phoSys, syst_shift);

	// inclusive category di-photon selection
	// FIXME pass smeared R9
	diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] ); ///??? serve ???
    	
	// Exclusive Modes
	int diphotonVBF_id = -1;
	int diphotonVHhad_id = -1;
	int diphotonVHlep_id = -1;
	VHmuevent = false;
	VHelevent = false;
	VBFevent = false;
	VHhadevent = false;
	
	// preselection 
	diphotonVBF_id = l.DiphotonMITPreSelection(leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] );

	diphoton_id = diphotonVBF_id; 
	if( diphoton_id == -1 ) { return false; }

	/// l.RescaleJetEnergy();
	/// postProcessJets(l,l.dipho_vtxind[diphoton_id]);
		
	TLorentzVector lead_p4, sublead_p4, Higgs;
	float lead_r9, sublead_r9;
	TVector3 * vtx;
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  
	if( Higgs.M() < massMin || Higgs.M() > massMax )  { return false; }
	
	// clean and sort jets
	std::vector<int> sorted_jets;
	for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) { 
	    TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
	    if( p4->DeltaR(lead_p4) > 0.5 && p4->DeltaR(sublead_p4) > 0.5 && p4->Pt() > 20. ) {
		sorted_jets.push_back(ijet);
	    }
	}
	std::sort(sorted_jets.begin(),sorted_jets.end(),
		  ClonesSorter<TLorentzVector,double,std::greater<double> >(l.jet_algoPF1_p4,&TLorentzVector::Pt));

	switchJetIdVertex( l, l.dipho_vtxind[diphoton_id] );

	/// select two highest pt jets passing loose PU jet Id
	int ijet1 = -1;
	int ijet2 = -1;
	int ijet3 = -1;
	for(size_t itjet=0; itjet<sorted_jets.size(); ++itjet ) {
	    int ijet = sorted_jets[itjet];
	    bool PUjetId = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kLoose);
	    if ( PUjetId && ijet1<0 ) ijet1 = ijet;
	    else if ( PUjetId && ijet2<0 ) ijet2 = ijet;
	    else if ( PUjetId && ijet3<0 ) ijet3 = ijet;
	    //else if ( ijet1!=-1 && ijet2!=-1 ) break;
	}

	if( ijet1 < 0  || ijet2 < 0) { return false; }

	// now dump variables in a flat tree
	TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet1);
	TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet2);
	TLorentzVector* jet3;
	if (ijet3>0) jet3 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet3);
	TLorentzVector dijet = (*jet1) + (*jet2);
	TLorentzVector diphoton = lead_p4+sublead_p4;
    

	//--- compute dphiJJ2 [http://arxiv.org/pdf/1001.3822.pdf]
	TLorentzVector* sumj1 = new TLorentzVector() ;
	TLorentzVector* sumj2 = new TLorentzVector();
	sumj1->SetPxPyPzE(0.,0.,0.,0.);
	sumj2->SetPxPyPzE(0.,0.,0.,0.);
	for(size_t itjet=0; itjet<sorted_jets.size(); ++itjet ) {
	    int ijet = sorted_jets[itjet];
	    bool PUjetId = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kLoose);
	    if (PUjetId==false) continue;
	    TLorentzVector* jet = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
	    if ( jet->Eta() < diphoton.Eta() ) (*sumj1) = (*sumj1) + (*jet);
	    if ( jet->Eta() > diphoton.Eta() ) (*sumj2) = (*sumj2) + (*jet);
	}


	if( dumpFlatTree ) {
	    tree_pho1pt        = lead_p4.Pt();
	    tree_pho2pt        = sublead_p4.Pt();
	    tree_diphopt       = diphoton.Pt();
	    tree_diphoM        = diphoton.M();
	    tree_diphoEta      = diphoton.Eta();
  
	    tree_jet1isMatched = l.jet_algoPF1_genMatched[ijet1];
	    tree_jet1genPt     = l.jet_algoPF1_genPt[ijet1];
	    tree_jet1genDr     = l.jet_algoPF1_genDr[ijet1];
	    tree_jet1PileupID  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet1], PileupJetIdentifier::kLoose);
	    tree_jet1pt        = jet1->Pt();
	    tree_jet1eta       = jet1->Eta();

	    tree_jet2isMatched = l.jet_algoPF1_genMatched[ijet2];
	    tree_jet2genPt     = l.jet_algoPF1_genPt[ijet2];
	    tree_jet2genDr     = l.jet_algoPF1_genDr[ijet2];
	    tree_jet2PileupID  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet2], PileupJetIdentifier::kLoose);
	    tree_jet2pt        = jet2->Pt();
	    tree_jet2eta       = jet2->Eta();

	    tree_zepp          = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta())); 
	    tree_mj1j2         = dijet.M();
	    tree_dphi          = fabs(diphoton.DeltaPhi(dijet));
	    tree_dijetEta      = dijet.Eta();
	    tree_dphiJJ        = fabs(jet1->DeltaPhi(*jet2));
	    tree_dphiJJ2       = fabs(sumj1->DeltaPhi(*sumj2));
	    
	    if (ijet3 > 0 ) 
		tree_deltaEta3 = fabs(jet3->Eta() - 0.5*(jet1->Eta() + jet2->Eta())); 
	    else 
		tree_deltaEta3 = -999; 
		
	    if ( cur_type < 0 ) 
		tree_isSignal  = true;
	    else 
		tree_isSignal  = false;

	    tree_mctype        = cur_type;

	    tree_entry++;

	    flatTree_->Fill();
	}


    }
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
