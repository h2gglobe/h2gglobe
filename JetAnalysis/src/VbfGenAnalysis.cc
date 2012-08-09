#include "../interface/VbfGenAnalysis.h"

#include "../interface/JetHandler.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#include "CMGTools/External/interface/PileupJetIdentifier.h"

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
VbfGenAnalysis::VbfGenAnalysis()  
{
    name_ = "VbfGenAnalysis";
}

// ----------------------------------------------------------------------------------------------------
VbfGenAnalysis::~VbfGenAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::Term(LoopAll& l) 
{
}

// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::Init(LoopAll& l) 
{
    StatAnalysis::Init(l);
}

// ----------------------------------------------------------------------------------------------------
bool VbfGenAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, 
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

    TLorentzVector lead_p4    = *((TLorentzVector*)l.gh_pho1_p4->At(0));
    TLorentzVector sublead_p4 = *((TLorentzVector*)l.gh_pho2_p4->At(0));
    TLorentzVector jet1       = *((TLorentzVector*)l.gh_vbfq1_p4->At(0));
    TLorentzVector jet2       = *((TLorentzVector*)l.gh_vbfq2_p4->At(0));

    if(jet1.Pt() < jet2.Pt())
      std::swap(jet1, jet2);

    if(lead_p4.Pt() < sublead_p4.Pt())
      std::swap(lead_p4, sublead_p4);

    TLorentzVector diphoton   = lead_p4 + sublead_p4;
    TLorentzVector dijet      = jet1 + jet2;

    myVBF_MVA = -2.;
    VBFevent = true;
    category = 1;

    myVBFLeadJPt= jet1.Pt();
    myVBFSubJPt = jet2.Pt();
    myVBF_Mjj   = dijet.M();
    myVBFdEta   = fabs(jet1.Eta() - jet2.Eta());
    myVBFZep    = fabs(diphoton.Eta() - 0.5*(jet1.Eta() + jet2.Eta()));
    myVBFdPhi   = fabs(diphoton.DeltaPhi(dijet));
    myVBF_Mgg   = diphoton.M();
    myVBFDiPhoPtOverM   = diphoton.Pt()   / myVBF_Mgg;
    myVBFLeadPhoPtOverM = lead_p4.Pt()    / myVBF_Mgg;
    myVBFSubPhoPtOverM  = sublead_p4.Pt() / myVBF_Mgg;
    myVBF_deltaPhiJJ = jet1.DeltaPhi(jet2);
    myVBF_deltaPhiGamGam = lead_p4.DeltaPhi(sublead_p4);
    myVBF_etaJJ = (jet1.Eta() + jet2.Eta())/2;

    TVector3 boost = diphoton.BoostVector();
    TLorentzVector jet1Boosted = jet1, jet2Boosted = jet2, leadingBoosted = lead_p4, subleadingBoosted = sublead_p4;
    jet1Boosted.Boost(-boost);
    jet2Boosted.Boost(-boost);
    leadingBoosted.Boost(-boost);
    subleadingBoosted.Boost(-boost);

    myVBF_thetaJ1 = leadingBoosted.Angle(jet1Boosted.Vect());
    myVBF_thetaJ2 = leadingBoosted.Angle(jet2Boosted.Vect());

    return (myVBFLeadPhoPtOverM > 0.5 && myVBFSubPhoPtOverM > 0.3) && (jet1.Pt() > 20 && jet2.Pt() > 20) && myVBF_Mjj > 50;


    //// // event selection
    //// if( ! skipSelection ) {
    //// 	
    //// 	// first apply corrections and smearing on the single photons 
    //// 	smeared_pho_energy.clear(); smeared_pho_energy.resize(l.pho_n,0.); 
    //// 	smeared_pho_r9.clear();     smeared_pho_r9.resize(l.pho_n,0.); 
    //// 	smeared_pho_weight.clear(); smeared_pho_weight.resize(l.pho_n,1.);
    //// 	applySinglePhotonSmearings(smeared_pho_energy, smeared_pho_r9, smeared_pho_weight, cur_type, l, energyCorrected, energyCorrectedError,
    //// 				   phoSys, syst_shift);
    //// 
    //// 	// inclusive category di-photon selection
    //// 	// FIXME pass smeared R9
    //// 	diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] ); ///??? serve ???
    //// 	
    //// 	// Exclusive Modes
    //// 	int diphotonVBF_id = -1;
    //// 	int diphotonVHhad_id = -1;
    //// 	int diphotonVHlep_id = -1;
    //// 	VHmuevent = false;
    //// 	VHelevent = false;
    //// 	VBFevent = false;
    //// 	VHhadevent = false;
    //// 	
    //// 	// preselection 
    //// 	diphotonVBF_id = l.DiphotonMITPreSelection(leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] );
    //// 
    //// 	diphoton_id = diphotonVBF_id; 
    //// 	if( diphoton_id == -1 ) { return false; }
    //// 
    //// 	/// l.RescaleJetEnergy();
    //// 	/// postProcessJets(l,l.dipho_vtxind[diphoton_id]);
    //// 		
    //// 	TLorentzVector lead_p4, sublead_p4, Higgs;
    //// 	float lead_r9, sublead_r9;
    //// 	TVector3 * vtx;
    //// 	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  
    //// 	if( Higgs.M() < massMin || Higgs.M() > massMax )  { return false; }
    //// 	
    //// 	// clean and sort jets
    //// 	std::vector<int> sorted_jets;
    //// 	for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) { 
    //// 	    TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
    //// 	    if( p4->DeltaR(lead_p4) > 0.5 && p4->DeltaR(sublead_p4) > 0.5 && p4->Pt() > 20. ) {
    //// 		sorted_jets.push_back(ijet);
    //// 	    }
    //// 	}
    //// 	std::sort(sorted_jets.begin(),sorted_jets.end(),
    //// 		  ClonesSorter<TLorentzVector,double,std::greater<double> >(l.jet_algoPF1_p4,&TLorentzVector::Pt));
    //// 
    //// 	switchJetIdVertex( l, l.dipho_vtxind[diphoton_id] );
    //// 
    //// 	/// select two highest pt jets passing loosePfId 
    //// 	int ijet1 = -1;
    //// 	int ijet2 = -1;
    //// 	for(size_t itjet=0; itjet<sorted_jets.size(); ++itjet ) {
    //// 	    int ijet = sorted_jets[itjet];
    //// 	    bool PUjetId = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kLoose);
    //// 	    if ( PUjetId && ijet1<0 ) ijet1 = ijet;
    //// 	    else if ( PUjetId && ijet2<0 ) ijet2 = ijet;
    //// 	    else if ( ijet1!=-1 && ijet2!=-1 ) break;
    //// 	}
    //// 
    //// 	if( ijet1 < 0  || ijet2 < 0) { return false; }
    //// 
    //// 	// now dump variables in a flat tree
    //// 	TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet1);
    //// 	TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet2);
    //// 	TLorentzVector dijet = (*jet1) + (*jet2);
    //// 	TLorentzVector diphoton = lead_p4+sublead_p4;
    //// 
    //// 	if( dumpFlatTree ) {
    //// 	    tree_pho1pt        = lead_p4.Pt();
    //// 	    tree_pho2pt        = sublead_p4.Pt();
    //// 	    tree_diphopt       = diphoton.Pt();
    //// 	    tree_diphoM        = diphoton.M();
    //// 	    tree_diphoEta      = diphoton.Eta();
    //// 
    //// 	    tree_jet1isMatched = l.jet_algoPF1_genMatched[ijet1];
    //// 	    tree_jet1genPt     = l.jet_algoPF1_genPt[ijet1];
    //// 	    tree_jet1genDr     = l.jet_algoPF1_genDr[ijet1];
    //// 	    tree_jet1PileupID  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet1], PileupJetIdentifier::kLoose);
    //// 	    tree_jet1pt        = jet1->Pt();
    //// 	    tree_jet1eta       = jet1->Eta();
    //// 
    //// 	    tree_jet2isMatched = l.jet_algoPF1_genMatched[ijet2];
    //// 	    tree_jet2genPt     = l.jet_algoPF1_genPt[ijet2];
    //// 	    tree_jet2genDr     = l.jet_algoPF1_genDr[ijet2];
    //// 	    tree_jet2PileupID  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet2], PileupJetIdentifier::kLoose);;
    //// 	    tree_jet2pt        = jet2->Pt();
    //// 	    tree_jet2eta       = jet2->Eta();
    //// 
    //// 	    tree_zepp          = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta())); 
    //// 	    tree_mj1j2         = dijet.M();
    //// 	    tree_dphi          = fabs(diphoton.DeltaPhi(dijet));
    //// 	    tree_dijetEta      = dijet.Eta();
    //// 
    //// 	    tree_mctype        = cur_type;
    //// 	    
    //// 	    if ( cur_type < 0 ) 
    //// 		tree_isSignal  = true;
    //// 	    else 
    //// 		tree_isSignal  = false;
    //// 
    //// 	    tree_entry++;
    //// 
    //// 	    flatTree_->Fill();
    //// 	}
    //// 
    //// 
    //// }
}


// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA, 
				    int category, float weight, bool isCorrectVertex, int diphoton_id) 
{
    if( category>0 && fillOptTree ) {
	l.FillTree("run",l.run);
	l.FillTree("lumis",l.lumis);
	l.FillTree("event",l.event);
	l.FillTree("mass",mass);
	l.FillTree("weight",weight);
	l.FillTree("category",category);
	l.FillTree("diphotonMVA",diphotonMVA);
	l.FillTree("vbfMVA",myVBF_MVA);
	l.FillTree("VBFevent", VBFevent);
	if( myVBF_MVA > -2. || VBFevent ) {
	    l.FillTree("deltaPhiJJ",myVBF_deltaPhiJJ);
	    l.FillTree("deltaPhiGamGam", myVBF_deltaPhiGamGam);
	    l.FillTree("etaJJ", myVBF_etaJJ);
	    l.FillTree("thetaJ1", myVBF_thetaJ1);
	    l.FillTree("thetaJ2", myVBF_thetaJ2);

	    l.FillTree("leadJPt", myVBFLeadJPt);
	    l.FillTree("subleadJPt", myVBFSubJPt);
	    l.FillTree("MJJ", myVBF_Mjj);
	    l.FillTree("deltaEtaJJ", myVBFdEta);
	    l.FillTree("Zep", myVBFZep);
	    l.FillTree("deltaPhiJJGamGam", myVBFdPhi);
	    l.FillTree("MGamGam", myVBF_Mgg);
	    l.FillTree("diphoPtOverM", myVBFDiPhoPtOverM);
	    l.FillTree("leadPtOverM", myVBFLeadPhoPtOverM);
	    l.FillTree("subleadPtOverM", myVBFSubPhoPtOverM);
	}
	l.FillTree("sampleType",cur_type);
	//// l.FillTree("isCorrectVertex",isCorrectVertex);
	//// l.FillTree("metTag",VHmetevent);
	//// l.FillTree("eleTag",VHelevent);
	//// l.FillTree("muTag",VHmuevent);
	
	TLorentzVector lead_p4, sublead_p4, Higgs;
	float lead_r9, sublead_r9;
	TVector3 * vtx;
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  
	l.FillTree("leadPt",(float)lead_p4.Pt());
	l.FillTree("subleadPt",(float)sublead_p4.Pt());
	l.FillTree("leadR9",lead_r9);
	l.FillTree("subleadR9",sublead_r9);
	l.FillTree("sigmaMrv",sigmaMrv);
	l.FillTree("sigmaMwv",sigmaMwv);

	vtxAna_.setPairID(diphoton_id);
	float vtxProb = vtxAna_.vertexProbability(l.vtx_std_evt_mva->at(diphoton_id), l.vtx_std_n);
	float altMass = 0.;
	if( l.vtx_std_n > 1 ) {
	    int altvtx = (*l.vtx_std_ranked_list)[diphoton_id][1];
	    altMass = ( l.get_pho_p4( l.dipho_leadind[diphoton_id], altvtx, &smeared_pho_energy[0]) + 
			l.get_pho_p4( l.dipho_subleadind[diphoton_id], altvtx, &smeared_pho_energy[0]) ).M();
	}
	l.FillTree("altMass",altMass);
	l.FillTree("vtxProb",vtxProb);
    }
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
