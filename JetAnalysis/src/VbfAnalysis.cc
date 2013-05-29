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
TreeVariables::TreeVariables() : 
    entry(0),		
    weight(0),		
    sampleweight(0),		
    pho1pt(0),		
    pho2pt(0),		
    diphopt(0),		
    diphoM(0),		
    diphoEta(999),	
    dijetEta(999),	
    jet1isMatched(false),	
    jet2isMatched(false),	
    jet1genPt(0),	
    jet2genPt(0),	
    jet1genDr(0),	
    jet2genDr(0),	
    jet1Pt(0),		
    jet2Pt(0),		
    jet1Eta(999),		
    jet2Eta(999),		
    zepp(999),		
    mj1j2(-999),		
    dphi(999),		
    dphiJJ(999),		
    dphiJJ2(999),		
    deltaEta3(999),	
    jet1PileupID(0),	
    jet2PileupID(0),	
    isSignal(false),	
    mctype(0),		
    diphomva(-999),	
    pho1Eta(999),		
    pho2Eta(999),		
    pho1r9(999),		
    pho2r9(999),		
    pho1idMva(-999),	
    pho2idMva(-999),	
    pho1sEoverE(999),	
    pho2sEoverE(999),	
    pho1sEoverEsmear(999),	
    pho2sEoverEsmear(999),	
    diphosMoverM(999),	
    diphosMoverMwv(999),	
    diphovtxProb(-1),
    jet1(-1),
    jet2(-1),
    jet3(-1) 
{}

// ----------------------------------------------------------------------------------------------------
VbfAnalysis::VbfAnalysis()  
{
    name_ = "VbfAnalysis";

    recomputeJetId = false;
    expoMatching   = false;
    dumpFlatTree   = false;
    requireTwoJets = true;
    
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
    if( dumpFlatTree ) {
	//// if( outputFile_ ) {
	////     outputFile_->cd();
	//// } else {
	////     l.outputFile->cd();
	//// }
	//// flatTree_->Write();
	//// if( outputFile_ ) {
	////     outputFile_->Close();
	//// }
    }
}

// ----------------------------------------------------------------------------------------------------
void VbfAnalysis::Init(LoopAll& l) 
{  
    if( l.runZeeValidation ) {
	leadEtCut = 15;
	subleadEtCut = 15;
	leadEtVBFCut = 15.;
	subleadEtVBFCut = 15.;
	massMin = 60, massMax = 120.;
	applyPtoverM = false;
    }
    

    MassFactorizedMvaAnalysis::Init(l);
    if( jetHandler_ == 0 ) {
    	jetHandler_ = new JetHandler(jetHandlerCfg, l);
    }

    if( dumpFlatTree ) {
	if( flatTree_ == 0 ) { 
	    //// outputFile_ = TFile::Open("vbfAnalysisTree_"+l.histFileName,"recreate");
	    //// flatTree_ = new TTree("flatTree","flatTree");
	    l.InitTrees("vbfAnalysis");
	    tree_.entry = 0;
	}
	
	l.BookExternalTreeBranch( "entry",         &tree_.entry, "vbfAnalysis" );         
	l.BookExternalTreeBranch( "weight",        &tree_.entry, "vbfAnalysis" );         
	l.BookExternalTreeBranch( "sampleweight",  &tree_.entry, "vbfAnalysis" );         
	l.BookExternalTreeBranch( "pho1pt",        &tree_.pho1pt, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "pho2pt",        &tree_.pho2pt, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "diphopt",       &tree_.diphopt, "vbfAnalysis" );       
	l.BookExternalTreeBranch( "diphoM",        &tree_.diphoM, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "diphoEta",      &tree_.diphoEta, "vbfAnalysis" );      
	l.BookExternalTreeBranch( "dijetEta",      &tree_.dijetEta, "vbfAnalysis" );      
	l.BookExternalTreeBranch( "jet1",          &tree_.jet1, "vbfAnalysis" ); 
	l.BookExternalTreeBranch( "jet2",          &tree_.jet2, "vbfAnalysis" ); 
	l.BookExternalTreeBranch( "jet3",          &tree_.jet3, "vbfAnalysis" ); 
	l.BookExternalTreeBranch( "jet1isMatched", &tree_.jet1isMatched, "vbfAnalysis" ); 
	l.BookExternalTreeBranch( "jet2isMatched", &tree_.jet2isMatched, "vbfAnalysis" ); 
	l.BookExternalTreeBranch( "jet1genPt",     &tree_.jet1genPt, "vbfAnalysis" );     
	l.BookExternalTreeBranch( "jet2genPt",     &tree_.jet2genPt, "vbfAnalysis" );     
	l.BookExternalTreeBranch( "jet1genDr",     &tree_.jet1genDr, "vbfAnalysis" );     
	l.BookExternalTreeBranch( "jet2genDr",     &tree_.jet2genDr, "vbfAnalysis" );     
	l.BookExternalTreeBranch( "jet1Pt",        &tree_.jet1Pt, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "jet2Pt",        &tree_.jet2Pt, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "jet1Eta",       &tree_.jet1Eta, "vbfAnalysis" );       
	l.BookExternalTreeBranch( "jet2Eta",       &tree_.jet2Eta, "vbfAnalysis" );       
	l.BookExternalTreeBranch( "zepp",          &tree_.zepp, "vbfAnalysis" );          
	l.BookExternalTreeBranch( "mj1j2",         &tree_.mj1j2, "vbfAnalysis" );         
	l.BookExternalTreeBranch( "dphi",          &tree_.dphi, "vbfAnalysis" );          
	l.BookExternalTreeBranch( "dphiJJ",        &tree_.dphiJJ, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "dphiJJ2",       &tree_.dphiJJ2, "vbfAnalysis" );       
	l.BookExternalTreeBranch( "deltaEta3",     &tree_.deltaEta3, "vbfAnalysis" );     
	l.BookExternalTreeBranch( "jet1PileupID",  &tree_.jet1PileupID, "vbfAnalysis" );  
	l.BookExternalTreeBranch( "jet2PileupID",  &tree_.jet2PileupID, "vbfAnalysis" );  
	l.BookExternalTreeBranch( "isSignal",      &tree_.isSignal, "vbfAnalysis" );      
	l.BookExternalTreeBranch( "mctype",        &tree_.mctype, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "diphomva",      &tree_.diphomva, "vbfAnalysis" );      
	l.BookExternalTreeBranch( "pho1Eta",       &tree_.pho1Eta, "vbfAnalysis" );       
	l.BookExternalTreeBranch( "pho2Eta",       &tree_.pho2Eta, "vbfAnalysis" );       
	l.BookExternalTreeBranch( "pho1r9",        &tree_.pho1r9, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "pho2r9",        &tree_.pho2r9, "vbfAnalysis" );        
	l.BookExternalTreeBranch( "pho1idMva",     &tree_.pho1idMva, "vbfAnalysis" );    
	l.BookExternalTreeBranch( "pho2idMva",     &tree_.pho2idMva, "vbfAnalysis" );    
	l.BookExternalTreeBranch( "pho1sEoverE",   &tree_.pho1sEoverE, "vbfAnalysis" );   
	l.BookExternalTreeBranch( "pho2sEoverE",   &tree_.pho2sEoverE, "vbfAnalysis" );   
	l.BookExternalTreeBranch( "pho1sEoverEsmear",   &tree_.pho1sEoverEsmear, "vbfAnalysis" );   
	l.BookExternalTreeBranch( "pho2sEoverEsmear",   &tree_.pho2sEoverEsmear, "vbfAnalysis" );   
	l.BookExternalTreeBranch( "diphosMoverM",  &tree_.diphosMoverM, "vbfAnalysis" );  
	l.BookExternalTreeBranch( "diphosMoverMwv",&tree_.diphosMoverMwv, "vbfAnalysis" );
	l.BookExternalTreeBranch( "diphovtxProb",  &tree_.diphovtxProb, "vbfAnalysis" );  
   }
}

// ----------------------------------------------------------------------------------------------------
void VbfAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    /// dumpFlatTree=true;
    /// flatTree_ = new TTree("flatTree","flatTree");
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
	
	//// // preselection 
	//// diphotonVBF_id = l.DiphotonMITPreSelection(leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] );
	//// diphoton_id = diphotonVBF_id; 
	
	if( diphoton_id == -1 ) { return false; }
	
	/// l.RescaleJetEnergy();
	/// postProcessJets(l,l.dipho_vtxind[diphoton_id]);
		
	TLorentzVector lead_p4, sublead_p4, Higgs;
	float lead_r9, sublead_r9;
	TVector3 * vtx;
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  
	if( Higgs.M() < massMin || Higgs.M() > massMax )  { return false; }
	
	
	// ---- evaluate dipho MVA
	diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );

        // Mass Resolution of the Event
        massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories,beamspotSigma);
        float vtx_mva  = l.vtx_std_evt_mva->at(diphoton_id);
        sigmaMrv = massResolutionCalculator->massResolutionEonly();
        sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
        float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
        // easy to calculate vertex probability from vtx mva output
        float vtxProb   = 1.-0.49*(vtx_mva+1.0); 

        float phoid_mvaout_lead = ( dataIs2011 ? 
				    l.photonIDMVA(diphoton_index.first,l.dipho_vtxind[diphoton_id],
						  lead_p4,bdtTrainingPhilosophy.c_str()) :
				    l.photonIDMVANew(diphoton_index.first,l.dipho_vtxind[diphoton_id],
						     lead_p4,bdtTrainingPhilosophy.c_str()) );
        float phoid_mvaout_sublead = ( dataIs2011 ? 
				       l.photonIDMVA(diphoton_index.second,l.dipho_vtxind[diphoton_id],
						     sublead_p4,bdtTrainingPhilosophy.c_str()) : 
				       l.photonIDMVANew(diphoton_index.second,l.dipho_vtxind[diphoton_id],
							sublead_p4,bdtTrainingPhilosophy.c_str()) );
	// apply di-photon level smearings and corrections
        int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
        if( cur_type != 0 && doMCSmearing ) {
	    applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, 
				   phoid_mvaout_lead,phoid_mvaout_sublead,
				   diPhoSys, syst_shift);
            isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
        }

	float diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id] ,
					      vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,
					      bdtTrainingPhilosophy.c_str(),
					      phoid_mvaout_lead,phoid_mvaout_sublead);
	
		
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
	
	TLorentzVector* jet1=(ijet1>=0?(TLorentzVector*)l.jet_algoPF1_p4->At(ijet1):0);
	TLorentzVector* jet2=(ijet2>=0?(TLorentzVector*)l.jet_algoPF1_p4->At(ijet2):0);
	TLorentzVector* jet3=(ijet3>=0?(TLorentzVector*)l.jet_algoPF1_p4->At(ijet3):0);
	TLorentzVector sumj1;
	TLorentzVector sumj2;
	TLorentzVector dijet;
	
	TLorentzVector diphoton = lead_p4+sublead_p4;

	if( ijet1 >= 0 && ijet2 >= 0) {

	    dijet = (*jet1) + (*jet2);
	    //--- compute dphiJJ2 [http://arxiv.org/pdf/1001.3822.pdf]
	    sumj1.SetPxPyPzE(0.,0.,0.,0.);
	    sumj2.SetPxPyPzE(0.,0.,0.,0.);
	    for(size_t itjet=0; itjet<sorted_jets.size(); ++itjet ) {
		int ijet = sorted_jets[itjet];
		bool PUjetId = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kLoose);
		if (PUjetId==false) continue;
		TLorentzVector* jet = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
		if ( jet->Eta() < diphoton.Eta() ) sumj1 = sumj1 + (*jet);
		if ( jet->Eta() > diphoton.Eta() ) sumj2 = sumj1 + (*jet);
	    }
	} else if( requireTwoJets ) {
	    return false;
	}

	// now dump variables in a flat tree
	if( dumpFlatTree ) {
	    tree_ = default_;
	    
	    tree_.weight        = evweight;
	    tree_.sampleweight  = l.weight;
	    tree_.pho1pt        = lead_p4.Pt();
	    tree_.pho2pt        = sublead_p4.Pt();
	    tree_.diphopt       = diphoton.Pt();
	    tree_.diphoM        = diphoton.M();
	    tree_.diphoEta      = diphoton.Eta();
	    tree_.diphomva      = diphobdt_output;
	    
	    tree_.jet1 = ijet1;
	    tree_.jet2 = ijet2;
	    tree_.jet3 = ijet3;
	    
	    if( ijet1 > 0 ) {
		tree_.jet1isMatched = l.jet_algoPF1_genMatched[ijet1];
		tree_.jet1genPt     = l.jet_algoPF1_genPt[ijet1];
		tree_.jet1genDr     = l.jet_algoPF1_genDr[ijet1];
		tree_.jet1PileupID  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet1], PileupJetIdentifier::kLoose);
		tree_.jet1Pt        = jet1->Pt();
		tree_.jet1Eta       = jet1->Eta();
	    }
	    if( ijet2 > 0 ) {
		tree_.jet2isMatched = l.jet_algoPF1_genMatched[ijet2];
		tree_.jet2genPt     = l.jet_algoPF1_genPt[ijet2];
		tree_.jet2genDr     = l.jet_algoPF1_genDr[ijet2];
		tree_.jet2PileupID  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet2], PileupJetIdentifier::kLoose);
		tree_.jet2Pt        = jet2->Pt();
		tree_.jet2Eta       = jet2->Eta();
	    }
	    
	    if( ijet1 > 0 && ijet2 > 0 ) {
		tree_.zepp          = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta())); 
		tree_.mj1j2         = dijet.M();
		tree_.dphi          = fabs(diphoton.DeltaPhi(dijet));
		tree_.dijetEta      = dijet.Eta();
		tree_.dphiJJ        = fabs(jet1->DeltaPhi(*jet2));
		tree_.dphiJJ2       = fabs(sumj1.DeltaPhi(sumj2));
	    }
	    
	    tree_.pho1r9         = lead_r9;
	    tree_.pho2r9         = sublead_r9;
	    tree_.pho1idMva     = phoid_mvaout_lead;
	    tree_.pho2idMva     = phoid_mvaout_sublead;
	    tree_.pho1sEoverE    = massResolutionCalculator->leadPhotonResolutionNoSmear();
	    tree_.pho2sEoverE    = massResolutionCalculator->subleadPhotonResolutionNoSmear();
	    tree_.pho1sEoverEsmear = massResolutionCalculator->leadPhotonResolution();   
	    tree_.pho2sEoverEsmear = massResolutionCalculator->subleadPhotonResolution();
	    tree_.diphosMoverM   = sigmaMrv;
	    tree_.diphosMoverMwv = sigmaMwv;
	    tree_.diphovtxProb   = vtxProb;
	    
	    if (ijet3 > 0 ) {
		tree_.deltaEta3 = fabs(jet3->Eta() - 0.5*(jet1->Eta() + jet2->Eta())); 
	    }
	    tree_.isSignal  = (cur_type < 0);
	    tree_.mctype    = cur_type;
	    tree_.entry++;

	    l.FillTreeContainer();

	    fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, diphoton_id,
			     category, isCorrectVertex, evweight, vtx, l, 0, -1, 0, -1 );
	}
	return false;
    }
    
    return false;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
