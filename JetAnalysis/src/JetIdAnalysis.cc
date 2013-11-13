#include "../interface/JetIdAnalysis.h"

#include "../interface/JetHandler.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#include "CMGTools/External/interface/PileupJetIdentifier.h"

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
JetIdAnalysis::JetIdAnalysis()  
{
    name_ = "JetIdAnalysis";

    recomputeJetId = false;
    expoMatching   = false;
    dumpFlatTree   = false;
    runZmumuValidation=false;
    saveAllJets=false;
    flatTree_ = 0;
    outputFile_ = 0;

   
}

// ----------------------------------------------------------------------------------------------------
JetIdAnalysis::~JetIdAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::Term(LoopAll& l) 
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
void JetIdAnalysis::Init(LoopAll& l) 
{
    if( l.runZeeValidation ) {
	leadEtCut = 15;
	subleadEtCut = 15;
	leadEtVBFCut = 15.;
	subleadEtVBFCut = 15.;
	massMin = 60, massMax = 120.;
	applyPtoverM = false;
    }
    StatAnalysis::Init(l);
    if( jetHandler_ == 0 ) {
    	jetHandler_ = new JetHandler(jetHandlerCfg, l);
    }

    if( dumpFlatTree ) {
	if( flatTree_ == 0 ) { 
	    outputFile_ = TFile::Open("jetid_"+l.histFileName,"recreate");
	    flatTree_ = new TTree("flatTree","flatTree");
	}
    	jetHandler_->bookFlatTree( flatTree_ );
	flatTree_->Branch( "ievent", &tree_ievent );
	flatTree_->Branch( "ijet", &tree_ijet );
	flatTree_->Branch( "isMatched", &tree_isMatched );
	flatTree_->Branch( "jetGenPt", &tree_genPt );
	flatTree_->Branch( "jetHenDr", &tree_genDr );
	flatTree_->Branch( "puweight", &tree_puweight );
	flatTree_->Branch( "npu", &tree_npu );
	flatTree_->Branch( "isData", &tree_isData );
	
	flatTree_->Branch( "njets", &tree_njets );
	flatTree_->Branch( "jetLooseID", &tree_jetLooseID );
	flatTree_->Branch( "simpleId",&tree_simpleId );
	flatTree_->Branch( "fullId",&tree_fullId );
	flatTree_->Branch( "cutbasedId",&tree_cutbasedId );
	flatTree_->Branch( "simpleDiscriminant",&tree_simpleDiscriminant );
	flatTree_->Branch( "fullDiscriminant",&tree_fullDiscriminant );
	flatTree_->Branch( "cutbasedDiscriminant",&tree_cutbasedDiscriminant );
	
	flatTree_->Branch( "dimuonMass",&tree_dimuonMass);
	flatTree_->Branch( "dimuonPt",&tree_dimuonPt);
	flatTree_->Branch( "dphiZJet",&tree_dphiZJet);

	
    }
}

// ----------------------------------------------------------------------------------------------------
bool JetIdAnalysis::SkimEvents(LoopAll& l, int jentry){
    
    if (runZmumuValidation) {
	return true;
    }
    else{
	SkimEvents(l,jentry);
    }
}
  
// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) {
  
    if (!runZmumuValidation) GetBranches(t,s);
   
}


// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    dumpFlatTree=true;
    flatTree_ = new TTree("flatTree","flatTree");
}

// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
}

// ----------------------------------------------------------------------------------------------------
bool JetIdAnalysis::SelectEventsReduction(LoopAll&, int)
{
    return true;
}

void switchJetIdVertex(LoopAll &l, int ivtx) ;

// ----------------------------------------------------------------------------------------------------
bool JetIdAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, 
				 int & category, int & diphoton_id,
				 bool & isCorrectVertex,
				 float &kinematic_bdtout,
				 bool isSyst, 
				 float syst_shift, bool skipSelection,
				 BaseGenLevelSmearer *genSys, BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys)
{
    assert( isSyst || ! skipSelection );

    int cur_type = l.itype[l.current];
    float sampleweight = l.sampleContainer[l.current_sample_index].weight();
    /// diphoton_id = -1;
    
    if (!runZmumuValidation){

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
	diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
	//// diphoton_id = l.DiphotonCiCSelection(l.phoNOCUTS, l.phoNOCUTS, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
	
	// N-1 plots
	if( ! isSyst ) {
	    int diphoton_nm1_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoNOCUTS, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
	    if(diphoton_nm1_id>-1) {
		float eventweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
		float myweight=1.;
		if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
		ClassicCatsNm1Plots(l, diphoton_nm1_id, &smeared_pho_energy[0], eventweight, myweight);
	    }
	}
	
	// Exclusive Modes
	int diphotonVBF_id = -1;
	int diphotonVHhad_id = -1;
	int diphotonVHlep_id = -1;
	VHmuevent = false;
	VHelevent = false;
	VBFevent = false;
	VHhadevent = false;
	
	diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,false, &smeared_pho_energy[0], true); 
	/// diphotonVBF_id = l.DiphotonCiCSelection(l.phoNOCUTS, l.phoNOCUTS, leadEtVBFCut, subleadEtVBFCut, 4,false, &smeared_pho_energy[0], true); 
	
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
	    if( p4->DeltaR(lead_p4) > 0.5 && p4->DeltaR(sublead_p4) > 0.5 ) {
		sorted_jets.push_back(ijet);
	    }
	}
	std::sort(sorted_jets.begin(),sorted_jets.end(),
		  ClonesSorter<TLorentzVector,double,std::greater<double> >(l.jet_algoPF1_p4,&TLorentzVector::Pt));
	
	///// if( recomputeJetId ) {
	/////     for(size_t itjet=0; itjet<sorted_jets.size(); ++itjet ) {
	///// 	jetHandler_->computeMva(sorted_jets[itjet],l.dipho_vtxind[diphoton_id]);
	/////     }
	///// }
	switchJetIdVertex( l, l.dipho_vtxind[diphoton_id] );

	/// Different jet ID options
	std::vector<std::vector<unsigned char> > jetids(12,std::vector<unsigned char>(l.jet_algoPF1_n,false)); 
	for(size_t itjet=0; itjet<sorted_jets.size(); ++itjet ) {
	    int & ijet = sorted_jets[itjet];
	    int ijetid = 0;
	    TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
	    
	    bool pfloose = (l.jet_algoPF1_pfloose[ijet] );
	    
	    bool sloose  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kLoose);
	    bool floose  = PileupJetIdentifier::passJetId(l.jet_algoPF1_full_wp_level[ijet], PileupJetIdentifier::kLoose);
	    bool cmedium = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kMedium);
	    bool smedium = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kMedium) ;
	    bool fmedium = PileupJetIdentifier::passJetId(l.jet_algoPF1_full_wp_level[ijet], PileupJetIdentifier::kMedium)   ;
	    
	    jetids[ijetid][ijet]   = true; // No selection
	    jetids[++ijetid][ijet] = pfloose;
	    jetids[++ijetid][ijet] = sloose ; 
	    jetids[++ijetid][ijet] = floose ; 
	    jetids[++ijetid][ijet] = cmedium; 
	    jetids[++ijetid][ijet] = smedium;
	    jetids[++ijetid][ijet] = fmedium;
	    jetids[++ijetid][ijet] = pfloose && sloose ; 
	    jetids[++ijetid][ijet] = pfloose && floose ; 
	    jetids[++ijetid][ijet] = pfloose && cmedium; 
	    jetids[++ijetid][ijet] = pfloose && smedium;
	    jetids[++ijetid][ijet] = pfloose && fmedium;
			    
	    /// Fill control plots
	    if( p4->Pt() > 20. && fabs(p4->Eta()) < 4.5 ) {
		bool isMatched, notMatched;
		if( expoMatching ) {
		    isMatched = 
			l.jet_algoPF1_genDr[ijet] < 0.1 + 0.3*exp(-0.05*(l.jet_algoPF1_genPt[ijet]-10.)) 
			& fabs(l.jet_algoPF1_genPt[ijet] - p4->Pt())/l.jet_algoPF1_genPt[ijet] < 0.5;
		    notMatched = ! isMatched || ( l.jet_algoPF1_genDr[ijet] > 0.3 && l.jet_algoPF1_genPt[ijet] < 10. );
		} else {
		    isMatched = l.jet_algoPF1_genMatched[ijet] && l.jet_algoPF1_genDr[ijet] < 0.3 && l.jet_algoPF1_genPt[ijet] > 10.;
		    notMatched = ! l.jet_algoPF1_genMatched[ijet] 
			|| ( l.jet_algoPF1_genDr[ijet] > 0.3 && l.jet_algoPF1_genPt[ijet] < 10. );
		}
		if( l.event>65200 && l.event<65300 ) {
		    cout << "event: " << l.event << " jet: " << itjet
			 << " pt: " << p4->Pt() << " pt_raw: " << p4->Pt() / l.jet_algoPF1_erescale[ijet] 
			 << " eta: " << p4->Eta()
			 << " isMatched: " << isMatched << " matchDr: " << l.jet_algoPF1_genDr[ijet] 
			 << " matchPt: " << l.jet_algoPF1_genPt[ijet] << std::endl  ;
		}
		for(; ijetid>=0; --ijetid) {
		    bool passId = (bool)jetids[ijetid][ijet];
		    if( ! passId ) { continue; }
		    bool kinOnly = ijetid >= 1; // Fill the jet ID distributions only for no-selection
		    fillJetIdPlots(l,ijet,ijetid,1.,"",kinOnly);
		    if( isMatched && itjet == 1 ) {
			fillJetIdPlots(l,ijet,ijetid,1.,"matched_",kinOnly);
		    } else if( notMatched ) {
			fillJetIdPlots(l,ijet,ijetid,1.,"unmatched_",kinOnly);
		    }
		}
		if( dumpFlatTree ) {
		    tree_ievent = l.event;
		    tree_ijet = itjet;
		    tree_isMatched = l.jet_algoPF1_genMatched[ijet];
		    tree_genPt = l.jet_algoPF1_genPt[ijet];
		    tree_genDr = l.jet_algoPF1_genDr[ijet];
		    tree_njets = l.jet_algoPF1_n;
		    tree_jetLooseID = pfloose;
		    if( isMatched && itjet == 1 && cur_type < 0 ) {
			jetHandler_->fillFromJet(ijet,l.dipho_vtxind[diphoton_id]);
		    } else if ( notMatched ) {
			jetHandler_->fillFromJet(ijet,l.dipho_vtxind[diphoton_id]);
		    }
		    flatTree_->Fill();
		}
	    }
	}
	if( diphotonVBF_id != -1 ) {
	    float eventweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
	    float myweight=1.;
	    if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
	    
	    l.FillHist2D("vtxid_vs_nvtx",0,l.vtx_std_n,l.dipho_vtxind[diphotonVBF_id],eventweight);
	    l.FillHist2D("vtxid_vs_rho",0,l.rho_algo1,l.dipho_vtxind[diphotonVBF_id],eventweight);
	    
	    /// Test effect of different jet IDs on the event selection
	    for(int ijetid=jetids.size()-1; ijetid>=0;--ijetid) {
		int ijet1, ijet2;
		VBFevent=VBFTag2012(ijet1, ijet2, 
				    l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight,(bool*)&(jetids[ijetid][0]));
		if( VBFevent ) {
		    diphoton_id = diphotonVBF_id;
		    l.FillHist2D("dijet_count_vs_nvtx",ijetid,l.vtx_std_n,1.,eventweight);
		    l.FillHist2D("dijet_count_vs_rho",ijetid,l.rho_algo1,1.,eventweight);
		    l.FillHist2D("vtxid_vs_nvtx",ijetid+1,l.vtx_std_n,l.dipho_vtxind[diphotonVBF_id],eventweight);
		    l.FillHist2D("vtxid_vs_rho",ijetid+1,l.rho_algo1,l.dipho_vtxind[diphotonVBF_id],eventweight);
		} else {
		    l.FillHist2D("dijet_count_vs_nvtx",ijetid,l.vtx_std_n,0.,eventweight);
		    l.FillHist2D("dijet_count_vs_rho",ijetid,l.rho_algo1,0.,eventweight);
		}
	    }
	}
    }

    // if we selected any di-photon, compute the Higgs candidate kinematics
    // and compute the event category 
    if (diphoton_id > -1 ) {
        diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );

        // bring all the weights together: lumi & Xsection, single gammas, pt kfactor
	evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
	if( ! isSyst ) {
	    l.countersred[diPhoCounter_]++;
	}
	
        TLorentzVector lead_p4, sublead_p4, Higgs;
        float lead_r9, sublead_r9;
        TVector3 * vtx;
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  
      
        // FIXME pass smeared R9
	category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),Higgs.Pt()/Higgs.M(),nEtaCategories,nR9Categories,nPtCategories,nPtOverMCategories);
	mass     = Higgs.M();

	// apply di-photon level smearings and corrections
        int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),Higgs.Pt()/Higgs.M(),nEtaCategories,nR9Categories,0,0);
        if( cur_type != 0 && doMCSmearing ) {
	    applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, zero_, zero_,
				   diPhoSys, syst_shift);
            isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
        }
        float ptHiggs = Higgs.Pt();
      
	// sanity check
        assert( evweight >= 0. ); 
	
	// fill control plots and counters
	if( ! isSyst ) {
	    l.FillCounter( "Accepted", weight );
	    l.FillCounter( "Smeared", evweight );
	    sumaccept += weight;
	    sumsmear += evweight;
	    fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, diphoton_id, category,
			     isCorrectVertex, evweight, vtx, l);
	}

	// see if the event falls into an exclusive category
	computeExclusiveCategory(l, category, diphoton_index, Higgs.Pt(), Higgs.M() );
  
	return true;
    }
    }

    //**** runZmumuValidation
    else {

	// check if there is a good Zmumu candidate
	bool selectEvent=false;
	int goodMuon1=-1, goodMuon2=-1;
	DiMuonSelection(l, goodMuon1, goodMuon2, selectEvent);
	if( ! selectEvent ) { return false; }
	
	TLorentzVector *lead_p4    = (TLorentzVector*) l.mu_glo_p4->At(goodMuon1);
	TLorentzVector *sublead_p4 = (TLorentzVector*) l.mu_glo_p4->At(goodMuon2);

	// clean and sort jets
	std::vector<int> sorted_jets;
	for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) { 
	    TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
	    if( p4->DeltaR(*lead_p4) > 0.5 && p4->DeltaR(*sublead_p4) > 0.5 ) {
		sorted_jets.push_back(ijet);
	    }
	}
	std::sort(sorted_jets.begin(),sorted_jets.end(),
		  ClonesSorter<TLorentzVector,double,std::greater<double> >(l.jet_algoPF1_p4,&TLorentzVector::Pt));


	// jet mc matching  
	if (cur_type !=0) 
	    l.doJetMatching(*l.jet_algoPF1_p4,*l.genjet_algo1_p4,l.jet_algoPF1_genMatched,l.jet_algoPF1_vbfMatched,l.jet_algoPF1_bgenMatched,l.jet_algoPF1_cgenMatched,l.jet_algoPF1_lgenMatched,l.jet_algoPF1_genPt,l.jet_algoPF1_genDr);


	// post process jets (recompute mvas and wp) - use vertex 0
	postProcessJets(l, 0) ;
	switchJetIdVertex( l, 0);		
	
	// loop over sorted and cleaned jets
	for(size_t itjet=0; itjet<sorted_jets.size(); ++itjet ) {
	    int & ijet = sorted_jets[itjet];
	    int ijetid = 0;
	    TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
	    
	    if( dumpFlatTree ) {
		tree_isData=1;
		tree_ievent = l.event;
		tree_ijet = itjet;
		tree_puweight = 1;
		if (cur_type != 0){
		    tree_isMatched = l.jet_algoPF1_genMatched[ijet];
		    tree_genPt = l.jet_algoPF1_genPt[ijet];
		    tree_genDr = l.jet_algoPF1_genDr[ijet];
		    tree_npu   = l.pu_n;
		    tree_puweight = getPuWeight( l.pu_n, cur_type, &(l.sampleContainer[l.current_sample_index]),0);
		    tree_isData=0;
		}
		tree_njets = l.jet_algoPF1_n;
		tree_jetLooseID = l.jet_algoPF1_pfloose[ijet];
		
		jetHandler_->fillFromJet(ijet,0);
		tree_simpleId   =  l.jet_algoPF1_simple_wp_level[ijet];
		tree_fullId     =  l.jet_algoPF1_full_wp_level[ijet];
		tree_cutbasedId =  l.jet_algoPF1_cutbased_wp_level[ijet];
		tree_simpleDiscriminant   =  l.jet_algoPF1_simple_mva[ijet];
		tree_fullDiscriminant     =  l.jet_algoPF1_full_mva[ijet];
		tree_cutbasedDiscriminant =  -999;
		
		tree_dimuonPt   = (*lead_p4 + *sublead_p4).Pt() ;
		tree_dimuonMass = (*lead_p4 + *sublead_p4).M() ;
		tree_dphiZJet   = (*lead_p4 + *sublead_p4).DeltaPhi(*p4);

		// save infos for leading jet
		if (itjet==0)
		    flatTree_->Fill();
		//if flag set to true, save also other jets
		if (itjet>0 && saveAllJets) 
		    flatTree_->Fill();
	    }
	}
	
	return true;
	
    }


    return false;
}


// lazy plot filling macro
#define FILL_JETID_PLOT(LABEL,VARIABLE,ICAT,IJET,WEI) \
    l.FillHist( Form("%sjet_" # VARIABLE, LABEL ), ICAT, l.jet_algoPF1_ ## VARIABLE [IJET], WEI )

// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::fillJetIdPlots(LoopAll & l, int ijet, int icat, float wei, const char * label, 
				   bool kinOnly)
{
    static float etabins[] = { 0., 2.5, 2.75, 4.0, 999. };
    static int netabins = sizeof(etabins) / sizeof(float);

    static float ptbins[] = { 20., 22., 25., 30., 35., 40., 50., 999999. };
    static int nptbins = sizeof(ptbins) / sizeof(float);
    
    TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
    l.FillHist2D(Form("%sjet_pt_eta",label),icat,p4->Pt(),p4->Eta(),wei);
    
    if( kinOnly ) { return; }
    assert( icat == 0 );
    int etacat = ( std::lower_bound( etabins, etabins + netabins, fabs(p4->Eta()) ) - etabins - 1 );
    int ptcat  = ( std::lower_bound( ptbins, ptbins + nptbins, fabs(p4->Pt()) )     - ptbins - 1 );

    // Jet ID variables pt inclusive
    FILL_JETID_PLOT(label,frac01         , etacat, ijet, wei ); 
    FILL_JETID_PLOT(label,frac02         , etacat, ijet, wei );
    FILL_JETID_PLOT(label,frac03         , etacat, ijet, wei );
    FILL_JETID_PLOT(label,frac04         , etacat, ijet, wei );
    FILL_JETID_PLOT(label,frac05         , etacat, ijet, wei );
    FILL_JETID_PLOT(label,beta 	         , etacat, ijet, wei );
    FILL_JETID_PLOT(label,betaStar       , etacat, ijet, wei );
    FILL_JETID_PLOT(label,betaStarClassic, etacat, ijet, wei );
    FILL_JETID_PLOT(label,dR2Mean        , etacat, ijet, wei );
    FILL_JETID_PLOT(label,nNeutrals      , etacat, ijet, wei );
    FILL_JETID_PLOT(label,nCharged       , etacat, ijet, wei );
    FILL_JETID_PLOT(label,simple_mva     , etacat, ijet, wei );
    FILL_JETID_PLOT(label,full_mva       , etacat, ijet, wei );

    // Jet ID variables in pt bins
    FILL_JETID_PLOT(label,frac01         , etacat+ (ptcat+1)*(netabins-1), ijet, wei ); 
    FILL_JETID_PLOT(label,frac02         , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,frac03         , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,frac04         , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,frac05         , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,beta 	         , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,betaStar       , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,betaStarClassic, etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,dR2Mean        , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,nNeutrals      , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,nCharged       , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,simple_mva     , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
    FILL_JETID_PLOT(label,full_mva       , etacat+ (ptcat+1)*(netabins-1), ijet, wei );
}

void switchJetIdVertex(LoopAll &l, int ivtx) 
{
    for(int ii=0; ii<l.jet_algoPF1_n; ++ii) {
	l.jet_algoPF1_beta[ii]              = (*l.jet_algoPF1_beta_ext)[ii][ivtx];
    	l.jet_algoPF1_betaStar[ii]          = (*l.jet_algoPF1_betaStar_ext)[ii][ivtx];
    	l.jet_algoPF1_betaStarClassic[ii]   = (*l.jet_algoPF1_betaStarClassic_ext)[ii][ivtx];
    	l.jet_algoPF1_simple_mva[ii]        = (*l.jet_algoPF1_simple_mva_ext)[ii][ivtx];
    	l.jet_algoPF1_full_mva[ii]          = (*l.jet_algoPF1_full_mva_ext)[ii][ivtx];
    	l.jet_algoPF1_simple_wp_level[ii]   = (*l.jet_algoPF1_simple_wp_level_ext)[ii][ivtx];
    	l.jet_algoPF1_full_wp_level[ii]     = (*l.jet_algoPF1_full_wp_level_ext)[ii][ivtx];
    	l.jet_algoPF1_cutbased_wp_level[ii] = (*l.jet_algoPF1_cutbased_wp_level_ext)[ii][ivtx];
    }
}




void JetIdAnalysis::DiMuonSelection(LoopAll & l, int& goodMuon1, int& goodMuon2, bool& isZcandidate)
{

    float ZMASS = 91.188;
    float minDeltaM = 9999;
    
    goodMuon1 = -1 ;
    goodMuon2 = -1 ;
    
    std::vector<int> goodMuonIndex;
    int nGoodMuons  = 0;
    
    //---preselect muons
    for ( unsigned int imu=0; imu< l.mu_glo_n; ++imu ) {
	// muon pt
	float muonpt = ((TLorentzVector*) l.mu_glo_p4->At(imu)) -> Pt();
	// require tight muon : tight ID and isolation
	if ( !(l.MuonTightID2012(imu, 0)) || !(l.MuonIsolation2012(imu,  muonpt))) continue;
	goodMuonIndex.push_back(imu);
	nGoodMuons++;
    }
    
    //  std::cout <<  goodMuonIndex.size() << "  " << nGoodMuons << std::endl;
    
    if (nGoodMuons > 1){
	for ( unsigned int i = 0; i < goodMuonIndex.size(); ++i ) {
	    for ( unsigned int j = i+1; j < goodMuonIndex.size(); ++j ) {
		int ii =  goodMuonIndex.at(i);
		int jj =  goodMuonIndex.at(j);
		
		TLorentzVector *mu1  = (TLorentzVector*) l.mu_glo_p4->At(ii);
		TLorentzVector *mu2  = (TLorentzVector*) l.mu_glo_p4->At(jj);
		float invmass = (*mu1+*mu2).M();
		if ( fabs(invmass-ZMASS) < 30 && fabs(invmass-ZMASS) < minDeltaM ) {
		    goodMuon1 = ii;
		    goodMuon2 = jj;
		    minDeltaM = fabs(invmass-ZMASS);
		}
	    }
	}
    }
    
    if (goodMuon1!=-1 && goodMuon2!=-1) isZcandidate = true;
    
}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
