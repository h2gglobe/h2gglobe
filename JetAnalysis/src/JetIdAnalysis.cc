#include "../interface/JetIdAnalysis.h"

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

    minBiasRefName = "aux/minBiasRef.root";
}

// ----------------------------------------------------------------------------------------------------
JetIdAnalysis::~JetIdAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::Term(LoopAll& l) 
{
    /// l.outputFile->cd();
}

// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::Init(LoopAll& l) 
{
     StatAnalysis::Init(l);
}

// ----------------------------------------------------------------------------------------------------
void JetIdAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
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

// ----------------------------------------------------------------------------------------------------
bool JetIdAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
			      bool & isCorrectVertex,
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
	diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
	
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
	
	l.RescaleJetEnergy();

	diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,false, &smeared_pho_energy[0], true); 
	
	TLorentzVector lead_p4, sublead_p4, Higgs;
	float lead_r9, sublead_r9;
	TVector3 * vtx;
	if( diphotonVBF_id != -1 ) { diphoton_id = diphotonVBF_id; }
	if( diphoton_id == -1 ) { return false; }
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  
	
	/// Different jet ID options
	if( l.dipho_vtxind[diphoton_id] == 0 && (*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1. ) { 
	    std::vector<std::vector<unsigned char> > jetids(7,std::vector<unsigned char>(l.jet_algoPF1_n,true)); 
	    for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) {
		int ijetid = 0;
		TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
		
		bool pfloose = (l.jet_algoPF1_pfloose[ijet]  && fabs(p4->Eta()) < 3.) || (l.jet_algoPF1_nNeutrals[ijet] > 1);
		
		bool sloose  = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kLoose);
		bool floose  = PileupJetIdentifier::passJetId(l.jet_algoPF1_full_wp_level[ijet], PileupJetIdentifier::kLoose);
		bool cmedium = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kMedium);
		bool smedium = PileupJetIdentifier::passJetId(l.jet_algoPF1_simple_wp_level[ijet], PileupJetIdentifier::kMedium) ;
		bool fmedium = PileupJetIdentifier::passJetId(l.jet_algoPF1_full_wp_level[ijet], PileupJetIdentifier::kMedium)   ;
		
		jetids[++ijetid][ijet] = pfloose;
		jetids[++ijetid][ijet] = pfloose && sloose ; 
		jetids[++ijetid][ijet] = pfloose && floose ; 
		jetids[++ijetid][ijet] = pfloose && cmedium; 
		jetids[++ijetid][ijet] = pfloose && smedium;
		jetids[++ijetid][ijet] = pfloose && fmedium;
		
		/// Fill control plots
		if( p4->Pt() > 20. && fabs(p4->Eta()) < 4.7 && p4->DeltaR(lead_p4) > 0.5 && p4->DeltaR(sublead_p4) > 0.5 ) {
		    bool isMatched = l.jet_algoPF1_genMatched[ijet] && l.jet_algoPF1_genDr[ijet] < 0.3 && l.jet_algoPF1_genPt[ijet] > 10.;
		    for(; ijetid>=0; --ijetid) {
			bool passId = (bool)jetids[ijetid][ijet];
			if( ! passId ) { continue; }
			bool kinOnly = ijetid >= 1; // Fill the jet ID distributions only for no-selection
			fillJetIdPlots(l,ijet,ijetid,1.,"",kinOnly);
			if( isMatched ) {
			    fillJetIdPlots(l,ijet,ijetid,1.,"matched_",kinOnly);
			} else {
			    fillJetIdPlots(l,ijet,ijetid,1.,"unmatched_",kinOnly);
			}
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
		    VBFevent=VBFTag2012(l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight,(bool*)&(jetids[ijetid][0]));
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
	category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
	mass     = Higgs.M();

	// apply di-photon level smearings and corrections
        int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
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
	    fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, category, isCorrectVertex, evweight, l );
	}

	// see if the event falls into an exclusive category
	computeExclusiveCategory(l, category, diphoton_index, Higgs.Pt() );
  
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
    static float etabins[] = { 0., 2.5, 2.75, 3.0, 999. };
    static int netabins = sizeof(etabins) / sizeof(float);

    static float ptbins[] = { 20., 30., 50., 999999. };
    static int nptbins = sizeof(ptbins) / sizeof(float);
    
    TLorentzVector * p4 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
    l.FillHist2D(Form("%sjet_pt_eta",label),icat,p4->Pt(),p4->Eta(),wei);
    
    if( kinOnly ) { return; }
    assert( icat == 0 );
    int etacat = ( std::lower_bound( etabins, etabins + netabins, fabs(p4->Eta()) ) - etabins - 1 );
    int ptcat  = ( std::lower_bound( ptbins, ptbins + nptbins, fabs(p4->Pt()) )     - ptbins - 1 );
    /// std::cout << fabs(p4->Eta()) << " " << p4->Pt() << " " << etacat << " " << ptcat << std::endl;
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
    // Jet ID variables pt bins
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

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
