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

    isLHE    = false;
    isVBFNLO = false;
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
    doSystematics = false;
    StatAnalysis::Init(l);
}

// ----------------------------------------------------------------------------------------------------
bool VbfGenAnalysis::SkimEvents(LoopAll& l, int)
{
    if( isLHE ) {
	if( l.gh_higgs_p4 == 0 ) {
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
	}
	SetNullHiggs(l);
	
	std::vector<int> gen_photons;
	std::vector<int> gen_quarks;
	for(int ii=0; ii<l.gp_n; ++ii) {
	    if( isVBFNLO ) {
		if( l.gp_status[ii] == 1 && l.gp_pdgid[ii] != 0 ) {
		    float eta = fabs( ((TLorentzVector *)l.gp_p4->At(ii))->Eta() );
		    if( l.gp_pdgid[ii] == 22 ) {
			if( eta < 2.5 ) {
			    gen_photons.push_back(ii);
			}
		    } else if ( abs(l.gp_pdgid[ii]) < 7 ) {
			if( eta < 4.7 ) {
			    gen_quarks.push_back(ii);
			}
		    }
		}
	    }
	}
	if( gen_photons.size() > 1 && gen_quarks.size() > 1 ) {
	    std::sort(gen_photons.begin(),gen_photons.end(),
		      ClonesSorter<TLorentzVector,double,std::greater<double> >(l.gp_p4,&TLorentzVector::Pt));
	    std::sort(gen_quarks.begin(),gen_quarks.end(),
		      ClonesSorter<TLorentzVector,double,std::greater<double> >(l.gp_p4,&TLorentzVector::Pt));

	    TLorentzVector & pho1 = *((TLorentzVector *)l.gp_p4->At(gen_photons[0]));
	    TLorentzVector & pho2 = *((TLorentzVector *)l.gp_p4->At(gen_photons[1]));

	    TLorentzVector & q1 = *((TLorentzVector *)l.gp_p4->At(gen_quarks[0]));
	    TLorentzVector & q2 = *((TLorentzVector *)l.gp_p4->At(gen_quarks[1]));
	    
	    *((TLorentzVector *)l.gh_pho1_p4->At(0)) =  pho1;
	    *((TLorentzVector *)l.gh_pho2_p4->At(0)) =  pho2;
	    
	    *((TLorentzVector*)l.gh_higgs_p4 ) = pho1 + pho2;
	    
	    *((TLorentzVector *)l.gh_vbfq1_p4->At(0)) = q1;
	    *((TLorentzVector *)l.gh_vbfq2_p4->At(0)) = q2;
	    l.gh_vbfq1_pdgid=l.gp_pdgid[gen_quarks[0]];
	    l.gh_vbfq2_pdgid=l.gp_pdgid[gen_quarks[1]];
	} else {
	    return false;
	}
    }

    return true;
}


// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s )
{

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
}

// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA, 
				    int category, float weight, bool isCorrectVertex, int diphoton_id) 
{
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
	/// l.FillTree("deltaPhiJJ",myVBF_deltaPhiJJ);
	/// l.FillTree("deltaPhiGamGam", myVBF_deltaPhiGamGam);
	/// l.FillTree("etaJJ", myVBF_etaJJ);
	/// l.FillTree("thetaJ1", myVBF_thetaJ1);
	/// l.FillTree("thetaJ2", myVBF_thetaJ2);
	
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
    float lead_r9 = 0., sublead_r9 = 0.;
    TVector3 * vtx;
    fillGenDiphoton(lead_p4, sublead_p4, Higgs,l);
    l.FillTree("leadPt",(float)lead_p4.Pt());
    l.FillTree("subleadPt",(float)sublead_p4.Pt());
    l.FillTree("leadR9",lead_r9);
    l.FillTree("subleadR9",sublead_r9);
    //// l.FillTree("sigmaMrv",sigmaMrv);
    //// l.FillTree("sigmaMwv",sigmaMwv);
    //// 
    //// vtxAna_.setPairID(diphoton_id);
    //// float vtxProb = vtxAna_.vertexProbability(l.vtx_std_evt_mva->at(diphoton_id), l.vtx_std_n);
    //// float altMass = 0.;
    //// if( l.vtx_std_n > 1 ) {
    //// 	int altvtx = (*l.vtx_std_ranked_list)[diphoton_id][1];
    //// 	altMass = ( l.get_pho_p4( l.dipho_leadind[diphoton_id], altvtx, &smeared_pho_energy[0]) + 
    //// 		    l.get_pho_p4( l.dipho_subleadind[diphoton_id], altvtx, &smeared_pho_energy[0]) ).M();
    //// }
    //// l.FillTree("altMass",altMass);
    //// l.FillTree("vtxProb",vtxProb);
    
    
    fillControlPlots(lead_p4, sublead_p4, Higgs, category, weight, l);
    
}

// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::fillGenDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs, LoopAll &l)
{
    lead_p4 = *((TLorentzVector *)l.gh_pho1_p4->At(0));
    sublead_p4 = *((TLorentzVector *)l.gh_pho2_p4->At(0));
    
    Higgs = lead_p4 + sublead_p4;
}


// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs, 				      
				      int category, float evweight, LoopAll & l )
{   
    // control plots 
    if( category>=0 ) { 
	fillControlPlots( lead_p4, sublead_p4, Higgs, 0, evweight, l ); 
    }
    float mass = Higgs.M();
    l.FillHist("all_mass",category+1, Higgs.M(), evweight);
    if( mass>=massMin && mass<=massMax  ) {
	l.FillHist("mass",category+1, Higgs.M(), evweight);
	l.FillHist("eta",category+1, Higgs.Eta(), evweight);
	l.FillHist("pt",category+1, Higgs.Pt(), evweight);
	
	l.FillHist("pho_pt",category+1,lead_p4.Pt(), evweight);
	l.FillHist("pho1_pt",category+1,lead_p4.Pt(), evweight);
	l.FillHist("pho_eta",category+1,lead_p4.Eta(), evweight);
	l.FillHist("pho1_eta",category+1,lead_p4.Eta(), evweight);
	
	l.FillHist("pho_pt",category+1,sublead_p4.Pt(), evweight);
	l.FillHist("pho2_pt",category+1,sublead_p4.Pt(), evweight);
	l.FillHist("pho_eta",category+1,sublead_p4.Eta(), evweight);
	l.FillHist("pho2_eta",category+1,sublead_p4.Eta(), evweight);
	
	if( mvaVbfSelection ) {
	    l.FillHist("vbf_mva",category+1,myVBF_MVA,evweight);
	    if (VBFevent){
		float myweight =  1;
		float sampleweight = l.sampleContainer[l.current_sample_index].weight;
		if(evweight*sampleweight!=0) myweight=evweight/sampleweight;
		l.FillCutPlots(category+1,1,"_sequential",evweight,myweight); 
	    }
	}
    }
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
