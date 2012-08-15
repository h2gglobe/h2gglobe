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

    fillGhBranches = false;
    useGenJets     = false;
    analyzeJetVariablesOnly = false;
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
bool VbfGenAnalysis::SkimEvents(LoopAll& l, int jentry)
{
  
   if( fillGhBranches ) {
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
	
	l.b_gp_n->GetEntry(jentry);
	l.b_gp_p4->GetEntry(jentry);
	l.b_gp_status->GetEntry(jentry);
	l.b_gp_pdgid->GetEntry(jentry);

	std::vector<int> gen_photons;
	std::vector<int> gen_quarks;
	for(int ii=0; ii<l.gp_n; ++ii) {
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


	if( gen_photons.size() > 1 ) {
	    std::sort(gen_photons.begin(),gen_photons.end(),
		      ClonesSorter<TLorentzVector,double,std::greater<double> >(l.gp_p4,&TLorentzVector::Pt));
	    
	    TLorentzVector & pho1 = *((TLorentzVector *)l.gp_p4->At(gen_photons[0]));
	    TLorentzVector & pho2 = *((TLorentzVector *)l.gp_p4->At(gen_photons[1]));
	    
	    *((TLorentzVector *)l.gh_pho1_p4->At(0)) =  pho1;
	    *((TLorentzVector *)l.gh_pho2_p4->At(0)) =  pho2;
	    
	    *((TLorentzVector*)l.gh_higgs_p4->At(0)) = pho1 + pho2;
	    
	} else {
	    /// Require at least two photons
	    if(!analyzeJetVariablesOnly)
		return false;
	    for(int ii=0; ii<l.gp_n; ++ii) {
		if( l.gp_pdgid[ii] == 25 ) {
		    *((TLorentzVector*)l.gh_higgs_p4->At(0)) = *((TLorentzVector*)l.gp_p4->At(ii)) ;
		}
	    }
	}


	if( useGenJets ) {
	    /// PUT gen jet selection here
	    // clean and sort jets
	    std::vector<int> sorted_jets;
	    for(int ijet=0; ijet<l.genjet_algo1_n; ++ijet) { 
		TLorentzVector * p4 = (TLorentzVector*)l.genjet_algo1_p4->At(ijet);
		if( p4->DeltaR(  *((TLorentzVector *)l.gp_p4->At(gen_photons[0]))   ) > 0.5 && p4->DeltaR( *((TLorentzVector *)l.gp_p4->At(gen_photons[1]))) > 0.5  ) {
		    sorted_jets.push_back(ijet);
		}
	    }
	    std::sort(sorted_jets.begin(),sorted_jets.end(),
		      ClonesSorter<TLorentzVector,double,std::greater<double> >(l.genjet_algo1_p4,&TLorentzVector::Pt));
	    
	    if ( sorted_jets.size() > 1){
		TLorentzVector & q1 = *((TLorentzVector *)l.genjet_algo1_p4->At(sorted_jets[0]));
		TLorentzVector & q2 = *((TLorentzVector *)l.genjet_algo1_p4->At(sorted_jets[1]));
	    	*((TLorentzVector *)l.gh_vbfq1_p4->At(0)) = q1;
		*((TLorentzVector *)l.gh_vbfq2_p4->At(0)) = q2;
	    }
	    else {
		/// Also require two jets
		return false;
	    }
	} else {
	    if( gen_quarks.size() > 1 ) {
	
		std::sort(gen_quarks.begin(),gen_quarks.end(),
			  ClonesSorter<TLorentzVector,double,std::greater<double> >(l.gp_p4,&TLorentzVector::Pt));
		
		TLorentzVector & q1 = *((TLorentzVector *)l.gp_p4->At(gen_quarks[0]));
		TLorentzVector & q2 = *((TLorentzVector *)l.gp_p4->At(gen_quarks[1]));
		
		*((TLorentzVector *)l.gh_vbfq1_p4->At(0)) = q1;
		*((TLorentzVector *)l.gh_vbfq2_p4->At(0)) = q2;
		l.gh_vbfq1_pdgid=l.gp_pdgid[gen_quarks[0]];
		l.gh_vbfq2_pdgid=l.gp_pdgid[gen_quarks[1]];
	    } else {
		/// Also require two jets
		return false;
	    }
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
 
    if (!analyzeJetVariablesOnly){

	// do gen-level dependent first (e.g. k-factor); only for signal
	genLevWeight=1.;
	if(cur_type!=0 ) {
	    applyGenLevelSmearings(genLevWeight,gP4,l.pu_n,cur_type,genSys,syst_shift);
	}
	evweight = genLevWeight;
	
	
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
	category = 6;
	
	myVBFLeadJPt= jet1.Pt();
	myVBFSubJPt = jet2.Pt();
	myVBF_Mjj   = dijet.M();
	myVBFdEta   = fabs(jet1.Eta() - jet2.Eta());
	myVBFZep    = fabs(diphoton.Eta() - 0.5*(jet1.Eta() + jet2.Eta()));
	myVBFdPhi   = fabs(diphoton.DeltaPhi(dijet));
	myVBF_Mgg   = diphoton.M();
	mass = myVBF_Mgg;
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
	

	/*cout << "M: " << mass << "; lead/M " << myVBFLeadPhoPtOverM << "; sublead/M " << myVBFSubPhoPtOverM;
	  cout << endl << "jet1 Pt " << jet1.Pt() << "; jet2 pt " << jet2.Pt() << "; Mjj " << myVBF_Mjj << endl;
	  cout << "retval " << ((myVBFLeadPhoPtOverM > 0.5 && myVBFSubPhoPtOverM > 0.3) && (jet1.Pt() > 20 && jet2.Pt() > 20) && myVBF_Mjj > 50) << endl;*/
	
	if((fabs(lead_p4.Eta()) < 2.5 && fabs(sublead_p4.Eta()) < 2.5) && (myVBFLeadPhoPtOverM > 0.3 && myVBFSubPhoPtOverM > 25./120. && sublead_p4.Pt() > 25.) && (jet1.Pt() > 20. && jet2.Pt() > 20.) && myVBF_Mjj > 100.)
	    {
		
		if(myVBFLeadPhoPtOverM > 0.5 && myVBFSubPhoPtOverM > 0.3 && myVBFLeadJPt > 30. && myVBFSubJPt > 20. && myVBF_Mjj > 250. && TMath::Abs(myVBFdEta) > 3. && TMath::Abs(myVBFZep) < 2.5 && TMath::Abs(myVBFdPhi) > 2.6)
		    {
			
			if(myVBFSubJPt > 30. && myVBF_Mjj > 500.)
			    category = 4;
			else
			    category = 5;
		    }
		return true;
	    }
	else
	    return false;
    }
    else { // for lhe files without higgs decay to photons
	
	TLorentzVector jet1       = *((TLorentzVector*)l.gh_vbfq1_p4->At(0));
	TLorentzVector jet2       = *((TLorentzVector*)l.gh_vbfq2_p4->At(0));
	
	if(jet1.Pt() < jet2.Pt())
	    std::swap(jet1, jet2);
	
	TLorentzVector diphoton   = *((TLorentzVector*)l.gh_higgs_p4->At(0));
	TLorentzVector dijet      = jet1 + jet2;
	
	myVBF_MVA = -2.;
	VBFevent = true;
	category = 6;
	
	myVBFLeadJPt= jet1.Pt();
	myVBFSubJPt = jet2.Pt();
	myVBF_Mjj   = dijet.M();
	myVBFdEta   = fabs(jet1.Eta() - jet2.Eta());
	myVBFZep    = fabs(diphoton.Eta() - 0.5*(jet1.Eta() + jet2.Eta()));
	myVBFdPhi   = fabs(diphoton.DeltaPhi(dijet));
	myVBF_Mgg   = diphoton.M();
	mass = myVBF_Mgg;
	myVBFDiPhoPtOverM   = diphoton.Pt()   / myVBF_Mgg;
	myVBF_deltaPhiJJ = jet1.DeltaPhi(jet2);
	myVBF_etaJJ = (jet1.Eta() + jet2.Eta())/2;


	if( jet1.Pt() > 20. && jet2.Pt() > 20. && myVBF_Mjj > 100.)
	    {
		if( myVBFLeadJPt > 30. && myVBFSubJPt > 20. && myVBF_Mjj > 250. && TMath::Abs(myVBFdEta) > 3. && TMath::Abs(myVBFZep) < 2.5 && TMath::Abs(myVBFdPhi) > 2.6)
		    {
			if(myVBFSubJPt > 30. && myVBF_Mjj > 500.)
			    category = 4;
			else
			    category = 5;
		    }
		return true;
	    }
	else
	    return false;
    }

    return true;
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
    
    if( VBFevent ) {
	
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
    
    
    VbfGenAnalysis::fillControlPlots(lead_p4, sublead_p4, Higgs, category, weight, l);
    
    
}

// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::fillGenDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs, LoopAll &l)
{
    lead_p4 = *((TLorentzVector *)l.gh_pho1_p4->At(0));
    sublead_p4 = *((TLorentzVector *)l.gh_pho2_p4->At(0));
    
    if(lead_p4.Pt() < sublead_p4.Pt())
      std::swap(lead_p4, sublead_p4);

    Higgs = lead_p4 + sublead_p4;
}


// ----------------------------------------------------------------------------------------------------
void VbfGenAnalysis::fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs, 				      
				      int category, float evweight, LoopAll & l )
{   
    // control plots 
    
    if( category>=0 ) { 
        
	fillControlPlots( lead_p4, sublead_p4, Higgs, -1, evweight, l );
         
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
	
    
	
        if (VBFevent){
	   float myweight =  1;
	   float sampleweight = l.sampleContainer[l.current_sample_index].weight;
	   if(evweight*sampleweight!=0) myweight=evweight/sampleweight;
	   
           l.FillCutPlots(category+1,1,"_sequential",evweight,myweight); 
	   	
        }
    }
    
        
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
