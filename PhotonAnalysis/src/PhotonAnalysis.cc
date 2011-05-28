#include "../interface/PhotonAnalysis.h"


#include "Sorters.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0

using namespace std;



// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::PhotonAnalysis()  : 
	runStatAnalysis(false), doTriggerSelection(false),
	name_("PhotonAnalysis"),
	vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{
}

// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::~PhotonAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Init(LoopAll& l) 
{
	if(PADEBUG) 
		cout << "InitRealPhotonAnalysis START"<<endl;

	if( vtxVarNames.empty() ) {
		vtxVarNames.push_back("ptbal"), vtxVarNames.push_back("ptasym"), vtxVarNames.push_back("logsumpt2");
	}
	
	// trigger
	triggerSelections.push_back(TriggerSelection(-1,147116));
	triggerSelections.back().addpath("HLT_DoublePhoton15_L1R");
	triggerSelections.push_back(TriggerSelection(146428,148058));
	triggerSelections.back().addpath("HLT_DoublePhoton17_L1R");
	triggerSelections.push_back(TriggerSelection(148822,149294));
	triggerSelections.back().addpath("HLT_DoublePhoton22_L1R_v1");
	triggerSelections.push_back(TriggerSelection(160000,-1));
	triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");

	trigCounter_ = l.countersred.size();
	l.countersred.resize(trigCounter_+1);

	// CiC initialization
	// FIXME should move this to GeneralFunctions
	const int phoNCUTS = LoopAll::phoNCUTS;
	const int phoCiC6NCATEGORIES = LoopAll::phoCiC6NCATEGORIES;
	const int phoCiC4NCATEGORIES = LoopAll::phoCiC4NCATEGORIES;
	const int phoNCUTLEVELS = LoopAll::phoNCUTLEVELS;

	for(int iLevel=0; iLevel<phoNCUTLEVELS; ++iLevel) {
		float cic6_cuts_lead[phoNCUTS][phoCiC6NCATEGORIES];
		float cic6_cuts_sublead[phoNCUTS][phoCiC6NCATEGORIES];
		float cic4_cuts_lead[phoNCUTS][phoCiC4NCATEGORIES];
		float cic4_cuts_sublead[phoNCUTS][phoCiC4NCATEGORIES];
		l.SetPhotonCutsInCategories((LoopAll::phoCiCIDLevel)iLevel, &cic6_cuts_lead[0][0], &cic6_cuts_sublead[0][0], &cic4_cuts_lead[0][0], &cic4_cuts_sublead[0][0] );
		
		float * cic6_cuts_arrays_lead[phoNCUTS] = {
			&l.cic6_cut_lead_isosumoet[0][0], &l.cic6_cut_lead_isosumoetbad[0][0], &l.cic6_cut_lead_trkisooet[0][0], &l.cic6_cut_lead_sieie[0][0],
			&l.cic6_cut_lead_hovere[0][0], &l.cic6_cut_lead_r9[0][0], &l.cic6_cut_lead_drtotk_25_99[0][0], &l.cic6_cut_lead_pixel[0][0] 
		};
		
		float * cic6_cuts_arrays_sublead[phoNCUTS] = {
			&l.cic6_cut_sublead_isosumoet[0][0], &l.cic6_cut_sublead_isosumoetbad[0][0], &l.cic6_cut_sublead_trkisooet[0][0], 
			&l.cic6_cut_sublead_sieie[0][0], &l.cic6_cut_sublead_hovere[0][0], &l.cic6_cut_sublead_r9[0][0],
			&l.cic6_cut_sublead_drtotk_25_99[0][0], &l.cic6_cut_sublead_pixel[0][0]
		};

		float * cic4_cuts_arrays_lead[phoNCUTS] = {
			&l.cic4_cut_lead_isosumoet[0][0], &l.cic4_cut_lead_isosumoetbad[0][0], &l.cic4_cut_lead_trkisooet[0][0], &l.cic4_cut_lead_sieie[0][0],
			&l.cic4_cut_lead_hovere[0][0], &l.cic4_cut_lead_r9[0][0], &l.cic4_cut_lead_drtotk_25_99[0][0], &l.cic4_cut_lead_pixel[0][0] 
		};
		
		float * cic4_cuts_arrays_sublead[phoNCUTS] = {
			&l.cic4_cut_sublead_isosumoet[0][0], &l.cic4_cut_sublead_isosumoetbad[0][0], &l.cic4_cut_sublead_trkisooet[0][0], 
			&l.cic4_cut_sublead_sieie[0][0], &l.cic4_cut_sublead_hovere[0][0], &l.cic4_cut_sublead_r9[0][0],
			&l.cic4_cut_sublead_drtotk_25_99[0][0], &l.cic4_cut_sublead_pixel[0][0]
		};

		for(int iCut=0; iCut<phoNCUTS; ++iCut) {
			for(int iCat=0; iCat<phoCiC6NCATEGORIES; ++iCat) {
				cic6_cuts_arrays_lead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_lead[iCut][iCat];
				cic6_cuts_arrays_sublead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_sublead[iCut][iCat];
			}
			for(int iCat=0; iCat<phoCiC4NCATEGORIES; ++iCat) {
				cic4_cuts_arrays_lead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_lead[iCut][iCat];
				cic4_cuts_arrays_sublead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_sublead[iCut][iCat];
			}
		}
	}

	if (l.typerun == 2 || l.typerun == 1) {
	}
  
	if (runStatAnalysis) {  
		/*
		  rooContainer->AddRealVar("mass",50.,250.);
		  rooContainer->AddRealVar("mu",-0.04,-1.,-0.001);
		  
		  // -------------------------------------//
		  std::vector<const char*> pars(2,"t");	 
		  pars[0] = "mass";
		  pars[1] = "mu";
		  // -------------------------------------//
		  rooContainer->AddGenericPdf("exp",
		  "exp((@0)*(@1))",pars,10);
		  
		  rooContainer->CreateDataSet("mass");
		*/
	}

    /* -------------------------------------------------------------------------------------------
    Pileup Reweighting
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
    ----------------------------------------------------------------------------------------------  */
    if (puHist != "") {
        if(PADEBUG) 
            cout << "Opening PU file"<<endl;
        TFile* puFile = TFile::Open( puHist );
        if (puFile) {
            cout<<"Reweighting events for pileup."<<endl;
            TH1D *histo; 
            puFile->GetObject("pileup",histo);
            weights = l.generate_flat10_weights(histo);
	    puFile->Close();
        }
        else {
            cout<<"Error opening " <<puHist<<" pileup reweighting histogram, using 1.0"<<endl; 
            weights.resize(50);
            for (int i=0; i<weights.size(); i++) weights[i] = 1.0;
        }
        if(PADEBUG) 
            cout << "Opening PU file END"<<endl;
    } 

	if(PADEBUG) 
		cout << "InitRealPhotonAnalysis END"<<endl;

	// FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
	if(PADEBUG) 
		cout << "Analysis START"<<endl;
	
    // From Here is the Standard Dec Review Selection/ gen Level studies
	std::vector<int> pho_loose;
	// For use of the preselected photons in Dec Review, empty the vector
	pho_presel.clear();

    //PU reweighting
    int n_pu = l.pu_n;
    float weight = l.sampleContainer[l.current_sample_index].weight;
    if (l.itype[l.current] !=0 && puHist != "") {
        if(n_pu<weights.size()){
             weight *= weights[n_pu];
        }    
        else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
             cout<<"n_pu ("<< n_pu<<") too big ("<<weights.size()<<"), event will not be reweighted for pileup"<<endl;
        }
    }

	if( pho_presel.empty() ) { 
		PreselectPhotons(l,jentry);
	}
	TVector3 *calopos;	
	TLorentzVector *p4;
	for (size_t i=0; i<pho_presel.size(); ++i) {
		int ipho = pho_presel[i];

		p4 = (TLorentzVector *) l.pho_p4->At(ipho);
		calopos  = (TVector3 *) l.pho_calopos->At(ipho);
		
		// FIXME vertex
		float pt  = p4->Pt(); 
		float eta = fabs(calopos->Eta());
		//PreSelection
		if ( (! l.pho_haspixseed[ipho])
		     && pt > 30. 
		     && l.pho_hoe[ipho] <  0.1
		     && l.pho_trksumptsolidconedr03[ipho] < 2*(3.5 + 0.001*pt)
		     && l.pho_ecalsumetconedr03[ipho] < 2*(4.2 + 0.006*pt)
		     && l.pho_hcalsumetconedr03[ipho] < 2*(2.2 + 0.0025*pt)
		     &&((eta < 1.4442) || ((eta > 1.566) && (eta < 2.5))) ) {
			pho_loose.push_back(ipho);
		}
	}
	
	//Event Selection
	int n_pho_loose = pho_loose.size();
	//FillHist("h_n_sel",n_preselected_pho);
	
	// Sort Photons into Pt Order
	std::sort(pho_presel.begin(),pho_presel.end(),
		  ClonesSorter<TLorentzVector,double,std::greater<double> >(l.pho_p4,&TLorentzVector::Pt));
	
        // Regular Event Selection begins here
	float best_mass = 0.;
	float best_pt = -1;
	int category = -1;
	float min_r9;
	float max_eta;
	
	if (n_pho_loose > 1 ){
		
		int leading, nleading;
				
		leading  = pho_loose[0];
		nleading = pho_loose[1];
		
		TLorentzVector * leading_p4 = (TLorentzVector *) l.pho_p4->At(leading);
		TLorentzVector * nleading_p4 = (TLorentzVector *) l.pho_p4->At(nleading);
		
		
		if (leading_p4->Pt() > 40.){
			
			// Determine the Category of the event
			// -> Which histogram is filled
			min_r9  = min(l.pho_r9[leading], l.pho_r9[nleading]);
			
			TVector3 * leading_calopos  = (TVector3 *) l.pho_calopos->At(leading);
			TVector3 * nleading_calopos  = (TVector3 *) l.pho_calopos->At(nleading);
			max_eta = max(fabs(leading_calopos->Eta())
			              ,fabs(nleading_calopos->Eta()));
			
			if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 1;
			if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 2;
			if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;
			if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 4;
			
			// -------------------------------------------------------
			TLorentzVector Higgs = *leading_p4 + *nleading_p4;
			
			float mass = Higgs.M();
			float h_pt = Higgs.Pt();
			if (mass > 100. && mass < 150.){
				//Good event, passes preselection and acceptance cuts
				
				int pass_selection[2];
				int pass_isolation[2];
				int in_iso_gap[2];

				//Now do selection on leading photon
				pass_selection[0] = l.pho_hoe[leading] < 0.02
					&& (((l.pho_sieie[leading] < 0.01)  && (fabs(leading_calopos->Eta()) < 1.4442)) 
					    || (( l.pho_sieie[leading] < 0.028)
						&& ((fabs(leading_calopos->Eta()) < 2.5) && (fabs(leading_calopos->Eta()) > 1.566))) );
				pass_isolation[0] =  l.pho_trksumptsolidconedr03[leading] < (1.5 + 0.001*leading_p4->Pt())		
					&& l.pho_ecalsumetconedr03[leading] < (2.0 + 0.006*leading_p4->Pt())
					&& l.pho_hcalsumetconedr03[leading] < (2.0 + 0.0025*leading_p4->Pt());
				
				//Selection on next to leading photon
				pass_selection[1] = l.pho_hoe[nleading] < 0.02
					&& (((l.pho_sieie[nleading] < 0.01)  && (fabs(nleading_calopos->Eta()) < 1.4442)) 
					    || (( l.pho_sieie[nleading] < 0.028)
						&& ((fabs(nleading_calopos->Eta()) < 2.5) && (fabs(nleading_calopos->Eta()) > 1.566))) );
				pass_isolation[1] =  l.pho_trksumptsolidconedr03[nleading] < (1.5 + 0.001*nleading_p4->Pt())
					&& (l.pho_ecalsumetconedr03[nleading] < (2.0 + 0.006*nleading_p4->Pt()))
					&& (l.pho_hcalsumetconedr03[nleading] < (2.0 + 0.0025*nleading_p4->Pt()));
				
				
				if (pass_selection[0] && pass_isolation[0] ){
					
					if (pass_selection[1] && pass_isolation[1]){
//						cout << "mass is " << mass << " and higgs pt is " << h_pt << endl;
						l.FillHist("pho_pt",category,leading_p4->Pt(),weight);
						l.FillHist("pho_pt",category,nleading_p4->Pt(),weight);
						best_mass = mass;
						best_pt   = h_pt;
					}
				}
			}
		}
	}
	
	if (best_mass > 100)
		l.FillHist("mass",0, best_mass, weight);
	l.FillHist("pt",0, best_pt, weight);
	if (category > -1){
		l.FillHist("mass",category, best_mass, weight);
		l.FillHist("pt",category, best_pt, weight);
	}
	if(PADEBUG) 
		cout<<"myFillHistRed END"<<endl;
	
	
	if( runStatAnalysis ) {
		StatAnalysis(l,jentry);
	}
	
	///// for (int i=0; i<l.pho_n; i++) {
	///// 	TLorentzVector *p4 = (TLorentzVector *) l.pho_p4->At(i);
	///// 	// FIXME move Fill* in BaseAnalysis
	///// 	l.FillHist("pho_pt", p4->Pt());
	///// }
  	///// 
	///// Int_t in_endcap = 0;
	///// Float_t best_mass = 0;
	///// for (int i=0; i<l.pho_n-1; i++) {
	///// 	TLorentzVector *pg1= (TLorentzVector *) l.pho_p4->At(i);
	///// 	if (fabs(pg1->Eta()) > 1.479)
	///// 		in_endcap = 1;
	///// 
	///// 	for (int j=i+1; j<l.pho_n; j++) {
	///// 		TLorentzVector *pg2= (TLorentzVector *) l.pho_p4->At(j);
	///// 		if (fabs(pg2->Eta()) > 1.479)
	///// 			in_endcap = 1;
	///// 		TLorentzVector higgs = (*pg1) + (*pg2);
	///// 		Float_t mass = higgs.M();
	///// 		if (mass > best_mass)
	///// 			best_mass = mass;
	///// 	}
	///// }     
	///// 
	///// if (best_mass != 0) 
	///// 	l.FillHist("invmass", best_mass);
  	///// 
	///// if(PADEBUG) 
	///// 	cout<<"myFillHist END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::StatAnalysis(LoopAll& l, Int_t jentry) 
{
}




// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
	vtxAna_.setBranchAdresses(t,"vtx_std_");
	vtxAna_.getBranches(t,"vtx_std_",s);
}


// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::PreselectPhotons(LoopAll& l, int jentry) 
{
	// Photon preselection
	pho_acc.clear();
	pho_presel.clear();
	pho_presel_lead.clear();
	pho_sc_et.clear();
	l.pho_matchingConv->clear();

	for(int ipho=0; ipho<l.pho_n; ++ipho) {

	  // match all photons in the original tree with the conversions from the merged collection and save the indices
	  int iConv  =l.matchPhotonToConversion(ipho);
	  if ( iConv>=0 )
		  (*l.pho_matchingConv).push_back(l.matchPhotonToConversion(ipho));
	  else
		  (*l.pho_matchingConv).push_back(-1);
	  


	  TLorentzVector * p4 = (TLorentzVector *) l.pho_p4->At(ipho);
	  float eta  = fabs(((TVector3 *) l.pho_calopos->At(ipho))->Eta());
	  // photon et wrt 0,0,0
	  float et = p4->Energy() / cosh(eta);
	  pho_sc_et.push_back(et);
	  /// std::cerr << " " << p4->Pt() << " " << et << " " << eta;
	  
	  if( et < presel_scet2 || (eta>1.4442 && eta<1.566) || eta>presel_maxeta ) { 
		  /// std::cerr << std::endl;
	    continue;  
	  }
	  /// std::cerr << "keeping " << ipho << std::endl;
	  pho_acc.push_back(ipho);
	  
	  bool isEB = l.pho_isEB[ipho];
	  float & ecaliso = isEB ? presel_ecaliso_eb : presel_ecaliso_ee;
	  float & sieie = isEB ? presel_ecaliso_eb : presel_ecaliso_ee;
	  if( l.pho_ecalsumetconedr03[ipho] >= ecaliso ||  l.pho_sieie[ipho] >= sieie || l.pho_hoe[ipho] >= presel_hoe ) {
	    continue;
	  }
	  
          
	  //FIXME trigger matching
	  pho_presel.push_back(ipho);
	} 

	// sort preslected photons by et
	std::sort(pho_acc.begin(),pho_acc.end(),
		  SimpleSorter<float,std::greater<float> >(&pho_sc_et[0]));
	std::sort(pho_presel.begin(),pho_presel.end(),
		  SimpleSorter<float,std::greater<float> >(&pho_sc_et[0]));

	if( pho_presel.size() > 1 ) {
		for(size_t ipho=0; ipho<pho_presel.size()-1; ++ipho ) {
			assert( pho_sc_et[pho_presel[ipho]] >= pho_sc_et[pho_presel[ipho+1]] );
		}
	}
	if( pho_acc.size()>1 ) {
		for(size_t ipho=0; ipho<pho_acc.size()-1; ++ipho ) {
			assert( pho_sc_et[pho_acc[ipho]] >= pho_sc_et[pho_acc[ipho+1]] );
		}
	}
	
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::FillReductionVariables(LoopAll& l, int jentry) 
{
	if(PADEBUG) 
		cout<<"myFillReduceVar START"<<endl;
  	
	PreselectPhotons(l,jentry);
		
	if(PADEBUG) 
		cout<<"myFillReduceVar END"<<endl;

}

// ----------------------------------------------------------------------------------------------------
bool PhotonAnalysis::SelectEventsReduction(LoopAll& l, int jentry) 
{

	if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
	// require at least two reconstructed photons to store the event
	if( pho_acc.size() < 2 || pho_sc_et[ pho_acc[0] ] < presel_scet1 ) { return false; }
	
	int ipho1 = pho_acc[0];
	int ipho2 = pho_acc[1];
	
	if( pho_presel.size() > 1 ) { 
		// use the first two preselected photons for the vertex algorithm
		ipho1 = pho_presel[0]; 
		ipho2 = pho_presel[1]; 
	} else if(pho_presel.size() > 0 ) {
		// if only one photon was preselected use the highest preselected and the higest non preselect photons 
		//    the event will be discarded by the analysis anyway
		ipho1 = pho_presel[0]; 
		ipho2 = pho_acc[0] == ipho1 ? pho_acc[1] : pho_acc[0]; 
	}
	assert( ipho1 != ipho2 );

	if(PADEBUG)        cout << " SelectEventsReduction going to fill photon info " << endl;
	PhotonInfo pho1=l.fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions);
        PhotonInfo pho2=l.fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions);
        if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;

        // run vertex analysis
	l.vertexAnalysis(vtxAna_, pho1, pho2 );
        // select vertxe
	*l.vtx_std_ranked_list = l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames);
	if( l.vtx_std_ranked_list->size() != 0 ) {  
		l.vtx_std_sel = (*l.vtx_std_ranked_list)[0];
	} else {
		l.vtx_std_sel = 0;
		std::cerr << "NO VERTEX SELECTED " << l.event << " " << l.run << " " << std::endl;
	}
	// update the photons' pt
	for(int ipho=0; ipho<l.pho_n; ++ipho) {
		l.set_pho_p4(ipho, l.vtx_std_sel);
	}
	
	// fill ID variables
	l.FillCICInputs();
	l.FillCIC();

	return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SkimEvents(LoopAll& l, int jentry)
{
	l.b_pho_n->GetEntry(jentry);
	if( l.pho_n < 2 ) {
		return false;
	}

	// do not run trigger selection on MC
	int filetype = l.itype[l.current];
	bool skipTrigger = ! doTriggerSelection || filetype != 0 || triggerSelections.empty();
	if( ! skipTrigger ) {
		// get the trigger selection for this run 
		l.b_run->GetEntry(jentry);
		std::vector<TriggerSelection>::iterator isel = find(triggerSelections.begin(), triggerSelections.end(), l.run );
		if(isel == triggerSelections.end() ) {
			std::cerr << "No trigger selection for run " << l.run << "defined" << std::endl;
			return true;
		}
		// get the trigegr data
		l.b_hlt1_bit->GetEntry(jentry);
		l.b_hlt_path_names_HLT1->GetEntry(jentry);
		if( !  isel->pass(*(l.hlt_path_names_HLT1),*(l.hlt1_bit)) ) {
			return false;
		}
		l.countersred[trigCounter_]++;
	}
		
	return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
	return true;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
	vtxAna_.branches(outputTree,"vtx_std_");	

	l.pho_matchingConv = new  std::vector<int>();
	l.Branch_pho_matchingConv(outputTree);


	l.vtx_std_ranked_list = new std::vector<int>();
	l.pho_tkiso_recvtx_030_002_0000_10_01 = new std::vector<std::vector<float> >();
	l.pho_cic6cutlevel_lead = new std::vector<Short_t>();
	l.pho_cic6passcuts_lead = new std::vector<std::vector<UInt_t> >();
	l.pho_cic6cutlevel_sublead = new std::vector<Short_t>();
	l.pho_cic6passcuts_sublead = new std::vector<std::vector<UInt_t> >();
	l.pho_cic4cutlevel_lead = new std::vector<Short_t>();
	l.pho_cic4passcuts_lead = new std::vector<std::vector<UInt_t> >();
	l.pho_cic4cutlevel_sublead = new std::vector<Short_t>();
	l.pho_cic4passcuts_sublead = new std::vector<std::vector<UInt_t> >();



	l.Branch_vtx_std_ranked_list(outputTree);
	l.Branch_vtx_std_sel(outputTree);
	l.Branch_pho_tkiso_recvtx_030_002_0000_10_01(outputTree);
	l.Branch_pho_tkiso_badvtx_040_002_0000_10_01(outputTree);
	l.Branch_pho_drtotk_25_99(outputTree);


  //l.Branch_dipho_n(outputTree);
  //l.Branch_dipho_leadind(outputTree);
  //l.Branch_dipho_subleadind(outputTree);
  //l.Branch_dipho_vtxind(outputTree);
  //l.Branch_dipho_leadet(outputTree);
  //l.Branch_dipho_subleadet(outputTree);
  //l.Branch_dipho_leadeta(outputTree);
  //l.Branch_dipho_subleadeta(outputTree);
  //l.Branch_dipho_leadci6cindex(outputTree);
  //l.Branch_dipho_subleadci6cindex(outputTree);
  //l.Branch_dipho_leadci4cindex(outputTree);
  //l.Branch_dipho_subleadci4cindex(outputTree);
  //l.Branch_dipho_mass(outputTree);
  //l.Branch_dipho_pt(outputTree);
  //l.Branch_dipho_eta(outputTree);
  //l.Branch_dipho_phi(outputTree);
  //l.Branch_dipho_cts(outputTree);

	l.Branch_pho_cic6cutlevel_lead( outputTree );
	l.Branch_pho_cic6passcuts_lead( outputTree );
	l.Branch_pho_cic6cutlevel_sublead( outputTree );
	l.Branch_pho_cic6passcuts_sublead( outputTree );
	l.Branch_pho_cic4cutlevel_lead( outputTree );
	l.Branch_pho_cic4passcuts_lead( outputTree );
	l.Branch_pho_cic4cutlevel_sublead( outputTree );
	l.Branch_pho_cic4passcuts_sublead( outputTree );
	
}

