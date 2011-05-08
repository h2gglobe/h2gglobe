#include "../interface/PhotonAnalysis.h"

#include "Sorters.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::PhotonAnalysis()  : 
	runStatAnalysis(false),
	name_("PhotonAnalysis"),
	vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{}

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
						cout << "mass is " << mass << " and higgs pt is " << h_pt << endl;
						l.FillHist("pho_pt",category,leading_p4->Pt());
						l.FillHist("pho_pt",category,nleading_p4->Pt());
						best_mass = mass;
						best_pt   = h_pt;
					}
				}
			}
		}
	}
	
	if (best_mass > 100)
		l.FillHist("mass",0, best_mass);
	l.FillHist("pt",0, best_pt);
	if (category > -1){
		l.FillHist("mass",category, best_mass);
		l.FillHist("pt",category, best_pt);
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
	for(int ipho=0; ipho<l.pho_n; ++ipho) {
		
		TLorentzVector * p4 = (TLorentzVector *) l.pho_p4->At(ipho);
		float eta  = fabs(((TVector3 *) l.pho_calopos->At(ipho))->Eta());
		// photon et wrt 0,0,0
		float et = p4->Energy() / cosh(eta);
		pho_sc_et.push_back(et);
		
		if( et < presel_scet2 || (eta>1.4442 && eta<1.566) || eta>presel_maxeta ) { 
			continue;  
		}
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

        // run vertex analysis
	l.vertexAnalysis(vtxAna_, ipho1, ipho2 );
	
	// fill track isolation
	l.fillTrackIsolation(tkIso_ptmin,tkIso_outerCone,tkIso_innerCone,tkIso_etaStripHalfW,tkIso_dzmax,tkIso_dxymax);

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

	l.pho_ntrkhgg = new std::vector<std::vector<int> >();
	l.pho_trksumpthgg = new std::vector<std::vector<float> >();
	outputTree->Branch("pho_trksumpthgg", &l.pho_trksumpthgg);
	outputTree->Branch("pho_ntrkhgg", &l.pho_ntrkhgg);
	
	outputTree->Print();
}

