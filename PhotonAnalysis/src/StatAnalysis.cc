#include "../interface/StatAnalysis.h"

#include "Sorters.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
StatAnalysis::StatAnalysis()  : 
	runStatAnalysis(false),
	name_("StatAnalysis"),
	vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{
}

// ----------------------------------------------------------------------------------------------------
StatAnalysis::~StatAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Term(LoopAll& l) 
{
      l.rooContainer->FitToData("exp","bkg_mass",95,105,145,200);
      l.rooContainer->FitToSystematicSet("exp","bkg_mass","e-scale");

      
      std::string outputfilename = (std::string) l.histFileName;  
      l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m100","bkg_mass");
      l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m110","bkg_mass");
      l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m120","bkg_mass");
      l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m130","bkg_mass");
      l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m140","bkg_mass");

}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Init(LoopAll& l) 
{
	if(PADEBUG) 
		cout << "InitRealStatAnalysis START"<<endl;

	// call the base class initializer
	PhotonAnalysis::Init(l);
	
	// next three lines should ideally be configured @ run-time 
	l.rooContainer->SetNCategories(4);
	std::vector<std::string> sys(1,"e-scale");
	l.rooContainer->MakeSystematicStudy(sys);
	// ----------------------------------------------------

	l.rooContainer->AddRealVar("data_mass" ,95.,155.);
	l.rooContainer->AddRealVar("bkg_mass" ,95.,155.);
	l.rooContainer->AddRealVar("sig_mass_m100" ,95.,155.);
	l.rooContainer->AddRealVar("sig_mass_m110" ,95.,155.);
	l.rooContainer->AddRealVar("sig_mass_m120" ,95.,155.);
	l.rooContainer->AddRealVar("sig_mass_m130" ,95.,155.);
	l.rooContainer->AddRealVar("sig_mass_m140" ,95.,155.);
	l.rooContainer->AddRealVar("mu",-0.04,-2.,-0.001);
		  
	// -------------------------------------//
	std::vector<std::string> pars(2,"t");	 
	pars[0] = "mass";
	pars[1] = "mu";
	// -------------------------------------//

	l.rooContainer->AddGenericPdf("exp",
	  "exp((@0)*(@1))",pars,1,1000);
		  
	l.rooContainer->CreateDataSet("data_mass",30);
	l.rooContainer->CreateDataSet("bkg_mass",30);

	l.rooContainer->CreateDataSet("sig_mass_m100",30);
	l.rooContainer->CreateDataSet("sig_mass_m110",30);
	l.rooContainer->CreateDataSet("sig_mass_m120",30);
	l.rooContainer->CreateDataSet("sig_mass_m130",30);
	l.rooContainer->CreateDataSet("sig_mass_m140",30);
	l.rooContainer->MakeSystematics("mass","e-scale");
	
	if(PADEBUG) 
		cout << "InitRealStatAnalysis END"<<endl;

	// FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
   if(PADEBUG) 
 	cout << "Analysis START"<<endl;
   
   int cur_type = l.itype[l.current];
   float weight = l.sampleContainer[l.current_sample_index].weight;

   std::pair<int,int> diphoton_index;
   int leadLevel=2, subLevel=2;

   float leadEtCut   = 40.0;
   float subleadEtCut = 35.0;

   diphoton_index = l.DiphotonCiCSelection((LoopAll::phoCiCIDLevel) leadLevel, (LoopAll::phoCiCIDLevel) subLevel, leadEtCut, subleadEtCut, false); 
 
  if (diphoton_index.first > -1 && diphoton_index.second > -1){

   	 TLorentzVector *lead_p4 = (TLorentzVector*)l.pho_p4->At(diphoton_index.first);
   	 TLorentzVector *sublead_p4 = (TLorentzVector*)l.pho_p4->At(diphoton_index.second);
  	 TLorentzVector Higgs = *lead_p4 + *sublead_p4; 	

	 float mass = Higgs.M();

         TVector3 *calopos1 = (TVector3*) l.pho_calopos->At(diphoton_index.first);
         TVector3 *calopos2 = (TVector3*) l.pho_calopos->At(diphoton_index.second);
         float min_r9  = min(l.pho_r9[diphoton_index.first],l.pho_r9[diphoton_index.second]);
	 float max_eta = max(fabs(calopos1->Eta()),fabs(calopos2->Eta()));

         int category;
	 if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 0;
	 if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 1;
	 if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 2;
	 if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;
 
           float sys_error = 0.05;
 
            std::vector<float> mass_errors;
            for (int sys=1;sys<31;sys++){
                 float incr = (float)sys/10;
                 TLorentzVector lead_err = (1.-incr*sys_error)*(*lead_p4);
                 TLorentzVector nlead_err = (1.-incr*sys_error)*(*sublead_p4);
                 mass_errors.push_back((lead_err+nlead_err).M());
            }
            for (int sys=1;sys<31;sys++){
                 float incr = (float)sys/10;
                 TLorentzVector lead_err = (1.+incr*sys_error)*(*lead_p4);
                 TLorentzVector nlead_err = (1.+incr*sys_error)*(*sublead_p4);
                 mass_errors.push_back((lead_err+nlead_err).M());
            }
 

        l.rooContainer->InputSystematicSet("bkg_mass","e-scale",category,mass_errors);
         
	 if (cur_type == 0){
	   l.rooContainer->InputDataPoint("data_mass",category,mass);
	}
	 else if (cur_type > 0)
	   l.rooContainer->InputDataPoint("bkg_mass",category,mass,weight);
	 else if (cur_type == -1 || cur_type == -2 || cur_type == -3)
	   l.rooContainer->InputDataPoint("sig_mass_m100",category,mass,weight);
	 else if (cur_type == -4 || cur_type == -5 || cur_type == -6)
	   l.rooContainer->InputDataPoint("sig_mass_m110",category,mass,weight);
	 else if (cur_type == -7 || cur_type == -8 || cur_type == -9)
	   l.rooContainer->InputDataPoint("sig_mass_m120",category,mass,weight);
	 else if (cur_type == -10 || cur_type == -11 || cur_type == -12)
	   l.rooContainer->InputDataPoint("sig_mass_m130",category,mass,weight);
	 else if (cur_type == -13 || cur_type == -14 || cur_type == -15)
	   l.rooContainer->InputDataPoint("sig_mass_m140",category,mass,weight);
	
	}
	
	if(PADEBUG) 
		cout<<"myFillHistRed END"<<endl;
}

// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
	vtxAna_.setBranchAdresses(t,"vtx_std_");
	vtxAna_.getBranches(t,"vtx_std_",s);
}


// ----------------------------------------------------------------------------------------------------
void StatAnalysis::PreselectPhotons(LoopAll& l, int jentry) 
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
bool StatAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
	if (l.itype[l.current] != 0) return true;

	if (l.run < 160404) return true;
  if (l.run == 160431 && l.lumis>=19 && l.lumis<=218) return true;
  if (l.run == 160577 && l.lumis>=254 && l.lumis<=306) return true;
  if (l.run == 160578 && l.lumis>=6 && l.lumis<=53) return true;
  if (l.run == 160578 && l.lumis>=274 && l.lumis<=400) return true;
  if (l.run == 160871 && l.lumis>=68 && l.lumis<=208) return true;
  if (l.run == 160872 && l.lumis>=1 && l.lumis<=9) return true;
  if (l.run == 160872 && l.lumis>=25 && l.lumis<=35) return true;
  if (l.run == 160872 && l.lumis>=38 && l.lumis<=55) return true;
  if (l.run == 160873 && l.lumis>=1 && l.lumis<=147) return true;
  if (l.run == 160874 && l.lumis>=1 && l.lumis<=51) return true;
  if (l.run == 160874 && l.lumis>=97 && l.lumis<=113) return true;
  if (l.run == 160939 && l.lumis>=1 && l.lumis<=123) return true;
  if (l.run == 160940 && l.lumis>=1 && l.lumis<=79) return true;
  if (l.run == 160942 && l.lumis>=1 && l.lumis<=12) return true;
  if (l.run == 160943 && l.lumis>=1 && l.lumis<=54) return true;
  if (l.run == 160955 && l.lumis>=1 && l.lumis<=130) return true;
  if (l.run == 160955 && l.lumis>=133 && l.lumis<=138) return true;
  if (l.run == 160955 && l.lumis>=140 && l.lumis<=151) return true;
  if (l.run == 160955 && l.lumis>=153 && l.lumis<=154) return true;
  if (l.run == 160955 && l.lumis>=156 && l.lumis<=172) return true;
  if (l.run == 160955 && l.lumis>=175 && l.lumis<=201) return true;
  if (l.run == 160955 && l.lumis>=204 && l.lumis<=206) return true;
  if (l.run == 160956 && l.lumis>=2 && l.lumis<=65) return true;
  if (l.run == 160957 && l.lumis>=1 && l.lumis<=953) return true;
  if (l.run == 160998 && l.lumis>=2 && l.lumis<=252) return true;
  if (l.run == 161008 && l.lumis>=1 && l.lumis<=77) return true;
  if (l.run == 161016 && l.lumis>=2 && l.lumis<=300) return true;
  if (l.run == 162803 && l.lumis>=60 && l.lumis<=124) return true;
  if (l.run == 162803 && l.lumis>=135 && l.lumis<=139) return true;
  if (l.run == 162808 && l.lumis>=1 && l.lumis<=51) return true;
  if (l.run == 162811 && l.lumis>=1 && l.lumis<=340) return true;
  if (l.run == 162822 && l.lumis>=73 && l.lumis<=307) return true;
  if (l.run == 162825 && l.lumis>=1 && l.lumis<=184) return true;
  if (l.run == 162826 && l.lumis>=1 && l.lumis<=24) return true;
  if (l.run == 162828 && l.lumis>=1 && l.lumis<=85) return true;
  if (l.run == 162909 && l.lumis>=54 && l.lumis<=290) return true;
  if (l.run == 163046 && l.lumis>=1 && l.lumis<=133) return true;
  if (l.run == 163046 && l.lumis>=135 && l.lumis<=238) return true;
  if (l.run == 163069 && l.lumis>=73 && l.lumis<=452) return true;
  if (l.run == 163069 && l.lumis>=468 && l.lumis<=633) return true;
  if (l.run == 163071 && l.lumis>=1 && l.lumis<=161) return true;
  if (l.run == 163078 && l.lumis>=1 && l.lumis<=23) return true;
  if (l.run == 163232 && l.lumis>=110 && l.lumis<=149) return true;
  if (l.run == 163233 && l.lumis>=1 && l.lumis<=283) return true;
  if (l.run == 163234 && l.lumis>=1 && l.lumis<=66) return true;
  if (l.run == 163235 && l.lumis>=1 && l.lumis<=461) return true;
  if (l.run == 163237 && l.lumis>=1 && l.lumis<=213) return true;
  if (l.run == 163238 && l.lumis>=9 && l.lumis<=15) return true;
  if (l.run == 163252 && l.lumis>=60 && l.lumis<=137) return true;
  if (l.run == 163255 && l.lumis>=1 && l.lumis<=359) return true;
  if (l.run == 163255 && l.lumis>=412 && l.lumis<=844) return true;
  if (l.run == 163255 && l.lumis>=846 && l.lumis<=846) return true;
  if (l.run == 163255 && l.lumis>=848 && l.lumis<=977) return true;
  if (l.run == 163261 && l.lumis>=1 && l.lumis<=3) return true;
  if (l.run == 163261 && l.lumis>=10 && l.lumis<=126) return true;
  if (l.run == 163270 && l.lumis>=1 && l.lumis<=76) return true;
  if (l.run == 163270 && l.lumis>=79 && l.lumis<=96) return true;
  if (l.run == 163270 && l.lumis>=99 && l.lumis<=475) return true;
  if (l.run == 163270 && l.lumis>=479 && l.lumis<=527) return true;
  if (l.run == 163270 && l.lumis>=529 && l.lumis<=685) return true;
  if (l.run == 163270 && l.lumis>=695 && l.lumis<=928) return true;
  if (l.run == 163286 && l.lumis>=112 && l.lumis<=401) return true;
  if (l.run == 163289 && l.lumis>=1 && l.lumis<=388) return true;
  if (l.run == 163296 && l.lumis>=59 && l.lumis<=230) return true;
  if (l.run == 163296 && l.lumis>=232 && l.lumis<=585) return true;
  if (l.run == 163297 && l.lumis>=1 && l.lumis<=219) return true;
  if (l.run == 163300 && l.lumis>=1 && l.lumis<=616) return true;
  if (l.run == 163301 && l.lumis>=1 && l.lumis<=192) return true;
  if (l.run == 163302 && l.lumis>=1 && l.lumis<=190) return true;
  if (l.run == 163332 && l.lumis>=43 && l.lumis<=118) return true;
  if (l.run == 163332 && l.lumis>=224 && l.lumis<=264) return true;
  if (l.run == 163332 && l.lumis>=266 && l.lumis<=599) return true;
  if (l.run == 163332 && l.lumis>=601 && l.lumis<=639) return true;
  if (l.run == 163332 && l.lumis>=641 && l.lumis<=801) return true;
  if (l.run == 163333 && l.lumis>=1 && l.lumis<=106) return true;
  if (l.run == 163334 && l.lumis>=1 && l.lumis<=35) return true;
  if (l.run == 163334 && l.lumis>=37 && l.lumis<=37) return true;
  if (l.run == 163334 && l.lumis>=156 && l.lumis<=556) return true;
  if (l.run == 163337 && l.lumis>=1 && l.lumis<=18) return true;
  if (l.run == 163337 && l.lumis>=27 && l.lumis<=201) return true;
  if (l.run == 163337 && l.lumis>=203 && l.lumis<=426) return true;
  if (l.run == 163337 && l.lumis>=434 && l.lumis<=461) return true;
  if (l.run == 163338 && l.lumis>=1 && l.lumis<=164) return true;
  if (l.run == 163339 && l.lumis>=1 && l.lumis<=172) return true;
  if (l.run == 163340 && l.lumis>=1 && l.lumis<=488) return true;
  if (l.run == 163358 && l.lumis>=39 && l.lumis<=63) return true;
  if (l.run == 163369 && l.lumis>=1 && l.lumis<=94) return true;
  if (l.run == 163370 && l.lumis>=1 && l.lumis<=147) return true;
  if (l.run == 163371 && l.lumis>=1 && l.lumis<=107) return true;
  if (l.run == 163371 && l.lumis>=148 && l.lumis<=363) return true;
  if (l.run == 163372 && l.lumis>=1 && l.lumis<=52) return true;
  if (l.run == 163374 && l.lumis>=1 && l.lumis<=599) return true;
  if (l.run == 163374 && l.lumis>=603 && l.lumis<=863) return true;
  if (l.run == 163375 && l.lumis>=1 && l.lumis<=10) return true;
  if (l.run == 163376 && l.lumis>=1 && l.lumis<=20) return true;
  if (l.run == 163376 && l.lumis>=22 && l.lumis<=246) return true;
  if (l.run == 163378 && l.lumis>=1 && l.lumis<=81) return true;
  if (l.run == 163378 && l.lumis>=89 && l.lumis<=272) return true;
  if (l.run == 163378 && l.lumis>=306 && l.lumis<=615) return true;
  if (l.run == 163385 && l.lumis>=52 && l.lumis<=240) return true;
  if (l.run == 163385 && l.lumis>=244 && l.lumis<=406) return true;
  if (l.run == 163387 && l.lumis>=1 && l.lumis<=256) return true;
  if (l.run == 163402 && l.lumis>=37 && l.lumis<=582) return true;
  if (l.run == 163402 && l.lumis>=586 && l.lumis<=801) return true;
  if (l.run == 163475 && l.lumis>=30 && l.lumis<=295) return true;
  if (l.run == 163476 && l.lumis>=1 && l.lumis<=94) return true;
  if (l.run == 163476 && l.lumis>=98 && l.lumis<=212) return true;
  if (l.run == 163478 && l.lumis>=1 && l.lumis<=70) return true;
  if (l.run == 163479 && l.lumis>=1 && l.lumis<=175) return true;
  if (l.run == 163480 && l.lumis>=1 && l.lumis<=92) return true;
  if (l.run == 163480 && l.lumis>=96 && l.lumis<=188) return true;
  if (l.run == 163480 && l.lumis>=190 && l.lumis<=191) return true;
  if (l.run == 163481 && l.lumis>=1 && l.lumis<=72) return true;
  if (l.run == 163481 && l.lumis>=74 && l.lumis<=77) return true;
  if (l.run == 163481 && l.lumis>=79 && l.lumis<=79) return true;
  if (l.run == 163482 && l.lumis>=1 && l.lumis<=27) return true;
  if (l.run == 163482 && l.lumis>=48 && l.lumis<=48) return true;
  if (l.run == 163483 && l.lumis>=1 && l.lumis<=57) return true;
  if (l.run == 163582 && l.lumis>=1 && l.lumis<=22) return true;
  if (l.run == 163583 && l.lumis>=1 && l.lumis<=63) return true;
  if (l.run == 163583 && l.lumis>=65 && l.lumis<=92) return true;
  if (l.run == 163583 && l.lumis>=96 && l.lumis<=155) return true;
  if (l.run == 163583 && l.lumis>=157 && l.lumis<=173) return true;
  if (l.run == 163583 && l.lumis>=175 && l.lumis<=219) return true;
  if (l.run == 163584 && l.lumis>=1 && l.lumis<=56) return true;
  if (l.run == 163585 && l.lumis>=1 && l.lumis<=32) return true;
  if (l.run == 163586 && l.lumis>=1 && l.lumis<=75) return true;
  if (l.run == 163587 && l.lumis>=1 && l.lumis<=52) return true;
  if (l.run == 163588 && l.lumis>=1 && l.lumis<=8) return true;
  if (l.run == 163588 && l.lumis>=10 && l.lumis<=446) return true;
  if (l.run == 163589 && l.lumis>=1 && l.lumis<=49) return true;
  if (l.run == 163589 && l.lumis>=51 && l.lumis<=160) return true;
  if (l.run == 163596 && l.lumis>=1 && l.lumis<=29) return true;
  if (l.run == 163630 && l.lumis>=76 && l.lumis<=164) return true;
  if (l.run == 163630 && l.lumis>=176 && l.lumis<=185) return true;
  if (l.run == 163655 && l.lumis>=15 && l.lumis<=23) return true;
  if (l.run == 163657 && l.lumis>=1 && l.lumis<=140) return true;
  if (l.run == 163658 && l.lumis>=1 && l.lumis<=3) return true;
  if (l.run == 163659 && l.lumis>=1 && l.lumis<=374) return true;
  if (l.run == 163659 && l.lumis>=376 && l.lumis<=650) return true;
  if (l.run == 163659 && l.lumis>=652 && l.lumis<=705) return true;
  if (l.run == 163660 && l.lumis>=1 && l.lumis<=74) return true;
  if (l.run == 163661 && l.lumis>=1 && l.lumis<=17) return true;
  if (l.run == 163662 && l.lumis>=1 && l.lumis<=154) return true;
  if (l.run == 163663 && l.lumis>=1 && l.lumis<=106) return true;
  if (l.run == 163663 && l.lumis>=109 && l.lumis<=246) return true;
  if (l.run == 163664 && l.lumis>=1 && l.lumis<=119) return true;
  if (l.run == 163664 && l.lumis>=121 && l.lumis<=178) return true;
  if (l.run == 163668 && l.lumis>=1 && l.lumis<=53) return true;
  if (l.run == 163668 && l.lumis>=57 && l.lumis<=136) return true;
  if (l.run == 163668 && l.lumis>=140 && l.lumis<=213) return true;
  if (l.run == 163738 && l.lumis>=34 && l.lumis<=311) return true;
  if (l.run == 163757 && l.lumis>=1 && l.lumis<=40) return true;
  if (l.run == 163758 && l.lumis>=1 && l.lumis<=17) return true;
  if (l.run == 163758 && l.lumis>=19 && l.lumis<=220) return true;
  if (l.run == 163758 && l.lumis>=222 && l.lumis<=224) return true;
  if (l.run == 163758 && l.lumis>=236 && l.lumis<=276) return true;
  if (l.run == 163758 && l.lumis>=283 && l.lumis<=374) return true;
  if (l.run == 163758 && l.lumis>=376 && l.lumis<=466) return true;
  if (l.run == 163758 && l.lumis>=468 && l.lumis<=591) return true;
  if (l.run == 163759 && l.lumis>=1 && l.lumis<=60) return true;
  if (l.run == 163759 && l.lumis>=62 && l.lumis<=72) return true;
  if (l.run == 163759 && l.lumis>=74 && l.lumis<=456) return true;
  if (l.run == 163759 && l.lumis>=458 && l.lumis<=461) return true;
  if (l.run == 163759 && l.lumis>=463 && l.lumis<=482) return true;
  if (l.run == 163759 && l.lumis>=504 && l.lumis<=510) return true;
  if (l.run == 163760 && l.lumis>=1 && l.lumis<=162) return true;
  if (l.run == 163760 && l.lumis>=165 && l.lumis<=340) return true;
  if (l.run == 163761 && l.lumis>=1 && l.lumis<=203) return true;
  if (l.run == 163763 && l.lumis>=1 && l.lumis<=79) return true;
  if (l.run == 163765 && l.lumis>=1 && l.lumis<=321) return true;
  if (l.run == 163795 && l.lumis>=10 && l.lumis<=34) return true;
  if (l.run == 163795 && l.lumis>=36 && l.lumis<=36) return true;
  if (l.run == 163795 && l.lumis>=38 && l.lumis<=43) return true;
  if (l.run == 163796 && l.lumis>=1 && l.lumis<=182) return true;
  if (l.run == 163817 && l.lumis>=50 && l.lumis<=140) return true;
  if (l.run == 163817 && l.lumis>=154 && l.lumis<=205) return true;
  if (l.run == 163817 && l.lumis>=216 && l.lumis<=295) return true;
  if (l.run == 163817 && l.lumis>=305 && l.lumis<=346) return true;
  if (l.run == 163817 && l.lumis>=358 && l.lumis<=457) return true;
  if (l.run == 163817 && l.lumis>=561 && l.lumis<=603) return true;
  if (l.run == 163817 && l.lumis>=618 && l.lumis<=966) return true;
  if (l.run == 163869 && l.lumis>=79 && l.lumis<=123) return true;
  return false;
}
