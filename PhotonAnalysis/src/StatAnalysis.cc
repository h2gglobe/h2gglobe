#include "../interface/StatAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
StatAnalysis::StatAnalysis()  : 
	name_("StatAnalysis"),
	vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{
}

// ----------------------------------------------------------------------------------------------------
StatAnalysis::~StatAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Term(LoopAll& l) 
{

        std::string outputfilename = (std::string) l.histFileName;
        // Make Fits to the data-sets and systematic sets
        l.rooContainer->FitToData("data_pol_model","data_mass");  // Fit to full range of dataset
  
        l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","data_pol_model");
        // mode 0 as above, 1 if want to bin in sub range from fit,

        // Write the data-card for the Combinations Code, needs the output filename, makes binned analysis DataCard
	// Assumes the signal datasets will be called signal_name+"_mXXX"
        l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","data_pol_model","data_mass",1,50,1); // 1 means systematics from the fit effect only the backgroundi. last digit mode = 1 means this is an internal constraint fit 
        l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","bkg_mass_rebinned");

	SaclayText.close();

//	kfacFile->Close();
//	PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Init(LoopAll& l) 
{
	if(PADEBUG) 
		cout << "InitRealStatAnalysis START"<<endl;
	
	// Saclay text file
		SaclayText.open("ascii_events_204pb.txt");
  	//
	// These parameters are set in the configuration file
 	std::cout
		<< "\n"
		<< "-------------------------------------------------------------------------------------- \n"
		<< "StatAnalysis " << "\n"
		<< "-------------------------------------------------------------------------------------- \n"
		<< "leadEtCut "<< leadEtCut << "\n"
		<< "subleadEtCut "<< subleadEtCut << "\n"
		<< "doTriggerSelection "<< doTriggerSelection << "\n"
		<< "nEtaCategories "<< nEtaCategories << "\n"
		<< "nR9Categories "<< nR9Categories << "\n"		
		<< "nPtCategories "<< nPtCategories << "\n"		
		<< "doEscaleSyst "<< doEscaleSyst << "\n"
		<< "doEresolSyst "<< doEresolSyst << "\n"
		<< "efficiencyFile " << efficiencyFile << "\n"
		<< "doPhotonIdEffSyst "<< doPhotonIdEffSyst << "\n"
		<< "doVtxEffSyst "<< doVtxEffSyst << "\n"
		<< "doTriggerEffSyst "<< doTriggerEffSyst << "\n"
		<< "-------------------------------------------------------------------------------------- \n"
		<< std::endl;

	// avoid recalculated the CIC ID every time
	l.runCiC = false;
	// call the base class initializer
	PhotonAnalysis::Init(l);
	
	diPhoCounter_ = l.countersred.size();
	l.countersred.resize(diPhoCounter_+1);

	// initialize the analysis variables
	nCategories_ = nEtaCategories;
	if( nR9Categories != 0 ) nCategories_ *= nR9Categories;
	if( nPtCategories != 0 ) nCategories_ *= nPtCategories;

	// Numbers from Riccardo: "slide 1" of https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/252/1.html
	// scale variations defined as: (deltaM_data - deltaM_MC) / M_Z   =>   scale factor >0 means energy should increased to MC
	float scale_offset_EBHighR9         =  0.49e-2;
	float scale_offset_EBLowR9          = -0.11e-2;
	float scale_offset_EEHighR9         = -0.57e-2;
	float scale_offset_EELowR9          = 0.39e-2;
	float scale_offset_error_EBHighR9   = 0.07e-2;
	float scale_offset_error_EBLowR9    = 0.05e-2;
	float scale_offset_error_EEHighR9   = 0.24e-2;
	float scale_offset_error_EELowR9    = 0.19e-2;
	float smearing_sigma_EBHighR9       = 0.89e-2;
	float smearing_sigma_EBLowR9        = 1.99e-2;
	float smearing_sigma_EEHighR9       = 4.09e-2;
	float smearing_sigma_EELowR9        = 2.46e-2;
	float smearing_sigma_error_EBHighR9 = 0.25e-2;
	float smearing_sigma_error_EBLowR9  = 0.12e-2;
	float smearing_sigma_error_EEHighR9 = 0.38e-2;
	float smearing_sigma_error_EELowR9  = 0.52e-2;
	



	eSmearPars.categoryType = "2CatR9_EBEE";
	eSmearPars.n_categories = 4;

	// E scale is shifted for data, NOT for MC 
	eSmearPars.scale_offset["EBHighR9"] = 0.;
	eSmearPars.scale_offset["EBLowR9"]  = 0.;
	eSmearPars.scale_offset["EEHighR9"] = 0.;
	eSmearPars.scale_offset["EELowR9"]  = 0.;
	// E scale systematics are applied to MC, NOT to data
	eSmearPars.scale_offset_error["EBHighR9"] = scale_offset_error_EBHighR9;
	eSmearPars.scale_offset_error["EBLowR9"]  = scale_offset_error_EBLowR9;
	eSmearPars.scale_offset_error["EEHighR9"] = scale_offset_error_EEHighR9;
	eSmearPars.scale_offset_error["EELowR9"]  = scale_offset_error_EELowR9;
	// E resolution smearing applied to MC 
	eSmearPars.smearing_sigma["EBHighR9"] = smearing_sigma_EBHighR9;
	eSmearPars.smearing_sigma["EBLowR9"]  = smearing_sigma_EBLowR9;
	eSmearPars.smearing_sigma["EEHighR9"] = smearing_sigma_EEHighR9;
	eSmearPars.smearing_sigma["EELowR9"]  = smearing_sigma_EELowR9;
	// E resolution systematics applied to MC 
	eSmearPars.smearing_sigma_error["EBHighR9"] = smearing_sigma_error_EBHighR9;
	eSmearPars.smearing_sigma_error["EBLowR9"]  = smearing_sigma_error_EBLowR9;
	eSmearPars.smearing_sigma_error["EEHighR9"] = smearing_sigma_error_EEHighR9;
	eSmearPars.smearing_sigma_error["EELowR9"]  = smearing_sigma_error_EELowR9;

	eSmearDataPars.categoryType = "2CatR9_EBEE";
	eSmearDataPars.n_categories = 4;

	// initialize smearer specific to energy shifts in DATA; use opposite of energy scale shift
	eSmearDataPars.scale_offset["EBHighR9"] = -1*scale_offset_EBHighR9;
	eSmearDataPars.scale_offset["EBLowR9"]  = -1*scale_offset_EBLowR9;
	eSmearDataPars.scale_offset["EEHighR9"] = -1*scale_offset_EEHighR9;
	eSmearDataPars.scale_offset["EELowR9"]  = -1*scale_offset_EELowR9;
	// no energy scale systematics applied to data
	eSmearDataPars.scale_offset_error["EBHighR9"] = 0.;
	eSmearDataPars.scale_offset_error["EBLowR9"]  = 0.;
	eSmearDataPars.scale_offset_error["EEHighR9"] = 0.;
	eSmearDataPars.scale_offset_error["EELowR9"]  = 0.;
	// E resolution smearing NOT applied to data 
	eSmearDataPars.smearing_sigma["EBHighR9"] = 0.;
	eSmearDataPars.smearing_sigma["EBLowR9"]  = 0.;
	eSmearDataPars.smearing_sigma["EEHighR9"] = 0.;
	eSmearDataPars.smearing_sigma["EELowR9"]  = 0.;
	// E resolution systematics NOT applied to data 
	eSmearDataPars.smearing_sigma_error["EBHighR9"] = 0.;
	eSmearDataPars.smearing_sigma_error["EBLowR9"]  = 0.;
	eSmearDataPars.smearing_sigma_error["EEHighR9"] = 0.;
	eSmearDataPars.smearing_sigma_error["EELowR9"]  = 0.;


	effSmearPars.categoryType = "2CatR9_EBEE";
	effSmearPars.n_categories = 4;
	effSmearPars.efficiency_file = efficiencyFile;

	diPhoEffSmearPars.n_categories = 8;
	diPhoEffSmearPars.efficiency_file = efficiencyFile;

	// energy scale systematics to MC
	eScaleSmearer = new EnergySmearer( eSmearPars );
	eScaleSmearer->name("E_scale");
	eScaleSmearer->doEnergy(true);
	eScaleSmearer->scaleOrSmear(true);
	photonSmearers_.push_back(eScaleSmearer);
	// energy scale corrections to Data
	eScaleDataSmearer = new EnergySmearer( eSmearDataPars );
	eScaleDataSmearer->name("E_scale_data");
	eScaleDataSmearer->doEnergy(true);
	eScaleDataSmearer->scaleOrSmear(true);
	//photonDataSmearers_.push_back(eScaleDataSmearer); // must not be included among MC smearers; will be singled out upon need // GF
	// energy resolution smearing
	eResolSmearer = new EnergySmearer( eSmearPars );
	eResolSmearer->name("E_res");
	eResolSmearer->doEnergy(false);
	eResolSmearer->scaleOrSmear(false);
	photonSmearers_.push_back(eResolSmearer);
	// photon ID efficiency 
	idEffSmearer = new EfficiencySmearer( effSmearPars );
	idEffSmearer->name("idEff");
	idEffSmearer->setEffName("ratioTP");
	idEffSmearer->init();
	photonSmearers_.push_back(idEffSmearer);
	// Vertex ID
	vtxEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );   // triplicate TF1's here
	vtxEffSmearer->name("vtxEff");
	vtxEffSmearer->setEffName("ratioVertex");
	vtxEffSmearer->doVtxEff(true);
	vtxEffSmearer->init();
 	diPhotonSmearers_.push_back(vtxEffSmearer);
	// trigger efficiency
	triggerEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
	triggerEffSmearer->name("triggerEff");
	triggerEffSmearer->setEffName("effL1HLT");
	triggerEffSmearer->doVtxEff(false);
	triggerEffSmearer->init();
 	diPhotonSmearers_.push_back(triggerEffSmearer);
	
        // Define the number of categories for the statistical analysis and
	// the systematic sets to be formed

	// FIXME move these params to config file
	l.rooContainer->SetNCategories(nCategories_);
	systRange  = 1.; // in units of sigma
	nSystSteps = 3;    

	if( doEscaleSyst ) {
		systPhotonSmearers_.push_back( eScaleSmearer );
		std::vector<std::string> sys(1,eScaleSmearer->name());
		std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
		l.rooContainer->MakeSystematicStudy(sys,sys_t);
	}
	if( doEresolSyst ) {
		systPhotonSmearers_.push_back( eResolSmearer );
		std::vector<std::string> sys(1,eResolSmearer->name());
		std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
		l.rooContainer->MakeSystematicStudy(sys,sys_t);
	}
	if( doPhotonIdEffSyst ) {
		systPhotonSmearers_.push_back( idEffSmearer );
		std::vector<std::string> sys(1,idEffSmearer->name());
		std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
		l.rooContainer->MakeSystematicStudy(sys,sys_t);
	}
	if( doVtxEffSyst ) {
	        systDiPhotonSmearers_.push_back( vtxEffSmearer );
		std::vector<std::string> sys(1,vtxEffSmearer->name());
		std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
		l.rooContainer->MakeSystematicStudy(sys,sys_t);
	}
	if( doTriggerEffSyst ) {
	        systDiPhotonSmearers_.push_back( triggerEffSmearer );
		std::vector<std::string> sys(1,triggerEffSmearer->name());
		std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
		l.rooContainer->MakeSystematicStudy(sys,sys_t);
	}
	// ----------------------------------------------------
	// ----------------------------------------------------
	// Global systematics - Lumi
	l.rooContainer->AddGlobalSystematic("lumi",1.04,1.00);
	// ----------------------------------------------------

	// Create observables for shape-analysis with ranges
	l.rooContainer->AddObservable("mass" ,100.,150.);

	// FIXME, get these numbers from the LoopAll or maybe sampleContainer?
	l.rooContainer->AddConstant("IntLumi",204.);
	l.rooContainer->AddConstant("XSBR_105",0.0387684+0.00262016+0.003037036);
	l.rooContainer->AddConstant("XSBR_110",0.0390848+0.00275406+0.002902204);
	l.rooContainer->AddConstant("XSBR_115",0.0386169+0.00283716+0.002717667);
	l.rooContainer->AddConstant("XSBR_120",0.0374175+0.00285525+0.002286);
	l.rooContainer->AddConstant("XSBR_130",0.0319112+0.00260804+0.0019327068);
	l.rooContainer->AddConstant("XSBR_140",0.0235322+0.00204088+0.0012874228);	

	l.rooContainer->AddRealVar("pol0",-0.01,-1.5,1.5);
	l.rooContainer->AddRealVar("pol1",-0.01,-1.5,1.5);
	std::vector<std::string> data_pol_pars(2,"p");	 
	data_pol_pars[0] = "pol0";
	data_pol_pars[1] = "pol1";
	l.rooContainer->AddGenericPdf("data_pol_model",
	  "0","mass",data_pol_pars,62);	// >= 61 means RooPolynomial of order >= 1
	// ------------------------------------------------------
        
        // Background - Exponential
	l.rooContainer->AddRealVar("mu",-0.04);
		  
	std::vector<std::string> pars(1,"t");	 
	pars[0] = "mu";

	l.rooContainer->AddGenericPdf("background_model",
	  "","mass",pars,1); // 1 for exonential, no need for formula 
        // -----------------------------------------------------
        // -----------------------------------------------------

	// Make some data sets from the observables to fill in the event loop		  
	// Binning is for histograms (will also produce unbinned data sets)
	l.rooContainer->CreateDataSet("mass","data_mass"    ,50); // (100,110,150) -> for a window, else full obs range is taken 
	l.rooContainer->CreateDataSet("mass","bkg_mass"     ,50);    	  	
	l.rooContainer->CreateDataSet("mass","sig_mass_m105",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m110",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m115",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m120",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m130",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m140",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_rv_m105",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_rv_m110",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_rv_m115",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_rv_m120",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_rv_m130",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_rv_m140",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_wv_m105",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_wv_m110",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_wv_m115",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_wv_m120",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_wv_m130",50);    
	l.rooContainer->CreateDataSet("mass","sig_mass_wv_m140",50);    

	// Make more data sets to represent systematic shitfs , 
	l.rooContainer->MakeSystematics("mass","sig_mass_m105",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m110",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m115",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m120",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m130",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m140",-1);	
	
     /* -----------------------------------------------------------------------------------------
     KFactors Reweighting
     ------------------------------------------------------------------------------------------- */
     cout << "Opening KFactor file"<<endl;
     if (kfacHist != ""){
	cout << "Opening KFactor file"<<endl;
        TFile* kfacFile = TFile::Open( kfacHist );
        if (kfacFile){
	TH1D *temp;
	temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",110));
	thm110 =(TH1D*) temp->Clone(Form("Hmass%d",110));
	temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",120));
	thm120 =(TH1D*) temp->Clone(Form("Hmass%d",120));
	temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",130));
	thm130 =(TH1D*) temp->Clone(Form("Hmass%d",130));
	temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",140));
	thm140 =(TH1D*) temp->Clone(Form("Hmass%d",140));

	}
     } 
	
	if(PADEBUG) 
		cout << "InitRealStatAnalysis END"<<endl;

	// FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
   if(PADEBUG) 
     cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
   
   int cur_type = l.itype[l.current];
   float weight = l.sampleContainer[l.current_sample_index].weight;
    //PU reweighting
   unsigned int n_pu = l.pu_n;
    if (l.itype[l.current] !=0 && puHist != "") {
        if(n_pu<weights.size()){
	    //cout << n_pu<<endl;
	    //std::cout << weights[n_pu] << std::endl;
	    weight *= weights[n_pu];   // put this back
        }    
        else{ //should not happen 
             cout<<"n_pu ("<< n_pu<<") too big ("<<weights.size()
             <<"), event will not be reweighted for pileup"<<endl;
        }
    }
 // ------------------------------------------------------------
   //PT-H K-factors
   double gPT = 0;
   if (cur_type<0){
    for (int gi=0;gi<l.gp_n;gi++){
	if (l.gp_pdgid[gi]==25){
	 TLorentzVector *gP4 = (TLorentzVector*)l.gp_p4->At(gi);
	 gPT = gP4->Pt();
	 break;
	}
    }
   }


   if (cur_type == -13)      weight*=GetDifferentialKfactor(gPT,105);
   else if (cur_type == -17) weight*=GetDifferentialKfactor(gPT,110);
   else if (cur_type == -21) weight*=GetDifferentialKfactor(gPT,115);
   else if (cur_type == -25) weight*=GetDifferentialKfactor(gPT,120);
   else if (cur_type == -29) weight*=GetDifferentialKfactor(gPT,130);
   else if (cur_type == -33) weight*=GetDifferentialKfactor(gPT,140);

   // ------------------------------------------------------------

   // smear all of the photons!
   std::pair<int,int> diphoton_index;
   int leadLevel=2, subLevel=2;

   TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.vtx_std_sel);

   // Nominal smearing
   TClonesArray smeared_pho_p4("TLorentzVector",l.pho_n);
   std::vector<float> smeared_pho_r9(l.pho_n,0.); 
   std::vector<float> smeared_pho_weight(l.pho_n,1.);
   
   for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
       std::vector<std::vector<bool> > p;
       PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), l.pho_isEB[ipho], l.pho_r9[ipho],
				   l.PhotonCiCSelectionLevel(ipho,p,l.phoSUPERTIGHT) );
       // std::cerr << "Smearing photon " << ipho << " " << phoInfo.energy() << " " << phoInfo.r9() << " " << phoInfo.iDet() << std::endl;
       float pweight = 1.;
       // smear MC. But apply energy shift to data 
       if( cur_type != 0 ) { // if it's MC
	 for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
	   float sweight = 1.;
	   // if (   (*si)->name()!=std::string("E_scale")  )     {        (*si)->smearPhoton(phoInfo,sweight,0.);	       } //GF
	   (*si)->smearPhoton(phoInfo,sweight,0.);	   
	   pweight *= sweight;
	 }
       } else 	 {          // if it's data
	 float sweight = 1.;
	 eScaleDataSmearer->smearPhoton(phoInfo,sweight,0.);
	 pweight *= sweight;
       }
       new( smeared_pho_p4[ipho] )  TLorentzVector(phoInfo.p4(vtx->X(), vtx->Y(), vtx->Z() ));
       smeared_pho_r9[ipho] = phoInfo.r9();
       smeared_pho_weight[ipho] = pweight;
   }
   
   // FIXME pass smeared R9
   diphoton_index = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,false, &smeared_pho_p4 ); 
   
   if (diphoton_index.first > -1 && diphoton_index.second > -1){
       float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second];
       
       l.countersred[diPhoCounter_]++;

       TLorentzVector *lead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.first);
       TLorentzVector *sublead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.second);
       TLorentzVector Higgs = *lead_p4 + *sublead_p4; 	
       
       bool CorrectVertex;
       // FIXME pass smeared R9
       int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
       int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
       if( cur_type != 0 ) {
	       float rewei=1.;
	       float pth = Higgs.Pt();
	       for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
		   (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
	       }
	       evweight *= rewei;
	       // WARNING the order of the smeares matter here
	       if( pth != Higgs.Pt() ) {
		   category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
	       }
	       CorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
       }
       float mass    = Higgs.M();
       float ptHiggs = Higgs.Pt();
      

	// Ouput Text File For Saclay Analysis
	if (cur_type == 0)
	  SaclayText << Form("typ=%d weight=%f category=%d mass=%f ptH=%f ptLEAD=%f ptSUBLEAD=%f, run=%d, lumis=%d, event=%d",cur_type,evweight,category,mass,ptHiggs,lead_p4->Pt(),sublead_p4->Pt())<<l.run<< l.lumis<<l.event <<endl;
	else
	  SaclayText << Form("typ=%d weight=%f category=%d mass=%f ptH=%f ptLEAD=%f ptSUBLEAD=%f",cur_type,evweight,category,mass,ptHiggs,lead_p4->Pt(),sublead_p4->Pt())<<endl;

	// --------------------------------------------------------------------------------------------- 
       if (cur_type == 0 ){
	   l.rooContainer->InputDataPoint("data_mass",category,mass);
       }
       if (cur_type > 0 && cur_type != 3 && cur_type != 4)
	   l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
       else if (cur_type == 3 || cur_type == 4)
	   l.rooContainer->InputDataPoint("zee_mass",category,mass,evweight);

       else if (cur_type == -13|| cur_type == -14 || cur_type == -15|| cur_type == -16){
	   l.rooContainer->InputDataPoint("sig_mass_m105",category,mass,evweight);
	   if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m105",category,mass,evweight);
	   else l.rooContainer->InputDataPoint("sig_mass_wv_m105",category,mass,evweight);
       }
       else if (cur_type == -17 || cur_type == -18 || cur_type == -19|| cur_type == -20){
	   l.rooContainer->InputDataPoint("sig_mass_m110",category,mass,evweight);
	   if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m110",category,mass,evweight);
	   else l.rooContainer->InputDataPoint("sig_mass_wv_m110",category,mass,evweight);
       }
       else if (cur_type == -21 || cur_type == -22 || cur_type == -23|| cur_type == -24){
	   l.rooContainer->InputDataPoint("sig_mass_m115",category,mass,evweight);
	   if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m115",category,mass,evweight);
	   else l.rooContainer->InputDataPoint("sig_mass_wv_m115",category,mass,evweight);
       }
       else if (cur_type == -25 || cur_type == -26 || cur_type == -27|| cur_type == -28){
	   l.rooContainer->InputDataPoint("sig_mass_m120",category,mass,evweight);
	   if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m120",category,mass,evweight);
	   else l.rooContainer->InputDataPoint("sig_mass_wv_m120",category,mass,evweight);
       }
       else if (cur_type == -29 || cur_type == -30 || cur_type == -31|| cur_type == -32){
	   l.rooContainer->InputDataPoint("sig_mass_m130",category,mass,evweight);
	   if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m130",category,mass,evweight);
	   else l.rooContainer->InputDataPoint("sig_mass_wv_m130",category,mass,evweight);
       }
       else if (cur_type == -33 || cur_type == -34 || cur_type == -35|| cur_type == -36){
	   l.rooContainer->InputDataPoint("sig_mass_m140",category,mass,evweight);
	   if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m140",category,mass,evweight);
	   else l.rooContainer->InputDataPoint("sig_mass_wv_m140",category,mass,evweight);
       }
       
   }
   
   
   // Systematics
   if( cur_type != 0 ) { 
       // FIXME smearers apply only to MC now
       
       // fill steps for syst uncertainty study
       float systStep = systRange / (float)nSystSteps;
       // di-photon smearers systematics
       if (diphoton_index.first > -1 && diphoton_index.second > -1){
	       
	       TLorentzVector *lead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.first);
	       TLorentzVector *sublead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.second);
	       
	       for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		       std::vector<double> mass_errors;
		       std::vector<double> weights;
		       std::vector<int> categories;
		       
		       for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systRange ) { 
			       if( syst_shift == 0. ) { continue; } // skip the central value
			       TLorentzVector Higgs = *lead_p4 + *sublead_p4; 	
			       
			       // restart with 'fresh' wait for this round of systematics
			       float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second];

			       // FIXME pass smeared R9 and di-photon
			       int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
       			       int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
			       for(std::vector<BaseDiPhotonSmearer *>::iterator sj=diPhotonSmearers_.begin(); sj!= diPhotonSmearers_.end(); ++sj ) {
				       float swei=1.;
				       float pth = Higgs.Pt();
				       if( *si == *sj ) { 
					 (*si)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), syst_shift );
				       } else { 
					   (*sj)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
				       }
				       evweight *= swei;
				       if( pth != Higgs.Pt() ) {
					   category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
				       }
			       }
			       
			       float mass = Higgs.M();
			       categories.push_back(category);
			       mass_errors.push_back(mass);
			       weights.push_back(evweight);

		       }
       			 if (cur_type == -13|| cur_type == -14 || cur_type == -15|| cur_type == -16)
			       l.rooContainer->InputSystematicSet("sig_mass_m105",(*si)->name(),categories,mass_errors,weights);
       			 else if (cur_type == -17|| cur_type == -18 || cur_type == -19|| cur_type == -20)
			       l.rooContainer->InputSystematicSet("sig_mass_m110",(*si)->name(),categories,mass_errors,weights);
       			 else if (cur_type == -21|| cur_type == -22 || cur_type == -23|| cur_type == -24)
			       l.rooContainer->InputSystematicSet("sig_mass_m115",(*si)->name(),categories,mass_errors,weights);
       			 else if (cur_type == -25|| cur_type == -26 || cur_type == -27|| cur_type == -28)
			       l.rooContainer->InputSystematicSet("sig_mass_m120",(*si)->name(),categories,mass_errors,weights);
       			 else if (cur_type == -29|| cur_type == -30 || cur_type == -31|| cur_type == -32)
			       l.rooContainer->InputSystematicSet("sig_mass_m130",(*si)->name(),categories,mass_errors,weights);
       			 else if (cur_type == -33|| cur_type == -34 || cur_type == -35|| cur_type == -36)
			       l.rooContainer->InputSystematicSet("sig_mass_m140",(*si)->name(),categories,mass_errors,weights);
	       }

       }
       
       // loop over the smearers included in the systematics study
       for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	   std::vector<double> mass_errors;
	   std::vector<double> weights;
	   std::vector<int> categories;
	   
	   // loop over syst shift
	   for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systRange ) { 
	       if( syst_shift == 0. ) { continue; } // skip the central value
	       // smear the photons 
	       for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		   std::vector<std::vector<bool> > p;
		   PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), l.pho_isEB[ipho], l.pho_r9[ipho],
					       l.PhotonCiCSelectionLevel(ipho,p,l.phoSUPERTIGHT));
		   float pweight = 1.;
		   for(std::vector<BaseSmearer *>::iterator  sj=photonSmearers_.begin(); sj!= photonSmearers_.end(); ++sj ) {
		       float sweight = 1.;
		       if( *si == *sj ) {
			   // move the SNsmearer under study by syst_shift
			   (*si)->smearPhoton(phoInfo,sweight,syst_shift);
		       } else {
			 // for the other use the nominal points
			   (*sj)->smearPhoton(phoInfo,sweight,0.);
		       }
		       pweight *= sweight;
		   }
		   *((TLorentzVector *)smeared_pho_p4.At(ipho)) = phoInfo.p4(vtx->X(), vtx->Y(), vtx->Z() );
		   smeared_pho_r9[ipho] = phoInfo.r9();
		   smeared_pho_weight[ipho] = pweight;
	       }
	       
	       // analyze the event
	       // FIXME pass smeared R9
	       diphoton_index = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,false, &smeared_pho_p4 ); 
	       
	       if (diphoton_index.first > -1 && diphoton_index.second > -1){
		   
		   float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second];
		   
		   TLorentzVector *lead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.first);
		   TLorentzVector *sublead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.second);
		   TLorentzVector Higgs = *lead_p4 + *sublead_p4; 	
		   
		   int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
       		   int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		   if( cur_type != 0 ) {
		       for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
			   float rewei=1.;
			   float pth = Higgs.Pt();
			   (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			   if( pth != Higgs.Pt() ) {
			       category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
			   }
			   evweight *= rewei;
		       }
		   }
		   float mass = Higgs.M();
		   
	       	    categories.push_back(category);
	            mass_errors.push_back(mass);
	            weights.push_back(evweight);

	       } else {
		   mass_errors.push_back(0.);   
		   weights.push_back(0.);   
		   categories.push_back(-1);
	       }
	       
	   }
	 if (cur_type == -13|| cur_type == -14 || cur_type == -15|| cur_type == -16)
	       l.rooContainer->InputSystematicSet("sig_mass_m105",(*si)->name(),categories,mass_errors,weights);
	 else if (cur_type == -17|| cur_type == -18 || cur_type == -19|| cur_type == -20)
	       l.rooContainer->InputSystematicSet("sig_mass_m110",(*si)->name(),categories,mass_errors,weights);
	 else if (cur_type == -21|| cur_type == -22 || cur_type == -23|| cur_type == -24)
	       l.rooContainer->InputSystematicSet("sig_mass_m115",(*si)->name(),categories,mass_errors,weights);
	 else if (cur_type == -25|| cur_type == -26 || cur_type == -27|| cur_type == -28)
	       l.rooContainer->InputSystematicSet("sig_mass_m120",(*si)->name(),categories,mass_errors,weights);
	 else if (cur_type == -29|| cur_type == -30 || cur_type == -31|| cur_type == -32)
	       l.rooContainer->InputSystematicSet("sig_mass_m130",(*si)->name(),categories,mass_errors,weights);
	 else if (cur_type == -33|| cur_type == -34 || cur_type == -35|| cur_type == -36)
	       l.rooContainer->InputSystematicSet("sig_mass_m140",(*si)->name(),categories,mass_errors,weights);
       
       }
       
       
   }
   smeared_pho_p4.Delete();
   
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
 return true;
}
// ----------------------------------------------------------------------------------------------------
double StatAnalysis::GetDifferentialKfactor(double gPT, int Mass)
{


	if (Mass <=110 ) return thm110->GetBinContent(thm110->FindFixBin(gPT));
	else if (Mass ==120 ) return thm120->GetBinContent(thm120->FindFixBin(gPT));
	else if (Mass ==130 ) return thm130->GetBinContent(thm130->FindFixBin(gPT));
	else if (Mass ==140 ) return thm140->GetBinContent(thm140->FindFixBin(gPT));
	else if (Mass ==115 ) return (0.5*thm110->GetBinContent(thm110->FindFixBin(gPT)) +0.5*thm120->GetBinContent(thm120->FindFixBin(gPT)));
	return 1.0;
/*
	int  genMasses[4] = {110,120,130,140};
	if (Mass<=genMasses[0] ) return kfactorHistograms[0]->GetBinContent(kfactorHistograms[0]->FindBin(gPT));
	else if (Mass<genMasses[nMasses-1]) {

		TH1D *hm1,*hm2;
		double m1=0,m2=0;
		for (int m=0;m<nMasses;m++){
                	if (Mass<genMasses[m+1]){
				hm1=kfactorHistograms[m];
				hm2=kfactorHistograms[m+1];
				m1 = genMasses[m];
				m2 = genMasses[m+1];
			//	cout << "Gen Mass: "<< Mass << " Using "<<m1<< " " << m2<< " Hist name check " << hm1->GetName()<<" " <<hm2->GetName()<<endl;
				break;
			}
		}
		if ((int)Mass == (int)m1 ){
			 //cout << "Found the appropriate historgam "<<hm1->GetName()<<endl;
			 return hm1->GetBinContent(hm1->FindBin(gPT));
		} else {

			TH1D *hm = (TH1D*) hm1->Clone("hm");
			double alpha = ((float) (Mass-m1))/(m2-m1); // make sure ms are not integers
			hm->Add(hm1,hm2,alpha,(1-alpha));
			return hm->GetBinContent(hm->GetBinContent(hm->FindBin(gPT)));
		}

	}
	else return kfactorHistograms[nMasses-1]->GetBinContent(kfactorHistograms[nMasses-1]->FindBin(gPT));
*/
}


