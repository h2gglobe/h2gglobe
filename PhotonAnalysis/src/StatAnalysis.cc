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
{}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Term(LoopAll& l) 
{
        // Make Fits to the data-sets and systematic sets
        l.rooContainer->FitToData("data_exp_model","data_mass");  // Fit to full range of dataset
//        l.rooContainer->FitToData("background_model","bkg_mass");    // Fit to full set

//        l.rooContainer->FitToData("signal_model","sig_mass_m110");		 // fit to full obervable range
//        l.rooContainer->FitToData("signal_model","sig_mass_m120");		 // fit to full obervable range
//        l.rooContainer->FitToData("signal_model","sig_mass_m130");		 // fit to full obervable range
//        l.rooContainer->FitToData("signal_model","sig_mass_m140");		 // fit to full obervable range
  
        // fit to the systematic shits
        // l.rooContainer->FitToSystematicSet("signal_model","sig_mass_m120","E_scale");
  
        // Can create binned models from the results of the fits, should be same bins as other 
        // binned models for when plugging into limit setting.
        //l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","background_model","mass",105); 	   // number of bins only -> full range
       // l.rooContainer->GenerateBinnedPdf("bkg_mass_narrow","background_model","mass",100,0,105,155); // 0, will take range from given range (if no range, will default to full obs range) 
        l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","data_exp_model","data_mass",1,100,0); // 1 means systematics from the fit effect only the background last digit mode = 0 
        // mode 0 as above, 1 if want to bin in sub range from fit,
  	//l.rooContainer->CombineBinnedDatasets("bkg_mass_narrow","zee_mass");

        // Write the data-card for the Combinations Code, needs the output filename, makes binned analysis DataCard
        std::string outputfilename = (std::string) l.histFileName;
        l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m110","bkg_mass_rebinned");
        l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m115","bkg_mass_rebinned");
        l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m120","bkg_mass_rebinned");
        l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m130","bkg_mass_rebinned");
        l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass_m140","bkg_mass_rebinned");

}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Init(LoopAll& l) 
{
	if(PADEBUG) 
		cout << "InitRealStatAnalysis START"<<endl;
	
	// call the base class initializer
	PhotonAnalysis::Init(l);

	// FIXME: move to config file 
	doEscaleSyst = true;
	doEresolSyst = true;
	doPhotonIdEffSyst = false;
	doVtxEffSyst = true;
	nCategories = 4;

	eSmearPars.categoryType = "2CatR9_EBEE";

	// Number's from Adi's presentation 23 / 5
	eSmearPars.n_categories = nCategories; 
	// FIXME: 
	// . shall apply E-scale correction to data?  
	// . double-check signs
	eSmearPars.scale_offset["EBHighR9"] = 4.6e-3;
	eSmearPars.scale_offset["EBLowR9"]  = -1.9e-3;
	eSmearPars.scale_offset["EEHighR9"] = -7.6e-3;
	eSmearPars.scale_offset["EELowR9"]  = -2.1e-3;

	eSmearPars.scale_offset_error["EBHighR9"] = 7e-4;
	eSmearPars.scale_offset_error["EBLowR9"]  = 6e-4;
	eSmearPars.scale_offset_error["EEHighR9"] = 1.8e-3;
	eSmearPars.scale_offset_error["EELowR9"]  = 2.5e-4;

	eSmearPars.smearing_sigma["EBHighR9"] = 9.6e-3;
	eSmearPars.smearing_sigma["EBLowR9"]  = 1.348e-2;
	eSmearPars.smearing_sigma["EEHighR9"] = 2.5e-2;
	eSmearPars.smearing_sigma["EELowR9"]  = 2.3e-2;

	eSmearPars.smearing_sigma_error["EBHighR9"] = 1.e-3;
	eSmearPars.smearing_sigma_error["EBLowR9"]  = 6e-4;
	eSmearPars.smearing_sigma_error["EEHighR9"] = 2e-3;
	eSmearPars.smearing_sigma_error["EELowR9"]  = 2e-3;

	eSmearPars.efficiency_file                = std::string("/afs/cern.ch/user/f/franzoni/public/globe/effReweight_tgraph.root");
	
	eScaleSmearer = new EnergySmearer( eSmearPars );
	eScaleSmearer->name("E_scale");
	eScaleSmearer->doEnergy(true);
	eScaleSmearer->scaleOrSmear(true);
	
	eResolSmearer = new EnergySmearer( eSmearPars );
	eResolSmearer->name("E_res");
	eResolSmearer->doEnergy(false);
	eResolSmearer->scaleOrSmear(false);

	idEffSmearer = new EnergySmearer( eSmearPars );
	idEffSmearer->name("idEff");
	idEffSmearer->doEfficiencies(true);
	idEffSmearer->setEffName("ratioTP");
	idEffSmearer->initEfficiency();

	vtxEffSmearer = new EnergySmearer( eSmearPars );   // triplicate TF1's here
	vtxEffSmearer->name("vtxEff");
	vtxEffSmearer->doEfficiencies(true);
	vtxEffSmearer->setEffName("ratioVertex");
	vtxEffSmearer->initEfficiency();

 
	photonSmearers_.push_back(eScaleSmearer);
	photonSmearers_.push_back(eResolSmearer);
	photonSmearers_.push_back(idEffSmearer);
	photonSmearers_.push_back(vtxEffSmearer);
	

/*
//	r = new TRandom3();
	if (icat == 0) {escale = 0.39; eres = 1.58;} 	// EB high r9
	if (icat == 1) {escale = 0.19; eres = 1.8;} 	// EB low r9
	if (icat == 2) {escale = 1.83; eres = 3.61;} 	// EE high r9
	if (icat == 3) {escale = 2.02; eres = 3.17;} 	// EE low r9

	energySmearerParams = EnergySmearer::energySmearingParameters;
	energySmearerParams.n_categories=4;

	scale_offset["EBHighR9"] = 0.39;
	scale_offset["EBLowR9"]  = 0.19;
	scale_offset["EEHighR9"] = 1.83;
	scale_offset["EELowR9"]  = 2.02;

	scale_offset["EBHighR9"] = 0.39;
	scale_offset["EBLowR9"]  = 0.19;
	scale_offset["EEHighR9"] = 1.83;
	scale_offset["EELowR9"]  = 2.02;

	energySmearer = new EnergySmearer()
*/
        // Define the number of categories for the statistical analysis and
	// the systematic sets to be formed

	// FIXME move these params to config file
	l.rooContainer->SetNCategories(nCategories);
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
	if( doPhotonIdEffSyst ) {//GF this is here as a test, for the moment - there are other efficiencies too (VTX, trigger etc)
		systPhotonSmearers_.push_back( idEffSmearer );
		std::vector<std::string> sys(1,idEffSmearer->name());
		std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
		l.rooContainer->MakeSystematicStudy(sys,sys_t);
	}//GF
	if( doVtxEffSyst ) {//GF this is here as a test, for the moment - there are other efficiencies too (VTX, trigger etc)
		systPhotonSmearers_.push_back( vtxEffSmearer );
		std::vector<std::string> sys(1,vtxEffSmearer->name());
		std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
		l.rooContainer->MakeSystematicStudy(sys,sys_t);
	}//GF
	// ----------------------------------------------------

	// Create observables for shape-analysis with ranges
	l.rooContainer->AddObservable("mass" ,100.,200.);

	// Create parameters and pdfs for signal/background

	// Data shape - Exponential + Zee tail
	l.rooContainer->AddRealVar("mu_data",-0.04,-2.,-0.0001);
	l.rooContainer->AddRealVar("mu_tail_data",-0.02,-0.5,-0.0001);

	std::vector<std::string> data_exp_pars(1,"e");	 
	data_exp_pars[0] = "mu_data";
	std::vector<std::string> data_tail_pars(1,"e");	 
	data_tail_pars[0] = "mu_tail_data";

	l.rooContainer->AddGenericPdf("data_exp_model",
	  "","mass",data_exp_pars,1);//,0.8,0.5,0.99); // 1 for exonential, no need for formula
 
	l.rooContainer->AddGenericPdf("data_tail_model",
	  "","mass",data_tail_pars,1,0.2,0.,0.5); // 1 for exonential, no need for formula 

        std::vector<std::string> components_data(2,"c");
        components_data[0] = "data_exp_model";
        components_data[1] = "data_tail_model";
        l.rooContainer->ComposePdf("data_model","data_exp_model+data_tail_model",components_data,false); // true means use extended pdfs
	// ------------------------------------------------------
        
        // Background - Exponential
	l.rooContainer->AddRealVar("mu",-0.04,-2.,-0.001);
		  
	std::vector<std::string> pars(1,"t");	 
	pars[0] = "mu";

	l.rooContainer->AddGenericPdf("background_model",
	  "","mass",pars,1); // 1 for exonential, no need for formula 
        // -----------------------------------------------------
	
	// Signal -- CB shape
        l.rooContainer->AddRealVar("mean",120,100,150);
        l.rooContainer->AddRealVar("sigma",1.,0.1,2.);
        l.rooContainer->AddRealVar("slope",1.,0.5,5.);
        l.rooContainer->AddRealVar("n",2.,0.5,10.);

	std::vector<std::string> cb_pars(4,"cb");
	cb_pars[0] = "mean";
	cb_pars[1] = "sigma";
	cb_pars[2] = "slope";
	cb_pars[3] = "n";

	l.rooContainer->AddGenericPdf("cb_shape",
	  "","mass",cb_pars,4,0.8,0.,0.99);  // 4 for CB shape, no need for formula

	// Signal -- Gaussian
        l.rooContainer->AddRealVar("gsigma",2.5,1.,5.);

	std::vector<std::string> gau_pars(2,"gau");
	gau_pars[0] = "mean";
	gau_pars[1] = "gsigma";

	l.rooContainer->AddGenericPdf("gau_shape",
	  "","mass",gau_pars,2,0.1,0.,0.5);  // 3 for Gaussian shape, no need for formula. 

        // Add the  CB+Gaussian
        std::vector<std::string> components(2,"c");
        components[0] = "cb_shape";
        components[1] = "gau_shape";
        l.rooContainer->ComposePdf("signal_model","cb_shape+gau_shape",components,false); // true means use extended pdfs
											  // eg for a signal+bkg model
        // -----------------------------------------------------

	// Make some data sets from the observables to fill in the event loop		  
	// Binning is for histograms (will also produce unbinned data sets)
	l.rooContainer->CreateDataSet("mass","data_mass"    ,100); // (100,110,150) -> for a window, else full obs range is taken 
	l.rooContainer->CreateDataSet("mass","bkg_mass"     ,100);    	  	
	l.rooContainer->CreateDataSet("mass","sig_mass_m110",100);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m115",100);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m120",100);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m130",100);    
	l.rooContainer->CreateDataSet("mass","sig_mass_m140",100);    

	// Make more data sets to represent systematic shitfs , 
	l.rooContainer->MakeSystematics("mass","sig_mass_m110",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m115",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m120",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m130",-1);	
	l.rooContainer->MakeSystematics("mass","sig_mass_m140",-1);	
	
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
    //PU reweighting
   unsigned int n_pu = l.pu_n;
    if (l.itype[l.current] !=0 && puHist != "") {
        if(n_pu<weights.size()){
	     //cout << n_pu<<endl;
	     //std::cout << weights[n_pu] << std::endl;
             //weight *= weights[n_pu];
        }    
        else{ //should not happen 
             cout<<"n_pu ("<< n_pu<<") too big ("<<weights.size()
             <<"), event will not be reweighted for pileup"<<endl;
        }
    }

   // smear all of the photons!

/*
   for (int ipho=0 ;ipho<l.pho_n;ipho++){
	
	double escale;
	double eres;

	int icat = l.PhotonCategory(ipho,2,2);
	if (icat == 0) {escale = 0.39; eres = 1.58;} 	// EB high r9
	if (icat == 1) {escale = 0.19; eres = 1.8;} 	// EB low r9
	if (icat == 2) {escale = 1.83; eres = 3.61;} 	// EE high r9
	if (icat == 3) {escale = 2.02; eres = 3.17;} 	// EE low r9

        TLorentzVector p4 = *((TLorentzVector*)l.pho_p4->At(ipho));
	double newE = p4.E()+r->Gaus(escale,eres);
	p4.SetE(newE);
	cout << newE<<endl;;
	*((TLorentzVector*)l.pho_p4->At(ipho))= p4;	

   }
*/
   std::pair<int,int> diphoton_index;
   int leadLevel=2, subLevel=2;

   float leadEtCut   = 40.0;
   float subleadEtCut = 30.0;

   TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.vtx_std_sel);

   // Nominal smearing
   TClonesArray smeared_pho_p4("TLorentzVector",l.pho_n);
   std::vector<float> smeared_pho_r9(l.pho_n,0.); 
   std::vector<float> smeared_pho_weight(l.pho_n,1.);
   
   for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
       PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), l.pho_isEB[ipho], l.pho_r9[ipho] );
       // std::cerr << "Smearing photon " << ipho << " " << phoInfo.energy() << " " << phoInfo.r9() << " " << phoInfo.iDet() << std::endl;
       float pweight = 1.;
       // smear only MC. 
       // FIXME: shall apply E-scale correction to data?  
       if( cur_type != 0 ) {
	   for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
	       float sweight = 1.;
	       (*si)->smearPhoton(phoInfo,weight,0.);
	       pweight *= sweight;
	   }
       }
       
       // std::cerr << "Smeared photon " << ipho << " " << phoInfo.energy() << " " << phoInfo.r9() << " " << phoInfo.iDet() << std::endl;
       new( smeared_pho_p4[ipho] )  TLorentzVector(phoInfo.p4(vtx->X(), vtx->Y(), vtx->Z() ));
       smeared_pho_r9[ipho] = phoInfo.r9();
       smeared_pho_weight[ipho] = pweight;
   }
   
   // FIXME pass smeared R9
   diphoton_index = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,false, &smeared_pho_p4 ); 
   
   if (diphoton_index.first > -1 && diphoton_index.second > -1){
       
       float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second];
       
       TLorentzVector *lead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.first);
       TLorentzVector *sublead_p4 = (TLorentzVector*)smeared_pho_p4.At(diphoton_index.second);
       TLorentzVector Higgs = *lead_p4 + *sublead_p4; 	
	   
       // FIXME add di-photon smearings (and systematics?)
       float mass = Higgs.M();
       
       // 4 categories used for the limit setting.... 2 r9, 2 Eta
       // FIXME pass smeared R9
       int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,2,2,0);
       
       if (cur_type == 0 ){
	   l.rooContainer->InputDataPoint("data_mass",category,mass,evweight);
       }
       if (cur_type > 0 && cur_type != 3 && cur_type != 4)
	   l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
       else if (cur_type == 3 || cur_type == 4)
	   l.rooContainer->InputDataPoint("zee_mass",category,mass,evweight);
       //else if (cur_type == -1 || cur_type == -2 || cur_type == -3)
       //  l.rooContainer->InputDataPoint("sig_mass_m100",category,mass,evweight);
       else if (cur_type == -1 || cur_type == -2 || cur_type == -3)
	   l.rooContainer->InputDataPoint("sig_mass_m110",category,mass,evweight);
       else if (cur_type == -4 || cur_type == -5 || cur_type == -6)
	   l.rooContainer->InputDataPoint("sig_mass_m115",category,mass,evweight);
       else if (cur_type == -7 || cur_type == -8 || cur_type == -9)
	   l.rooContainer->InputDataPoint("sig_mass_m120",category,mass,evweight);
       else if (cur_type == -10 || cur_type == -11 || cur_type == -12)
	   l.rooContainer->InputDataPoint("sig_mass_m130",category,mass,evweight);
       else if (cur_type == -13 || cur_type == -14 || cur_type == -15)
	   l.rooContainer->InputDataPoint("sig_mass_m140",category,mass,evweight);
       
   }
   
   
   // Systematics
   
   if( cur_type != 0 ) {
       // FIXME smearers apply only to MC now
       
       
       // fill steps for syst uncertainty study
       float systStep = systRange / (float)nSystSteps;
       
       // loop over the smearers included in the systematics study
       for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	   std::vector<double> mass_errors;
	   std::vector<double> weights;
	   std::vector<int> categories;
	   
	   // loop over syst shift
	   for(float syst_shift=(-1*systRange); syst_shift<=systRange; syst_shift+=systRange ) { 
	       // skip the central value
	       if( syst_shift == 0. ) { continue; } 
	       
	       // smear the photons 
	       for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		   
		   PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), l.pho_isEB[ipho], l.pho_r9[ipho] );
		   float pweight = 1.;
		   
		   // move the smearer under study by syst_shift
		   (*si)->smearPhoton(phoInfo,weight,syst_shift);
		       
		   // for the other use the nominal points
		   for(std::vector<BaseSmearer *>::iterator  sj=photonSmearers_.begin(); sj!= photonSmearers_.end(); ++sj ) {
		       if( *si == *sj ) { continue; }
		       float sweight = 1.;
		       (*sj)->smearPhoton(phoInfo,weight,0.);
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
		   
		   // FIXME di-photon smearings
		   float mass = Higgs.M();
		   
		   // 4 categories used for the limit setting.... 2 r9, 2 Eta
		   // FIXME pass smeared R9
		   int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,2,2,0);
		   
		   categories.push_back(category);
		   mass_errors.push_back(mass);
		   weights.push_back(evweight);
		   
	       } else {
		   mass_errors.push_back(0.);   
		   weights.push_back(0.);   
		   categories.push_back(-1);
	       }
	       
	   }
       	   if (cur_type == -1 || cur_type == -2 || cur_type == -3)
	     l.rooContainer->InputSystematicSet("sig_mass_m110",(*si)->name(),categories,mass_errors,weights);
           else if (cur_type == -4 || cur_type == -5 || cur_type == -6)
	     l.rooContainer->InputSystematicSet("sig_mass_m115",(*si)->name(),categories,mass_errors,weights);
           else if (cur_type == -7 || cur_type == -8 || cur_type == -9)
	     l.rooContainer->InputSystematicSet("sig_mass_m120",(*si)->name(),categories,mass_errors,weights);
           else if (cur_type == -10 || cur_type == -11 || cur_type == -12)
	     l.rooContainer->InputSystematicSet("sig_mass_m130",(*si)->name(),categories,mass_errors,weights);
           else if (cur_type == -13 || cur_type == -14 || cur_type == -15)
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
