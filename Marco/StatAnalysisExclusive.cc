#include "../interface/StatAnalysisExclusive.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
StatAnalysisExclusive::StatAnalysisExclusive()  : 
    name_("StatAnalysisExclusive"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
}

// ----------------------------------------------------------------------------------------------------
StatAnalysisExclusive::~StatAnalysisExclusive() 
{
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::Term(LoopAll& l) 
{

    std::string outputfilename = (std::string) l.histFileName;
    // Make Fits to the data-sets and systematic sets
    l.rooContainer->FitToData("data_pol_model","data_mass");  // Fit to full range of dataset
  
//    l.rooContainer->WriteSpecificCategoryDataCards(outputfilename,"data_mass","sig_mass","data_pol_model");
//    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","data_pol_model");
    // mode 0 as above, 1 if want to bin in sub range from fit,

    // Write the data-card for the Combinations Code, needs the output filename, makes binned analysis DataCard
    // Assumes the signal datasets will be called signal_name+"_mXXX"
//    l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","data_pol_model","data_mass",1,50,1); // 1 means systematics from the fit effect only the backgroundi. last digit mode = 1 means this is an internal constraint fit 
//    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","bkg_mass_rebinned");

    eventListText.close();

    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;

//	kfacFile->Close();
//	PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::Init(LoopAll& l) 
{
    if(PADEBUG) 
	cout << "InitRealStatAnalysisExclusive START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;
    
    std::string outputfilename = (std::string) l.histFileName;
    eventListText.open(Form("%s_ascii_events.txt",outputfilename.c_str()));
    FillSignalLabelMap();
    //
    // These parameters are set in the configuration file
    std::cout
	<< "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "StatAnalysisExclusive " << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "leadEtCut "<< leadEtCut << "\n"
	<< "subleadEtCut "<< subleadEtCut << "\n"
	<< "doTriggerSelection "<< doTriggerSelection << "\n"
	<< "nEtaCategories "<< nEtaCategories << "\n"
	<< "nR9Categories "<< nR9Categories << "\n"		
	<< "nPtCategories "<< nPtCategories << "\n"		
	<< "doEscaleSyst "<< doEscaleSyst << "\n"
	<< "doEresolSyst "<< doEresolSyst << "\n"
	<< "doEcorrectionSyst "<< doEcorrectionSyst << "\n"
	<< "efficiencyFile " << efficiencyFile << "\n"
	<< "doPhotonIdEffSyst "<< doPhotonIdEffSyst << "\n"
	<< "doR9Syst "<< doR9Syst << "\n"
	<< "doVtxEffSyst "<< doVtxEffSyst << "\n"
	<< "doTriggerEffSyst "<< doTriggerEffSyst << "\n"
	<< "doKFactorSyst "<< doKFactorSyst << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< std::endl;

    // avoid recalculated the CIC ID every time
    l.runCiC = reRunCiC;
    // call the base class initializer
    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
	l.histoContainer[ind].setScale(1.);
    }
    
    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nCategories_ *= nR9Categories;
    if( nPtCategories != 0 ) nCategories_ *= nPtCategories;

    // CP

    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    
    int nVBFCategories  = ((int)includeVBF)*nVBFEtaCategories;
    int nVHadCategories = (int)includeVHad;
   
    nCategories_+=(nVBFCategories+nVHadCategories);

    l.SetCutVariables("cut_AllLeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("cut_AllSubJPt",        &myAllSubJPt);
    l.SetCutVariables("cut_AllLeadJEta",      &myAllLeadJEta);
    l.SetCutVariables("cut_AllSubJEta",       &myAllSubJEta);
    l.SetCutVariables("cut_All_Mjj",          &myAll_Mjj);
    l.SetCutVariables("cut_All_dEta",         &myAlldEta);
    l.SetCutVariables("cut_All_Zep",          &myAllZep);
    l.SetCutVariables("cut_All_dPhi",         &myAlldPhi);
    l.SetCutVariables("cut_All_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("cut_All_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("cut_All_Mgg_100-160",  &myAll_Mgg);
    l.SetCutVariables("cut_All_Mgg2_100-160", &myAll_Mgg);
    l.SetCutVariables("cut_AllPtHiggs",       &myAllPtHiggs);

    l.SetCutVariables("cut_VBFLeadJPt",       &myVBFLeadJPt);
    l.SetCutVariables("cut_VBFSubJPt",        &myVBFSubJPt);
    l.SetCutVariables("cut_VBF_Mjj",          &myVBF_Mjj);
    l.SetCutVariables("cut_VBF_dEta",         &myVBFdEta);
    l.SetCutVariables("cut_VBF_Zep",          &myVBFZep);
    l.SetCutVariables("cut_VBF_dPhi",         &myVBFdPhi);
    l.SetCutVariables("cut_VBF_Mgg0",         &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg2",         &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg4",         &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg10",        &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg4_100-160",        &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg2_100-160",        &myVBF_Mgg);
    
    l.SetCutVariables("cut_VHadLeadJPt",      &myVHadLeadJPt);
    l.SetCutVariables("cut_VHadSubJPt",       &myVHadSubJPt);
    l.SetCutVariables("cut_VHad_Mjj",         &myVHad_Mjj);
    l.SetCutVariables("cut_VHad_dEta",        &myVHaddEta);
    l.SetCutVariables("cut_VHad_Zep",         &myVHadZep);
    l.SetCutVariables("cut_VHad_dPhi",        &myVHaddPhi);
    l.SetCutVariables("cut_VHad_Mgg0",        &myVHad_Mgg);
    l.SetCutVariables("cut_VHad_Mgg2",        &myVHad_Mgg);
    l.SetCutVariables("cut_VHad_Mgg4",        &myVHad_Mgg);
    l.SetCutVariables("cut_VHad_Mgg10",        &myVHad_Mgg);
    l.SetCutVariables("cut_VHad_Mgg2_100-160",        &myVHad_Mgg);
    l.SetCutVariables("cut_VHad_Mgg4_100-160",        &myVHad_Mgg);

    // CP


    effSmearPars.categoryType = "2CatR9_EBEE";
    effSmearPars.n_categories = 4;
    effSmearPars.efficiency_file = efficiencyFile;

    diPhoEffSmearPars.n_categories = 8;
    diPhoEffSmearPars.efficiency_file = efficiencyFile;

    if( doEcorrectionSmear ) {
        // instance of this smearer done in PhotonAnalysis
        photonSmearers_.push_back(eCorrSmearer);
    }
    if( doEscaleSmear ) {
        photonSmearers_.push_back(eScaleSmearer);
    }
    if( doEresolSmear ) {
	// energy resolution smearing
	std::cerr << __LINE__ << std::endl; 
	eResolSmearer = new EnergySmearer( eSmearPars );
	eResolSmearer->name("E_res");
	eResolSmearer->doEnergy(false);
	eResolSmearer->scaleOrSmear(false);
	photonSmearers_.push_back(eResolSmearer);
    }
    if( doPhotonIdEffSmear ) {
	// photon ID efficiency 
	std::cerr << __LINE__ << std::endl; 
	idEffSmearer = new EfficiencySmearer( effSmearPars );
	idEffSmearer->name("idEff");
	idEffSmearer->setEffName("ratioTP");
	idEffSmearer->init();
	idEffSmearer->doPhoId(true);
	photonSmearers_.push_back(idEffSmearer);
    }
    if( doR9Smear ) {
	// R9 re-weighting
	r9Smearer = new EfficiencySmearer( effSmearPars );
	r9Smearer->name("r9Eff");
	r9Smearer->setEffName("ratioR9");
	r9Smearer->init();
	r9Smearer->doR9(true);
	photonSmearers_.push_back(r9Smearer);
    }
    if( doVtxEffSmear ) {
	// Vertex ID
	std::cerr << __LINE__ << std::endl; 
	vtxEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );   // triplicate TF1's here
	vtxEffSmearer->name("vtxEff");
	vtxEffSmearer->setEffName("ratioVertex");
	vtxEffSmearer->doVtxEff(true);
	vtxEffSmearer->init();
	diPhotonSmearers_.push_back(vtxEffSmearer);
    }
    if( doTriggerEffSmear ) {
	// trigger efficiency
	std::cerr << __LINE__ << std::endl; 
	triggerEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
	triggerEffSmearer->name("triggerEff");
	triggerEffSmearer->setEffName("effL1HLT");
	triggerEffSmearer->doVtxEff(false);
	triggerEffSmearer->init();
	diPhotonSmearers_.push_back(triggerEffSmearer);
    }
    if(doKFactorSmear) {
	// kFactor efficiency
	std::cerr << __LINE__ << std::endl; 
	kFactorSmearer = new KFactorSmearer( kfacHist );
	kFactorSmearer->name("kFactor");
	kFactorSmearer->init();
	genLevelSmearers_.push_back(kFactorSmearer);
    }

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed

    // FIXME move these params to config file
    l.rooContainer->SetNCategories(nCategories_);
    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;
    // RooContainer does not support steps different from 1 sigma
    //assert( ((float)nSystSteps) == systRange );
    if( doEcorrectionSmear && doEcorrectionSyst ) {
        // instance of this smearer done in PhotonAnalysis
        systPhotonSmearers_.push_back(eCorrSmearer);
	std::vector<std::string> sys(1,eCorrSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEscaleSmear && doEscaleSyst ) {
	systPhotonSmearers_.push_back( eScaleSmearer );
	std::vector<std::string> sys(1,eScaleSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
	systPhotonSmearers_.push_back( eResolSmearer );
	std::vector<std::string> sys(1,eResolSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonIdEffSmear && doPhotonIdEffSyst ) {
	systPhotonSmearers_.push_back( idEffSmearer );
	std::vector<std::string> sys(1,idEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doR9Smear && doR9Syst ) {
	systPhotonSmearers_.push_back( r9Smearer );
	std::vector<std::string> sys(1,r9Smearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doVtxEffSmear && doVtxEffSyst ) {
	systDiPhotonSmearers_.push_back( vtxEffSmearer );
	std::vector<std::string> sys(1,vtxEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doTriggerEffSmear && doTriggerEffSyst ) {
	systDiPhotonSmearers_.push_back( triggerEffSmearer );
	std::vector<std::string> sys(1,triggerEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if(doKFactorSmear && doKFactorSyst) {
	systGenLevelSmearers_.push_back(kFactorSmearer);
	std::vector<std::string> sys(1,kFactorSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
	
    // ----------------------------------------------------
    // ----------------------------------------------------
    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.06,1.00);
    // ----------------------------------------------------

    // Create observables for shape-analysis with ranges
    // l.rooContainer->AddObservable("mass" ,100.,150.);
    l.rooContainer->AddObservable("CMS_hgg_mass" ,massMin,massMax);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);

    // SM Model
    l.rooContainer->AddConstant("XSBR_ggh_150",0.01428);
    l.rooContainer->AddConstant("XSBR_vbf_150",0.001308);
    l.rooContainer->AddConstant("XSBR_wzh_150",0.000641);
    l.rooContainer->AddConstant("XSBR_tth_150",0.000066);
    l.rooContainer->AddConstant("XSBR_ggh_145",0.018820);
    l.rooContainer->AddConstant("XSBR_vbf_145",0.001676);
    l.rooContainer->AddConstant("XSBR_wzh_145",0.000891);
    l.rooContainer->AddConstant("XSBR_tth_145",0.000090);
    l.rooContainer->AddConstant("XSBR_ggh_140",0.0234109);
    l.rooContainer->AddConstant("XSBR_vbf_140",0.00203036);
    l.rooContainer->AddConstant("XSBR_wzh_140",0.001163597);
    l.rooContainer->AddConstant("XSBR_tth_140",0.000117189);
    l.rooContainer->AddConstant("XSBR_ggh_135",0.0278604);
    l.rooContainer->AddConstant("XSBR_vbf_135",0.002343);
    l.rooContainer->AddConstant("XSBR_wzh_135",0.001457559);
    l.rooContainer->AddConstant("XSBR_tth_135",0.000145053);
    l.rooContainer->AddConstant("XSBR_ggh_130",0.0319112);
    l.rooContainer->AddConstant("XSBR_vbf_130",0.00260804);
    l.rooContainer->AddConstant("XSBR_wzh_130",0.001759636);
    l.rooContainer->AddConstant("XSBR_tth_130",0.000173070);
    l.rooContainer->AddConstant("XSBR_ggh_125",0.0350599);
    l.rooContainer->AddConstant("XSBR_vbf_125",0.00277319);
    l.rooContainer->AddConstant("XSBR_wzh_125",0.002035123);
    l.rooContainer->AddConstant("XSBR_tth_125",0.000197718);
    l.rooContainer->AddConstant("XSBR_ggh_120",0.0374175);
    l.rooContainer->AddConstant("XSBR_vbf_120",0.00285525);
    l.rooContainer->AddConstant("XSBR_wzh_120",0.002285775);
    l.rooContainer->AddConstant("XSBR_tth_120",0.00021951);
    l.rooContainer->AddConstant("XSBR_ggh_123",0.0360696);
    l.rooContainer->AddConstant("XSBR_vbf_123",0.00281352);
    l.rooContainer->AddConstant("XSBR_wzh_123",0.00213681);
    l.rooContainer->AddConstant("XSBR_tth_123",0.00020663);
    l.rooContainer->AddConstant("XSBR_ggh_121",0.0369736);
    l.rooContainer->AddConstant("XSBR_vbf_121",0.00284082);
    l.rooContainer->AddConstant("XSBR_wzh_121",0.00223491);
    l.rooContainer->AddConstant("XSBR_tth_121",0.00021510);
    l.rooContainer->AddConstant("XSBR_ggh_115",0.0386169);
    l.rooContainer->AddConstant("XSBR_vbf_115",0.00283716);
    l.rooContainer->AddConstant("XSBR_wzh_115",0.002482089);
    l.rooContainer->AddConstant("XSBR_tth_115",0.000235578);
    l.rooContainer->AddConstant("XSBR_ggh_110",0.0390848);
    l.rooContainer->AddConstant("XSBR_vbf_110",0.00275406);
    l.rooContainer->AddConstant("XSBR_wzh_110",0.002654575);
    l.rooContainer->AddConstant("XSBR_tth_110",0.000247629);
    l.rooContainer->AddConstant("XSBR_ggh_105",0.0387684);
    l.rooContainer->AddConstant("XSBR_vbf_105",0.00262016);
    l.rooContainer->AddConstant("XSBR_wzh_105",0.002781962);
    l.rooContainer->AddConstant("XSBR_tth_105",0.000255074);

    // FF model	
    l.rooContainer->AddConstant("ff_XSBR_vbf_150",0.00259659);
    l.rooContainer->AddConstant("ff_XSBR_wzh_150",0.00127278);
    l.rooContainer->AddConstant("ff_XSBR_vbf_145",0.00387544);
    l.rooContainer->AddConstant("ff_XSBR_wzh_145",0.00205969);
    l.rooContainer->AddConstant("ff_XSBR_vbf_140",0.00565976);
    l.rooContainer->AddConstant("ff_XSBR_wzh_140",0.003243602);
    l.rooContainer->AddConstant("ff_XSBR_vbf_135",0.00825);
    l.rooContainer->AddConstant("ff_XSBR_wzh_135",0.00513225);
    l.rooContainer->AddConstant("ff_XSBR_vbf_130",0.0122324);
    l.rooContainer->AddConstant("ff_XSBR_wzh_130",0.00825316);
    l.rooContainer->AddConstant("ff_XSBR_vbf_125",0.0186494);
    l.rooContainer->AddConstant("ff_XSBR_wzh_125",0.01368598);
    l.rooContainer->AddConstant("ff_XSBR_vbf_123",0.022212);
    l.rooContainer->AddConstant("ff_XSBR_wzh_123",0.0168696);
    l.rooContainer->AddConstant("ff_XSBR_vbf_121",0.0266484);
    l.rooContainer->AddConstant("ff_XSBR_wzh_121",0.0209646);
    l.rooContainer->AddConstant("ff_XSBR_vbf_120",0.0293139);
    l.rooContainer->AddConstant("ff_XSBR_wzh_120",0.02346729);
    l.rooContainer->AddConstant("ff_XSBR_vbf_115",0.0482184);
    l.rooContainer->AddConstant("ff_XSBR_wzh_115",0.04218386);
    l.rooContainer->AddConstant("ff_XSBR_vbf_110",0.083181);
    l.rooContainer->AddConstant("ff_XSBR_wzh_110",0.08017625);
    l.rooContainer->AddConstant("ff_XSBR_vbf_105",0.151616);
    l.rooContainer->AddConstant("ff_XSBR_wzh_105",0.1609787);

    // Background modeling 
    l.rooContainer->AddRealVar("pol0",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol1",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol2",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modpol0","@0*@0","pol0");
    l.rooContainer->AddFormulaVar("modpol1","@0*@0","pol1");
    l.rooContainer->AddFormulaVar("modpol2","@0*@0","pol2");

    l.rooContainer->AddRealVar("expol0",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("expol1",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modexpol0","@0*@0","expol0");
    l.rooContainer->AddFormulaVar("modexpol1","@0*@0","expol1");
    
    // Generic PDF ok in the std analysis but excluisve channels need different models CP
    //l.rooContainer->AddGenericPdf("data_pol_model",
	  //"0","CMS_hgg_mass",data_pol_pars,73);	// >= 71 means RooBernstein of order >= 1
        

    int cats_with_std[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int cats_with_excl[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    for(int i=0; i<nEtaCategories*nR9Categories*nPtCategories; i++){
      cats_with_std[i]=1;
      cats_with_excl[i]=0;
    }


    std::vector<std::string> data_pol_pars(3,"p");	 
    data_pol_pars[0] = "modpol0";
    data_pol_pars[1] = "modpol1";
    data_pol_pars[2] = "modpol2";
    l.rooContainer->AddSpecificCategoryPdf(cats_with_std,"data_pol_model",
	  "0","CMS_hgg_mass",data_pol_pars,73);	// >= 71 means RooBernstein of order >= 1
    

    std::vector<std::string> data_excl_pars(2,"p");	 
    data_excl_pars[0] = "modexpol0";
    data_excl_pars[1] = "modexpol1";
    l.rooContainer->AddSpecificCategoryPdf(cats_with_excl, "data_pol_model",
	  "0","CMS_hgg_mass",data_excl_pars,72);	// >= 71 means RooBernstein of order >= 1

    // CP

    // -----------------------------------------------------
    // Make some data sets from the observables to fill in the event loop		  
    // Binning is for histograms (will also produce unbinned data sets)

    l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass"    ,nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass"     ,nDataBins);    	  	

    // Create Signal DataSets:
    for (int sig=105;sig<=150;sig+=5){
	// Needed to use S4 for the GGH 145 Signal which has the BUG so no 145 sample
	if (sig==145) continue;
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_ggh_mass_m%d",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_vbf_mass_m%d",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_wzh_mass_m%d",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_tth_mass_m%d",sig),nDataBins);   
 
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_ggh_mass_m%d_rv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_vbf_mass_m%d_rv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_wzh_mass_m%d_rv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_tth_mass_m%d_rv",sig),nDataBins);    
                                                                              
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_ggh_mass_m%d_wv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_vbf_mass_m%d_wv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_wzh_mass_m%d_wv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_tth_mass_m%d_wv",sig),nDataBins);    
    }

    // Also create the 121 and 123 test points
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m121",nDataBins);   
 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m121_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m121_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m121_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m121_rv",nDataBins);    
                                                                      
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m121_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m121_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m121_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m121_wv",nDataBins);    
	
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m123",nDataBins);   
 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m123_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m123_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m123_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m123_rv",nDataBins);    
                                                                      
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m123_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m123_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m123_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m123_wv",nDataBins);    


    // Make more datasets representing Systematic Shifts of various quantities

    for (int sig=105;sig<=150;sig+=5){
	if (sig==145) continue;
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_ggh_mass_m%d",sig),-1);	
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_vbf_mass_m%d",sig),-1);	
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_wzh_mass_m%d",sig),-1);	
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_tth_mass_m%d",sig),-1);	
    }

    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_ggh_mass_m121",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_vbf_mass_m121",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_wzh_mass_m121",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_tth_mass_m121",-1);
	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_ggh_mass_m123",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_vbf_mass_m123",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_wzh_mass_m123",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_tth_mass_m123",-1);	

    // Make sure the Map is filled
    FillSignalLabelMap();

    if(PADEBUG) 
	cout << "InitRealStatAnalysisExclusive END"<<endl;
	
    // FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::Analysis(LoopAll& l, Int_t jentry) 
{
    if(PADEBUG) 
	cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
   
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;
    l.FillCounter( "Processed", 1. );
    assert( weight > 0. );  
    l.FillCounter( "XSWeighted", weight );
    nevents+=1.;

    //PU reweighting
    unsigned int n_pu = l.pu_n;
    if ( cur_type !=0 && puHist != "") {
	bool hasSpecificWeight = weights.find( cur_type ) != weights.end() ; 
	if( cur_type < 0 && !hasSpecificWeight && jentry == 1 ) {
	    std::cerr  << "WARNING no pu weights specific for sample " << cur_type << std::endl;
	}
	std::vector<double> & puweights = hasSpecificWeight ? weights[ cur_type ] : weights[0]; 
	if(n_pu<puweights.size()){
	    weight *= puweights[n_pu];
	    sumwei+=puweights[n_pu]; 
	}    
	else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	    cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< l.itype[l.current]<<"], event will not be reweighted for pileup"<<endl;
	}
    }
    
    assert( weight >= 0. );  
    l.FillCounter( "PUWeighted", weight );
    
    if( jentry % 10000 ==  0 ) {
	    std::cout << " nevents " <<  nevents << " sumpuweights " << sumwei << " ratio " << sumwei / nevents 
		      << " equiv events " << sumev << " accepted " << sumaccept << " smeared " << sumsmear << " "  
		      <<  sumaccept / sumev << " " << sumsmear / sumaccept
		      << std::endl;
    }
    // ------------------------------------------------------------
    //PT-H K-factors
    double gPT = 0;
    TLorentzVector gP4(0,0,0,0);
    if (cur_type<0){            // if background sample, gP4 remains 4vect(0)
	for (int gi=0;gi<l.gp_n;gi++){
	    if (l.gp_pdgid[gi]==25){
		gP4 = *((TLorentzVector*)l.gp_p4->At(gi));
		gPT = gP4.Pt();
		break;
	    }
	}
    }

    // ------------------------------------------------------------

    // smear all of the photons!
    std::pair<int,int> diphoton_index;
   
    // do gen-level dependent first (e.g. k-factor); only for signal
    double genLevWeight=1; 
    if(cur_type!=0){
	for(std::vector<BaseGenLevelSmearer*>::iterator si=genLevelSmearers_.begin(); si!=genLevelSmearers_.end(); si++){
	    float genWeight=1;
	    (*si)->smearEvent( genWeight,gP4, l.pu_n, cur_type, 0. );
	    if( genWeight < 0. ) {
		std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		assert(0);
	    }
	    genLevWeight*=genWeight;
	}
    }

    // Nominal smearing
    std::vector<float> smeared_pho_energy(l.pho_n,0.); 
    std::vector<float> smeared_pho_r9(l.pho_n,0.); 
    std::vector<float> smeared_pho_weight(l.pho_n,1.);
   
    for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
	std::vector<std::vector<bool> > p;
	PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
				    // *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
				    ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
				    energyCorrected[ipho],
				    l.pho_isEB[ipho], l.pho_r9[ipho],
				    l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_),
				    (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
				    );
	float pweight = 1.;
	// smear MC. But apply energy shift to data 
	if( cur_type != 0 && doMCSmearing ) { // if it's MC
	    for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
		float sweight = 1.;
		(*si)->smearPhoton(phoInfo,sweight,l.run,0.);	   
		if( sweight < 0. ) {
			std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
			assert(0);
		}
		pweight *= sweight;
	    }
	} else if( cur_type == 0 ) {          // if it's data
	    float sweight = 1.;
	    if( doEcorrectionSmear )  { 
	      eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.); 
	    }
	    eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
	    pweight *= sweight;
	}
	smeared_pho_energy[ipho] = phoInfo.energy();
	smeared_pho_r9[ipho] = phoInfo.r9();
	smeared_pho_weight[ipho] = pweight;
    }
   
    sumev += weight;
    // FIXME pass smeared R9
    int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
    // CP
    int diphotonVBF_id = -1;
    int diphotonVHad_id = -1;
    bool VBFevent = false;
    bool VHadevent = false;
    std::pair<int,int> highestPtJets(-1,-1);
    if((includeVBF || includeVHad)&&l.jet_algoPF1_n>1) {
      RescaleJetEnergy(l);
      if(includeVBF) diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
      if(includeVHad) diphotonVHad_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHadCut, subleadEtVHadCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 

      if(diphotonVBF_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );

        bool VBFpresel = (highestPtJets.first!=-1)&&(highestPtJets.first!=-2);
  
        if(VBFpresel){
          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
	        
          TLorentzVector diphoton = lead_p4+sublead_p4;
          
          myAllLeadJPt = jet1->Pt();
          myAllSubJPt = jet2->Pt();
          myAllLeadJEta = jet1->Eta();
          myAllSubJEta = jet2->Eta();
          myAll_Mjj = dijet.M();
          myAlldEta = fabs(jet1->Eta() - jet2->Eta());
          myAllZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myAlldPhi = fabs(diphoton.DeltaPhi(dijet));
          myAll_Mgg =diphoton.M();
          myAllPtHiggs =diphoton.Pt();

          myVBFLeadJPt = jet1->Pt();
          myVBFSubJPt = jet2->Pt();
          myVBF_Mjj = dijet.M();
          myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
          myVBFZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVBFdPhi = fabs(diphoton.DeltaPhi(dijet));
          myVBF_Mgg =diphoton.M();

	  //should be done for both earlier
	  diphoton_index = std::make_pair( l.dipho_leadind[diphotonVBF_id],  l.dipho_subleadind[diphotonVBF_id] );
	  float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
	  float myweight=1.;
	  if(evweight*weight!=0) myweight=evweight/weight;
	  
          l.ApplyCutsFill(0,3,evweight, myweight);

          VBFevent = l.ApplyCutsFill(0,1,evweight, myweight);
        }
      }

      if(diphotonVHad_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
	      TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );

        bool VHadpresel = (highestPtJets.first!=-1)&&(highestPtJets.first!=-2);
  
        if(VHadpresel){
          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
	        
          TLorentzVector diphoton = lead_p4+sublead_p4;
          
          myVHadLeadJPt = jet1->Pt();
          myVHadSubJPt = jet2->Pt();
          myVHad_Mjj = dijet.M();
          myVHaddEta = fabs(jet1->Eta() - jet2->Eta());
          myVHadZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVHaddPhi = fabs(diphoton.DeltaPhi(dijet));
          myVHad_Mgg =diphoton.M();

	        float evweight = weight * smeared_pho_weight[l.dipho_leadind[diphotonVHad_id]] * smeared_pho_weight[l.dipho_vtxind[diphotonVHad_id]] * genLevWeight;
	        float myweight=1.;
	        if(evweight*weight!=0) myweight=evweight/weight;
	  
          VHadevent = l.ApplyCutsFill(0,2,evweight, myweight);
        }
      }
    }
    
    if(VBFevent) diphoton_id = diphotonVBF_id;
    else if(VHadevent) diphoton_id = diphotonVHad_id;
    
    // CP

    /// std::cerr << "Selected pair " << l.dipho_n << " " << diphoton_id << std::endl;
    if (diphoton_id > -1 ) {

	diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
    	// bring all the weights together: lumi & Xsection, single gammas, pt kfactor
	float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

	l.countersred[diPhoCounter_]++;

	TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
	float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
	TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

	bool CorrectVertex;
	
  
  // FIXME pass smeared R9
	int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
	int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
	if( cur_type != 0 && doMCSmearing ) {
	    float pth = Higgs.Pt();
	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
		float rewei=1.;
		(*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
		if( rewei < 0. ) {
		    std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		    assert(0);
		}
		evweight *= rewei;
	    }
	    CorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
	}
	float mass    = Higgs.M();
	float ptHiggs = Higgs.Pt();
      
	assert( evweight >= 0. ); 

	l.FillCounter( "Accepted", weight );
	l.FillCounter( "Smeared", evweight );
	sumaccept += weight;
 	sumsmear += evweight;

  if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
  else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;



	// control plots 
	l.FillHist("all_mass",0, Higgs.M(), evweight);
	l.FillHist("all_mass",category+1, Higgs.M(), evweight);
	if( mass>=massMin && mass<=massMax  ) {
		l.FillHist("mass",0, Higgs.M(), evweight);
		l.FillHist("pt",0, Higgs.Pt(), evweight);
		l.FillHist("eta",0, Higgs.Eta(), evweight);
		
		l.FillHist("pho_pt",0,lead_p4.Pt(), evweight);
		l.FillHist("pho1_pt",0,lead_p4.Pt(), evweight);
		l.FillHist("pho_eta",0,lead_p4.Eta(), evweight);
		l.FillHist("pho1_eta",0,lead_p4.Eta(), evweight);
		l.FillHist("pho_r9",0, lead_r9, evweight);
		l.FillHist("pho1_r9",0, lead_r9, evweight);
		
		l.FillHist("pho_pt",0,sublead_p4.Pt(), evweight);
		l.FillHist("pho2_pt",0,sublead_p4.Pt(), evweight);
		l.FillHist("pho_eta",0,sublead_p4.Eta(), evweight);
		l.FillHist("pho2_eta",0,sublead_p4.Eta(), evweight);
		l.FillHist("pho_r9",0, sublead_r9, evweight);
		l.FillHist("pho1_r9",0, sublead_r9, evweight);
		
		l.FillHist("mass",category+1, Higgs.M(), evweight);
		l.FillHist("pt",category+1, Higgs.Pt(), evweight);
		l.FillHist("eta",category+1, Higgs.Eta(), evweight);
		
		l.FillHist("pho_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho1_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,lead_p4.Eta(), evweight);
		l.FillHist("pho1_eta",category+1,lead_p4.Eta(), evweight);
		l.FillHist("pho_r9",category+1, lead_r9, evweight);
		l.FillHist("pho1_r9",category+1, lead_r9, evweight);
		
		l.FillHist("pho_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho2_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,sublead_p4.Eta(), evweight);
		l.FillHist("pho2_eta",category+1,sublead_p4.Eta(), evweight);
		l.FillHist("pho_r9",category+1, sublead_r9, evweight);
		l.FillHist("pho1_r9",category+1, sublead_r9, evweight);
		
		l.FillHist("pho_n",category+1,l.pho_n, evweight);
	}

	if (cur_type==0){
	  eventListText << setprecision(4) <<"Type = "<< cur_type <<  "Run = " << l.run << "  LS = " << l.lumis << "  Event = " << l.event << "  SelVtx = " << l.vtx_std_sel << "  CAT4 = " << category % 4 << "  ggM = " << mass << " gg_Pt =  " << ptHiggs;
	  eventListText << endl;
	}
       
	// --------------------------------------------------------------------------------------------- 
	if (cur_type == 0 ){
	    l.rooContainer->InputDataPoint("data_mass",category,mass);
	}
	if (cur_type > 0 && cur_type != 3 && cur_type != 4)
	    l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
	else if (cur_type < 0){
	    l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type),category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_rv",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_wv",category,mass,evweight);
	}
       
    }
   
   
    // Systematics
    if( cur_type != 0 && doMCSmearing ) { 
	// fill steps for syst uncertainty study
	float systStep = systRange / (float)nSystSteps;
	// di-photon smearers systematics
	if (diphoton_id > -1 ) {
	       
	    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
	 
	    for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
		std::vector<double> mass_errors;
		std::vector<double> weights;
		std::vector<int>    categories;
	   
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	     
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
		    double genLevWeightSyst=1; 
	     
		    for(std::vector<BaseGenLevelSmearer *>::iterator sj=genLevelSmearers_.begin(); sj!= genLevelSmearers_.end(); ++sj ) {
			float swei=1.;
			if( *si == *sj ) { 
			    (*si)->smearEvent(swei, gP4, l.pu_n, cur_type, syst_shift );
			} else {
			    (*sj)->smearEvent(swei, gP4, l.pu_n, cur_type, 0. );
			}
			genLevWeightSyst *= swei;
		    }
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeightSyst;
	     
		    float mass = Higgs.M();
        
        if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
        else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;
		    categories.push_back(category);
		    mass_errors.push_back(mass);
		    weights.push_back(evweight);
		}// end loop on systematics steps
	   
		if (cur_type < 0){
		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
		}
	    }// end loop on smearers 
		 

	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		std::vector<double> mass_errors;
		std::vector<double> weights;
		std::vector<int> categories;
		       
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
			       
		    // restart with 'fresh' wait for this round of systematics
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
			       
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
		    }
		    float mass = Higgs.M();
        if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
        else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;
		    categories.push_back(category);
		    mass_errors.push_back(mass);
		    weights.push_back(evweight);
		}
		if (cur_type < 0){
		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
		}
	    }

	}
       
	// loop over the smearers included in the systematics study
	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	    std::vector<double> mass_errors;
	    std::vector<double> weights;
	    std::vector<int> categories;
	   
	    // loop over syst shift
	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		if( syst_shift == 0. ) { continue; } // skip the central value
		// smear the photons 
		for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		    std::vector<std::vector<bool> > p;
		    //std::cout << "GF check: " <<  l.pho_residCorrEnergy[ipho] << "  " << l.pho_residCorrResn[ipho] << std::endl;
		    PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
						/// *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
						((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
						energyCorrected[ipho],
						l.pho_isEB[ipho], l.pho_r9[ipho],
						l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_));
		   
		    float pweight = 1.;
		    for(std::vector<BaseSmearer *>::iterator  sj=photonSmearers_.begin(); sj!= photonSmearers_.end(); ++sj ) {
			float sweight = 1.;
			if( *si == *sj ) {
			    // move the smearer under study by syst_shift
			    (*si)->smearPhoton(phoInfo,sweight,l.run,syst_shift);
			} else {
			    // for the other use the nominal points
			    (*sj)->smearPhoton(phoInfo,sweight,l.run,0.);
			}
			pweight *= sweight;
		    }
		    smeared_pho_energy[ipho] = phoInfo.energy();
		    smeared_pho_r9[ipho] = phoInfo.r9();
		    smeared_pho_weight[ipho] = pweight;
		}
	       
		// analyze the event
		// FIXME pass smeared R9
		int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
    int diphotonVBF_id = -1;
    int diphotonVHad_id = -1;
    bool VBFevent = false;
    bool VHadevent = false;
    std::pair<int,int> highestPtJets(-1,-1);
    if((includeVBF || includeVHad)&&l.jet_algoPF1_n>1) {
      if(includeVBF) diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
      if(includeVHad) diphotonVHad_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHadCut, subleadEtVHadCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 

      if(diphotonVBF_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
	      TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );

        bool VBFpresel = (highestPtJets.first!=-1)&&(highestPtJets.first!=-2);
  
        if(VBFpresel){
          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
	        
          TLorentzVector diphoton = lead_p4+sublead_p4;
          
          myVBFLeadJPt = jet1->Pt();
          myVBFSubJPt = jet2->Pt();
          myVBF_Mjj = dijet.M();
          myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
          myVBFZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVBFdPhi = fabs(diphoton.DeltaPhi(dijet));
          myVBF_Mgg =diphoton.M();

          VBFevent = l.ApplyCuts(0,1);
        }
      }
      
      if(diphotonVHad_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
	      TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );

        bool VHadpresel = (highestPtJets.first!=-1)&&(highestPtJets.first!=-2);
  
        if(VHadpresel){
          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
	        
          TLorentzVector diphoton = lead_p4+sublead_p4;
          
          myVHadLeadJPt = jet1->Pt();
          myVHadSubJPt = jet2->Pt();
          myVHad_Mjj = dijet.M();
          myVHaddEta = fabs(jet1->Eta() - jet2->Eta());
          myVHadZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVHaddPhi = fabs(diphoton.DeltaPhi(dijet));
          myVHad_Mgg =diphoton.M();

          VHadevent = l.ApplyCuts(0,2);
        }
      }
    }
    
    if(VBFevent) diphoton_id = diphotonVBF_id;
    else if(VHadevent) diphoton_id = diphotonVHad_id;
	       
		if (diphoton_id > -1 ) {
		   
		    diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] *genLevWeight;
		   
		    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
		   
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
		    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		    if( cur_type != 0 && doMCSmearing ) {
			for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
			    float rewei=1.;
			    float pth = Higgs.Pt();
			    (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			    evweight *= rewei;
			}
		    }
		    float mass = Higgs.M();
		   
              if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
              else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;
	       	    categories.push_back(category);
	            mass_errors.push_back(mass);
	            weights.push_back(evweight);

		} else {
		    mass_errors.push_back(0.);   
		    weights.push_back(0.);   
		    categories.push_back(-1);
		}
	       
	    }
	    if (cur_type < 0){
	       l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
	    }
       
	}
       
       
    }
   
    if(PADEBUG) 
	cout<<"myFillHistRed END"<<endl;
}

// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool StatAnalysisExclusive::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}
// ----------------------------------------------------------------------------------------------------
double StatAnalysisExclusive::GetDifferentialKfactor(double gPT, int Mass)
{

/*  
    if (Mass <=110 ) return thm110->GetBinContent(thm110->FindFixBin(gPT));
    else if (Mass ==120 ) return thm120->GetBinContent(thm120->FindFixBin(gPT));
    else if (Mass ==130 ) return thm130->GetBinContent(thm130->FindFixBin(gPT));
    else if (Mass ==140 ) return thm140->GetBinContent(thm140->FindFixBin(gPT));
    else if (Mass ==115 ) return (0.5*thm110->GetBinContent(thm110->FindFixBin(gPT)) +0.5*thm120->GetBinContent(thm120->FindFixBin(gPT)));
*/
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

void StatAnalysisExclusive::FillSignalLabelMap(){

	// Basically A Map of the ID (type) to the signal's name which can be filled Now:
  signalLabels[-57]="ggh_mass_m123";
  signalLabels[-58]="vbf_mass_m123";
  signalLabels[-60]="wzh_mass_m123";
  signalLabels[-59]="tth_mass_m123";
  signalLabels[-53]="ggh_mass_m121";
  signalLabels[-54]="vbf_mass_m121";
  signalLabels[-56]="wzh_mass_m121";
  signalLabels[-55]="tth_mass_m121";
  signalLabels[-65]="ggh_mass_m160";
  signalLabels[-66]="vbf_mass_m160";
  signalLabels[-68]="wzh_mass_m160";
  signalLabels[-67]="tth_mass_m160";
  signalLabels[-61]="ggh_mass_m155";
  signalLabels[-62]="vbf_mass_m155";
  signalLabels[-64]="wzh_mass_m155";
  signalLabels[-63]="tth_mass_m155";
  signalLabels[-49]="ggh_mass_m150";
  signalLabels[-50]="vbf_mass_m150";
  signalLabels[-52]="wzh_mass_m150";
  signalLabels[-51]="tth_mass_m150";
  signalLabels[-45]="ggh_mass_m145";
  signalLabels[-46]="vbf_mass_m145";
  signalLabels[-48]="wzh_mass_m145";
  signalLabels[-47]="tth_mass_m145";
  signalLabels[-33]="ggh_mass_m140";
  signalLabels[-34]="vbf_mass_m140";
  signalLabels[-36]="wzh_mass_m140";
  signalLabels[-35]="tth_mass_m140";
  signalLabels[-41]="ggh_mass_m135";
  signalLabels[-42]="vbf_mass_m135";
  signalLabels[-44]="wzh_mass_m135";
  signalLabels[-43]="tth_mass_m135";
  signalLabels[-29]="ggh_mass_m130";
  signalLabels[-30]="vbf_mass_m130";
  signalLabels[-32]="wzh_mass_m130";
  signalLabels[-31]="tth_mass_m130";
  signalLabels[-37]="ggh_mass_m125";
  signalLabels[-38]="vbf_mass_m125";
  signalLabels[-40]="wzh_mass_m125";
  signalLabels[-39]="tth_mass_m125";
  signalLabels[-25]="ggh_mass_m120";
  signalLabels[-26]="vbf_mass_m120";
  signalLabels[-28]="wzh_mass_m120";
  signalLabels[-27]="tth_mass_m120";
  signalLabels[-21]="ggh_mass_m115";
  signalLabels[-22]="vbf_mass_m115";
  signalLabels[-24]="wzh_mass_m115";
  signalLabels[-23]="tth_mass_m115";
  signalLabels[-17]="ggh_mass_m110";
  signalLabels[-18]="vbf_mass_m110";
  signalLabels[-19]="wzh_mass_m110";
  signalLabels[-20]="tth_mass_m110";
  signalLabels[-13]="ggh_mass_m105";
  signalLabels[-14]="vbf_mass_m105";
  signalLabels[-16]="wzh_mass_m105";
  signalLabels[-15]="tth_mass_m105";
  signalLabels[-69]="ggh_mass_m100";
  signalLabels[-70]="vbf_mass_m100";
  signalLabels[-72]="wzh_mass_m100";
  signalLabels[-71]="tth_mass_m100";
}

std::string StatAnalysisExclusive::GetSignalLabel(int id){
	
	// For the lazy man, can return a memeber of the map rather than doing it yourself
	std::map<int,std::string>::iterator it = signalLabels.find(id);

	if (it!=signalLabels.end()){
		return it->second;
		
	} else { 

		std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
		return "NULL";
	}
	
}

std::pair<int, int> StatAnalysisExclusive::Select2HighestPtJets(LoopAll& l, TLorentzVector& leadpho, TLorentzVector& subleadpho, float jtLMinPt, float jtTMinPt){

  std::pair<int, int> myJets(-1,-1);
  std::pair<int, int> fail(-1,-1);

  float dr2pho = 0.5;
  float dr2jet = 0.5;

  TLorentzVector* j1p4;
  TLorentzVector* j2p4;
  float j1pt=-1;
  float j2pt=-1;

  // select ighest pt jets
  // veto jets close to photons or each other
  for(int j1_i=0; j1_i<l.jet_algoPF1_n; j1_i++){
    j1p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(j1_i);
    if(fabs(j1p4->Eta()) > 4.7) continue;
    if(j1p4->DeltaR(leadpho) < dr2pho) continue;
    if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
    j1pt=j1p4->Pt();
    if(j1pt<jtTMinPt) continue;
    for(int j2_i=j1_i+1; j2_i<l.jet_algoPF1_n; j2_i++){
      j2p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(j2_i);
      if(fabs(j2p4->Eta()) > 4.7) continue;
      if(j2p4->DeltaR(leadpho) < dr2pho) continue;
      if(j2p4->DeltaR(subleadpho) < dr2pho) continue;
      if(j2p4->DeltaR(*j1p4) < dr2jet) continue;
      j2pt=j2p4->Pt();
      
      if(j2pt<jtTMinPt) continue;
      if(std::max(j1pt,j2pt)<jtLMinPt) continue;

      if(j1pt>j2pt){
        jtLMinPt=j1pt;
        jtTMinPt=j2pt;

        myJets.first = j1_i;
        myJets.second = j2_i;
      } else {
        jtLMinPt=j2pt;
        jtTMinPt=j1pt;

        myJets.first = j2_i;
        myJets.second = j1_i;
      }
    }
  }

  if(myJets.first==-1) return fail;


  return myJets;
}



int  StatAnalysisExclusive::RescaleJetEnergy(LoopAll& l) {
  for (int i = 0; i<l.jet_algoPF1_n; i++) {
    TLorentzVector * thisjet = (TLorentzVector *) l.jet_algoPF1_p4->At(i);
    *thisjet*=l.jet_algoPF1_erescale[i];
  }
  return 1;
}
