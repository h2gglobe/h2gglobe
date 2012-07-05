#include "../interface/StatAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>
#include <stdio.h>

#define PADEBUG 0

using namespace std;

void dumpPhoton(std::ostream & eventListText, int lab, 
		LoopAll & l, int ipho, int ivtx, TLorentzVector & phop4, float * pho_energy_array);
void dumpJet(std::ostream & eventListText, int lab, LoopAll & l, int ijet);

    ofstream met_sync;
    ofstream lep_sync;
    
// ----------------------------------------------------------------------------------------------------
StatAnalysis::StatAnalysis()  : 
    name_("StatAnalysis")
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
    doSystematics = true;   
    dataIs2011 = false;
    nVBFDijetJetCategories=2;
    scaleClusterShapes = true;
    dumpAscii = false;
    dumpMcAscii = false;
    unblind = false;
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
    std::string postfix=(dataIs2011?"":"_8TeV");
    l.rooContainer->FitToData("data_pol_model"+postfix,"data_mass");  // Fit to full range of dataset
  
    //    l.rooContainer->WriteSpecificCategoryDataCards(outputfilename,"data_mass","sig_mass","data_pol_model");
    //    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","data_pol_model");
    // mode 0 as above, 1 if want to bin in sub range from fit,

    // Write the data-card for the Combinations Code, needs the output filename, makes binned analysis DataCard
    // Assumes the signal datasets will be called signal_name+"_mXXX"
    //    l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","data_pol_model","data_mass",1,50,1); // 1 means systematics from the fit effect only the backgroundi. last digit mode = 1 means this is an internal constraint fit 
    //    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","bkg_mass_rebinned");

    eventListText.close();
    lep_sync.close();

    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;

    met_sync.close();

    //  kfacFile->Close();
    //  PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
        cout << "InitRealStatAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;
    
    met_sync.open ("met_sync.txt");
    
    std::string outputfilename = (std::string) l.histFileName;
    eventListText.open(Form("%s",l.outputTextFileName.c_str()));
    lep_sync.open ("lep_sync.txt");
    //eventListText.open(Form("%s_ascii_events.txt",outputfilename.c_str()));
    FillSignalLabelMap(l);
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
        << "doEcorrectionSyst "<< doEcorrectionSyst << "\n"
        << "efficiencyFile " << efficiencyFile << "\n"
        << "doPhotonIdEffSyst "<< doPhotonIdEffSyst << "\n"
        << "doR9Syst "<< doR9Syst << "\n"
        << "doVtxEffSyst "<< doVtxEffSyst << "\n"
        << "doTriggerEffSyst "<< doTriggerEffSyst << "\n"
        << "doKFactorSyst "<< doKFactorSyst << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << std::endl;

    // call the base class initializer
    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
        l.histoContainer[ind].setScale(1.);
    }
    
    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nInclusiveCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nInclusiveCategories_ *= nR9Categories;
    if( nPtCategories != 0 ) nInclusiveCategories_ *= nPtCategories;

    // mva removed cp march 8
    //if( useMVA ) nInclusiveCategories_ = nDiphoEventClasses;

    // CP

    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    
    int nVBFCategories   = ((int)includeVBF)*nVBFEtaCategories*nVBFDijetJetCategories;
    int nVHhadCategories = ((int)includeVHhad)*nVHhadEtaCategories;
    int nVHlepCategories = (int)includeVHlep * 2;
    int nVHmetCategories = (int)includeVHmet;  //met at analysis step
    
    nCategories_=(nInclusiveCategories_+nVBFCategories+nVHhadCategories+nVHlepCategories+nVHmetCategories);  //met at analysis step
//    nCategories_=(nInclusiveCategories_+nVBFCategories+nVHhadCategories+nVHlepCategories);
    

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
    if(doInterferenceSmear) {
        // interference efficiency
        std::cerr << __LINE__ << std::endl; 
        interferenceSmearer = new InterferenceSmearer(2.5e-2,0.);
        genLevelSmearers_.push_back(interferenceSmearer);
    }

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed

    // FIXME move these params to config file
    l.rooContainer->SetNCategories(nCategories_);
    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;

    if( doEcorrectionSmear && doEcorrectionSyst ) {
        // instance of this smearer done in PhotonAnalysis
        systPhotonSmearers_.push_back(eCorrSmearer);
        std::vector<std::string> sys(1,eCorrSmearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEscaleSmear && doEscaleSyst ) {
        systPhotonSmearers_.push_back( eScaleSmearer );
        std::vector<std::string> sys(1,eScaleSmearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
        systPhotonSmearers_.push_back( eResolSmearer );
        std::vector<std::string> sys(1,eResolSmearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonIdEffSmear && doPhotonIdEffSyst ) {
        systPhotonSmearers_.push_back( idEffSmearer );
        std::vector<std::string> sys(1,idEffSmearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doR9Smear && doR9Syst ) {
        systPhotonSmearers_.push_back( r9Smearer );
        std::vector<std::string> sys(1,r9Smearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doVtxEffSmear && doVtxEffSyst ) {
        systDiPhotonSmearers_.push_back( vtxEffSmearer );
        std::vector<std::string> sys(1,vtxEffSmearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doTriggerEffSmear && doTriggerEffSyst ) {
        systDiPhotonSmearers_.push_back( triggerEffSmearer );
        std::vector<std::string> sys(1,triggerEffSmearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if(doKFactorSmear && doKFactorSyst) {
        systGenLevelSmearers_.push_back(kFactorSmearer);
        std::vector<std::string> sys(1,kFactorSmearer->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    
    // ----------------------------------------------------
    // ----------------------------------------------------
    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.045,1.00);
    // ----------------------------------------------------

    // Create observables for shape-analysis with ranges
    // l.rooContainer->AddObservable("mass" ,100.,150.);
    l.rooContainer->AddObservable("CMS_hgg_mass" ,massMin,massMax);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);

    // SM Model
    l.rooContainer->AddConstant("XSBR_tth_155",0.00004370);
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
    std::string postfix=(dataIs2011?"":"_8TeV");
    /////// l.rooContainer->AddRealVar("CMS_hgg_pol0"+postfix,-0.1,-1.0,1.0);
    /////// l.rooContainer->AddRealVar("CMS_hgg_pol1"+postfix,-0.1,-1.0,1.0);
    /////// l.rooContainer->AddRealVar("CMS_hgg_pol2"+postfix,-0.1,-1.0,1.0);
    /////// l.rooContainer->AddRealVar("CMS_hgg_pol3"+postfix,-0.01,-1.0,1.0);
    /////// l.rooContainer->AddRealVar("CMS_hgg_pol4"+postfix,-0.01,-1.0,1.0);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modpol0"+postfix,"@0*@0","CMS_hgg_pol0"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modpol1"+postfix,"@0*@0","CMS_hgg_pol1"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modpol2"+postfix,"@0*@0","CMS_hgg_pol2"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modpol3"+postfix,"@0*@0","CMS_hgg_pol3"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modpol4"+postfix,"@0*@0","CMS_hgg_pol4"+postfix);
    /////// 
    /////// l.rooContainer->AddRealVar("CMS_hgg_quartic0"+postfix,-0.1,-1.0,1.0);
    /////// l.rooContainer->AddRealVar("CMS_hgg_quartic1"+postfix,-0.1,-1.0,1.0);
    /////// l.rooContainer->AddRealVar("CMS_hgg_quartic2"+postfix,-0.1,-1.0,1.0);
    /////// l.rooContainer->AddRealVar("CMS_hgg_quartic3"+postfix,-0.01,-1.0,1.0);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modquartic0"+postfix,"@0*@0","CMS_hgg_quartic0"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modquartic1"+postfix,"@0*@0","CMS_hgg_quartic1"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modquartic2"+postfix,"@0*@0","CMS_hgg_quartic2"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modquartic3"+postfix,"@0*@0","CMS_hgg_quartic3"+postfix);
    /////// 
    /////// l.rooContainer->AddRealVar("CMS_hgg_quad0"+postfix,-0.1,-1.5,1.5);
    /////// l.rooContainer->AddRealVar("CMS_hgg_quad1"+postfix,-0.01,-1.5,1.5);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modquad0"+postfix,"@0*@0","CMS_hgg_quad0"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modquad1"+postfix,"@0*@0","CMS_hgg_quad1"+postfix);
    /////// 
    /////// l.rooContainer->AddRealVar("CMS_hgg_cubic0"+postfix,-0.1,-1.5,1.5);
    /////// l.rooContainer->AddRealVar("CMS_hgg_cubic1"+postfix,-0.1,-1.5,1.5);
    /////// l.rooContainer->AddRealVar("CMS_hgg_cubic2"+postfix,-0.01,-1.5,1.5);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modcubic0"+postfix,"@0*@0","CMS_hgg_cubic0"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modcubic1"+postfix,"@0*@0","CMS_hgg_cubic1"+postfix);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modcubic2"+postfix,"@0*@0","CMS_hgg_cubic2"+postfix);
    /////// 
    /////// l.rooContainer->AddRealVar("CMS_hgg_lin0"+postfix,-0.01,-1.5,1.5);
    /////// l.rooContainer->AddFormulaVar("CMS_hgg_modlin0"+postfix,"@0*@0","CMS_hgg_lin0"+postfix);
    /////// 
    /////// // Generic PDF ok in the std analysis but excluisve channels need different models CP
    /////// //l.rooContainer->AddGenericPdf("data_pol_model",
    /////// //"0","CMS_hgg_mass",data_pol_pars,73); // >= 71 means RooBernstein of order >= 1
    //////// int cats_with_std[]     = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //////// int cats_with_lin[]     = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //////// int cats_with_quad[]    = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //////// int cats_with_cubic[]   = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //////// int cats_with_quartic[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    //////////// std::cout << "Number of categories: " << nCategories_ << std::endl;
    //////////// for(int i=0; i<nCategories_; i++){
    ////////////     if(i<nInclusiveCategories_) {
    ////////////         cats_with_std[i]=1;
    ////////////     } else if(i<nInclusiveCategories_+nVBFCategories){
    ////////////         /// cats_with_quad[i]=1;
    ////////////         cats_with_cubic[i]=1;
    ////////////     } else if(i<nInclusiveCategories_+nVBFCategories+nVHhadCategories){
    ////////////         cats_with_quad[i]=1;
    ////////////     } else if(i<nInclusiveCategories_+nVBFCategories+nVHhadCategories+nVHlepCategories){
    ////////////         cats_with_lin[i]=1;
    ////////////     } else {
    ////////////         cats_with_cubic[i]=1;
    ////////////     }
    //////////// }
    /////// 
    /////// std::vector<std::string> data_pol_pars(5,"p");   
    /////// data_pol_pars[0] = "CMS_hgg_modpol0"+postfix;
    /////// data_pol_pars[1] = "CMS_hgg_modpol1"+postfix;
    /////// data_pol_pars[2] = "CMS_hgg_modpol2"+postfix;
    /////// data_pol_pars[3] = "CMS_hgg_modpol3"+postfix;
    /////// data_pol_pars[4] = "CMS_hgg_modpol4"+postfix;
    /////// 
    /////// l.rooContainer->AddSpecificCategoryPdf(cats_with_std,"data_pol_model"+postfix,
    ///////                                        "0","CMS_hgg_mass",data_pol_pars,75);    // >= 71 means RooBernstein of order >= 1
    /////// 
    /////// std::vector<std::string> data_quartic_pars(4,"p");   
    /////// data_quartic_pars[0] = "CMS_hgg_modquartic0"+postfix;
    /////// data_quartic_pars[1] = "CMS_hgg_modquartic1"+postfix;
    /////// data_quartic_pars[2] = "CMS_hgg_modquartic2"+postfix;
    /////// data_quartic_pars[3] = "CMS_hgg_modquartic3"+postfix;
    /////// 
    /////// l.rooContainer->AddSpecificCategoryPdf(cats_with_quartic,"data_pol_model"+postfix,
    ///////                                        "0","CMS_hgg_mass",data_quartic_pars,74);    // >= 71 means RooBernstein of order >= 1
    /////// 
    /////// std::vector<std::string> data_cubic_pars(3,"p");     
    /////// data_cubic_pars[0] = "CMS_hgg_modcubic0"+postfix;
    /////// data_cubic_pars[1] = "CMS_hgg_modcubic1"+postfix;
    /////// data_cubic_pars[2] = "CMS_hgg_modcubic2"+postfix;
    /////// 
    /////// l.rooContainer->AddSpecificCategoryPdf(cats_with_cubic, "data_pol_model"+postfix,
    ///////                                        "0","CMS_hgg_mass",data_cubic_pars,73);  // >= 71 means RooBernstein of order >= 1
    /////// 
    /////// 
    /////// std::vector<std::string> data_quad_pars(2,"p");  
    /////// data_quad_pars[0] = "CMS_hgg_modquad0"+postfix;
    /////// data_quad_pars[1] = "CMS_hgg_modquad1"+postfix;
    /////// 
    /////// l.rooContainer->AddSpecificCategoryPdf(cats_with_quad, "data_pol_model"+postfix,
    ///////                                        "0","CMS_hgg_mass",data_quad_pars,72);   // >= 71 means RooBernstein of order >= 1
    /////// 
    /////// 
    /////// std::vector<std::string> data_lin_pars(1,"p");   
    /////// data_lin_pars[0] = "CMS_hgg_modlin0"+postfix;
    /////// 
    /////// l.rooContainer->AddSpecificCategoryPdf(cats_with_lin, "data_pol_model"+postfix,
    /////// "0","CMS_hgg_mass",data_lin_pars,71);    // >= 71 means RooBernstein of order >= 1
    // CP

    // -----------------------------------------------------
    // Configurable background model
    // if no configuration was given, set some defaults
    if( bkgPolOrderByCat.empty() ) {
	for(int i=0; i<nCategories_; i++){
	    if(i<nInclusiveCategories_) {
		bkgPolOrderByCat.push_back(5);
	    } else if(i<nInclusiveCategories_+nVBFCategories){
		bkgPolOrderByCat.push_back(3);
	    } else if(i<nInclusiveCategories_+nVBFCategories+nVHhadCategories){
		bkgPolOrderByCat.push_back(2);
	    } else if(i<nInclusiveCategories_+nVBFCategories+nVHhadCategories+nVHlepCategories){
		bkgPolOrderByCat.push_back(1);
	    }
	}
    }
    // build the model
    buildBkgModel(l, postfix);

    // -----------------------------------------------------
    // Make some data sets from the observables to fill in the event loop         
    // Binning is for histograms (will also produce unbinned data sets)
    l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass"    ,nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass"     ,nDataBins);            

    // Create Signal DataSets:
    for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
	////// for (int sig=105;sig<=150;sig+=5){
	//////     // Needed to use S4 for the GGH 145 Signal which has the BUG so no 145 sample
        ////// if ( sig==145 // && dataIs2011 
	//////     ) continue;
	int sig = sigPointsToBook[isig];
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

    // Make more datasets representing Systematic Shifts of various quantities
    for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
	////// for (int sig=105;sig<=150;sig+=5){
	//////     // Needed to use S4 for the GGH 145 Signal which has the BUG so no 145 sample
        ////// if ( sig==145 // && dataIs2011 
	//////     ) continue;
	int sig = sigPointsToBook[isig];
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_ggh_mass_m%d",sig),-1);    
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_vbf_mass_m%d",sig),-1);    
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_wzh_mass_m%d",sig),-1);    
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_tth_mass_m%d",sig),-1);    
    }

    // Make sure the Map is filled
    FillSignalLabelMap(l);

    if(PADEBUG) 
        cout << "InitRealStatAnalysis END"<<endl;
    
    // FIXME book of additional variables
}


// ----------------------------------------------------------------------------------------------------
void StatAnalysis::buildBkgModel(LoopAll& l, const std::string & postfix) 
{ 

    // sanity check
    assert( bkgPolOrderByCat.size() == nCategories_ );

    l.rooContainer->AddRealVar("CMS_hgg_pol6_0"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol6_1"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol6_2"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol6_3"+postfix,-0.01,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol6_4"+postfix,-0.01,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol6_5"+postfix,-0.01,-1.0,1.0);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol6_0"+postfix,"@0*@0","CMS_hgg_pol6_0"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol6_1"+postfix,"@0*@0","CMS_hgg_pol6_1"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol6_2"+postfix,"@0*@0","CMS_hgg_pol6_2"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol6_3"+postfix,"@0*@0","CMS_hgg_pol6_3"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol6_4"+postfix,"@0*@0","CMS_hgg_pol6_4"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol6_5"+postfix,"@0*@0","CMS_hgg_pol6_4"+postfix);

    l.rooContainer->AddRealVar("CMS_hgg_pol5_0"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol5_1"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol5_2"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol5_3"+postfix,-0.01,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol5_4"+postfix,-0.01,-1.0,1.0);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol5_0"+postfix,"@0*@0","CMS_hgg_pol5_0"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol5_1"+postfix,"@0*@0","CMS_hgg_pol5_1"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol5_2"+postfix,"@0*@0","CMS_hgg_pol5_2"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol5_3"+postfix,"@0*@0","CMS_hgg_pol5_3"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol5_4"+postfix,"@0*@0","CMS_hgg_pol5_4"+postfix);

    l.rooContainer->AddRealVar("CMS_hgg_quartic0"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_quartic1"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_quartic2"+postfix,-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_quartic3"+postfix,-0.01,-1.0,1.0);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic0"+postfix,"@0*@0","CMS_hgg_quartic0"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic1"+postfix,"@0*@0","CMS_hgg_quartic1"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic2"+postfix,"@0*@0","CMS_hgg_quartic2"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic3"+postfix,"@0*@0","CMS_hgg_quartic3"+postfix);

    l.rooContainer->AddRealVar("CMS_hgg_quad0"+postfix,-0.1,-1.5,1.5);
    l.rooContainer->AddRealVar("CMS_hgg_quad1"+postfix,-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquad0"+postfix,"@0*@0","CMS_hgg_quad0"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquad1"+postfix,"@0*@0","CMS_hgg_quad1"+postfix);
    
    l.rooContainer->AddRealVar("CMS_hgg_cubic0"+postfix,-0.1,-1.5,1.5);
    l.rooContainer->AddRealVar("CMS_hgg_cubic1"+postfix,-0.1,-1.5,1.5);
    l.rooContainer->AddRealVar("CMS_hgg_cubic2"+postfix,-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("CMS_hgg_modcubic0"+postfix,"@0*@0","CMS_hgg_cubic0"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modcubic1"+postfix,"@0*@0","CMS_hgg_cubic1"+postfix);
    l.rooContainer->AddFormulaVar("CMS_hgg_modcubic2"+postfix,"@0*@0","CMS_hgg_cubic2"+postfix);
    
    l.rooContainer->AddRealVar("CMS_hgg_lin0"+postfix,-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("CMS_hgg_modlin0"+postfix,"@0*@0","CMS_hgg_lin0"+postfix);

    // prefix for models parameters
    std::map<int,std::string> parnames;
    parnames[1] = "modlin";
    parnames[2] = "modquad";
    parnames[3] = "modcubic";
    parnames[4] = "modquartic";
    parnames[5] = "modpol5_";
    parnames[6] = "modpol6_";
    
    // map order to categories flags + parameters names
    std::map<int, std::pair<std::vector<int>, std::vector<std::string> > > catmodels;
    // fill the map
    for(int icat=0; icat<nCategories_; ++icat) {
	// get the poly order for this category
	int catmodel = bkgPolOrderByCat[icat];
	std::vector<int> & catflags = catmodels[catmodel].first;
	std::vector<std::string> & catpars = catmodels[catmodel].second;
	// if this is the first time we find this order, build the parameters
	if( catflags.empty() ) {
	    assert( catpars.empty() );
	    // by default no category has the new model
	    catflags.resize(nCategories_, 0);
	    std::string & parname = parnames[catmodel];
	    for(int iorder = 0; iorder<catmodel; ++iorder) {
		catpars.push_back( Form( "CMS_hgg_%s%d%s", parname.c_str(), iorder, +postfix.c_str() ) );
	    }
	} else {
	    assert( catflags.size() == nCategories_ && catpars.size() == catmodel );
	}
	// chose category order
	catflags[icat] = 1;
    }

    // now loop over the models and allocate the pdfs
    /// for(size_t imodel=0; imodel<catmodels.size(); ++imodel ) {
    for(std::map<int, std::pair<std::vector<int>, std::vector<std::string> > >::iterator modit = catmodels.begin();
	modit!=catmodels.end(); ++modit ) {
	std::vector<int> & catflags = modit->second.first;
	std::vector<std::string> & catpars = modit->second.second;
	
	l.rooContainer->AddSpecificCategoryPdf(&catflags[0],"data_pol_model"+postfix,
					       "0","CMS_hgg_mass",catpars,70+catpars.size()); 
	// >= 71 means RooBernstein of order >= 1
    }
}

// ----------------------------------------------------------------------------------------------------
bool StatAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
    if(PADEBUG) 
        cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
   
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;
    float sampleweight = l.sampleContainer[l.current_sample_index].weight;

    // Set reRunCiC Only if this is an MC event since scaling of R9 and Energy isn't done at reduction
    if (cur_type==0) {
        l.runCiC=reRunCiCForData;
    } else {
        l.runCiC = true;
    } 
    if (l.runZeeValidation) l.runCiC=true;

    // make sure that rho is properly set
    if( l.version >= 13 && forcedRho < 0. ) {
	l.rho = l.rho_algo1;
    }

    l.FillCounter( "Processed", 1. );
    assert( weight > 0. );  
    l.FillCounter( "XSWeighted", weight );
    nevents+=1.;

    //PU reweighting
    double pileupWeight=getPuWeight( l.pu_n, cur_type, &(l.sampleContainer[l.current_sample_index]), jentry == 1);
    sumwei +=pileupWeight;
    weight *= pileupWeight;
    sumev  += weight;
    
    assert( weight >= 0. );  
    l.FillCounter( "PUWeighted", weight );
    
    if( jentry % 1000 ==  0 ) {
        std::cout << " nevents " <<  nevents << " sumpuweights " << sumwei << " ratio " << sumwei / nevents 
                  << " equiv events " << sumev << " accepted " << sumaccept << " smeared " << sumsmear << " "  
                  <<  sumaccept / sumev << " " << sumsmear / sumaccept
                  << std::endl;
    }
    // ------------------------------------------------------------
    //PT-H K-factors
    double gPT = 0;
    TLorentzVector gP4(0,0,0,0);
    if (cur_type<0){
	gP4 = l.GetHiggs();
	gPT = gP4.Pt();
    }

    //Calculate cluster shape variables prior to shape rescaling
    for (int ipho=0;ipho<l.pho_n;ipho++){
	l.pho_s4ratio[ipho]  = l.pho_e2x2[ipho]/l.pho_e5x5[ipho];
	float rr2=l.pho_eseffsixix[ipho]*l.pho_eseffsixix[ipho]+l.pho_eseffsiyiy[ipho]*l.pho_eseffsiyiy[ipho];
	l.pho_ESEffSigmaRR[ipho] = 0.0; 
	if(rr2>0. && rr2<999999.) {
	    l.pho_ESEffSigmaRR[ipho] = sqrt(rr2);
	}
    }

    // Data driven MC corrections to cluster shape variables and energy resolution estimate
    if (cur_type !=0 && scaleClusterShapes ){
        rescaleClusterVariables(l);
    }
    if( reRunVtx ) {
	reVertex(l);
    }

    // Re-apply JEC and / or recompute JetID
    if(includeVBF || includeVHhad) { postProcessJets(l); }
    
    // Analyse the event assuming nominal values of corrections and smearings
    float mass, evweight, diphotonMVA;
    int diphoton_id, category;
    bool isCorrectVertex;

    if( AnalyseEvent(l,jentry, weight, gP4, mass,  evweight, category, diphoton_id, isCorrectVertex,diphotonMVA) ) {
	// feed the event to the RooContainer 
	FillRooContainer(l, cur_type, mass, diphotonMVA, category, evweight, isCorrectVertex);
    }
    
    // Systematics uncertaities for the binned model
    // We re-analyse the event several times for different values of corrections and smearings
    if( cur_type < 0 && doMCSmearing && doSystematics ) { 
        
        // fill steps for syst uncertainty study
        float systStep = systRange / (float)nSystSteps;
	
	float syst_mass, syst_weight, syst_diphotonMVA;
	int syst_category;
 	std::vector<double> mass_errors;
 	std::vector<double> mva_errors;
	std::vector<double> weights;
	std::vector<int>    categories;
	
        if (diphoton_id > -1 ) {
     
	    // gen-level systematics, i.e. ggH k-factor for the moment
            for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
		mass_errors.clear(), weights.clear(), categories.clear(), mva_errors.clear();
		
                for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
                    if( syst_shift == 0. ) { continue; } // skip the central value
		    syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
		    
		    // re-analyse the event without redoing the event selection as we use nominal values for the single photon
		    // corrections and smearings
		    AnalyseEvent(l, jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,syst_diphotonMVA,
				 true, syst_shift, true, *si, 0, 0 );
		    
		    AccumulateSyst( cur_type, syst_mass, syst_diphotonMVA, syst_category, syst_weight,
				    mass_errors, mva_errors, categories, weights);
		}
		
		FillRooContainerSyst(l, (*si)->name(), cur_type, mass_errors, mva_errors, categories, weights);
	    }
	    
	    // di-photon systematics: vertex efficiency and trigger 
	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		mass_errors.clear(), weights.clear(), categories.clear(), mva_errors.clear();
		
                for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
                    if( syst_shift == 0. ) { continue; } // skip the central value
		    syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
		    
		    // re-analyse the event without redoing the event selection as we use nominal values for the single photon
		    // corrections and smearings
		    AnalyseEvent(l,jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,syst_diphotonMVA,
				 true, syst_shift, true,  0, 0, *si );
		    
		    AccumulateSyst( cur_type, syst_mass, syst_diphotonMVA, syst_category, syst_weight,
                                    mass_errors, mva_errors, categories, weights);
		}

		FillRooContainerSyst(l, (*si)->name(), cur_type, mass_errors, mva_errors, categories, weights);
	    }
	}
	
	int diphoton_id_syst;	
	// single photon level systematics: several
	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	    mass_errors.clear(), weights.clear(), categories.clear(), mva_errors.clear();
	    
	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		if( syst_shift == 0. ) { continue; } // skip the central value
		syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
		
		// re-analyse the event redoing the event selection this time
		AnalyseEvent(l,jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id_syst, isCorrectVertex,syst_diphotonMVA,
			     true, syst_shift, false,  0, *si, 0 );
		
		AccumulateSyst( cur_type, syst_mass, syst_diphotonMVA, syst_category, syst_weight,
				mass_errors, mva_errors, categories, weights);
	    }
	    
	    FillRooContainerSyst(l, (*si)->name(), cur_type, mass_errors, mva_errors, categories, weights);
	}
    }

    if(PADEBUG) 
        cout<<"myFillHistRed END"<<endl;

    return (diphoton_id > -1);
}

// ----------------------------------------------------------------------------------------------------
bool StatAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4,
				float & mass, float & evweight, int & category, int & diphoton_id, bool & isCorrectVertex,
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
    int ijet1=0, ijet2=0;
   
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

	// Fill CiC efficiency plots for ggH, mH=124
	if (cur_type==-73) fillSignalEfficiencyPlots(weight, l);

	// inclusive category di-photon selection
	// FIXME pass smeared R9
	diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0], false, false, cicCutLevels ); 
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
	int diphotonVHmet_id = -1; //met at analysis step
	VHmuevent = false;
	VHelevent = false;
	VBFevent = false;
	VHhadevent = false;
	VHmetevent = false; //met at analysis step
	
	// lepton tag
	if(includeVHlep){
	    diphotonVHlep_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHlepCut, subleadEtVHlepCut, 4, false, &smeared_pho_energy[0], true, true );
	    //Add tighter cut on dr to tk
	    if(dataIs2011){
            if(l.pho_drtotk_25_99[l.dipho_leadind[diphotonVHlep_id]] < 1 || l.pho_drtotk_25_99[l.dipho_subleadind[diphotonVHlep_id]] < 1) diphotonVHlep_id = -1;
	        VHmuevent=MuonTag2011(l, diphotonVHlep_id, &smeared_pho_energy[0]);
	        VHelevent=ElectronTag2011(l, diphotonVHlep_id, &smeared_pho_energy[0]);
	    } else {
	        VHmuevent=MuonTag2012(l, diphotonVHlep_id, &smeared_pho_energy[0],lep_sync);
	        VHelevent=ElectronTag2012(l, diphotonVHlep_id, &smeared_pho_energy[0],lep_sync);
        }
    }
	
	//Met tag //met at analysis step
	if(includeVHmet){
	    diphotonVHmet_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHmetCut, subleadEtVHmetCut, 4, false, &smeared_pho_energy[0], true);
	    if(diphotonVHmet_id>-1) {
            TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHmet_id], l.dipho_vtxind[diphoton_id] , &smeared_pho_energy[0]);
            TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHmet_id], l.dipho_vtxind[diphoton_id] , &smeared_pho_energy[0]);	    
            TLorentzVector TwoPhoton_Vector = (lead_p4) + (sublead_p4);  
            float m_gamgam = TwoPhoton_Vector.M();
                
                PhotonAnalysis::MetCorrections2012_Simple( l, lead_p4 , sublead_p4 );
                
		if (m_gamgam>100 && m_gamgam<180 && ( lead_p4.Pt() > 45*m_gamgam/120.) && sublead_p4.Pt() > 25. && l.shiftscaleMET_pt>70) {
		  met_sync << " run: " << l.run
		    << "\tevent: " << l.event
		    << "\tleadPt: " << lead_p4.Pt()
		    << "\tsubleadPt: " << sublead_p4.Pt()
		    << "\tdiphomass: " << m_gamgam
		    << "\traw_met: " << l.met_pfmet
		    << "\traw_met_phi: " << l.met_phi_pfmet
		    << "\tshifted_met: " << l.shiftMET_pt
		    << "\tcorrected_met: " << l.shiftscaleMET_pt
		    << "\tcorrected_met_phi: " << l.shiftscaleMET_phi
		    << "\tjet_algoPF1_n: " << l.jet_algoPF1_n
		    << endl;
	        }		
	    VHmetevent=METTag2012(l, diphotonVHmet_id, &smeared_pho_energy[0]);
	   }
	}
	
	// VBF+hadronic VH
	if((includeVBF || includeVHhad)&&l.jet_algoPF1_n>1 && !isSyst /*avoid rescale > once*/) {
	    l.RescaleJetEnergy();
	}

	if( !( VHmuevent || VHelevent ) && (includeVBF || includeVHhad)) {
	    if(includeVBF) {
		diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,false, &smeared_pho_energy[0], true); 
		
		float eventweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
		float myweight=1.;
		if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
		
		VBFevent= ( dataIs2011 ? 
			    VBFTag2011(l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight) :
			    VBFTag2012(ijet1, ijet2, l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight) )
		    ;
	    }
	    if(includeVHhad) {
		diphotonVHhad_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHhadCut, subleadEtVHhadCut, 4,false, &smeared_pho_energy[0], true); 
		
		float eventweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
		float myweight=1.;
		if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
		
		VHhadevent = VHhadronicTag2011(l, diphotonVHhad_id, &smeared_pho_energy[0], true, eventweight, myweight);
	    }
	}
	
	// priority of analysis:  lepton tag, vbf, VH hadronic
	if(includeVHlep && (VHelevent || VHmuevent)) {
	    diphoton_id = diphotonVHlep_id;
	} else if(includeVBF&&VBFevent) {
	    diphoton_id = diphotonVBF_id;	
	} else if(includeVHmet&&VHmetevent) {
	    diphoton_id = diphotonVHmet_id;
	} else if(includeVHhad&&VHhadevent) {
	    diphoton_id = diphotonVHhad_id;
	}
	// End exclusive mode selection
    }
    
    //// std::cout << isSyst << " " << diphoton_id << " " << sumaccept << std::endl;
    
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
        bool defaultvtx=false;
        if(( (includeVHlep && (VHelevent || VHmuevent))) && !(includeVBF&&VBFevent) ) defaultvtx=true;
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id, defaultvtx);  
      
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
	    fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, diphoton_id,  
			     category, isCorrectVertex, evweight, l );
	}

	// see if the event falls into an exclusive category
	computeExclusiveCategory(l, category, diphoton_index, Higgs.Pt() );
  
        if (dumpAscii && !isSyst && (cur_type==0||dumpMcAscii) && mass>=massMin && mass<=massMax ) {
	    
	    if( unblind ) {
		eventListText << "run:" << l.run 
			      << "\tlumi:" << l.lumis 
			      << "\tevent:" << l.event 
			      << "\tcat:" << category
			      << "\tscEta1:" << ((TVector3*)l.sc_xyz->At(l.pho_scind[l.dipho_leadind[diphoton_id]]))->Eta()
			      << "\tscEta2:" << ((TVector3*)l.sc_xyz->At(l.pho_scind[l.dipho_subleadind[diphoton_id]]))->Eta()
			      << "\tmass:" << Higgs.M()
			      << std::endl;
	    } else { 
		eventListText << "run:" << l.run 
			      << "\tlumi:" << l.lumis 
			      << "\tevent:" << l.event 
			      << "\trho:" << l.rho_algo1 
			      << "\tenergy1:" << lead_p4.Energy() 
			      << "\tenergy2:" << sublead_p4.Energy()
			      << "\tscEta1:" << ((TVector3*)l.sc_xyz->At(l.pho_scind[l.dipho_leadind[diphoton_id]]))->Eta()
			      << "\tscEta2:" << ((TVector3*)l.sc_xyz->At(l.pho_scind[l.dipho_subleadind[diphoton_id]]))->Eta()
			      << "\tr91:" << l.pho_r9[l.dipho_leadind[diphoton_id]]
			      << "\tr92:" << l.pho_r9[l.dipho_subleadind[diphoton_id]]
	    	;
	    }
	    vtxAna_.setPairID(diphoton_id);
	    std::vector<int> & vtxlist = l.vtx_std_ranked_list->at(diphoton_id);
	    for(size_t ii=0; ii<3; ++ii ) {
		eventListText << "\tvertexId"<< ii+1 <<":" << (ii < vtxlist.size() ? vtxlist[ii] : -1);
	    }
	    for(size_t ii=0; ii<3; ++ii ) {
		eventListText << "\tvertexMva"<< ii+1 <<":" << (ii < vtxlist.size() ? vtxAna_.mva(vtxlist[ii]) : -2.);
	    }
	    eventListText << "\tptbal:"   << vtxAna_.ptbal(0)
	    		  << "\tptasym:"  << vtxAna_.ptasym(0)
	    		  << "\tlogspt2:" << vtxAna_.logsumpt2(0)
	    		  << "\tp2conv:"  << vtxAna_.pulltoconv(0)
	    	;
	    dumpPhoton(eventListText,1,l,l.dipho_leadind[diphoton_id],l.dipho_vtxind[diphoton_id],lead_p4,&smeared_pho_energy[0]);
	    dumpPhoton(eventListText,2,l,l.dipho_subleadind[diphoton_id],l.dipho_vtxind[diphoton_id],sublead_p4,&smeared_pho_energy[0]);
	    if( VBFevent ) {
	    	eventListText << "\tnvtx:" << l.vtx_std_n 
	    		      << "\tjetPt1:"  << ( (TLorentzVector*)l.jet_algoPF1_p4->At(ijet1) )->Pt()
	    		      << "\tjetPt2:"  << ( (TLorentzVector*)l.jet_algoPF1_p4->At(ijet2) )->Pt()
	    		      << "\tjetEta1:" << ( (TLorentzVector*)l.jet_algoPF1_p4->At(ijet1) )->Eta()
	    		      << "\tjetEta2:" << ( (TLorentzVector*)l.jet_algoPF1_p4->At(ijet2) )->Eta()
	    	    ;
	    	dumpJet(eventListText,1,l,ijet1);
	    	dumpJet(eventListText,2,l,ijet2);
	    }
	    eventListText << std::endl;
	}
	
	return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA, 
				    int category, float weight, bool isCorrectVertex) 
{


    if (cur_type == 0 ){
            l.rooContainer->InputDataPoint("data_mass",category,mass);
    }
    if (cur_type > 0 && cur_type != 3 && cur_type != 4) {
            l.rooContainer->InputDataPoint("bkg_mass",category,mass,weight);
	}
    else if (cur_type < 0){
            l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type),category,mass,weight);
            if (isCorrectVertex) l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_rv",category,mass,weight);
            else l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_wv",category,mass,weight);
    }
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::AccumulateSyst(int cur_type, float mass, float diphotonMVA, 
				  int category, float weight,
				  std::vector<double> & mass_errors,
				  std::vector<double> & mva_errors,
				  std::vector<int>    & categories,
				  std::vector<double> & weights)
{
    categories.push_back(category);
    mass_errors.push_back(mass);
    weights.push_back(weight);
}


// ----------------------------------------------------------------------------------------------------
void StatAnalysis::FillRooContainerSyst(LoopAll& l, const std::string &name, int cur_type,
			  std::vector<double> & mass_errors, std::vector<double> & mva_errors,
			  std::vector<int>    & categories, std::vector<double> & weights) 
{	
    if (cur_type < 0){
	// feed the modified signal model to the RooContainer
	l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),name,categories,mass_errors,weights);
    }
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::computeExclusiveCategory(LoopAll & l, int & category, std::pair<int,int> diphoton_index, float pt)
{
    if(VBFevent)        {
    category=nInclusiveCategories_ + 
	    l.DiphotonCategory(diphoton_index.first,diphoton_index.second,pt,nVBFEtaCategories,1,1) 
	    + nVBFEtaCategories*l.DijetSubCategory(myVBF_Mjj,myVBFLeadJPt,myVBFSubJPt,nVBFDijetJetCategories)
	    ;
    } else if(VHhadevent) { category=nInclusiveCategories_ + ( (int)includeVBF )*nVBFEtaCategories
	    + l.DiphotonCategory(diphoton_index.first,diphoton_index.second,pt,nVHhadEtaCategories,1,1); 
    } else if(VHmuevent) {
	category=nInclusiveCategories_ + ( (int)includeVBF )*nVBFEtaCategories + ( (int)includeVHhad )*nVHhadEtaCategories;  
    } else if(VHelevent) { 
	category=nInclusiveCategories_ + ( (int)includeVBF )*nVBFEtaCategories + ( (int)includeVHhad )*nVHhadEtaCategories + (int)includeVHlep;
    } else if(VHmetevent) { 
	category=nInclusiveCategories_ + ( (int)includeVBF )*nVBFEtaCategories + ( (int)includeVHhad )*nVHhadEtaCategories + (int)includeVHlep + (int)includeVHmet;  //met at analysis step
    }
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs, 
				    float lead_r9, float sublead_r9, int diphoton_id,
				    int category, bool isCorrectVertex, float evweight, LoopAll & l )
{       
    // control plots 
    if( category>=0 ) { 
	fillControlPlots( lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, diphoton_id, -1, isCorrectVertex, evweight, l ); 
    }
    float mass = Higgs.M();
    l.FillHist("all_mass",category+1, Higgs.M(), evweight);
    if( mass>=massMin && mass<=massMax  ) {
	l.FillHist("mass",category+1, Higgs.M(), evweight);
	l.FillHist("eta",category+1, Higgs.Eta(), evweight);
	l.FillHist("pt",category+1, Higgs.Pt(), evweight);
	if( isCorrectVertex ) { l.FillHist("pt_rv",category+1, Higgs.Pt(), evweight); }
	l.FillHist("nvtx",category+1, l.vtx_std_n, evweight);
        if( isCorrectVertex ) { l.FillHist("nvtx_rv",category+1, l.vtx_std_n, evweight); }
	
	vtxAna_.setPairID(diphoton_id);
	float vtxProb = vtxAna_.vertexProbability(l.vtx_std_evt_mva->at(diphoton_id), l.vtx_std_n);
	l.FillHist2D("probmva_pt",category+1, Higgs.Pt(), l.vtx_std_evt_mva->at(diphoton_id), evweight);
	l.FillHist2D("probmva_nvtx",category+1, l.vtx_std_n, l.vtx_std_evt_mva->at(diphoton_id), evweight);
	if( isCorrectVertex ) { l.FillHist2D("probmva_rv_nvtx",category+1, l.vtx_std_n, l.vtx_std_evt_mva->at(diphoton_id), evweight); }
	l.FillHist2D("vtxprob_pt",category+1, Higgs.Pt(), vtxProb, evweight);
	l.FillHist2D("vtxprob_nvtx",category+1, l.vtx_std_n, vtxProb, evweight);
	std::vector<int> & vtxlist = l.vtx_std_ranked_list->at(diphoton_id);
	size_t maxv = std::min(vtxlist.size(),(size_t)5);
	for(size_t ivtx=0; ivtx<maxv; ++ivtx) {
	    int vtxid = vtxlist.at(ivtx);
	    l.FillHist(Form("vtx_mva_%d",ivtx),category+1,vtxAna_.mva(ivtx),evweight);
	    if( ivtx > 0 ) {
		l.FillHist(Form("vtx_dz_%d",ivtx),category+1,
			   vtxAna_.vertexz(ivtx)-vtxAna_.vertexz(l.dipho_vtxind[diphoton_id]),evweight);
	    }
	}
	l.FillHist("vtx_nconv",vtxAna_.nconv(0));

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
	l.FillHist("pho2_r9",category+1, sublead_r9, evweight);

	l.FillHist("pho_n",category+1,l.pho_n, evweight);
	
	l.FillHist("pho_rawe",category+1,l.sc_raw[l.pho_scind[l.dipho_leadind[diphoton_id]]], evweight);
	l.FillHist("pho_rawe",category+1,l.sc_raw[l.pho_scind[l.dipho_subleadind[diphoton_id]]], evweight);
    }
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::fillSignalEfficiencyPlots(float weight, LoopAll & l)
{
    //Fill histograms to use as denominator (kinematic pre-selection only) and numerator (selection applied)
    //for single photon ID efficiency calculation.
    int diphoton_id_kinpresel = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,-1.,applyPtoverM, &smeared_pho_energy[0], true ); 
    if (diphoton_id_kinpresel>-1) {

	TLorentzVector lead_p4, sublead_p4, Higgs;
	float lead_r9, sublead_r9;
	TVector3 * vtx;
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id_kinpresel);  

	int ivtx = l.dipho_vtxind[diphoton_id_kinpresel];
	int lead = l.dipho_leadind[diphoton_id_kinpresel];
	int sublead = l.dipho_subleadind[diphoton_id_kinpresel];
	int leadpho_category = l.PhotonCategory(lead, 2, 2);
	int subleadpho_category = l.PhotonCategory(sublead, 2, 2);
	float leadEta = ((TVector3 *)l.sc_xyz->At(l.pho_scind[lead]))->Eta();
	float subleadEta = ((TVector3 *)l.sc_xyz->At(l.pho_scind[sublead]))->Eta();

	float evweight = weight * smeared_pho_weight[lead] * smeared_pho_weight[sublead] * genLevWeight;

	//Fill eta and pt distributions after pre-selection only (efficiency denominator)
	l.FillHist("pho1_pt_presel",0,lead_p4.Pt(), evweight);
	l.FillHist("pho2_pt_presel",0,sublead_p4.Pt(), evweight);
	l.FillHist("pho1_eta_presel",0,leadEta, evweight);
	l.FillHist("pho2_eta_presel",0,subleadEta, evweight);

	l.FillHist("pho1_pt_presel",leadpho_category+1,lead_p4.Pt(), evweight);
	l.FillHist("pho2_pt_presel",subleadpho_category+1,sublead_p4.Pt(), evweight);
	l.FillHist("pho1_eta_presel",leadpho_category+1,leadEta, evweight);
	l.FillHist("pho2_eta_presel",subleadpho_category+1,subleadEta, evweight);

	//Apply single photon CiC selection and fill eta and pt distributions (efficiency numerator)
	std::vector<std::vector<bool> > ph_passcut;
	if( l.PhotonCiCSelectionLevel(lead, ivtx, ph_passcut, 4, 0, &smeared_pho_energy[0]) >=  (LoopAll::phoCiCIDLevel) l.phoSUPERTIGHT) {
	    l.FillHist("pho1_pt_sel",0,lead_p4.Pt(), evweight);
	    l.FillHist("pho1_eta_sel",0,leadEta, evweight);
	    l.FillHist("pho1_pt_sel",leadpho_category+1,lead_p4.Pt(), evweight);
	    l.FillHist("pho1_eta_sel",leadpho_category+1,leadEta, evweight);
	}
	if( l.PhotonCiCSelectionLevel(sublead, ivtx, ph_passcut, 4, 1, &smeared_pho_energy[0]) >=  (LoopAll::phoCiCIDLevel) l.phoSUPERTIGHT ) {
	    l.FillHist("pho2_pt_sel",0,sublead_p4.Pt(), evweight);
	    l.FillHist("pho2_eta_sel",0,subleadEta, evweight);
	    l.FillHist("pho2_pt_sel",subleadpho_category+1,sublead_p4.Pt(), evweight);
	    l.FillHist("pho2_eta_sel",subleadpho_category+1,subleadEta, evweight);
	}
    }
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool StatAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}
// ----------------------------------------------------------------------------------------------------
double StatAnalysis::GetDifferentialKfactor(double gPT, int Mass)
{
    return 1.0;
}

void StatAnalysis::FillSignalLabelMap(LoopAll & l)
{
    std::map<int,std::pair<TString,double > > & signalMap = l.signalNormalizer->SignalType();
    
    for( std::map<int,std::pair<TString,double > >::iterator it=signalMap.begin();
	 it!=signalMap.end(); ++it ) {
	signalLabels[it->first] = it->second.first+Form("_mass_m%1.0f", it->second.second);
    }
    
    /////////// // Basically A Map of the ID (type) to the signal's name which can be filled Now:
    /////////// signalLabels[-57]="ggh_mass_m123";
    /////////// signalLabels[-58]="vbf_mass_m123";
    /////////// signalLabels[-60]="wzh_mass_m123";
    /////////// signalLabels[-59]="tth_mass_m123";
    /////////// signalLabels[-53]="ggh_mass_m121";
    /////////// signalLabels[-54]="vbf_mass_m121";
    /////////// signalLabels[-56]="wzh_mass_m121";
    /////////// signalLabels[-55]="tth_mass_m121";
    /////////// signalLabels[-65]="ggh_mass_m160";
    /////////// signalLabels[-66]="vbf_mass_m160";
    /////////// signalLabels[-68]="wzh_mass_m160";
    /////////// signalLabels[-67]="tth_mass_m160";
    /////////// signalLabels[-61]="ggh_mass_m155";
    /////////// signalLabels[-62]="vbf_mass_m155";
    /////////// signalLabels[-64]="wzh_mass_m155";
    /////////// signalLabels[-63]="tth_mass_m155";
    /////////// signalLabels[-49]="ggh_mass_m150";
    /////////// signalLabels[-50]="vbf_mass_m150";
    /////////// signalLabels[-52]="wzh_mass_m150";
    /////////// signalLabels[-51]="tth_mass_m150";
    /////////// signalLabels[-45]="ggh_mass_m145";
    /////////// signalLabels[-46]="vbf_mass_m145";
    /////////// signalLabels[-48]="wzh_mass_m145";
    /////////// signalLabels[-47]="tth_mass_m145";
    /////////// signalLabels[-33]="ggh_mass_m140";
    /////////// signalLabels[-34]="vbf_mass_m140";
    /////////// signalLabels[-36]="wzh_mass_m140";
    /////////// signalLabels[-35]="tth_mass_m140";
    /////////// signalLabels[-41]="ggh_mass_m135";
    /////////// signalLabels[-42]="vbf_mass_m135";
    /////////// signalLabels[-44]="wzh_mass_m135";
    /////////// signalLabels[-43]="tth_mass_m135";
    /////////// signalLabels[-29]="ggh_mass_m130";
    /////////// signalLabels[-30]="vbf_mass_m130";
    /////////// signalLabels[-32]="wzh_mass_m130";
    /////////// signalLabels[-31]="tth_mass_m130";
    /////////// signalLabels[-37]="ggh_mass_m125";
    /////////// signalLabels[-38]="vbf_mass_m125";
    /////////// signalLabels[-40]="wzh_mass_m125";
    /////////// signalLabels[-39]="tth_mass_m125";
    /////////// signalLabels[-25]="ggh_mass_m120";
    /////////// signalLabels[-26]="vbf_mass_m120";
    /////////// signalLabels[-28]="wzh_mass_m120";
    /////////// signalLabels[-27]="tth_mass_m120";
    /////////// signalLabels[-21]="ggh_mass_m115";
    /////////// signalLabels[-22]="vbf_mass_m115";
    /////////// signalLabels[-24]="wzh_mass_m115";
    /////////// signalLabels[-23]="tth_mass_m115";
    /////////// signalLabels[-17]="ggh_mass_m110";
    /////////// signalLabels[-18]="vbf_mass_m110";
    /////////// signalLabels[-19]="wzh_mass_m110";
    /////////// signalLabels[-20]="tth_mass_m110";
    /////////// signalLabels[-13]="ggh_mass_m105";
    /////////// signalLabels[-14]="vbf_mass_m105";
    /////////// signalLabels[-16]="wzh_mass_m105";
    /////////// signalLabels[-15]="tth_mass_m105";
    /////////// signalLabels[-69]="ggh_mass_m100";
    /////////// signalLabels[-70]="vbf_mass_m100";
    /////////// signalLabels[-72]="wzh_mass_m100";
    /////////// signalLabels[-71]="tth_mass_m100";
}

std::string StatAnalysis::GetSignalLabel(int id){
    
    // For the lazy man, can return a memeber of the map rather than doing it yourself
    std::map<int,std::string>::iterator it = signalLabels.find(id);

    if (it!=signalLabels.end()){
        return it->second;
        
    } else { 

        std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
        return "NULL";
    }
    
}

void StatAnalysis::rescaleClusterVariables(LoopAll &l){

    // Data-driven MC scalings 
    for (int ipho=0;ipho<l.pho_n;ipho++){

	if (dataIs2011) {

	    if( scaleR9Only ) {
		double R9_rescale = (l.pho_isEB[ipho]) ? 1.0048 : 1.00492 ;
		l.pho_r9[ipho]*=R9_rescale;
	    } else {
		l.pho_r9[ipho]*=1.0035;
		if (l.pho_isEB[ipho]){ l.pho_sieie[ipho] = (0.87*l.pho_sieie[ipho]) + 0.0011 ;}
	    	else {l.pho_sieie[ipho]*=0.99;}
		l.sc_seta[l.pho_scind[ipho]]*=0.99;  
		l.sc_sphi[l.pho_scind[ipho]]*=0.99;  
	    }

	} else {
	    //2012 rescaling from here https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/752/1/1/2/1/3.html

	    if (l.pho_isEB[ipho]) {
		l.pho_r9[ipho] = 1.0045*l.pho_r9[ipho] + 0.0010;
	    } else {
		l.pho_r9[ipho] = 1.0086*l.pho_r9[ipho] - 0.0007;
	    }
	    if( !scaleR9Only ) {
		if (l.pho_isEB[ipho]) {
		    l.pho_s4ratio[ipho] = 1.01894*l.pho_s4ratio[ipho] - 0.01034;
		    l.pho_sieie[ipho] = 0.891832*l.pho_sieie[ipho] + 0.0009133;
		    l.pho_etawidth[ipho] =  1.04302*l.pho_etawidth[ipho] - 0.000618;
		    l.sc_sphi[l.pho_scind[ipho]] =  1.00002*l.sc_sphi[l.pho_scind[ipho]] - 0.000371;
		} else {
		    l.pho_s4ratio[ipho] = 1.04969*l.pho_s4ratio[ipho] - 0.03642;
		    l.pho_sieie[ipho] = 0.99470*l.pho_sieie[ipho] + 0.00003;
		    l.pho_etawidth[ipho] =  0.903254*l.pho_etawidth[ipho] + 0.001346;
		    l.sc_sphi[l.pho_scind[ipho]] =  0.99992*l.sc_sphi[l.pho_scind[ipho]] - 0.00000048;
		    //Agreement not to rescale ES shape (https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/789/1/1/1/1/1/1/2/1/1.html)
		    //if (l.pho_ESEffSigmaRR[ipho]>0) l.pho_ESEffSigmaRR[ipho] = 1.00023*l.pho_ESEffSigmaRR[ipho] + 0.0913;
		}
	    }
	}
    // Scale DYJets sample for now
    /*
    if (l.itype[l.current]==6){
    if (l.pho_isEB[ipho]) {
        energyCorrectedError[ipho] = 1.02693*energyCorrectedError[ipho]-0.0042793;
    } else {
        energyCorrectedError[ipho] = 1.01372*energyCorrectedError[ipho]+0.000156943;
    }
    }
    */
    }
}

void StatAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    if( doEresolSmear ) {
	eResolSmearer->resetRandom();
    }
}

void dumpJet(std::ostream & eventListText, int lab, LoopAll & l, int ijet)
{
    eventListText << std::setprecision(4) << std::scientific
		  << "\tjec"      << lab << ":" << l.jet_algoPF1_erescale[ijet]
		  << "\tbetaStar" << lab << ":" << l.jet_algoPF1_betaStarClassic[ijet]
		  << "\tRMS"      << lab << ":" << l.jet_algoPF1_dR2Mean[ijet]
	;
}

void dumpPhoton(std::ostream & eventListText, int lab, 
		LoopAll & l, int ipho, int ivtx, TLorentzVector & phop4, float * pho_energy_array)
{
    float val_tkisobad = -99;
    for(int iv=0; iv < l.vtx_std_n; iv++) {
	if((*l.pho_pfiso_mycharged04)[ipho][iv] > val_tkisobad) {
	    val_tkisobad = (*l.pho_pfiso_mycharged04)[ipho][iv];
	}
    }
    TLorentzVector phop4_badvtx = l.get_pho_p4( ipho, l.pho_tkiso_badvtx_id[ipho], pho_energy_array  );

    float val_tkiso        = (*l.pho_pfiso_mycharged03)[ipho][ivtx];
    float val_ecaliso      = l.pho_pfiso_myphoton03[ipho];
    float val_ecalisobad   = l.pho_pfiso_myphoton04[ipho];
    float val_sieie        = l.pho_sieie[ipho];
    float val_hoe          = l.pho_hoe[ipho];
    float val_r9           = l.pho_r9[ipho];
    float val_conv         = l.pho_isconv[ipho];
    
    float rhofacbad=0.23, rhofac=0.09;
    
    float val_isosumoet    = (val_tkiso    + val_ecaliso    - l.rho_algo1 * rhofac )   * 50. / phop4.Et();
    float val_isosumoetbad = (val_tkisobad + val_ecalisobad - l.rho_algo1 * rhofacbad) * 50. / phop4_badvtx.Et();
    
    // tracker isolation cone energy divided by Et
    float val_trkisooet    = (val_tkiso) * 50. / phop4.Pt();

    eventListText << std::setprecision(4) << std::scientific
		  << "\tchIso03"  << lab << ":" << val_tkiso
		  << "\tphoIso03" << lab << ":" << val_ecaliso 
		  << "\tchIso04"  << lab << ":" << val_tkisobad
		  << "\tphoIso04" << lab << ":" << val_ecalisobad 
		  << "\tsIeIe"    << lab << ":" << val_sieie
		  << "\thoe"      << lab << ":" << val_hoe
		  << "\tecalIso"  << lab << ":" << l.pho_ecalsumetconedr03[ipho]
		  << "\thcalIso"  << lab << ":" << l.pho_hcalsumetconedr03[ipho]
		  << "\ttrkIso"   << lab << ":" << l.pho_trksumpthollowconedr03[ipho]
		  << "\tchIso02"  << lab << ":" << (*l.pho_pfiso_mycharged02)[ipho][ivtx]
		  << "\teleVeto"  << lab << ":" << !val_conv
	;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
