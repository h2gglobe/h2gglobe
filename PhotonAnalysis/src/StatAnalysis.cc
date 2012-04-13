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

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
    doSystematics = true;    
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
  
    //    l.rooContainer->WriteSpecificCategoryDataCards(outputfilename,"data_mass","sig_mass","data_pol_model");
    //    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","data_pol_model");
    // mode 0 as above, 1 if want to bin in sub range from fit,

    // Write the data-card for the Combinations Code, needs the output filename, makes binned analysis DataCard
    // Assumes the signal datasets will be called signal_name+"_mXXX"
    //    l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","data_pol_model","data_mass",1,50,1); // 1 means systematics from the fit effect only the backgroundi. last digit mode = 1 means this is an internal constraint fit 
    //    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","bkg_mass_rebinned");

    eventListText.close();

    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;

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
    
    std::string outputfilename = (std::string) l.histFileName;
    eventListText.open(Form("%s",l.outputTextFileName.c_str()));
    //eventListText.open(Form("%s_ascii_events.txt",outputfilename.c_str()));
    FillSignalLabelMap();
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

    // avoid recalculated the CIC ID every time
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
    
    int nVBFCategories  = ((int)includeVBF)*nVBFEtaCategories;
    int nVHhadCategories = ((int)includeVHhad)*nVHhadEtaCategories;
    int nVHlepCategories = (int)includeVHlep * 2;
   
    nCategories_=(nInclusiveCategories_+nVBFCategories+nVHhadCategories+nVHlepCategories);

    

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
    l.rooContainer->AddRealVar("CMS_hgg_pol0",-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol1",-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol2",-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol3",-0.01,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_pol4",-0.01,-1.0,1.0);
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol0","@0*@0","CMS_hgg_pol0");
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol1","@0*@0","CMS_hgg_pol1");
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol2","@0*@0","CMS_hgg_pol2");
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol3","@0*@0","CMS_hgg_pol3");
    l.rooContainer->AddFormulaVar("CMS_hgg_modpol4","@0*@0","CMS_hgg_pol4");

    l.rooContainer->AddRealVar("CMS_hgg_quartic0",-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_quartic1",-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_quartic2",-0.1,-1.0,1.0);
    l.rooContainer->AddRealVar("CMS_hgg_quartic3",-0.01,-1.0,1.0);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic0","@0*@0","CMS_hgg_quartic0");
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic1","@0*@0","CMS_hgg_quartic1");
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic2","@0*@0","CMS_hgg_quartic2");
    l.rooContainer->AddFormulaVar("CMS_hgg_modquartic3","@0*@0","CMS_hgg_quartic3");

    l.rooContainer->AddRealVar("CMS_hgg_quad0",-0.1,-1.5,1.5);
    l.rooContainer->AddRealVar("CMS_hgg_quad1",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("CMS_hgg_modquad0","@0*@0","CMS_hgg_quad0");
    l.rooContainer->AddFormulaVar("CMS_hgg_modquad1","@0*@0","CMS_hgg_quad1");
    
    l.rooContainer->AddRealVar("CMS_hgg_cubic0",-0.1,-1.5,1.5);
    l.rooContainer->AddRealVar("CMS_hgg_cubic1",-0.1,-1.5,1.5);
    l.rooContainer->AddRealVar("CMS_hgg_cubic2",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("CMS_hgg_modcubic0","@0*@0","CMS_hgg_cubic0");
    l.rooContainer->AddFormulaVar("CMS_hgg_modcubic1","@0*@0","CMS_hgg_cubic1");
    l.rooContainer->AddFormulaVar("CMS_hgg_modcubic2","@0*@0","CMS_hgg_cubic2");
    
    l.rooContainer->AddRealVar("CMS_hgg_lin0",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("CMS_hgg_modlin0","@0*@0","CMS_hgg_lin0");
    
    // Generic PDF ok in the std analysis but excluisve channels need different models CP
    //l.rooContainer->AddGenericPdf("data_pol_model",
    //"0","CMS_hgg_mass",data_pol_pars,73); // >= 71 means RooBernstein of order >= 1
    int cats_with_std[]     = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int cats_with_lin[]     = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int cats_with_quad[]    = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int cats_with_cubic[]   = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int cats_with_quartic[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    for(int i=0; i<nCategories_; i++){
        if(i<nInclusiveCategories_) {
            cats_with_std[i]=1;
            // mva removed cp march 8
            //if(useMVA) {
            //  if(i<5) {
            //    cats_with_std[i]=1;
            //  } else {
            //    cats_with_quartic[i]=1;
            //  }
            //} else {
            //  cats_with_std[i]=1;
            //}
        } else if(i<nInclusiveCategories_+((int)includeVBF)*nVBFEtaCategories){
            cats_with_quad[i]=1;
        } else if(i<nInclusiveCategories_+((int)includeVBF)*nVBFEtaCategories+((int)includeVHhad)*nVHhadEtaCategories){
            cats_with_quad[i]=1;
        } else {
            cats_with_lin[i]=1;
        }  
    }



    std::vector<std::string> data_pol_pars(5,"p");   
    data_pol_pars[0] = "CMS_hgg_modpol0";
    data_pol_pars[1] = "CMS_hgg_modpol1";
    data_pol_pars[2] = "CMS_hgg_modpol2";
    data_pol_pars[3] = "CMS_hgg_modpol3";
    data_pol_pars[4] = "CMS_hgg_modpol4";
    
    l.rooContainer->AddSpecificCategoryPdf(cats_with_std,"data_pol_model",
                                           "0","CMS_hgg_mass",data_pol_pars,75);    // >= 71 means RooBernstein of order >= 1
    
    std::vector<std::string> data_quartic_pars(4,"p");   
    data_quartic_pars[0] = "CMS_hgg_modquartic0";
    data_quartic_pars[1] = "CMS_hgg_modquartic1";
    data_quartic_pars[2] = "CMS_hgg_modquartic2";
    data_quartic_pars[3] = "CMS_hgg_modquartic3";
    
    l.rooContainer->AddSpecificCategoryPdf(cats_with_quartic,"data_pol_model",
                                           "0","CMS_hgg_mass",data_quartic_pars,74);    // >= 71 means RooBernstein of order >= 1

    std::vector<std::string> data_cubic_pars(3,"p");     
    data_cubic_pars[0] = "CMS_hgg_modcubic0";
    data_cubic_pars[1] = "CMS_hgg_modcubic1";
    data_cubic_pars[2] = "CMS_hgg_modcubic2";
    
    l.rooContainer->AddSpecificCategoryPdf(cats_with_cubic, "data_pol_model",
                                           "0","CMS_hgg_mass",data_cubic_pars,73);  // >= 71 means RooBernstein of order >= 1
    

    std::vector<std::string> data_quad_pars(2,"p");  
    data_quad_pars[0] = "CMS_hgg_modquad0";
    data_quad_pars[1] = "CMS_hgg_modquad1";
    
    l.rooContainer->AddSpecificCategoryPdf(cats_with_quad, "data_pol_model",
                                           "0","CMS_hgg_mass",data_quad_pars,72);   // >= 71 means RooBernstein of order >= 1
    
    
    std::vector<std::string> data_lin_pars(1,"p");   
    data_lin_pars[0] = "CMS_hgg_modlin0";
    
    l.rooContainer->AddSpecificCategoryPdf(cats_with_lin, "data_pol_model",
                                           "0","CMS_hgg_mass",data_lin_pars,71);    // >= 71 means RooBernstein of order >= 1


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
    // additional TTH155 needed for interpolation -> sigh!

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
    float sampleweight = l.sampleContainer[l.current_sample_index].weight;

    // Set reRunCiC Only if this is an MC event since scaling of R9 and Energy isn't done at reduction
    if (cur_type==0) {
        l.runCiC=reRunCiCForData;
    } else {
        l.runCiC = true;
    } 

    l.FillCounter( "Processed", 1. );
    assert( weight > 0. );  
    l.FillCounter( "XSWeighted", weight );
    nevents+=1.;

    //PU reweighting
    double pileupWeight=getPuWeight( l.pu_n, cur_type, jentry == 1);
    sumwei +=pileupWeight;
    weight *= pileupWeight;
    sumev  += weight;
    
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
    if (cur_type<0){
	gP4 = l.GetHiggs();
	gPT = gP4.Pt();
    }

    // TEMPORARY FIX -------------------------------------------------------------------------------------------------------//
    // Scale all the r9 of the photons in the MC
    // For now we just let it use the index but we specifically Change the r9 in the branch AFTER Energy regression smearing
    // Ideally we want to pass a smeared r9 too and apply after energy corrections, currently the smeared_pho_r9 isnt used!
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    if (cur_type !=0){
        for (int ipho=0;ipho<l.pho_n;ipho++){
            double R9_rescale = (l.pho_isEB[ipho]) ? 1.0048 : 1.00492 ;
            l.pho_r9[ipho]*=R9_rescale;
        }
    }
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//

    // Analyse the event assuming nominal values of corrections and smearings
    float mass, evweight;
    int diphoton_id, category;
    bool isCorrectVertex;
    if( AnalyseEvent(l,jentry, weight, gP4, mass,  evweight, category, diphoton_id, isCorrectVertex) ) {
	// feed the event to the RooContainer 
        if (cur_type == 0 ){
            l.rooContainer->InputDataPoint("data_mass",category,mass);
        }
        if (cur_type > 0 && cur_type != 3 && cur_type != 4) {
            l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
	}
        else if (cur_type < 0){
            l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type),category,mass,evweight);
            if (isCorrectVertex) l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_rv",category,mass,evweight);
            else l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_wv",category,mass,evweight);
        }
    }
    
    // Systematics uncertaities for the binned model
    // We re-analyse the event several times for different values of corrections and smearings
    if( cur_type < 0 && doMCSmearing && doSystematics ) { 
        
        // fill steps for syst uncertainty study
        float systStep = systRange / (float)nSystSteps;
	
	float syst_mass, syst_weight;
	int syst_category;
 	std::vector<double> mass_errors;
	std::vector<double> weights;
	std::vector<int>    categories;
	
        if (diphoton_id > -1 ) {
     
	    // gen-level systematics, i.e. ggH k-factor for the moment
            for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
		mass_errors.clear(), weights.clear(), categories.clear();
		
                for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
                    if( syst_shift == 0. ) { continue; } // skip the central value
		    syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
		    
		    // re-analyse the event without redoing the event selection as we use nominal values for the single photon
		    // corrections and smearings
		    AnalyseEvent(l, jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,
				 true, syst_shift, true, *si, 0, 0 );
		    		    
		    categories.push_back(syst_category);
		    mass_errors.push_back(syst_mass);
                    weights.push_back(syst_weight);
		}
		if (cur_type < 0){
		    // feed the modified signal model to the RooContainer
		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
		}
	    }
	    
	    // di-photon systematics: vertex efficiency and trigger 
	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		mass_errors.clear(), weights.clear(), categories.clear();
		
                for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
                    if( syst_shift == 0. ) { continue; } // skip the central value
		    syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
		    
		    // re-analyse the event without redoing the event selection as we use nominal values for the single photon
		    // corrections and smearings
		    AnalyseEvent(l,jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,
				 true, syst_shift, true,  0, 0, *si );
		    
		    categories.push_back(syst_category);
		    mass_errors.push_back(syst_mass);
                    weights.push_back(syst_weight);
		}
		if (cur_type < 0){
		    // feed the modified signal model to the RooContainer
		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
		}
	    }
	}
	
	// single photon level systematics: several
	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	    mass_errors.clear(), weights.clear(), categories.clear();
	    
	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		if( syst_shift == 0. ) { continue; } // skip the central value
		syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
		
		// re-analyse the event redoing the event selection this time
		AnalyseEvent(l,jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,
			     true, syst_shift, false,  0, *si, 0 );
		
		categories.push_back(syst_category);
		mass_errors.push_back(syst_mass);
		weights.push_back(syst_weight);
	    }
	    if (cur_type < 0){
		// feed the modified signal model to the RooContainer
		l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
	    }
	}
    }

    if(PADEBUG) 
        cout<<"myFillHistRed END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
bool StatAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4,
				float & mass, float & evweight, int & category, int & diphoton_id, bool & isCorrectVertex,
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
	
	// lepton tag
	if(includeVHlep){
	    diphotonVHlep_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHlepCut, subleadEtVHlepCut, 4, false, &smeared_pho_energy[0], true );
	    //Add tighter cut on dr to tk
	    if(l.pho_drtotk_25_99[l.dipho_leadind[diphotonVHlep_id]] < 1 || l.pho_drtotk_25_99[l.dipho_subleadind[diphotonVHlep_id]] < 1) diphotonVHlep_id = -1;
	    VHmuevent=MuonTag2011(l, diphotonVHlep_id, &smeared_pho_energy[0]);
	    VHelevent=ElectronTag2011(l, diphotonVHlep_id, &smeared_pho_energy[0]);
	}
	
	// VBF+hadronic VH
	if((includeVBF || includeVHhad)&&l.jet_algoPF1_n>1) {
	    l.RescaleJetEnergy();
	}

	if( !( VHmuevent || VHelevent ) && (includeVBF || includeVHhad)) {
	    if(includeVBF) {
		diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,false, &smeared_pho_energy[0], true); 
		
		float eventweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
		float myweight=1.;
		if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
		
		VBFevent=VBFTag2011(l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight);
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
	} else if(includeVHhad&&VHhadevent) {
	    diphoton_id = diphotonVHhad_id;
	}
	// End exclusive mode selection
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
	    fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, category, evweight, l );
	}

	// see if the event falls into an exclusive category
	computeExclusiveCategory(l, category, diphoton_index, Higgs.Pt() );
  
        if (!isSyst && (cur_type==0||PADEBUG) ) {
            eventListText << "Type = "<< cur_type <<  " Run = " << l.run << "  LS = " << l.lumis << "  Event = " <<  l.event << "  ggM = " << mass 
			  <<  " CAT " << category << " Vertex = " <<  l.dipho_vtxind[diphoton_id];
            eventListText << endl;
        }
	
	return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::computeExclusiveCategory(LoopAll & l, int & category, std::pair<int,int> diphoton_index, float pt)
{
    if(VBFevent)        category=nInclusiveCategories_ + l.DiphotonCategory(diphoton_index.first,diphoton_index.second,pt,nVBFEtaCategories,1,1);
    else if(VHhadevent) category=nInclusiveCategories_ + ( (int)includeVBF )*nVBFEtaCategories
	+ l.DiphotonCategory(diphoton_index.first,diphoton_index.second,pt,nVHhadEtaCategories,1,1);
    else if(VHmuevent) category=nInclusiveCategories_ + ( (int)includeVBF )*nVBFEtaCategories + ( (int)includeVHhad )*nVHhadEtaCategories;
    else if(VHelevent) category=nInclusiveCategories_ + ( (int)includeVBF )*nVBFEtaCategories + ( (int)includeVHhad )*nVHhadEtaCategories + (int)includeVHlep;
    
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs, float lead_r9, float sublead_r9, 
				    int category, float evweight, LoopAll & l )
{       
    // control plots 
    float mass = Higgs.M();
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

void StatAnalysis::FillSignalLabelMap()
{
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


void StatAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    eResolSmearer->resetRandom();
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
