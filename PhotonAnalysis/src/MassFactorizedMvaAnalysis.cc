#include "../interface/MassFactorizedMvaAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
MassFactorizedMvaAnalysis::MassFactorizedMvaAnalysis()  : 
    name_("MassFactorizedMvaAnalysis"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
}

// ----------------------------------------------------------------------------------------------------
MassFactorizedMvaAnalysis::~MassFactorizedMvaAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void MassFactorizedMvaAnalysis::Term(LoopAll& l) 
{

    if (! l.is_subjob){ // no need to waste time when running a subjob
        std::string outputfilename = (std::string) l.histFileName;
        l.rooContainer->FitToData("data_pol_model","data_mass");  // Fit to full range of dataset
    }

    eventListText.close();
    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;

}

// ----------------------------------------------------------------------------------------------------
void MassFactorizedMvaAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
        cout << "InitRealMassFactorizedMvaAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;

    //  std::string outputfilename = (std::string) l.histFileName;
    eventListText.open(Form("%s",l.outputTextFileName.c_str()));
    FillSignalLabelMap();
    //
    // These parameters are set in the configuration file
    std::cout
        << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << "MassFactorizedMvaAnalysis " << "\n"
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

    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
        l.histoContainer[ind].setScale(1.);
    }

    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    nInclusiveCategories_ = 4;

  
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
        eResolSmearer->doEnergy(false); // allows for future reweighting also
        eResolSmearer->scaleOrSmear(false);
        photonSmearers_.push_back(eResolSmearer);
    }
    if( doRegressionSmear ) {
        // energy regression. smearing
        std::cerr << __LINE__ << std::endl; 
        eRegressionSmearer = new EnergySmearer( eSmearPars );
        eRegressionSmearer->name("regSig");
        eRegressionSmearer->doEnergy(false);// allows for future reweighting also
        eRegressionSmearer->doRegressionSigma(true);
        photonSmearers_.push_back(eRegressionSmearer);
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
    if( doPhotonMvaIdSmear ) {
        // trigger efficiency
        std::cerr << __LINE__ << std::endl; 
        photonMvaIdSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
        photonMvaIdSmearer->name("phoIdMva");
        photonMvaIdSmearer->setEffName("effL1HLT");
        photonMvaIdSmearer->doVtxEff(false);
        photonMvaIdSmearer->doMvaIdEff(true);
        photonMvaIdSmearer->init();
        diPhotonSmearers_.push_back(photonMvaIdSmearer);
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
    if (bdtTrainingPhilosophy == "UCSD") {
        l.rooContainer->SetNCategories(8);
    } else if (bdtTrainingPhilosophy == "MIT") {
        if (includeVBF)  l.rooContainer->SetNCategories(6);
        else l.rooContainer->SetNCategories(4);
    }

    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;
    l.rooContainer->SaveRooDataHists(true);
    l.rooContainer->Verbose(false);

    if( doEcorrectionSmear && doEcorrectionSyst ) {
        // instance of this smearer done in PhotonAnalysis
        systPhotonSmearers_.push_back(eCorrSmearer);
        std::vector<std::string> sys(1,eCorrSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEscaleSmear && doEscaleSyst ) {
        systPhotonSmearers_.push_back( eScaleSmearer );
        std::vector<std::string> sys(1,eScaleSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
        systPhotonSmearers_.push_back( eResolSmearer );
        std::vector<std::string> sys(1,eResolSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doRegressionSmear && doRegressionSyst ) {
        systPhotonSmearers_.push_back( eRegressionSmearer );
        std::vector<std::string> sys(1,eRegressionSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonIdEffSmear && doPhotonIdEffSyst ) {
        systPhotonSmearers_.push_back( idEffSmearer );
        std::vector<std::string> sys(1,idEffSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doR9Smear && doR9Syst ) {
        systPhotonSmearers_.push_back( r9Smearer );
        std::vector<std::string> sys(1,r9Smearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doVtxEffSmear && doVtxEffSyst ) {
        systDiPhotonSmearers_.push_back( vtxEffSmearer );
        std::vector<std::string> sys(1,vtxEffSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doTriggerEffSmear && doTriggerEffSyst ) {
        systDiPhotonSmearers_.push_back( triggerEffSmearer );
        std::vector<std::string> sys(1,triggerEffSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonMvaIdSmear && doPhotonMvaIdSyst ) {
        systDiPhotonSmearers_.push_back( photonMvaIdSmearer );
        std::vector<std::string> sys(1,photonMvaIdSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if(doKFactorSmear && doKFactorSyst) {
        systGenLevelSmearers_.push_back(kFactorSmearer);
        std::vector<std::string> sys(1,kFactorSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }

    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.045,1.00);
    // ----------------------------------------------------

    // Create observables for shape-analysis with ranges
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

    l.rooContainer->AddRealVar("pol0",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol1",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol2",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol3",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol4",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modpol0","@0*@0","pol0");
    l.rooContainer->AddFormulaVar("modpol1","@0*@0","pol1");
    l.rooContainer->AddFormulaVar("modpol2","@0*@0","pol2");
    l.rooContainer->AddFormulaVar("modpol3","@0*@0","pol3");
    l.rooContainer->AddFormulaVar("modpol4","@0*@0","pol4");

    if (bdtTrainingPhilosophy=="UCSD"){
        // UCSD BDT Categories

        std::vector<std::string> data_pol5_pars(5,"p");   
        data_pol5_pars[0] = "modpol0";
        data_pol5_pars[1] = "modpol1";
        data_pol5_pars[2] = "modpol2";
        data_pol5_pars[3] = "modpol3";
        data_pol5_pars[4] = "modpol4";
        l.rooContainer->AddGenericPdf("data_pol_model","0","CMS_hgg_mass",data_pol5_pars,75); // >= 71 means RooBernstein of order >= 1

    } else if (bdtTrainingPhilosophy=="MIT"){
        // MIT BDT Categories
        // Background modeling - Separate Polynomial models for the different categories in MVA
        /*
          int poly3cats[5] = {1,1,0,0,0};
          int poly5cats[5] = {0,0,1,1,1};

          std::vector<std::string> data_pol5_pars(5,"p");   
          data_pol5_pars[0] = "modpol0";
          data_pol5_pars[1] = "modpol1";
          data_pol5_pars[2] = "modpol2";
          data_pol5_pars[3] = "modpol3";
          data_pol5_pars[4] = "modpol4";
          l.rooContainer->AddSpecificCategoryPdf(poly5cats,"data_pol_model",
          "0","CMS_hgg_mass",data_pol5_pars,75);  // >= 71 means RooBernstein of order >= 1
          std::vector<std::string> data_pol3_pars(3,"p");   
          data_pol3_pars[0] = "modpol0";
          data_pol3_pars[1] = "modpol1";
          data_pol3_pars[2] = "modpol2";
          l.rooContainer->AddSpecificCategoryPdf(poly3cats,"data_pol_model",
          "0","CMS_hgg_mass",data_pol3_pars,73);    // >= 71 means RooBernstein of order >= 1
        */
        if (includeVBF){
            int poly3cats[6] = {0,0,0,0,1,1};
            int poly4cats[6] = {1,0,0,0,0,0};
            int poly5cats[6] = {0,1,1,1,0,0};

            std::vector<std::string> data_pol5_pars(5,"p");   
            data_pol5_pars[0] = "modpol0";
            data_pol5_pars[1] = "modpol1";
            data_pol5_pars[2] = "modpol2";
            data_pol5_pars[3] = "modpol3";
            data_pol5_pars[4] = "modpol4";
            l.rooContainer->AddSpecificCategoryPdf(poly5cats,"data_pol_model",
                                                   "0","CMS_hgg_mass",data_pol5_pars,75);  // >= 71 means RooBernstein of order >= 1

            std::vector<std::string> data_pol4_pars(4,"p");   
            data_pol4_pars[0] = "modpol0";
            data_pol4_pars[1] = "modpol1";
            data_pol4_pars[2] = "modpol2";
            data_pol4_pars[3] = "modpol3";
            l.rooContainer->AddSpecificCategoryPdf(poly4cats,"data_pol_model",
                                                   "0","CMS_hgg_mass",data_pol4_pars,74);    // >= 71 means RooBernstein of order >= 1

            std::vector<std::string> data_pol3_pars(3,"p");   
            data_pol3_pars[0] = "modpol0";
            data_pol3_pars[1] = "modpol1";
            data_pol3_pars[2] = "modpol2";
            l.rooContainer->AddSpecificCategoryPdf(poly3cats,"data_pol_model",
                                                   "0","CMS_hgg_mass",data_pol3_pars,73);    // >= 71 means RooBernstein of order >= 1
        }
        else{
            int poly2cats[4] = {1,0,0,0};
            int poly4cats[4] = {0,1,0,0};
            int poly5cats[4] = {0,0,1,1};

            std::vector<std::string> data_pol5_pars(5,"p");   
            data_pol5_pars[0] = "modpol0";
            data_pol5_pars[1] = "modpol1";
            data_pol5_pars[2] = "modpol2";
            data_pol5_pars[3] = "modpol3";
            data_pol5_pars[4] = "modpol4";
            l.rooContainer->AddSpecificCategoryPdf(poly5cats,"data_pol_model",
                                                   "0","CMS_hgg_mass",data_pol5_pars,75);  // >= 71 means RooBernstein of order >= 1

            std::vector<std::string> data_pol4_pars(4,"p");   
            data_pol4_pars[0] = "modpol0";
            data_pol4_pars[1] = "modpol1";
            data_pol4_pars[2] = "modpol2";
            data_pol4_pars[3] = "modpol3";
            l.rooContainer->AddSpecificCategoryPdf(poly4cats,"data_pol_model",
                                                   "0","CMS_hgg_mass",data_pol4_pars,74);    // >= 71 means RooBernstein of order >= 1

            std::vector<std::string> data_pol2_pars(2,"p");   
            data_pol2_pars[0] = "modpol0";
            data_pol2_pars[1] = "modpol1";
            l.rooContainer->AddSpecificCategoryPdf(poly2cats,"data_pol_model",
                                                   "0","CMS_hgg_mass",data_pol2_pars,72);    // >= 71 means RooBernstein of order >= 1
            // -----------------------------------------------------
        }
    }

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
    /*
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

    */

    // Make more datasets representing Systematic Shifts of various quantities

    for (int sig=105;sig<=150;sig+=5){
        if (sig==145) continue;
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_ggh_mass_m%d",sig),-1);  
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_vbf_mass_m%d",sig),-1);  
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_wzh_mass_m%d",sig),-1);  
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_tth_mass_m%d",sig),-1);  
    }

    /*
      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_ggh_mass_m121",-1);  
      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_vbf_mass_m121",-1);  
      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_wzh_mass_m121",-1);  
      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_tth_mass_m121",-1);

      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_ggh_mass_m123",-1);  
      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_vbf_mass_m123",-1);  
      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_wzh_mass_m123",-1);  
      l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_tth_mass_m123",-1);  
    */

    // Make sure the Map is filled
    FillSignalLabelMap();

    // Initialize all MVA ---------------------------------------------------//
    l.SetAllMVA();
    // UCSD
    l.tmvaReaderID_UCSD->BookMVA("Gradient"      ,photonLevelMvaUCSD.c_str()  );
    l.tmvaReader_dipho_UCSD->BookMVA("Gradient"  ,eventLevelMvaUCSD.c_str()   );
    // New ID MVA
    l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevelNewIDMVA_EB.c_str());
    l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevelNewIDMVA_EE.c_str());
    // MIT 
    l.tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost",photonLevelMvaMIT_EB.c_str());
    l.tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost",photonLevelMvaMIT_EE.c_str());
    l.tmvaReader_dipho_MIT->BookMVA("Gradient"   ,eventLevelMvaMIT.c_str()    );
    // ----------------------------------------------------------------------//

    if(PADEBUG) 
        cout << "InitRealMassFactorizedMvaAnalysis END"<<endl;

    // FIXME book of additional variables
}

////////// // ----------------------------------------------------------------------------------------------------
////////// void MassFactorizedMvaAnalysis::Analysis(LoopAll& l, Int_t jentry) 
////////// {
//////////     if(PADEBUG) 
//////////         cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
////////// 
////////// 
//////////     int cur_type = l.itype[l.current];
//////////     float weight = l.sampleContainer[l.current_sample_index].weight;
//////////     float newweight = l.sampleContainer[l.current_sample_index].weight;
//////////     double pileupWeight=1.; 
//////////   
//////////     l.FillCounter( "Processed", 1. );
//////////     assert( weight > 0. );  
//////////     l.FillCounter( "XSWeighted", weight );
//////////     nevents+=1.;
////////// 
//////////     // Set reRunCiC Only if this is an MC event since scaling of R9 and Energy isn't done at reduction
//////////     if (cur_type==0) {
//////////         l.runCiC=reRunCiCForData;
//////////     } else {
//////////         l.runCiC = true;
//////////     }
//////////     if (l.runZeeValidation) l.runCiC=true;
////////// 
//////////     // -----------------------------------------------------------------------------------------------
////////// 
//////////     //PU reweighting
//////////     unsigned int n_pu = l.pu_n;
//////////     if ( cur_type !=0 && puHist != "" && cur_type < 100 ) {
//////////         bool hasSpecificWeight = weights.find( cur_type ) != weights.end() ; 
//////////         if( cur_type < 0 && !hasSpecificWeight && jentry == 1 ) {
//////////             std::cerr  << "WARNING no pu weights specific for sample " << cur_type << std::endl;
//////////         }
//////////         std::vector<double> & puweights = hasSpecificWeight ? weights[ cur_type ] : weights[0]; 
//////////         if(n_pu<puweights.size()){
//////////             weight *= puweights[n_pu];
//////////             sumwei+=puweights[n_pu]; 
//////////             pileupWeight *= puweights[n_pu];
//////////         }    
//////////         else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
//////////             cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< l.itype[l.current]<<"], event will not be reweighted for pileup"<<endl;
//////////         }
//////////     }
////////// 
//////////     assert( weight >= 0. );  
//////////     l.FillCounter( "PUWeighted", weight );
////////// 
//////////     if( jentry % 10000 ==  0 ) {
//////////         std::cout << " nevents " <<  nevents << " sumpuweights " << sumwei << " ratio " << sumwei / nevents 
//////////                   << " equiv events " << sumev << " accepted " << sumaccept << " smeared " << sumsmear << " "  
//////////                   <<  sumaccept / sumev << " " << sumsmear / sumaccept
//////////                   << std::endl;
//////////     }
//////////     // ------------------------------------------------------------
//////////     //PT-H K-factors
//////////     double gPT = 0;
//////////     TLorentzVector gP4(0,0,0,0);
//////////     if (cur_type<0){
////////// 	    gP4 = l.GetHiggs();
////////// 	    gPT = gP4.Pt();
//////////     }
//////////     // ------------------------------------------------------------
////////// 
//////////     // TEMPORARY FIX -------------------------------------------------------------------------------------------------------//
//////////     // Scale all the r9 of the photons (and also some other variables) in the MC for better agreement
//////////     // For now we just let it use the index but we specifically Change the r9 in the branch AFTER Energy regression smearing
//////////     // Ideally we want to pass a smeared r9 too and apply after energy corrections, currently the smeared_pho_r9 isnt used!
//////////     // ---------------------------------------------------------------------------------------------------------------------//
//////////     // ---------------------------------------------------------------------------------------------------------------------//
//////////     // ---------------------------------------------------------------------------------------------------------------------//
////////// 
//////////     if (cur_type !=0){
//////////         for (int ipho=0;ipho<l.pho_n;ipho++){
//////////             l.pho_isEB[ipho]=fabs((*((TVector3*)l.sc_xyz->At(l.pho_scind[ipho]))).Eta())<1.5;
//////////             //double R9_rescale = (l.pho_isEB[ipho]) ? 1.0048 : 1.00492 ;
//////////             //l.pho_r9[ipho]*=R9_rescale;
//////////             l.pho_r9[ipho]*=1.0035;
//////////             if (l.pho_isEB[ipho]){ l.pho_sieie[ipho] = (0.87*l.pho_sieie[ipho]) + 0.0011 ;}
//////////             else {l.pho_sieie[ipho]*=0.99;}
//////////             l.sc_seta[l.pho_scind[ipho]]*=0.99;  
//////////             l.sc_sphi[l.pho_scind[ipho]]*=0.99;  
//////////             energyCorrectedError[ipho] *=(l.pho_isEB[ipho]) ? 1.07 : 1.045 ;
//////////         }
//////////     }
//////////     // ---------------------------------------------------------------------------------------------------------------------//
//////////     // ---------------------------------------------------------------------------------------------------------------------//
//////////     // ---------------------------------------------------------------------------------------------------------------------//
//////////     // ---------------------------------------------------------------------------------------------------------------------//
////////// 
//////////     // Analyse the event assuming nominal values of corrections and smearings
//////////     float mass, evweight;
//////////     int diphoton_id, category;
//////////     bool isCorrectVertex;
//////////     if( AnalyseEvent(l,jentry, weight, gP4, mass,  evweight, category, diphoton_id, isCorrectVertex,zero_) ) {
////////// 	// feed the event to the RooContainer 
//////////         if (cur_type == 0 ){
//////////             l.rooContainer->InputDataPoint("data_mass",category,mass);
//////////         }
//////////         if (cur_type > 0 && cur_type != 3 && cur_type != 4) {
//////////             l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
////////// 	}
//////////         else if (cur_type < 0){
//////////             l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type),category,mass,evweight);
//////////         }
//////////     }
////////// 
//////////     // Systematics uncertaities for the binned model
//////////     // We re-analyse the event several times for different values of corrections and smearings
//////////     if( cur_type < 0 && doMCSmearing && doSystematics ) { 
//////////         
//////////         // fill steps for syst uncertainty study
//////////         float systStep = systRange / (float)nSystSteps;
////////// 	
////////// 	    float syst_mass, syst_weight;
////////// 	    int syst_category;
//////////  	    std::vector<double> mass_errors;
////////// 	    std::vector<double> weights;
////////// 	    std::vector<int>    categories;
////////// 	
//////////         if (diphoton_id > -1 ) {
//////////      
////////// 	    // gen-level systematics, i.e. ggH k-factor for the moment
//////////             for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
////////// 		mass_errors.clear(), weights.clear(), categories.clear();
////////// 		
//////////                 for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
//////////                     if( syst_shift == 0. ) { continue; } // skip the central value
////////// 		    syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
////////// 		    
////////// 		    // re-analyse the event without redoing the event selection as we use nominal values for the single photon
////////// 		    // corrections and smearings
////////// 		    AnalyseEvent(l, jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,zero_,
////////// 				 true, syst_shift, true, *si, 0, 0 );
////////// 		    		    
////////// 		    categories.push_back(syst_category);
////////// 		    mass_errors.push_back(syst_mass);
//////////                     weights.push_back(syst_weight);
////////// 		}
////////// 		if (cur_type < 0){
////////// 		    // feed the modified signal model to the RooContainer
////////// 		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
////////// 		}
////////// 	    }
////////// 	    
////////// 	    // di-photon systematics: vertex efficiency and trigger 
////////// 	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
////////// 		mass_errors.clear(), weights.clear(), categories.clear();
////////// 		
//////////                 for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
//////////                     if( syst_shift == 0. ) { continue; } // skip the central value
////////// 		    syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
////////// 		    
////////// 		    // re-analyse the event without redoing the event selection as we use nominal values for the single photon
////////// 		    // corrections and smearings
////////// 		    AnalyseEvent(l,jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,zero_,
////////// 				 true, syst_shift, true,  0, 0, *si );
////////// 		    
////////// 		    categories.push_back(syst_category);
////////// 		    mass_errors.push_back(syst_mass);
//////////                     weights.push_back(syst_weight);
////////// 		}
////////// 		if (cur_type < 0){
////////// 		    // feed the modified signal model to the RooContainer
////////// 		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
////////// 		}
////////// 	    }
////////// 	}
////////// 	
////////// 	// single photon level systematics: several
////////// 	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
////////// 	    mass_errors.clear(), weights.clear(), categories.clear();
////////// 	    
////////// 	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
////////// 		if( syst_shift == 0. ) { continue; } // skip the central value
////////// 		syst_mass     =  0., syst_category = -1, syst_weight   =  0.;
////////// 		
////////// 		// re-analyse the event redoing the event selection this time
////////// 		AnalyseEvent(l,jentry, weight, gP4, syst_mass,  syst_weight, syst_category, diphoton_id, isCorrectVertex,zero_,
////////// 			     true, syst_shift, false,  0, *si, 0 );
////////// 		
////////// 		categories.push_back(syst_category);
////////// 		mass_errors.push_back(syst_mass);
////////// 		weights.push_back(syst_weight);
////////// 	    }
////////// 	    if (cur_type < 0){
////////// 		// feed the modified signal model to the RooContainer
////////// 		l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
////////// 	    }
////////// 	}
//////////     }
////////// 
//////////     if(PADEBUG) 
//////////         cout<<"myFillHistRed END"<<endl;
////////// }

bool MassFactorizedMvaAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4,
				float & mass, float & evweight, int & category, int & diphoton_id, bool & isCorrectVertex, float &kinematic_bdtout,
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
	
    if (!skipSelection){
        // first apply corrections and smearing on the single photons 
        smeared_pho_energy.clear(); smeared_pho_energy.resize(l.pho_n,0.); 
        smeared_pho_r9.clear();     smeared_pho_r9.resize(l.pho_n,0.); 
        smeared_pho_weight.clear(); smeared_pho_weight.resize(l.pho_n,1.);
        applySinglePhotonSmearings(smeared_pho_energy, smeared_pho_r9, smeared_pho_weight, cur_type, l, energyCorrected, energyCorrectedError,
                       phoSys, syst_shift);

        // inclusive category di-photon selection
        // FIXME pass smeared R9
        diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] ); 
            
        // Exclusive Modes
        int diphotonVBF_id = -1;
        VBFevent = false;
        
        // VBF
        if((includeVBF || includeVHhad)&&l.jet_algoPF1_n>1 && !isSyst /*avoid rescale > once*/) {
            l.RescaleJetEnergy();
        }

        if(includeVBF) {
            diphotonVBF_id = l.DiphotonMITPreSelection(leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] );
            float eventweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
            float myweight=1.;
            if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
            
	    VBFevent= ( dataIs2011 ? 
			VBFTag2011(l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight) :
			VBFTag2012(l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight) )
		    ;
        }
        
        if(includeVBF&&VBFevent) {
            diphoton_id = diphotonVBF_id;
        }
    }
    // if we selected any di-photon, compute the Higgs candidate kinematics
    // and compute the event category
    if (PADEBUG) std::cout << "Found a Diphoton , diphoton ID " <<diphoton_id << std::endl; 
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
	    mass     = Higgs.M();
        float ptHiggs = Higgs.Pt();

	// For Zee validation, reweight MC pT distribution to match data 
	if( l.runZeeValidation && cur_type != 0) {
	  if (zeePtBinLowEdge.size() != zeePtWeight.size()) {
	    std::cout << "Array size mismatch: zeePtBinLowEdge[" << zeePtBinLowEdge.size()
		      << "], zeePtWeight[" << zeePtWeight.size() << "]" <<diphoton_id << std::endl;
	  }
	  for (int i=0; i<zeePtBinLowEdge.size(); i++) {
	    float zeePtBinHighEdge = 999.;
	    if (i<zeePtBinLowEdge.size()-1) zeePtBinHighEdge = zeePtBinLowEdge[i+1];
	    if (ptHiggs>zeePtBinLowEdge[i] && ptHiggs<zeePtBinHighEdge) {
	      evweight *= zeePtWeight[i];
	      //cout << ptHiggs << " " << zeePtWeight[i] << endl;
	      break;
	    }
	  }
	}

        // Mass Resolution of the Event
        massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
        float vtx_mva  = l.vtx_std_evt_mva->at(diphoton_id);
        float sigmaMrv = massResolutionCalculator->massResolutionEonly();
        float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
        float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
        // easy to calculate vertex probability from vtx mva output
        float vtxProb   = 1.-0.49*(vtx_mva+1.0); /// should better use this: vtxAna_.vertexProbability(vtx_mva); PM

        float phoid_mvaout_lead = l.photonIDMVANew(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4,bdtTrainingPhilosophy.c_str());
        float phoid_mvaout_sublead = l.photonIDMVANew(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4,bdtTrainingPhilosophy.c_str());

	    // apply di-photon level smearings and corrections
        int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
        if( cur_type != 0 && doMCSmearing ) {
	    applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, phoid_mvaout_lead,phoid_mvaout_sublead,
				   diPhoSys, syst_shift);
            isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
        }
                           
        // Must be calculated after photon id has potentially been smeared
        float diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id] ,vtxProb,lead_p4,sublead_p4 ,sigmaMrv,sigmaMwv,sigmaMeonly ,bdtTrainingPhilosophy.c_str() ,phoid_mvaout_lead,phoid_mvaout_sublead);
        kinematic_bdtout = diphobdt_output;
        bool isEBEB  = (lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);
        category = GetBDTBoundaryCategory(diphobdt_output,isEBEB,VBFevent);
        if (diphobdt_output>=0.05) computeExclusiveCategory(l,category,diphoton_index,Higgs.Pt()); 
 
        if (PADEBUG) std::cout << " Diphoton Category " <<category <<std::endl;
    	// sanity check
        assert( evweight >= 0. ); 
	
	    // fill control plots and counters
	    if( ! isSyst ) {
	        l.FillCounter( "Accepted", weight );
	        l.FillCounter( "Smeared", evweight );
	        sumaccept += weight;
	        sumsmear += evweight;
	        if (l.runZeeValidation) {
		  fillZeeControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, phoid_mvaout_lead, phoid_mvaout_sublead, diphobdt_output, sigmaMrv, sigmaMwv, vtxProb, diphoton_id, category, evweight, l );
		} else {
		  fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9,  diphoton_id, 
				   category, isCorrectVertex, evweight, l );
		}
	    }
	
        //if (cur_type==0 && mass >= 100. && mass < 180. && !isSyst /*should never be if running data anyway*/){
        if (mass >= 100. && mass < 180. && !isSyst){
        //if (1>0){
            eventListText <<"Type="<< cur_type 
            <<" Run=" << l.run 
            <<" LS=" << l.lumis 
            <<" Event=" << l.event 
            <<" BDTCAT=" << category 
            <<" ggM=" << mass 
            <<" gg_Pt=" << ptHiggs 
            <<" LeadPhotonPhoid=" <<phoid_mvaout_lead 
            <<" SubleadPhotonPhoid=" <<phoid_mvaout_sublead 
            <<" diphotonBDT=" << diphobdt_output 
            <<" photon1Eta=" << lead_p4.Eta() 
            <<" photon2Eta="<<sublead_p4.Eta() 
            <<" sigmaMrv="<<sigmaMrv 
            <<" sigmaMwv=" << sigmaMwv 
            <<" photon1Pt="<<lead_p4.Pt()
            <<" photon2Pt="<<sublead_p4.Pt() 
            <<" vtxProb="<<vtxProb 
            <<" cosDphi="<<TMath::Cos(lead_p4.Phi() - sublead_p4.Phi()) 
            <<" r9_1=" <<lead_r9 
            <<" r9_2=" <<sublead_r9 
            <<" E1="<<lead_p4.E()
            <<" E2="<<sublead_p4.E()
            <<" reresraw1="<<energyCorrectedError[diphoton_index.first]/energyCorrected[diphoton_index.first] 
            <<" reresraw2="<<energyCorrectedError[diphoton_index.second]/energyCorrected[diphoton_index.second]
            <<" reresraw1JBD="<<l.pho_regr_energyerr[diphoton_index.first]/l.pho_regr_energy[diphoton_index.first] 
            <<" reresraw2JBD="<<l.pho_regr_energyerr[diphoton_index.second]/l.pho_regr_energy[diphoton_index.second]
            <<" etcorecal1="<<l.pho_ecalsumetconedr03[diphoton_index.first] - 0.012*lead_p4.Et()
            <<" etcorecal2="<<l.pho_ecalsumetconedr03[diphoton_index.second] - 0.012*sublead_p4.Et()
            <<" etcorhcal1="<<l.pho_hcalsumetconedr03[diphoton_index.first] - 0.005*lead_p4.Et()
            <<" etcorhcal2="<<l.pho_hcalsumetconedr03[diphoton_index.second] - 0.005*sublead_p4.Et()
            <<" etcortrkiso1="<<l.pho_trksumpthollowconedr03[diphoton_index.first] - 0.002*lead_p4.Et()
            <<" etcortrkiso2="<<l.pho_trksumpthollowconedr03[diphoton_index.second] - 0.002*sublead_p4.Et()
            <<" trkhollowcone1="<<l.pho_trksumpthollowconedr03[diphoton_index.first] 
            <<" trkhollowcone2="<<l.pho_trksumpthollowconedr03[diphoton_index.second]
            <<" pucorrhcalecal1="<<(l.pho_ecalsumetconedr03[diphoton_index.first]+l.pho_hcalsumetconedr03[diphoton_index.first]-l.rho*0.17)
            <<" pucorrhcalecal2="<<(l.pho_ecalsumetconedr03[diphoton_index.second]+l.pho_hcalsumetconedr03[diphoton_index.second]-l.rho*0.17)
            <<" cictrkiso1="<<(*l.pho_tkiso_recvtx_030_002_0000_10_01)[diphoton_index.first][l.dipho_vtxind[diphoton_id]] 
            <<" cictrkiso2="<<(*l.pho_tkiso_recvtx_030_002_0000_10_01)[diphoton_index.second][l.dipho_vtxind[diphoton_id]] 
            <<" specialphoton1="<<photonInfoCollection[diphoton_index.first].isSphericalPhoton()  
            <<" specialphoton2="<<photonInfoCollection[diphoton_index.second].isSphericalPhoton()
            <<" specialphotongen1="<<l.CheckSphericalPhoton(diphoton_index.first) 
            <<" specialphotongen2="<<l.CheckSphericalPhoton(diphoton_index.second)
            <<" hovere1="<<l.pho_hoe[diphoton_index.first]
            <<" hovere2="<<l.pho_hoe[diphoton_index.second]
            <<" sigietaieta1="<<l.pho_sieie[diphoton_index.first]
            <<" sigietaieta2="<<l.pho_sieie[diphoton_index.second]
            <<" VBFevent="<<VBFevent
            <<" FileName="<<l.files[l.current];
            eventListText << endl;
        }
	    return true;
    }
    return false;
}

// ----------------------------------------------------------------------------------------------------
void MassFactorizedMvaAnalysis::fillZeeControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
	                        const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
				float phoid_mvaout_lead, float phoid_mvaout_sublead, 
                                float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
				int diphoton_id, int category, float evweight, LoopAll & l )
{

  float mass = Higgs.M();
  float ptHiggs = Higgs.Pt();

  float cos_dphi = TMath::Cos(lead_p4.Phi()-sublead_p4.Phi());
  float pho1_sigmaE = energyCorrectedError[l.dipho_leadind[diphoton_id]];
  float pho2_sigmaE = energyCorrectedError[l.dipho_subleadind[diphoton_id]];

  l.FillHist("mass",0, mass, evweight);
  l.FillHist("mass",category+1, mass, evweight);
  if (ptHiggs<20.) {
    l.FillHist("mass_pt0to20",0, mass, evweight);
  } else if (ptHiggs<40.) {
    l.FillHist("mass_pt20to40",0, mass, evweight);
  } else if (ptHiggs<60.) {
    l.FillHist("mass_pt40to60",0, mass, evweight);
  } else if (ptHiggs<100.) {
    l.FillHist("mass_pt60to100",0, mass, evweight);
  } else {
    l.FillHist("mass_pt100up",0, mass, evweight);
  }

  if (mass>=100. && mass<=180.) l.FillHist("bdtout_m100to180",0,diphobdt_output,evweight);

  if( mass>=60. && mass<=120.  ) {

    l.FillHist("bdtout",0,diphobdt_output,evweight);
    if (fabs(lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442)){
      l.FillHist("bdtoutEB",0,diphobdt_output,evweight);
    } else if (fabs(lead_p4.Eta() > 1.566 ) && fabs(sublead_p4.Eta()>1.566)){
      l.FillHist("bdtoutEE",0,diphobdt_output,evweight);
    } else {
      l.FillHist("bdtoutEBEE",0,diphobdt_output,evweight);
    }

    l.FillHist("nvtx",0, l.vtx_std_n, evweight);
    l.FillHist("nvtx",category+1, l.vtx_std_n, evweight);
    l.FillHist("pt",0, ptHiggs, evweight);
    l.FillHist("pt",category+1, ptHiggs, evweight);

    if (fabs(lead_p4.Eta()) < 1.4442) {
      l.FillHist("pho1_sigmaEOverE_EB",0,pho1_sigmaE/lead_p4.E(), evweight);
      l.FillHist("pho1_sigmaEOverE_EB_up",0,(pho1_sigmaE/1.07*1.17)/lead_p4.E(), evweight);
      l.FillHist("pho1_sigmaEOverE_EB_down",0,(pho1_sigmaE/1.07*0.97)/lead_p4.E(), evweight);
    } else {
      l.FillHist("pho1_sigmaEOverE_EE",0,pho1_sigmaE/lead_p4.E(), evweight);
      l.FillHist("pho1_sigmaEOverE_EE_up",0,(pho1_sigmaE/1.045*1.145)/lead_p4.E(), evweight);
      l.FillHist("pho1_sigmaEOverE_EE_down",0,(pho1_sigmaE/1.045*0.945)/lead_p4.E(), evweight);
    }

    if (fabs(sublead_p4.Eta()) < 1.4442) {
      l.FillHist("pho2_sigmaEOverE_EB",0,pho2_sigmaE/sublead_p4.E(), evweight);
      l.FillHist("pho2_sigmaEOverE_EB_up",0,(pho2_sigmaE/1.07*1.17)/sublead_p4.E(), evweight);
      l.FillHist("pho2_sigmaEOverE_EB_down",0,(pho2_sigmaE/1.07*0.97)/sublead_p4.E(), evweight);
    } else {
      l.FillHist("pho2_sigmaEOverE_EE",0,pho2_sigmaE/sublead_p4.E(), evweight);
      l.FillHist("pho2_sigmaEOverE_EE_up",0,(pho2_sigmaE/1.045*1.145)/sublead_p4.E(), evweight);
      l.FillHist("pho2_sigmaEOverE_EE_down",0,(pho2_sigmaE/1.045*0.945)/sublead_p4.E(), evweight);
    }

    if (fabs(lead_p4.Eta()) < 1.4442) {
      l.FillHist("pho1_phoidMva_EB",0,phoid_mvaout_lead, evweight);
      l.FillHist("pho1_phoidMva_EB_up",0,(phoid_mvaout_lead+0.025), evweight);
      l.FillHist("pho1_phoidMva_EB_down",0,(phoid_mvaout_lead-0.025), evweight);
    } else {
      l.FillHist("pho1_phoidMva_EE",0,phoid_mvaout_lead, evweight);
      l.FillHist("pho1_phoidMva_EE_up",0,(phoid_mvaout_lead+0.025), evweight);
      l.FillHist("pho1_phoidMva_EE_down",0,(phoid_mvaout_lead-0.025), evweight);
    }

    if (fabs(sublead_p4.Eta()) < 1.4442) {
      l.FillHist("pho2_phoidMva_EB",0,phoid_mvaout_sublead, evweight);
      l.FillHist("pho2_phoidMva_EB_up",0,(phoid_mvaout_sublead+0.025), evweight);
      l.FillHist("pho2_phoidMva_EB_down",0,(phoid_mvaout_sublead-0.025), evweight);
    } else {
      l.FillHist("pho2_phoidMva_EE",0,phoid_mvaout_sublead, evweight);
      l.FillHist("pho2_phoidMva_EE_up",0,(phoid_mvaout_sublead+0.025), evweight);
      l.FillHist("pho2_phoidMva_EE_down",0,(phoid_mvaout_sublead-0.025), evweight);
    }

    if (diphobdt_output > 0.05) {

      l.FillHist("pho1_phoidMva",0,phoid_mvaout_lead, evweight);
      l.FillHist("pho1_phoidMva",category+1,phoid_mvaout_lead, evweight);
      l.FillHist("pho2_phoidMva",0,phoid_mvaout_sublead, evweight);
      l.FillHist("pho2_phoidMva",category+1,phoid_mvaout_sublead, evweight);
      l.FillHist("sigmaMOverM",0,sigmaMrv/mass, evweight);
      l.FillHist("sigmaMOverM",category+1,sigmaMrv/mass, evweight);
      l.FillHist("sigmaMOverM_wrongVtx",0,sigmaMwv/mass, evweight);
      l.FillHist("sigmaMOverM_wrongVtx",category+1,sigmaMwv/mass, evweight);
      l.FillHist("vtxProb",0,vtxProb, evweight);
      l.FillHist("vtxProb",category+1,vtxProb, evweight);
      l.FillHist("pho1_ptOverM",0,lead_p4.Pt()/mass, evweight);
      l.FillHist("pho1_ptOverM",category+1,lead_p4.Pt()/mass, evweight);
      l.FillHist("pho2_ptOverM",0,sublead_p4.Pt()/mass, evweight);
      l.FillHist("pho2_ptOverM",category+1,sublead_p4.Pt()/mass, evweight);
      l.FillHist("pho1_eta",0,lead_p4.Eta(), evweight);
      l.FillHist("pho1_eta",category+1,lead_p4.Eta(), evweight);
      l.FillHist("pho2_eta",0,sublead_p4.Eta(), evweight);
      l.FillHist("pho2_eta",category+1,sublead_p4.Eta(), evweight);
      l.FillHist("cosDeltaPhi",0,cos_dphi, evweight);
      l.FillHist("cosDeltaPhi",category+1,cos_dphi, evweight);

      l.FillHist("r9",0,lead_r9, evweight);
      l.FillHist("r9",0,sublead_r9, evweight);
      l.FillHist("r9",category+1,lead_r9, evweight);
      l.FillHist("r9",category+1,sublead_r9, evweight);

      if (fabs(lead_p4.Eta()) < 1.4442 && fabs(sublead_p4.Eta())<1.4442) {
	l.FillHist("sigmaMOverM_EB",0,sigmaMrv/mass, evweight);
	l.FillHist("sigmaMOverM_wrongVtx_EB",0,sigmaMwv/mass, evweight);
	if (lead_r9>0.93 && sublead_r9>0.93) {
	  l.FillHist("sigmaMOverM_EB",1,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EB",1,sigmaMwv/mass, evweight);
	} else if (lead_r9<0.93 && sublead_r9<0.93) {
	  l.FillHist("sigmaMOverM_EB",3,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EB",3,sigmaMwv/mass, evweight);
	} else {
	  l.FillHist("sigmaMOverM_EB",2,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EB",2,sigmaMwv/mass, evweight);
	}
      } else if (fabs(lead_p4.Eta()) > 1.566 && fabs(sublead_p4.Eta())>1.566) {
	l.FillHist("sigmaMOverM_EE",0,sigmaMrv/mass, evweight);
	l.FillHist("sigmaMOverM_wrongVtx_EE",0,sigmaMwv/mass, evweight);
	if (lead_r9>0.93 && sublead_r9>0.93) {
	  l.FillHist("sigmaMOverM_EE",1,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EE",1,sigmaMwv/mass, evweight);
	} else if (lead_r9<0.93 && sublead_r9<0.93) {
	  l.FillHist("sigmaMOverM_EE",3,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EE",3,sigmaMwv/mass, evweight);
	} else {
	  l.FillHist("sigmaMOverM_EE",2,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EE",2,sigmaMwv/mass, evweight);
	}
      } else {
	l.FillHist("sigmaMOverM_EBEE",0,sigmaMrv/mass, evweight);
	l.FillHist("sigmaMOverM_wrongVtx_EBEE",0,sigmaMwv/mass, evweight);
	if (lead_r9>0.93 && sublead_r9>0.93) {
	  l.FillHist("sigmaMOverM_EBEE",1,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EBEE",1,sigmaMwv/mass, evweight);
	} else if (lead_r9<0.93 && sublead_r9<0.93) {
	  l.FillHist("sigmaMOverM_EBEE",3,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EBEE",3,sigmaMwv/mass, evweight);
	} else {
	  l.FillHist("sigmaMOverM_EBEE",2,sigmaMrv/mass, evweight);
	  l.FillHist("sigmaMOverM_wrongVtx_EBEE",2,sigmaMwv/mass, evweight);
	}
      }

    }

  }

}


// ----------------------------------------------------------------------------------------------------
int MassFactorizedMvaAnalysis::GetBDTBoundaryCategory(float bdtout, bool isEB, bool VBFevent){

    if (bdtTrainingPhilosophy=="UCSD"){
        if (isEB) { // 6 Categories for the EB-EB 
            if (bdtout < -0.30) return 5;
            if (bdtout >=-0.30 && bdtout < 0.00) return 4;
            if (bdtout >= 0.00 && bdtout < 0.30) return 3;
            if (bdtout >= 0.30 && bdtout < 0.60) return 2;
            if (bdtout >= 0.60 && bdtout < 0.70) return 1;
            if (bdtout >= 0.70) return 0;
        }
        else {// 2 Categories for the EB/EE 
            if (bdtout <  0.1) return 7;
            if (bdtout >= 0.1) return 6;
        }

    } else if (bdtTrainingPhilosophy=="MIT"){
        /*
          if (bdtout >=-0.50 && bdtout < 0.3) return 4;
          //     if (bdtout < 0.3) return 4;
          if (bdtout >= 0.3 && bdtout < 0.65) return 3;
          if (bdtout >= 0.65 && bdtout < 0.84) return 2;
          if (bdtout >= 0.84 && bdtout < 0.90) return 1;
          if (bdtout >= 0.90) return 0;
          return -1;
        */
        //       if (bdtout >=-0.50 && bdtout < 0.05) return 4;
        if (VBFevent) {
            //if (bdtout >= 0.05 && VBFevent==1) return 4;
            //if (bdtout >= 0.05 && VBFevent==2) return 5;
            if (bdtout >= 0.05) return 4;
            else return -1;
        } else {
            if (bdtout >= 0.05 && bdtout < 0.545) return 3;
            if (bdtout >= 0.545 && bdtout < 0.74) return 2;
            if (bdtout >= 0.74 && bdtout < 0.89) return 1;
            if (bdtout >= 0.89) return 0;
            return -1;
        }

    } else std::cerr << "No BDT Philosophy known - " << bdtTrainingPhilosophy << std::endl;
}

void MassFactorizedMvaAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    eResolSmearer->resetRandom();
}



// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
