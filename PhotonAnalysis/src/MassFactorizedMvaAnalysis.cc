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

  std::string outputfilename = (std::string) l.histFileName;
  l.rooContainer->FitToData("data_pol_model","data_mass");  // Fit to full range of dataset

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

  // Define the number of categories for the statistical analysis and
  // the systematic sets to be formed

  // FIXME move these params to config file
  if (bdtTrainingPhilosophy == "UCSD") {
    l.rooContainer->SetNCategories(8);
  } else if (bdtTrainingPhilosophy == "MIT") {
    if (includeVBF)  l.rooContainer->SetNCategories(5);
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
      int poly3cats[5] = {0,0,0,0,1};
      int poly4cats[5] = {1,0,0,0,0};
      int poly5cats[5] = {0,1,1,1,0};

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
  // MIT 
  l.tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost",photonLevelMvaMIT_EB.c_str());
  l.tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost",photonLevelMvaMIT_EE.c_str());
  l.tmvaReader_dipho_MIT->BookMVA("Gradient"   ,eventLevelMvaMIT.c_str()    );
  // ----------------------------------------------------------------------//

  if(PADEBUG) 
    cout << "InitRealMassFactorizedMvaAnalysis END"<<endl;

  // FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void MassFactorizedMvaAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
  if(PADEBUG) 
    cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;


  int cur_type = l.itype[l.current];
  float weight = l.sampleContainer[l.current_sample_index].weight;
  l.FillCounter( "Processed", 1. );
  assert( weight > 0. );  
  l.FillCounter( "XSWeighted", weight );
  nevents+=1.;

  // Set reRunCiC Only if this is an MC event since scaling of R9 and Energy isn't done at reduction
  if (cur_type==0) {
    l.runCiC=reRunCiCForData;
  } else {
    l.runCiC = true;
  }

  // -----------------------------------------------------------------------------------------------

  // Re-do Vertexing
  
   //PhotonAnalysis::FillReductionVariables(l, jentry);
  // PhotonAnalysis::SelectEventsReduction(l, jentry);

  //PU reweighting
  unsigned int n_pu = l.pu_n;
  if ( cur_type !=0 && puHist != "" && cur_type < 100 ) {
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

  // TEMPORARY FIX -------------------------------------------------------------------------------------------------------//
  // Scale all the r9 of the photons (and also some other variables) in the MC for better agreement
  // For now we just let it use the index but we specifically Change the r9 in the branch AFTER Energy regression smearing
  // Ideally we want to pass a smeared r9 too and apply after energy corrections, currently the smeared_pho_r9 isnt used!
  // ---------------------------------------------------------------------------------------------------------------------//
  // ---------------------------------------------------------------------------------------------------------------------//
  // ---------------------------------------------------------------------------------------------------------------------//

  if (cur_type !=0){
    for (int ipho=0;ipho<l.pho_n;ipho++){
      l.pho_isEB[ipho]=fabs((*((TVector3*)l.sc_xyz->At(l.pho_scind[ipho]))).Eta())<1.5;
      //double R9_rescale = (l.pho_isEB[ipho]) ? 1.0048 : 1.00492 ;
      //l.pho_r9[ipho]*=R9_rescale;
      l.pho_r9[ipho]*=1.0035;
      if (l.pho_isEB[ipho]){ l.pho_sieie[ipho] = (0.87*l.pho_sieie[ipho]) + 0.0011 ;}
      else {l.pho_sieie[ipho]*=0.99;}
      l.sc_seta[l.pho_scind[ipho]]*=0.99;  
      l.sc_sphi[l.pho_scind[ipho]]*=0.99;  
      energyCorrectedError[ipho] *=(l.pho_isEB[ipho]) ? 1.07 : 1.045 ;
    }
  }
  // ---------------------------------------------------------------------------------------------------------------------//
  // ---------------------------------------------------------------------------------------------------------------------//
  // ---------------------------------------------------------------------------------------------------------------------//
  // ---------------------------------------------------------------------------------------------------------------------//

  photonInfoCollection.clear(); // this is a member of PhotonAnalysis, just clear it
  for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
    std::vector<std::vector<bool> > p;
    PhotonReducedInfo phoInfo (// *((TVector3*)l.pho_calopos->At(ipho)), 
        *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
        ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
        energyCorrected[ipho],
        l.pho_isEB[ipho], l.pho_r9[ipho],
        l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_),
        (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
        );
    if (l.CheckSphericalPhoton(ipho)) phoInfo.setSphericalPhoton(true);
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
    photonInfoCollection.push_back(phoInfo);
  }

  sumev += weight;

  // Get Ready for VBF Tagging
  bool VBFevent = false;
  double leadEtCutVBF = 55.;
  int diphoton_id=-1;
    // VBF-TAGGING -------------------------------------------------------------------------- //
    // CP // NW Use Same pre-selected events but tag as VBF if pass Jet Requirements

    if(includeVBF) {
      PhotonAnalysis::RescaleJetEnergy(l);
      int diphoton_id_vbf = l.DiphotonMITPreSelection(leadEtCutVBF,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 

      if (diphoton_id_vbf > -1){
        diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id_vbf],  l.dipho_subleadind[diphoton_id_vbf] );

        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
        TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
        float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id_vbf]];
        float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id_vbf]];
        TLorentzVector Higgs = lead_p4 + sublead_p4;   
      // JET MET Corrections
      float jet1ptcut =0.0;
      float jet2ptcut =0.0;
      bool crosscheck = false;
      std::pair<int,int> highestPtJets(-1,-1);

      highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );
      bool VBFpresel = (highestPtJets.first!=-1)&&(highestPtJets.second!=-1);

      // Make sure VBF event is set to false before checking for the Jets
      VBFevent = false; 
      if(VBFpresel){

        TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
        TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
        TLorentzVector dijet = (*jet1) + (*jet2);

        myAllLeadJPt = jet1->Pt();
        myAllSubJPt = jet2->Pt();
        myAllLeadJEta = jet1->Eta();
        myAllSubJEta = jet2->Eta();
        myAll_Mjj = dijet.M();
        myAlldEta = fabs(jet1->Eta() - jet2->Eta());
        myAllZep  = fabs(Higgs.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
        myAlldPhi = fabs(Higgs.DeltaPhi(dijet));
        myAll_Mgg =Higgs.M();
        myAllPtHiggs =Higgs.Pt();

        myVBFLeadJPt = jet1->Pt();
        myVBFSubJPt = jet2->Pt();
        myVBF_Mjj = dijet.M();
        myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
        myVBFZep  = fabs(Higgs.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
        myVBFdPhi = fabs(Higgs.DeltaPhi(dijet));
        myVBF_Mgg =Higgs.M();


        // Cannot Get Apply cuts to work -> Need to discuss with C.Palmer, for now, jyst apply the cuts
        //l.ApplyCutsFill(0,3,evweight, myweight);
        //VBFevent = l.ApplyCutsFill(0,5,evweight, myweight);
        //VBFevent = l.ApplyCutsFill(0,1,evweight, myweight);
        if (myVBFLeadJPt  >30.){
          if (myVBFSubJPt >20.){
            if (myVBF_Mjj   >350.){
              if (myVBFdEta   >3.5){
                if (myVBFdPhi   >2.6){
                  if (myVBFZep    <2.5){
                    VBFevent=true;
			  diphoton_id = diphoton_id_vbf;
                  }
                }
              }   
            }
          }
        }
        /*  
            if(VBFevent && cur_type==-18) 
            std::cout << setprecision(4) <<  "Run = " << l.run << "  LS = " << l.lumis <<
            "  Event = " << l.event << "  SelVtx = " << l.dipho_vtxind[diphoton_id] 
            << "  CAT4 = " << CAT4 << "  ggM = " << myVBF_Mgg << " ggPt =  " << myAllPtHiggs 
            << "  jetEta1 = " << jet1->Eta() << "  jetEta2 = " << jet2->Eta()
            << "  jetPhi1 = " << jet1->Phi() << "  jetPhi2 = " << jet2->Phi()
            <<  "  jetEt1 = " << jet1->Et() << "  jetEt2 = "  << jet2->Et()
            << " Mjj " << myVBF_Mjj
            << " dEtajj " << myVBFdEta 
            << " Zeppenfeld " << myVBFZep
            << " dPhijjgg " << myVBFdPhi << " VBF itype " <<cur_type << std::endl;
         */
      }
    }
  }
    // CP // NW VBF Tagging
    // --------------------- END VBF-TAGGING --------------------------------------------------------//

  // Divergence from StatAnalysis Here! Apply loose pre-selection to select photons
  // FIXME pass smeared R9
  if (diphoton_id<0 ){ // then the VBF selection failed at some point and so we try for inclusive
  if (bdtTrainingPhilosophy=="MIT"){
    diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 
  } else if (bdtTrainingPhilosophy=="UCSD"){
    diphoton_id = l.DiphotonCiCSelection(l.phoLOOSE, l.phoLOOSE, leadEtCut, subleadEtCut, nPhotonCategories_,applyPtoverM, &smeared_pho_energy[0] ); 
  }
  }

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

    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);

    bool CorrectVertex;
    // FIXME pass smeared R9
    if( cur_type != 0 && doMCSmearing && cur_type < 100) { 
      float pth = Higgs.Pt();
      for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
        float rewei=1.;
        (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), zero_ ,zero_,0.);
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

    // Mass Resolution of the Event
    //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
    massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);

    // Make sure we know about the additional smearing category
//    massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
//    massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

    float vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
//    float sigmaMrv = massResolutionCalculator->massResolutionCorrVtx();
    float sigmaMrv = massResolutionCalculator->massResolutionEonly();
    float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
    float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
    // easy to calculate vertex probability from vtx mva output
    float vtxProb   = 1.-0.49*(vtx_mva+1.0);

    float diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second
        ,l.dipho_vtxind[diphoton_id]
        ,vtxProb,lead_p4,sublead_p4
        ,sigmaMrv,sigmaMwv,sigmaMeonly
        ,bdtTrainingPhilosophy.c_str());
    float phoid_mvaout_lead = l.photonIDMVA(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4 
        ,bdtTrainingPhilosophy.c_str());
    float phoid_mvaout_sublead = l.photonIDMVA(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4 
        ,bdtTrainingPhilosophy.c_str());

    bool isEBEB  = (lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);
    int category = GetBDTBoundaryCategory(diphobdt_output, isEBEB,VBFevent);

    assert( evweight >= 0. ); 

    l.FillCounter( "Accepted", weight );
    l.FillCounter( "Smeared", evweight );
    sumaccept += weight;
    sumsmear += evweight;

    // control plots 
    l.FillHist("all_mass",0, Higgs.M(), evweight);

    float rhofac=0.17;
    float rhofacbad=0.52;
    if( mass>=massMin && mass<=massMax  ) {

      l.FillHist("bdtout",0,diphobdt_output,evweight);
      if (fabs(lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442)){
        l.FillHist("bdtoutEB",0,diphobdt_output,evweight);
      } else {
        l.FillHist("bdtoutEE",0,diphobdt_output,evweight);
      }
      l.FillHist("phoid_mvaout_lead",0,phoid_mvaout_lead,evweight);
      l.FillHist("phoid_mvaout_sublead",0,phoid_mvaout_sublead,evweight);
      // Photon ID Input 
      if (l.pho_isEB) {
        l.FillHist("hoe",0,l.pho_hoe[diphoton_index.first],evweight);
        //l.FillHist("hoe",0,l.pho_hoe[diphoton_index.second],evweight);
        l.FillHist("sieie",0,l.pho_sieie[diphoton_index.first],evweight);
        //l.FillHist("sieie",0,l.pho_sieie[diphoton_index.second],evweight);
        l.FillHist("tiso1",0,(*(l.pho_tkiso_recvtx_030_002_0000_10_01))[diphoton_index.first][l.dipho_vtxind[diphoton_id]] + l.pho_ecalsumetconedr03[diphoton_index.first]+l.pho_hcalsumetconedr04[diphoton_index.first]-l.rho*rhofac,evweight);
        l.FillHist("tiso2",0,l.pho_tkiso_badvtx_040_002_0000_10_01[diphoton_index.first] +l.pho_ecalsumetconedr03[diphoton_index.first]+l.pho_hcalsumetconedr04[diphoton_index.first]-l.rho*rhofacbad,evweight);
        l.FillHist("tiso3",0,(*(l.pho_tkiso_recvtx_030_002_0000_10_01))[diphoton_index.first][l.dipho_vtxind[diphoton_id]],evweight);
        l.FillHist("sieip",0,l.pho_sieip[diphoton_index.first],evweight);
        l.FillHist("sipip",0,TMath::Sqrt(l.pho_sipip[diphoton_index.first]),evweight);
      }
      /*  
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

          l.FillHist("mass",selectioncategory+1, Higgs.M(), evweight);
          l.FillHist("pt",selectioncategory+1, Higgs.Pt(), evweight);
          l.FillHist("eta",selectioncategory+1, Higgs.Eta(), evweight);

          l.FillHist("pho_pt",selectioncategory+1,lead_p4.Pt(), evweight);
          l.FillHist("pho1_pt",selectioncategory+1,lead_p4.Pt(), evweight);
          l.FillHist("pho_eta",selectioncategory+1,lead_p4.Eta(), evweight);
          l.FillHist("pho1_eta",selectioncategory+1,lead_p4.Eta(), evweight);
          l.FillHist("pho_r9",selectioncategory+1, lead_r9, evweight);
          l.FillHist("pho1_r9",selectioncategory+1, lead_r9, evweight);

          l.FillHist("pho_pt",selectioncategory+1,sublead_p4.Pt(), evweight);
          l.FillHist("pho2_pt",selectioncategory+1,sublead_p4.Pt(), evweight);
          l.FillHist("pho_eta",selectioncategory+1,sublead_p4.Eta(), evweight);
          l.FillHist("pho2_eta",selectioncategory+1,sublead_p4.Eta(), evweight);
          l.FillHist("pho_r9",selectioncategory+1, sublead_r9, evweight);
          l.FillHist("pho1_r9",selectioncategory+1, sublead_r9, evweight);

          l.FillHist("pho_n",selectioncategory+1,l.pho_n, evweight);
       */
    }

      if (cur_type==0 && mass >= 100. && mass < 180.){
        eventListText <<"Type="<< cur_type <<  " Run=" << l.run << "  LS=" << l.lumis << "  Event=" << l.event << " BDTCAT=" << category << " ggM=" << mass << " gg_Pt=" << ptHiggs << " LeadPhotonPhoid=" <<phoid_mvaout_lead << " SubleadPhotonPhoid=" <<phoid_mvaout_sublead << " diphotonBDT=" << diphobdt_output << " photon1Eta=" << lead_p4.Eta() <<" photon2Eta="<<sublead_p4.Eta() << " sigmaMrv="<<sigmaMrv << " sigmaMwv=" << sigmaMwv << " photon1Pt="<<lead_p4.Pt()<<" photon2Pt="<<sublead_p4.Pt() << " vtxProb="<<vtxProb <<" cosDphi="<<TMath::Cos(lead_p4.Phi() - sublead_p4.Phi()) << " r9_1=" <<lead_r9 <<" r9_2=" <<sublead_r9 
<<" E1="<<lead_p4.E()<<" E2="<<sublead_p4.E()
<<" reresraw1="<<energyCorrectedError[diphoton_index.first]/energyCorrected[diphoton_index.first] <<" reresraw2="<<energyCorrectedError[diphoton_index.second]/energyCorrected[diphoton_index.second]
<<" reresraw1JBD="<<l.pho_regr_energyerr[diphoton_index.first]/l.pho_regr_energy[diphoton_index.first] <<" reresraw2JBD="<<l.pho_regr_energyerr[diphoton_index.second]/l.pho_regr_energy[diphoton_index.second]
<<" etcorecal1="<<l.pho_ecalsumetconedr03[diphoton_index.first] - 0.012*lead_p4.Et()
<<" etcorecal2="<<l.pho_ecalsumetconedr03[diphoton_index.second] - 0.012*sublead_p4.Et()
<<" etcorhcal1="<<l.pho_hcalsumetconedr03[diphoton_index.first] - 0.005*lead_p4.Et()
<<" etcorhcal2="<<l.pho_hcalsumetconedr03[diphoton_index.second] - 0.005*sublead_p4.Et()
<<" etcortrkiso1="<<l.pho_trksumpthollowconedr03[diphoton_index.first] - 0.002*lead_p4.Et()
<<" etcortrkiso2="<<l.pho_trksumpthollowconedr03[diphoton_index.second] - 0.002*sublead_p4.Et()
<<" etcortrkisoabs1="<<l.pho_trksumpthollowconedr03[diphoton_index.first] 
<<" etcortrkisoabs2="<<l.pho_trksumpthollowconedr03[diphoton_index.second]
<<" specialphoton1="<<photonInfoCollection[diphoton_index.first].isSphericalPhoton()	
<<" specialphoton2="<<photonInfoCollection[diphoton_index.second].isSphericalPhoton()
<<" specialphotongen1="<<l.CheckSphericalPhoton(diphoton_index.first)	
<<" specialphotongen2="<<l.CheckSphericalPhoton(diphoton_index.second)


<<" FileName="<<l.histFileName;
        eventListText << endl;

/*
	std::cout << "Event="<<l.event<< " selvtx="<< l.dipho_vtxind[diphoton_id] << " vtxmva=" << l.vtx_std_evt_mva->at(diphoton_id)<<std::endl;
	vtxAna_.setPairID(diphoton_index.first,diphoton_index.second);
	std::cout << " ptbal " << vtxAna_.ptbal(l.dipho_vtxind[diphoton_id]) << std::endl;
	std::cout << " ptasym " << vtxAna_.ptasym(l.dipho_vtxind[diphoton_id]) << std::endl;
	std::cout << " logsumpt2 " << vtxAna_.logsumpt2(l.dipho_vtxind[diphoton_id]) << std::endl;
	std::cout << " mva " << vtxAna_.mva(l.dipho_vtxind[diphoton_id]) << std::endl;
	
	std::cout << "The Rest -- "<<std::endl;
	for (int vid=0;vid<l.vtx_std_n;vid++){
	std::cout << " VTX -- "<< vid << std::endl;
	std::cout << " ptbal " << vtxAna_.ptbal(vid) << std::endl;
	std::cout << " ptasym " << vtxAna_.ptasym(vid) << std::endl;
	std::cout << " logsumpt2 " << vtxAna_.logsumpt2(vid) << std::endl;
	std::cout << " nconv " << vtxAna_.nconv(vid) << std::endl;
	std::cout << " mva " << vtxAna_.mva(vid) << std::endl;
	std::cout << " sumpt " << vtxAna_.sumpt(vid) << std::endl;
	std::cout << " limpulltoconv " << vtxAna_.limpulltoconv(vid) << std::endl;


	}

       for (int cv=0;cv<l.conv_n;cv++){
	std::cout << " conv_n" << cv <<std::endl;
	std::cout << " conv_p4.x" << ((TLorentzVector*)l.conv_p4->At(cv))->X() <<std::endl;
	std::cout << " conv_p4.y" << ((TLorentzVector*)l.conv_p4->At(cv))->Y() <<std::endl;
	std::cout << " conv_p4.z" << ((TLorentzVector*)l.conv_p4->At(cv))->Z() <<std::endl;
	std::cout << " conv_ntracks" << l.conv_ntracks[cv] <<std::endl;
	std::cout << " conv_eoverp" << l.conv_eoverp[cv] <<std::endl;
	std::cout << " conv_chi2_probability" << l.conv_chi2_probability[cv] <<std::endl;
	
	}

	std::cout << "Photon E1 =" <<lead_p4.E() <<std::endl;	
	std::cout << "Photon E2 =" <<sublead_p4.E() <<std::endl;	
	std::cout << "regr sigma1 =" << l.pho_regr_energyerr[diphoton_index.first]<<std::endl;	
	std::cout << "regr sigma2 =" << l.pho_regr_energyerr[diphoton_index.second]<<std::endl;	
	

*/
      }

    // --------------------------------------------------------------------------------------------- 
    if (cur_type == 0 ){
      l.rooContainer->InputDataPoint("data_mass",category,mass);
    }
    if (cur_type > 0 && cur_type != 3 && cur_type != 4){
      l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
    }
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
          int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);

          float mass = Higgs.M();
          float ptHiggs = Higgs.Pt();

          // Mass Resolution of the Event
          //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
    	  //massResolutionCalculator->Setup(l,photonLeadInfo,photonSubleadInfo,diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
    	  massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
    	  // Make sure we know about the additional smearing category
    	  //massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
    	  //massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

          float vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
          //float sigmaMrv = massResolutionCalculator->massResolutionCorrVtx();
          float sigmaMrv = massResolutionCalculator->massResolutionEonly();
          float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
          float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
          // easy to calculate vertex probability from vtx mva output
          float vtxProb   = 1.-0.49*(vtx_mva+1.0);

          float diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second
              ,l.dipho_vtxind[diphoton_id]
              ,vtxProb,lead_p4,sublead_p4
              ,sigmaMrv,sigmaMwv,sigmaMeonly
              ,bdtTrainingPhilosophy.c_str());

          bool isEBEB  = (lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);
          int category = GetBDTBoundaryCategory(diphobdt_output,isEBEB,VBFevent);

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

          int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);

	  float photon_idMVA1=l.photonIDMVA(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4,"MIT");
	  float photon_idMVA2=l.photonIDMVA(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4,"MIT");

          bool isEBEB  = (lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);
          for(std::vector<BaseDiPhotonSmearer *>::iterator sj=diPhotonSmearers_.begin(); sj!= diPhotonSmearers_.end(); ++sj ) {
            float swei=1.;
            float pth = Higgs.Pt();
            if( *si == *sj ) { 
              (*si)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), photon_idMVA1,photon_idMVA2,syst_shift);
            } else { 
              (*sj)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)) ,photon_idMVA1,photon_idMVA2,0.);
            }
            evweight *= swei;
          }
          float mass = Higgs.M();
          float ptHiggs = Higgs.Pt();

          // Mass Resolution of the Event
          //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
    	  //massResolutionCalculator->Setup(l,photonLeadInfo,photonSubleadInfo,diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
    	  massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
    	   // Make sure we know about the additional smearing category
    	  //massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
    	  //massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

          float vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
          //float sigmaMrv = massResolutionCalculator->massResolutionCorrVtx();
          float sigmaMrv = massResolutionCalculator->massResolutionEonly();
          float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
          float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
          // easy to calculate vertex probability from vtx mva output
          float vtxProb   = 1.-0.49*(vtx_mva+1.0);

          float diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second
              ,l.dipho_vtxind[diphoton_id]
              ,vtxProb,lead_p4,sublead_p4
              ,sigmaMrv,sigmaMwv,sigmaMeonly
              ,bdtTrainingPhilosophy.c_str()
	      ,photon_idMVA1,photon_idMVA2);

          int category = GetBDTBoundaryCategory(diphobdt_output,isEBEB,VBFevent);

          categories.push_back(category);
          mass_errors.push_back(mass);
          weights.push_back(evweight);
        }
        if (cur_type < 0){
          l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
        }
      }

    }


    // loop over the smearers included in the systematics study -> These can alter the selection so need to re-analyze
    for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
      std::vector<double> mass_errors;
      std::vector<double> weights;
      std::vector<int> categories;

      // loop over syst shift
      for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
        if( syst_shift == 0. ) { continue; } // skip the central value
        // smear the photons 
  	photonInfoCollection.clear(); // this is a member of PhotonAnalysis, just clear it
        for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
          std::vector<std::vector<bool> > p;
          //std::cout << "GF check: " <<  l.pho_residCorrEnergy[ipho] << "  " << l.pho_residCorrResn[ipho] << std::endl;
          PhotonReducedInfo phoInfo (// *((TVector3*)l.pho_calopos->At(ipho)), 
              *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
              ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
              energyCorrected[ipho],
              l.pho_isEB[ipho], l.pho_r9[ipho],
              l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_));
    	  if (l.CheckSphericalPhoton(ipho)) phoInfo.setSphericalPhoton(true);

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
	  photonInfoCollection.push_back(phoInfo);
        }

        // analyze the event
        // FIXME pass smeared R9
        // Get Ready for VBF Tagging
        bool VBFevent = false;

        int diphoton_id=-1;

          // VBF-TAGGING -------------------------------------------------------------------------- //
          // CP // NW Use Same pre-selected events but tag as VBF if pass Jet Requirements
          // Reset VBF event to False
          VBFevent = false;
          if(includeVBF) {
          int diphoton_id_vbf = l.DiphotonMITPreSelection(leadEtCutVBF,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 
	  if (diphoton_id_vbf >-1){
            diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id_vbf],  l.dipho_subleadind[diphoton_id_vbf] );

            TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
            TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
            float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id_vbf]];
            float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id_vbf]];
            TLorentzVector Higgs = lead_p4 + sublead_p4;   
            // JET MET Corrections // No need to reset the pointers again
            //PhotonAnalysis::RescaleJetEnergy(l);
            float jet1ptcut =0.0;
            float jet2ptcut =0.0;
            bool crosscheck = false;
            std::pair<int,int> highestPtJets(-1,-1);

            highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );
            bool VBFpresel = (highestPtJets.first!=-1)&&(highestPtJets.second!=-1);

            if(VBFpresel){

              TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
              TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
              TLorentzVector dijet = (*jet1) + (*jet2);

              myAllLeadJPt = jet1->Pt();
              myAllSubJPt = jet2->Pt();
              myAllLeadJEta = jet1->Eta();
              myAllSubJEta = jet2->Eta();
              myAll_Mjj = dijet.M();
              myAlldEta = fabs(jet1->Eta() - jet2->Eta());
              myAllZep  = fabs(Higgs.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
              myAlldPhi = fabs(Higgs.DeltaPhi(dijet));
              myAll_Mgg =Higgs.M();
              myAllPtHiggs =Higgs.Pt();

              myVBFLeadJPt = jet1->Pt();
              myVBFSubJPt = jet2->Pt();
              myVBF_Mjj = dijet.M();
              myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
              myVBFZep  = fabs(Higgs.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
              myVBFdPhi = fabs(Higgs.DeltaPhi(dijet));
              myVBF_Mgg =Higgs.M();


              // Cannot Get Apply cuts to work -> Need to discuss with C.Palmer, for now, jyst apply the cuts
              //l.ApplyCutsFill(0,3,evweight, myweight);
              //VBFevent = l.ApplyCutsFill(0,5,evweight, myweight);
              //VBFevent = l.ApplyCutsFill(0,1,evweight, myweight);
              if (myVBFLeadJPt  >30.){
                if (myVBFSubJPt >20.){
                  if (myVBF_Mjj   >350.){
                    if (myVBFdEta   >3.5){
                      if (myVBFdPhi   >2.6){
                        if (myVBFZep    <2.5){
                          VBFevent=true;
			  diphoton_id = diphoton_id_vbf;
                        }
                      }
                    }   
                  }
                }
              }

              /*
                 if(VBFevent) 
                 std::cout << setprecision(4) <<  "Run = " << l.run << "  LS = " << l.lumis <<
                 "  Event = " << l.event << "  SelVtx = " << l.dipho_vtxind[diphoton_id] 
                 << "  CAT4 = " << CAT4 << "  ggM = " << myVBF_Mgg << " ggPt =  " << myAllPtHiggs 
                 << "  jetEta1 = " << jet1->Eta() << "  jetEta2 = " << jet2->Eta()
                 << "  jetPhi1 = " << jet1->Phi() << "  jetPhi2 = " << jet2->Phi()
                 <<  "  jetEt1 = " << jet1->Et() << "  jetEt2 = "  << jet2->Et()
                 << " Mjj " << myVBF_Mjj
                 << " dEtajj " << myVBFdEta 
                 << " Zeppenfeld " << myVBFZep
                 << " dPhijjgg " << myVBFdPhi << " VBF itype " <<cur_type << std::endl;
               */
            }
          }
	}
          // CP // NW VBF Tagging
          // --------------------- END VBF-TAGGING --------------------------------------------------------//
	
	if (diphoton_id < 0){ // failed to find a VBF 
        if (bdtTrainingPhilosophy=="MIT"){
          diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 
        } else if (bdtTrainingPhilosophy=="UCSD"){
          diphoton_id = l.DiphotonCiCSelection(l.phoLOOSE, l.phoLOOSE, leadEtCut, subleadEtCut, nPhotonCategories_,applyPtoverM, &smeared_pho_energy[0] ); 
        }
	}

        if (diphoton_id > -1 ) {

          diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
          float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] *genLevWeight;

          TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
          TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
          TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
          TLorentzVector Higgs = lead_p4 + sublead_p4; 


          int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
          if( cur_type != 0 && doMCSmearing ) {
            for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
              float rewei=1.;
              float pth = Higgs.Pt();
              (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)),zero_,zero_,0.);
              evweight *= rewei;
            }
          }
          float mass = Higgs.M();
          float ptHiggs = Higgs.Pt();

          // Mass Resolution of the Event
          //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
    	   massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
    	  // Make sure we know about the additional smearing category
    	  //massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
   	  //massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

          float vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
          //float sigmaMrv = massResolutionCalculator->massResolutionCorrVtx();
          float sigmaMrv = massResolutionCalculator->massResolutionEonly();
          float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
          float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
          // easy to calculate vertex probability from vtx mva output
          float vtxProb   = 1.-0.49*(vtx_mva+1.0);

          float diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second
              ,l.dipho_vtxind[diphoton_id]
              ,vtxProb,lead_p4,sublead_p4
              ,sigmaMrv,sigmaMwv,sigmaMeonly
              ,bdtTrainingPhilosophy.c_str());

          bool isEBEB  = (lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);

          // Recalculate BDT since kinematics could change now.  

          int category = GetBDTBoundaryCategory(diphobdt_output,isEBEB,VBFevent);

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
void MassFactorizedMvaAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
  vtxAna_.setBranchAdresses(t,"vtx_std_");
  vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool MassFactorizedMvaAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
  return true;
}
// ----------------------------------------------------------------------------------------------------
double MassFactorizedMvaAnalysis::GetDifferentialKfactor(double gPT, int Mass)
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
  //  cout << "Gen Mass: "<< Mass << " Using "<<m1<< " " << m2<< " Hist name check " << hm1->GetName()<<" " <<hm2->GetName()<<endl;
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

void MassFactorizedMvaAnalysis::FillSignalLabelMap(){

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
      if (bdtout >= 0.05) return 4;
      else return -1;
    } else {
      if (bdtout >= 0.05 && bdtout < 0.55) return 3;
      if (bdtout >= 0.55 && bdtout < 0.72) return 2;
      if (bdtout >= 0.72 && bdtout < 0.89) return 1;
      if (bdtout >= 0.89) return 0;
      return -1;
    }

  } else std::cerr << "No BDT Philosophy known - " << bdtTrainingPhilosophy << std::endl;
}

std::string MassFactorizedMvaAnalysis::GetSignalLabel(int id){

  // For the lazy man, can return a memeber of the map rather than doing it yourself
  std::map<int,std::string>::iterator it = signalLabels.find(id);

  if (it!=signalLabels.end()){
    return it->second;

  } else { 

    std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
    return "NULL";
  }

}

void MassFactorizedMvaAnalysis::ResetAnalysis(){
  // Reset Random Variable on the EnergyResolution Smearer
  eResolSmearer->resetRandom();
}
