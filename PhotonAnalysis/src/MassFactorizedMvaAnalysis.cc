#include "../interface/MassFactorizedMvaAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
MassFactorizedMvaAnalysis::MassFactorizedMvaAnalysis()  : 
    name_("MassFactorizedMvaAnalysis")
                            //vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;

    forceStdPlotsOnZee = false;

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
        // if (dataIs2011) l.rooContainer->FitToData("data_pol_model_7TeV","data_mass");
        if (dataIs2011) l.rooContainer->FitToData("data_pol_model","data_mass");
        else l.rooContainer->FitToData("data_pol_model_8TeV","data_mass");  // Fit to full range of dataset
    }

    eventListText.close();
    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;
    
    // default categories: Jan16
    bdtCategoryBoundaries.push_back(-0.05);
    bdtCategoryBoundaries.push_back(0.5);
    bdtCategoryBoundaries.push_back(0.71);
    bdtCategoryBoundaries.push_back(0.88);
    bdtCategoryBoundaries.push_back(1.);
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
    FillSignalLabelMap(l);
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
        << "doPdfSmearerSyst "<< doPdfWeightSyst << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << std::endl;

    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
        l.histoContainer[ind].setScale(1.);
    }

    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

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
	setupEscaleSmearer();
    }
    if( doEresolSmear ) {
        // energy resolution smearing
    setupEresolSmearer();
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
    if(doPdfWeightSmear) {
        // PdfWeights efficiency (For now only consider QCD Scale Uncertainty 
        std::cerr << __LINE__ << std::endl; 
        pdfWeightSmearer = new PdfWeightSmearer( pdfWeightHist,"up","down");
        pdfWeightSmearer->name("pdfWeight");
        pdfWeightSmearer->init();
        genLevelSmearers_.push_back(pdfWeightSmearer);
    }
    if(doInterferenceSmear) {
        // interference efficiency
        std::cerr << __LINE__ << std::endl; 
        interferenceSmearer = new InterferenceSmearer(2.5e-2,0.);
        genLevelSmearers_.push_back(interferenceSmearer);
    }

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed
    // initialize the analysis variables
    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;

    std::sort(bdtCategoryBoundaries.begin(),bdtCategoryBoundaries.end(), std::greater<float>() );
    nInclusiveCategories_ = bdtCategoryBoundaries.size()-1;

    nVBFCategories   = ((int)includeVBF)*( (mvaVbfSelection && !multiclassVbfSelection) ? mvaVbfCatBoundaries.size()-1 : nVBFEtaCategories*nVBFDijetJetCategories );
    if(includeVHlep){
        nVHlepCategories = nElectronCategories + nMuonCategories;
    }
    if(includeVHmet){
    nVHmetCategories = nMetCategories;
    }
    
    std::sort(mvaVbfCatBoundaries.begin(),mvaVbfCatBoundaries.end(), std::greater<float>() );
    if (multiclassVbfSelection) {
        std::vector<int> vsize;
        vsize.push_back((int)multiclassVbfCatBoundaries0.size());
        vsize.push_back((int)multiclassVbfCatBoundaries1.size());
        vsize.push_back((int)multiclassVbfCatBoundaries2.size());
        std::sort(vsize.begin(),vsize.end(), std::greater<int>());
        // sanity check: there sould be at least 2 vectors with size==2
        if (vsize[0]<2 || vsize[1]<2 ){
            std::cout << "Not enough category boundaries:" << std::endl;
            std::cout << "multiclassVbfCatBoundaries0 size = " << multiclassVbfCatBoundaries0.size() << endl;
            std::cout << "multiclassVbfCatBoundaries1 size = " << multiclassVbfCatBoundaries1.size() << endl;
            std::cout << "multiclassVbfCatBoundaries2 size = " << multiclassVbfCatBoundaries2.size() << endl;
            assert( 0 );
        }
        nVBFCategories   = vsize[0]-1;
        cout << "@@@@@@@@@@@@@@@@@ 	nVBFCategories = " << 	nVBFCategories << endl;
        std::sort(multiclassVbfCatBoundaries0.begin(),multiclassVbfCatBoundaries0.end(), std::greater<float>() );
        std::sort(multiclassVbfCatBoundaries1.begin(),multiclassVbfCatBoundaries1.end(), std::greater<float>() );
        std::sort(multiclassVbfCatBoundaries2.begin(),multiclassVbfCatBoundaries2.end(), std::greater<float>() );
    }

    nCategories_=(nInclusiveCategories_+nVBFCategories+nVHlepCategories+nVHmetCategories);

    if (bdtTrainingPhilosophy == "UCSD") {
        l.rooContainer->SetNCategories(8);
    } else if (bdtTrainingPhilosophy == "MIT") {
        l.rooContainer->SetNCategories(nCategories_);
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
	setupEscaleSyst(l);
        //// systPhotonSmearers_.push_back( eScaleSmearer );
        //// std::vector<std::string> sys(1,eScaleSmearer->name());
        //// std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        //// l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
	setupEresolSyst(l);
        //// systPhotonSmearers_.push_back( eResolSmearer );
        //// std::vector<std::string> sys(1,eResolSmearer->name());
        //// std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        //// l.rooContainer->MakeSystematicStudy(sys,sys_t);
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
    if(doPdfWeightSmear && doPdfWeightSyst) {
        systGenLevelSmearers_.push_back(pdfWeightSmearer);
        std::vector<std::string> sys(1,pdfWeightSmearer->name());
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
    for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
        int sig = sigPointsToBook[isig];
        l.rooContainer->AddConstant(Form("XSBR_ggh_%d",sig),l.signalNormalizer->GetXsection(double(sig),"ggh")*l.signalNormalizer->GetBR(double(sig)));
    }
    /*
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
    */
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

    l.rooContainer->AddRealVar("pol0_8TeV",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol1_8TeV",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol2_8TeV",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol3_8TeV",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol4_8TeV",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol5_8TeV",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modpol0_8TeV","@0*@0","pol0_8TeV");
    l.rooContainer->AddFormulaVar("modpol1_8TeV","@0*@0","pol1_8TeV");
    l.rooContainer->AddFormulaVar("modpol2_8TeV","@0*@0","pol2_8TeV");
    l.rooContainer->AddFormulaVar("modpol3_8TeV","@0*@0","pol3_8TeV");
    l.rooContainer->AddFormulaVar("modpol4_8TeV","@0*@0","pol4_8TeV");
    l.rooContainer->AddFormulaVar("modpol5_8TeV","@0*@0","pol5_8TeV");

    if (bdtTrainingPhilosophy=="UCSD"){
        // UCSD BDT Categories

        std::vector<std::string> data_pol5_pars(6,"p");   
        data_pol5_pars[0] = "modpol0_8TeV";
        data_pol5_pars[1] = "modpol1_8TeV";
        data_pol5_pars[2] = "modpol2_8TeV";
        data_pol5_pars[3] = "modpol3_8TeV";
        data_pol5_pars[4] = "modpol4_8TeV";
        data_pol5_pars[5] = "modpol5_8TeV";
        l.rooContainer->AddGenericPdf("data_pol_model_8TeV","0","CMS_hgg_mass",data_pol5_pars,76); // >= 71 means RooBernstein of order >= 1

    } else if (bdtTrainingPhilosophy=="MIT"){
    
    // -----------------------------------------------------
    // Configurable background model
    // if no configuration was given, set some defaults
    std::string postfix=(dataIs2011?"":"_8TeV");
    if( bkgPolOrderByCat.empty() ) {
        for(int i=0; i<nCategories_; i++){
        if(i<1) {
            bkgPolOrderByCat.push_back(4);
        } else if(i<nInclusiveCategories_) {
            bkgPolOrderByCat.push_back(5);
        } else if(i<nInclusiveCategories_+nVBFCategories){
            bkgPolOrderByCat.push_back(3);
        } else if(i<nInclusiveCategories_+nVBFCategories+nVHlepCategories){
            bkgPolOrderByCat.push_back(-1);
                } else if(i<nInclusiveCategories_+nVBFCategories+nVHlepCategories+nVHmetCategories){
                    bkgPolOrderByCat.push_back(3);
                }
        }
    }
    // build the model
    buildBkgModel(l, postfix);
    }
    
    // Make some data sets from the observables to fill in the event loop      
    // Binning is for histograms (will also produce unbinned data sets)
    l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass"    ,nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass"     ,nDataBins);          
    
    // Create Signal DataSets:
    for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
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
    int sig = sigPointsToBook[isig];
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_ggh_mass_m%d",sig),-1);    
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_vbf_mass_m%d",sig),-1);    
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_wzh_mass_m%d",sig),-1);    
        l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_tth_mass_m%d",sig),-1);    
    }

    // Make sure the Map is filled
    FillSignalLabelMap(l);

    // Initialize all MVA ---------------------------------------------------//
    l.SetAllMVA();
    // UCSD
    l.tmvaReaderID_UCSD->BookMVA("Gradient"      ,photonLevelMvaUCSD.c_str()  );
    l.tmvaReader_dipho_UCSD->BookMVA("Gradient"  ,eventLevelMvaUCSD.c_str()   );
    // New ID MVA
    if( photonLevelNewIDMVA_EB != "" && photonLevelNewIDMVA_EE != "" ) {
    l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevelNewIDMVA_EB.c_str());
    l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevelNewIDMVA_EE.c_str());
    } else { 
    assert( dataIs2011 );
    }
    // MIT 
    if( photonLevelMvaMIT_EB != "" && photonLevelMvaMIT_EE != "" ) {
    l.tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost",photonLevelMvaMIT_EB.c_str());
    l.tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost",photonLevelMvaMIT_EE.c_str());
    } else {
    assert( ! dataIs2011 );
    }
    l.tmvaReader_dipho_MIT->BookMVA("Gradient"   ,eventLevelMvaMIT.c_str()    );
    // ----------------------------------------------------------------------//

    if(PADEBUG) 
        cout << "InitRealMassFactorizedMvaAnalysis END"<<endl;

    // FIXME book of additional variables
}



///// float MassFactorizedMvaAnalysis::getDiphoMva(LoopAll & l, int diphotonId, bool smear, float syst_shift) 
///// {
/////     massResolutionCalculator->Setup(l,&photonInfoCollection[l.dipho_leadind[diphotonId]],&photonInfoCollection[l.dipho_subleadind[diphotonId]],diphotonId,
///// 				    massResoPars,nR9Categories,nEtaCategories,beamspotSigma);
/////     float vtx_mva  = l.vtx_std_evt_mva->at(diphotonId);
/////     sigmaMrv = massResolutionCalculator->massResolutionEonly();
/////     sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
/////     float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
/////     // easy to calculate vertex probability from vtx mva output
/////     float vtxProb   = 1.-0.49*(vtx_mva+1.0); /// should better use this: vtxAna_.setPairID(diphotonId); vtxAna_.vertexProbability(vtx_mva); PM
/////     
/////     float phoid_mvaout_lead = ( dataIs2011 ? 
///// 				l.photonIDMVA(l.dipho_leadind[diphotonId],l.dipho_vtxind[diphotonId],
///// 					      lead_p4,bdtTrainingPhilosophy.c_str()) :
///// 				l.photonIDMVANew(l.dipho_leadind[diphotonId],l.dipho_vtxind[diphotonId],
///// 						 lead_p4,bdtTrainingPhilosophy.c_str()));
/////     float phoid_mvaout_sublead = ( dataIs2011 ? 
///// 				   l.photonIDMVA(l.dipho_subleadind[diphotonId],l.dipho_vtxind[diphotonId],
///// 						 sublead_p4,bdtTrainingPhilosophy.c_str()) : 
///// 				   l.photonIDMVANew(l.dipho_subleadind[diphotonId],l.dipho_vtxind[diphotonId],
///// 						    sublead_p4,bdtTrainingPhilosophy.c_str()));
/////     // apply di-photon level smearings and corrections
/////     int selectioncategory = l.DiphotonCategory(l.dipho_leadind[diphotonId],l.dipho_subleadind[diphotonId],Higgs.Pt(),nEtaCategories,nR9Categories,0);
/////     if( smear && ur_type != 0 ) {
///// 	applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, 
///// 			       phoid_mvaout_lead,phoid_mvaout_sublead,
///// 			       diPhoSys, syst_shift);
/////     }
/////     
/////     return l.diphotonMVA(l.dipho_leadind[diphotonId],l.dipho_subleadind[diphotonId],l.dipho_vtxind[diphotonId] ,
///// 			 vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,
///// 			 bdtTrainingPhilosophy.c_str(),
///// 			 phoid_mvaout_lead,phoid_mvaout_sublead);
///// }

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
        
    int mu_ind=-1;
    int el_ind=-1;

    int muVtx=-1;
    int elVtx=-1;
    
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
	vbfIjet1=-1, vbfIjet2=-1;
        VBFevent = false;
        
        int diphotonVHlep_id = -1;
        VHmuevent = false;
        VHmuevent_cat=0;
        VHelevent = false;
        VHelevent_cat=0;
        
        int diphotonVHmet_id = -1; 
        VHmetevent = false; 
        VHmetevent_cat=0; 

        int diphotonVHhad_id = -1;
        VHhadevent = false;
       

        if(includeVHlep){
            float eventweight = weight * genLevWeight;
            //float eventweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
            float myweight=1.;
            if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
            VHmuevent=MuonTag2012B(l, diphotonVHlep_id, mu_ind, muVtx, VHmuevent_cat, &smeared_pho_energy[0], lep_sync, true, phoidMvaCut, eventweight,  smeared_pho_weight, !isSyst);
            //ZWithFakeGammaCS(l, &smeared_pho_energy[0]);
            ElectronStudies2012B(l, &smeared_pho_energy[0], true,  phoidMvaCut, eventweight, myweight, jentry);
            int diphotonVH_ele_id=-1;
            VHelevent=ElectronTag2012B(l, diphotonVH_ele_id, el_ind, elVtx, VHelevent_cat, &smeared_pho_energy[0], lep_sync, true, phoidMvaCut, eventweight,  smeared_pho_weight, !isSyst);
            // FIXME  need to un-flag events failing the diphoton mva cut.
            
            if(!VHmuevent && VHelevent){
                diphotonVHlep_id=diphotonVH_ele_id;
            }
        }

        if(includeVHmet && !dataIs2011) {
	    //	    std::cout << "+++PFMET UNCORR " << l.met_pfmet << std::endl;
            if(!isSyst) VHmetevent=METTag2012B(l, diphotonVHmet_id, VHmetevent_cat, &smeared_pho_energy[0], met_sync, true, phoidMvaCut, false); 
            if(isSyst)  VHmetevent=METTag2012B(l, diphotonVHmet_id, VHmetevent_cat, &smeared_pho_energy[0], met_sync, true, phoidMvaCut, true); 
            // FIXME  need to un-flag events failing the diphoton mva cut.
        }

        // VBF
        if((includeVBF || includeVHhad)&&l.jet_algoPF1_n>1 && !isSyst /*avoid rescale > once*/) {
            l.RescaleJetEnergy();
        }

        if(includeVBF) {   
            diphotonVBF_id = l.DiphotonMITPreSelection(leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] );
            float eventweight = weight * smeared_pho_weight[l.dipho_leadind[diphotonVBF_id]] * smeared_pho_weight[l.dipho_subleadind[diphotonVBF_id]] * genLevWeight;
            float myweight=1.;
            if(eventweight*sampleweight!=0) myweight=eventweight/sampleweight;
            
            VBFevent= ( dataIs2011 ? 
                VBFTag2011(l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight) :
                VBFTag2012(vbfIjet1, vbfIjet2, l, diphotonVBF_id, &smeared_pho_energy[0], true, eventweight, myweight) );

        }
    
        if(includeVHlep&&VHmuevent){
            diphoton_id = diphotonVHlep_id;
        } else if (includeVHlep&&VHelevent){
            diphoton_id = diphotonVHlep_id;
        } else if(includeVBF&&VBFevent) {
            diphoton_id = diphotonVBF_id;
        } else if(includeVHmet&&VHmetevent) {
            diphoton_id = diphotonVHmet_id;
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
        
        // apply beamspot reweighting if necessary
        if(reweighBeamspot && cur_type!=0) {
            evweight*=BeamspotReweight(vtx->Z(),((TVector3*)l.gv_pos->At(0))->Z());
        }
    
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
        massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,
					massResoPars,nR9Categories,nEtaCategories,beamspotSigma);
        float vtx_mva  = l.vtx_std_evt_mva->at(diphoton_id);
        sigmaMrv = massResolutionCalculator->massResolutionEonly();
        sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
        float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
        // easy to calculate vertex probability from vtx mva output
        float vtxProb   = 1.-0.49*(vtx_mva+1.0); /// should better use this: vtxAna_.setPairID(diphoton_id); vtxAna_.vertexProbability(vtx_mva); PM

        float phoid_mvaout_lead = ( dataIs2011 ? 
                    l.photonIDMVA(diphoton_index.first,l.dipho_vtxind[diphoton_id],
                          lead_p4,bdtTrainingPhilosophy.c_str()) :
                    l.photonIDMVANew(diphoton_index.first,l.dipho_vtxind[diphoton_id],
                             lead_p4,bdtTrainingPhilosophy.c_str()) );
        float phoid_mvaout_sublead = ( dataIs2011 ? 
                       l.photonIDMVA(diphoton_index.second,l.dipho_vtxind[diphoton_id],
                             sublead_p4,bdtTrainingPhilosophy.c_str()) : 
                       l.photonIDMVANew(diphoton_index.second,l.dipho_vtxind[diphoton_id],
					sublead_p4,bdtTrainingPhilosophy.c_str()));
	// apply di-photon level smearings and corrections
        int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
        if( cur_type != 0 && doMCSmearing ) {
	    applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, 
				   phoid_mvaout_lead,phoid_mvaout_sublead,
				   diPhoSys, syst_shift);
            isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
        }
	
        // Must be calculated after photon id has potentially been smeared
        //fillTrainTree(l,diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id] ,
	/// vtxProb,lead_p4,sublead_p4 ,sigmaMrv,sigmaMwv,sigmaMeonly ,bdtTrainingPhilosophy.c_str() ,phoid_mvaout_lead,phoid_mvaout_sublead);
        float diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id] ,
                          vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,
                          bdtTrainingPhilosophy.c_str(),
                          phoid_mvaout_lead,phoid_mvaout_sublead);
        kinematic_bdtout = diphobdt_output;

        bool isEBEB  = fabs(lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);
        category = GetBDTBoundaryCategory(diphobdt_output,isEBEB,VBFevent);
        if (diphobdt_output>=bdtCategoryBoundaries.back()) { 
            computeExclusiveCategory(l,category,diphoton_index,Higgs.Pt()); 
        }

        if (fillOptTree) {
            std::string name;
            if (genSys != 0)
            name = genSys->name();
            if (phoSys != 0)
                name = phoSys->name();
            if (diPhoSys != 0)
                name = diPhoSys->name();

            if (!isSyst)
            fillOpTree(l, lead_p4, sublead_p4, vtxProb, diphoton_index, diphoton_id, phoid_mvaout_lead, phoid_mvaout_sublead, weight, 
                    mass, sigmaMrv, sigmaMwv, Higgs, diphobdt_output, category, VBFevent, myVBF_Mjj, myVBFLeadJPt, 
                    myVBFSubJPt, nVBFDijetJetCategories, isSyst, "no-syst");
            else
            fillOpTree(l, lead_p4, sublead_p4, vtxProb, diphoton_index, diphoton_id, phoid_mvaout_lead, phoid_mvaout_sublead, weight, 
                    mass, sigmaMrv, sigmaMwv, Higgs, diphobdt_output, category, VBFevent, myVBF_Mjj, myVBFLeadJPt, 
                    myVBFSubJPt, nVBFDijetJetCategories, isSyst, name);

        }

        if (PADEBUG) std::cout << " Diphoton Category " <<category <<std::endl;
        // sanity check
        assert( evweight >= 0. ); 

        // dump BS trees in requested
        if (!isSyst && cur_type!=0 && saveBSTrees_) saveBSTrees(l,evweight,category,Higgs, vtx, (TVector3*)l.gv_pos->At(0),diphobdt_output);

        // save trees for unbinned datacards
        if (!isSyst && cur_type<0 && saveDatacardTrees_) saveMassFacDatCardTree(l,cur_type,category, evweight, diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,bdtTrainingPhilosophy.c_str(),phoid_mvaout_lead,phoid_mvaout_sublead);

        // save trees for IC spin analysis
        if (!isSyst && saveSpinTrees_) saveSpinTree(l,category,evweight,Higgs,lead_p4,sublead_p4,diphoton_index.first,diphoton_index.second,diphobdt_output,sigmaMrv,sigmaMwv,massResolutionCalculator->leadPhotonResolution(),massResolutionCalculator->leadPhotonResolutionNoSmear(),massResolutionCalculator->subleadPhotonResolution(),massResolutionCalculator->subleadPhotonResolutionNoSmear(),vtxProb,phoid_mvaout_lead,phoid_mvaout_sublead);






        // fill control plots and counters
        if( ! isSyst ) {
            
            l.FillCounter( "Accepted", weight );
            l.FillCounter( "Smeared", evweight );
            sumaccept += weight;
            sumsmear += evweight;
            if (l.runZeeValidation && !forceStdPlotsOnZee) {
                fillZeeControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, phoid_mvaout_lead, phoid_mvaout_sublead, diphobdt_output, sigmaMrv, sigmaMwv, vtxProb, diphoton_id, category, selectioncategory, evweight, l );
            } else {
	            if (category>-1) {
	    	        fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9,  diphoton_id, 
	    			     category, isCorrectVertex, evweight, vtx, l, muVtx, mu_ind, elVtx, el_ind, diphobdt_output ); 
                } else {
	    	        fillControlPlots(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9,  diphoton_id, 
	    			     -10, isCorrectVertex, evweight, vtx, l, muVtx, mu_ind, elVtx, el_ind, diphobdt_output ); 
                }
            }
            if (fillEscaleTrees) fillEscaleTree(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, phoid_mvaout_lead, phoid_mvaout_sublead, diphobdt_output, sigmaMrv, sigmaMwv, vtxProb, diphoton_id, category, selectioncategory, evweight, l );
        }
    
        if (dumpAscii && category > -1 && mass >= 100. && mass < 180. && !isSyst){
            // New ascii event list for syncrhonizing MVA Preselection + Diphoton MVA
            eventListText <<"type:"<< cur_type 
              << "    run:"   <<  l.run
              << "    lumi:"  <<  l.lumis
              << "    event:" <<  l.event
              << "    rho:" <<    l.rho
        // Preselection Lead
              << "    r9_1:"  <<  lead_r9
              << "    sceta_1:"   << (photonInfoCollection[diphoton_index.first]).caloPosition().Eta() 
              << "    hoe_1:" <<  l.pho_hoe[diphoton_index.first]
              << "    sigieie_1:" <<  l.pho_sieie[diphoton_index.first]
              << "    ecaliso_1:" <<  l.pho_ecalsumetconedr03[diphoton_index.first] - 0.012*lead_p4.Et()
              << "    hcaliso_1:" <<  l.pho_hcalsumetconedr03[diphoton_index.first] - 0.005*lead_p4.Et()
              << "    trckiso_1:" <<  l.pho_trksumpthollowconedr03[diphoton_index.first] - 0.002*lead_p4.Et()
              << "    chpfiso2_1:" <<  (*l.pho_pfiso_mycharged02)[diphoton_index.first][l.dipho_vtxind[diphoton_id]] 
              << "    chpfiso3_1:" <<  (*l.pho_pfiso_mycharged03)[diphoton_index.first][l.dipho_vtxind[diphoton_id]] 
        // Preselection SubLead
              << "    r9_2:"  <<  sublead_r9
              << "    sceta_2:"   << (photonInfoCollection[diphoton_index.second]).caloPosition().Eta() 
              << "    hoe_2:" <<  l.pho_hoe[diphoton_index.second]
              << "    sigieie_2:" <<  l.pho_sieie[diphoton_index.second]
              << "    ecaliso_2:" <<  l.pho_ecalsumetconedr03[diphoton_index.second] - 0.012*sublead_p4.Et()
              << "    hcaliso_2:" <<  l.pho_hcalsumetconedr03[diphoton_index.second] - 0.005*sublead_p4.Et()
              << "    trckiso_2:" <<  l.pho_trksumpthollowconedr03[diphoton_index.second] - 0.002*sublead_p4.Et()
              << "    chpfiso2_2:" <<  (*l.pho_pfiso_mycharged02)[diphoton_index.second][l.dipho_vtxind[diphoton_id]] 
              << "    chpfiso3_2:" <<  (*l.pho_pfiso_mycharged03)[diphoton_index.second][l.dipho_vtxind[diphoton_id]] 
        // Diphoton MVA inputs
              << "    ptH:"  <<  ptHiggs 
              << "    phoid_1:"   <<  phoid_mvaout_lead 
              << "    phoid_2:"   <<  phoid_mvaout_sublead
              << "    phoeta_1:"  <<  lead_p4.Eta() 
              << "    phoeta_2:"  <<  sublead_p4.Eta() 
              << "    sigmrv:"    <<  sigmaMrv 
              << "    bsw:"       <<  beamspotWidth 
              << "    sigmwv:"    <<  sigmaMwv 
              << "    pt_1/m:"      <<  lead_p4.Pt()/mass
              << "    pt_2/m:"      <<  sublead_p4.Pt()/mass
              << "    vtxprob:"   <<  vtxProb 
              << "    cosdphi:"   <<  TMath::Cos(lead_p4.Phi() - sublead_p4.Phi()) 
              << "    diphoBDT:"  <<  diphobdt_output 
        // Extra
              << "    mgg:"       <<  mass 
              << "    e_1:"       <<  lead_p4.E()
              << "    e_2:"       <<  sublead_p4.E()
              << "    eerr_1:"    <<  massResolutionCalculator->leadPhotonResolutionNoSmear()
              << "    eerr_2:"    <<  massResolutionCalculator->subleadPhotonResolutionNoSmear()
              << "    eerrsmeared_1:" << massResolutionCalculator->leadPhotonResolution()   
              << "    eerrsmeared_2:" << massResolutionCalculator->subleadPhotonResolution()
              << "    vbfevent:"  <<  VBFevent
              << "    muontag:"   <<  VHmuevent
              << "    electrontag:"<< VHelevent
              << "    mettag:"    <<  VHmetevent
              << "    evcat:"     <<  category
              << "    FileName:"  <<  l.files[l.current]
        // VBF MVA
              << "    vbfmva: "   <<  myVBF_MVA;
        // Vertex MVA
            vtxAna_.setPairID(diphoton_id);
            std::vector<int> & vtxlist = l.vtx_std_ranked_list->at(diphoton_id);
            // Conversions
            PhotonInfo p1 = l.fillPhotonInfos(l.dipho_leadind[diphoton_id], vtxAlgoParams.useAllConversions, 0); // WARNING using default photon energy: it's ok because we only re-do th$
            PhotonInfo p2 = l.fillPhotonInfos(l.dipho_subleadind[diphoton_id], vtxAlgoParams.useAllConversions, 0); // WARNING using default photon energy: it's ok because we only re-do$
            int convindex1 = l.matchPhotonToConversion(diphoton_index.first,vtxAlgoParams.useAllConversions);
            int convindex2 = l.matchPhotonToConversion(diphoton_index.second,vtxAlgoParams.useAllConversions);

            for(size_t ii=0; ii<3; ++ii ) {
                eventListText << "\tvertexId"<< ii+1 <<":" << (ii < vtxlist.size() ? vtxlist[ii] : -1);
            }
            for(size_t ii=0; ii<3; ++ii ) {
                eventListText << "\tvertexMva"<< ii+1 <<":" << (ii < vtxlist.size() ? vtxAna_.mva(vtxlist[ii]) : -2.);
            }
            for(size_t ii=1; ii<3; ++ii ) {
                eventListText << "\tvertexdeltaz"<< ii+1 <<":" << (ii < vtxlist.size() ? vtxAna_.vertexz(vtxlist[ii])-vtxAna_.vertexz(vtxlist[0]) : -999.);
            }
            eventListText << "\tptbal:"   << vtxAna_.ptbal(vtxlist[0])
                          << "\tptasym:"  << vtxAna_.ptasym(vtxlist[0])
                          << "\tlogspt2:" << vtxAna_.logsumpt2(vtxlist[0])
                          << "\tp2conv:"  << vtxAna_.pulltoconv(vtxlist[0])
                          << "\tnconv:"   << vtxAna_.nconv(vtxlist[0]);

        //Photon IDMVA inputs
            double pfchargedisobad03=0.;
            for(int ivtx=0; ivtx<l.vtx_std_n; ivtx++) {
                pfchargedisobad03 = l.pho_pfiso_mycharged03->at(diphoton_index.first).at(ivtx) > pfchargedisobad03 ? l.pho_pfiso_mycharged03->at(diphoton_index.first).at(ivtx) : pfchargedisobad03;
            }

            eventListText << "\tscetawidth_1: " << l.pho_etawidth[diphoton_index.first]
                          << "\tscphiwidth_1: " << l.sc_sphi[diphoton_index.first]
                          << "\tsieip_1: " << l.pho_sieip[diphoton_index.first]
                          << "\tbc_e2x2_1: " << l.pho_e2x2[diphoton_index.first]
                          << "\tpho_e5x5_1: " << l.bc_s25[l.sc_bcseedind[l.pho_scind[diphoton_index.first]]]
                          << "\ts4ratio_1: " << l.pho_s4ratio[diphoton_index.first]
                          << "\tpfphotoniso03_1: " << l.pho_pfiso_myphoton03[diphoton_index.first]
                          << "\tpfchargedisogood03_1: " << l.pho_pfiso_mycharged03->at(diphoton_index.first).at(vtxlist[0])
                          << "\tpfchargedisobad03_1: " << pfchargedisobad03
                          << "\teseffsigmarr_1: " << l.pho_ESEffSigmaRR[diphoton_index.first];
            pfchargedisobad03=0.;
            for(int ivtx=0; ivtx<l.vtx_std_n; ivtx++) {
                pfchargedisobad03 = l.pho_pfiso_mycharged03->at(diphoton_index.second).at(ivtx) > pfchargedisobad03 ? l.pho_pfiso_mycharged03->at(diphoton_index.second).at(ivtx) : pfchargedisobad03;
            }

            eventListText << "\tscetawidth_2: " << l.pho_etawidth[diphoton_index.second]
                          << "\tscphiwidth_2: " << l.sc_sphi[diphoton_index.second]
                          << "\tsieip_2: " << l.pho_sieip[diphoton_index.second]
                          << "\tbc_e2x2_2: " << l.pho_e2x2[diphoton_index.second]
                          << "\tpho_e5x5_2: " << l.bc_s25[l.sc_bcseedind[l.pho_scind[diphoton_index.second]]]
                          << "\ts4ratio_2: " << l.pho_s4ratio[diphoton_index.second]
                          << "\tpfphotoniso03_2: " << l.pho_pfiso_myphoton03[diphoton_index.second]
                          << "\tpfchargedisogood03_2: " << l.pho_pfiso_mycharged03->at(diphoton_index.second).at(vtxlist[0])
                          << "\tpfchargedisobad03_2: " << pfchargedisobad03
                          << "\teseffsigmarr_2: " << l.pho_ESEffSigmaRR[diphoton_index.second];


            if (convindex1!=-1) {
                eventListText 
                << "    convVtxZ1:"  <<  vtxAna_.vtxZFromConv(p1)
                << "    convVtxdZ1:"  <<  vtxAna_.vtxZFromConv(p1)-vtxAna_.vertexz(vtxlist[0])
                << "    convRes1:"   << vtxAna_.vtxdZFromConv(p1) 
                << "    convChiProb1:"  <<  l.conv_chi2_probability[convindex1]
                << "    convNtrk1:"  <<  l.conv_ntracks[convindex1]
                << "    convindex1:"  <<  convindex1
                << "    convvtxZ1:" << ((TVector3*) l.conv_vtx->At(convindex1))->Z()
                << "    convvtxR1:" << ((TVector3*) l.conv_vtx->At(convindex1))->Perp()
                << "    convrefittedPt1:" << ((TVector3*) l.conv_refitted_momentum->At(convindex1))->Pt();
            } else {
                eventListText 
                << "    convVtxZ1:"  <<  -999
                << "    convVtxdZ1:"  <<  -999
                << "    convRes1:"    <<  -999
                << "    convChiProb1:"  <<  -999
                << "    convNtrk1:"  <<  -999
                << "    convindex1:"  <<  -999
                << "    convvtxZ1:" << -999
                << "    convvtxR1:" << -999
                << "    convrefittedPt1:" << -999;
            }
            if (convindex2!=-1) {
                eventListText 
                << "    convVtxZ2:"  <<  vtxAna_.vtxZFromConv(p2)
                << "    convVtxdZ2:"  <<  vtxAna_.vtxZFromConv(p2)-vtxAna_.vertexz(vtxlist[0])
                << "    convRes2:"   << vtxAna_.vtxdZFromConv(p2)
                << "    convChiProb2:"  <<  l.conv_chi2_probability[convindex2]
                << "    convNtrk2:"  <<  l.conv_ntracks[convindex2]
                << "    convindex2:"  <<  convindex2
                << "    convvtxZ2:" << ((TVector3*) l.conv_vtx->At(convindex2))->Z()
                << "    convvtxR2:" << ((TVector3*) l.conv_vtx->At(convindex2))->Perp()
                << "    convrefittedPt2:" << ((TVector3*) l.conv_refitted_momentum->At(convindex2))->Pt();
            } else {
                eventListText 
                << "    convVtxZ2:"  <<  -999
                << "    convVtxdZ2:"  <<  -999
                << "    convRes2:"    <<  -999
                << "    convChiProb2:"  <<  -999
                << "    convNtrk2:"  <<  -999
                << "    convindex2:"  <<  -999
                << "    convvtxZ2:" << -999
                << "    convvtxR2:" << -999
                << "    convrefittedPt2:" << -999;
            }
            
            if(VHelevent){
                TLorentzVector* myel = (TLorentzVector*) l.el_std_p4->At(el_ind);
                TLorentzVector* myelsc = (TLorentzVector*) l.el_std_sc->At(el_ind);
                float thiseta = fabs(myelsc->Eta());
                double Aeff=0.;
                if(thiseta<1.0)                   Aeff=0.135;
                if(thiseta>=1.0 && thiseta<1.479) Aeff=0.168;
                if(thiseta>=1.479 && thiseta<2.0) Aeff=0.068;
                if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.116;
                if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.162;
                if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.241;
                if(thiseta>=2.4)                  Aeff=0.23;
                float thisiso=l.el_std_pfiso_charged[el_ind]+std::max(l.el_std_pfiso_neutral[el_ind]+l.el_std_pfiso_photon[el_ind]-l.rho*Aeff,0.);
    
                TLorentzVector elpho1=*myel + lead_p4;
                TLorentzVector elpho2=*myel + sublead_p4;

                eventListText 
                    << "    elind:"<<       el_ind
                    << "    elpt:"<<        myel->Pt()
                    << "    eleta:"<<       myel->Eta()
                    << "    elsceta:"<<     myelsc->Eta()
                    << "    elmva:"<<       l.el_std_mva_nontrig[el_ind]
                    << "    eliso:"<<       thisiso
                    << "    elisoopt:"<<    thisiso/myel->Pt()
                    << "    elAeff:"<<      Aeff
                    << "    eldO:"<<        fabs(l.el_std_D0Vtx[el_ind][elVtx])
                    << "    eldZ:"<<        fabs(l.el_std_DZVtx[el_ind][elVtx])
                    << "    elmishits:"<<   l.el_std_hp_expin[el_ind]
                    << "    elconv:"<<      l.el_std_conv[el_ind]
                    << "    eldr1:"<<       myel->DeltaR(lead_p4)
                    << "    eldr2:"<<       myel->DeltaR(sublead_p4)
                    << "    elmeg1:"<<      elpho1.M()
                    << "    elmeg2:"<<      elpho2.M();

            } else {
                eventListText 
                    << "    elind:"<<       -1
                    << "    elpt:"<<        -1
                    << "    eleta:"<<       -1
                    << "    elsceta:"<<     -1
                    << "    elmva:"<<       -1
                    << "    eliso:"<<       -1
                    << "    elisoopt:"<<    -1
                    << "    elAeff:"<<      -1
                    << "    eldO:"<<        -1
                    << "    eldZ:"<<        -1
                    << "    elmishits:"<<   -1
                    << "    elconv:"<<      -1
                    << "    eldr1:"<<       -1
                    << "    eldr2:"<<       -1
                    << "    elmeg1:"<<      -1
                    << "    elmeg2:"<<      -1;
            }

            if(VHmetevent){
                TLorentzVector myMet = l.METCorrection2012B(lead_p4, sublead_p4); 
                float corrMet    = myMet.Pt();
                float corrMetPhi = myMet.Phi();

                eventListText 
                    << "    metuncor:"<<        l.met_pfmet
                    << "    metphiuncor:"<<     l.met_phi_pfmet
                    << "    metcor:"<<          corrMet
                    << "    metphicor:"<<       corrMetPhi;
            } else {
                eventListText 
                    << "    metuncor:"<<        -1
                    << "    metphiuncor:"<<     -1
                    << "    metcor:"<<          -1
                    << "    metphicor:"<<       -1;
            }

            eventListText << endl;
        }
	return (l.runZeeValidation || (saveSpinTrees_ && mass>=massMin && mass<=massMax) || (category >= 0 && mass>=massMin && mass<=massMax));
    }
    return false;
}


// ----------------------------------------------------------------------------------------------------
void MassFactorizedMvaAnalysis::fillEscaleTree(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
                           const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
					       float phoid_mvaout_lead, float phoid_mvaout_sublead, 
					       float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
					       int diphoton_id, int category, int selectioncategory, float evweight, LoopAll & l )
{

    int lead = l.dipho_leadind[diphoton_id];
    int sublead = l.dipho_subleadind[diphoton_id];

    float mass = Higgs.M();
    float ptHiggs = Higgs.Pt();

    float cos_dphi = TMath::Cos(lead_p4.Phi()-sublead_p4.Phi());
    float pho1_sigmaE = energyCorrectedError[lead];
    float pho2_sigmaE = energyCorrectedError[sublead];

    float pho1_e5x5  = l.pho_e5x5[l.pho_scind[lead]];
    float pho2_e5x5  = l.pho_e5x5[l.pho_scind[sublead]];

    float pho1_energy_noregr = ((TLorentzVector*)l.pho_p4->At(lead))->Energy();
    float pho2_energy_noregr = ((TLorentzVector*)l.pho_p4->At(sublead))->Energy();

    float pho1_scEnergy  = ((TLorentzVector *)l.sc_p4->At(l.pho_scind[lead]))->Energy();
    float pho1_scEraw  = l.sc_raw[l.pho_scind[lead]];
    float pho1_scEpresh  = l.sc_pre[l.pho_scind[lead]];
    float pho1_sceta  = ((TVector3 *)l.sc_xyz->At(l.pho_scind[lead]))->Eta();
    float pho1_scphi  = ((TVector3 *)l.sc_xyz->At(l.pho_scind[lead]))->Phi();

    float pho2_scEnergy  = ((TLorentzVector *)l.sc_p4->At(l.pho_scind[sublead]))->Energy();
    float pho2_scEraw  = l.sc_raw[l.pho_scind[sublead]];
    float pho2_scEpresh  = l.sc_pre[l.pho_scind[sublead]];
    float pho2_sceta  = ((TVector3 *)l.sc_xyz->At(l.pho_scind[sublead]))->Eta();
    float pho2_scphi  = ((TVector3 *)l.sc_xyz->At(l.pho_scind[sublead]))->Phi();

    float el1_eta=-99.;
    float el1_phi=-99.;
    float el2_eta=-99.;
    float el2_phi=-99.;
    for (int iel=0;iel<l.el_std_n;iel++){
	if (l.el_std_scind[iel] == l.pho_scind[lead]) {
	    el1_eta  = ((TLorentzVector *)l.el_std_p4->At(iel))->Eta();
	    el1_phi  = ((TLorentzVector *)l.el_std_p4->At(iel))->Phi();
	}
	if (l.el_std_scind[iel] == l.pho_scind[sublead]) {
	    el2_eta  = ((TLorentzVector *)l.el_std_p4->At(iel))->Eta();
	    el2_phi  = ((TLorentzVector *)l.el_std_p4->At(iel))->Phi();
	}
    }

    TVector3* vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

    l.FillTree("run",l.run);
    l.FillTree("event",l.event);
    l.FillTree("lumi",l.lumis);
    l.FillTree("weight", evweight);
    l.FillTree("category",category);
    l.FillTree("category_baseline",selectioncategory);

    l.FillTree("pho1_energy",lead_p4.E());
    l.FillTree("pho1_energy_regr",l.pho_regr_energy[lead]);
    l.FillTree("pho1_energy_noregr",pho1_energy_noregr);
    l.FillTree("pho1_pt",lead_p4.Pt());
    l.FillTree("pho1_eta",lead_p4.Eta());
    l.FillTree("pho1_phi",lead_p4.Phi());
    l.FillTree("pho1_r9",lead_r9);
    l.FillTree("pho1_phoidMva",phoid_mvaout_lead);
    l.FillTree("pho1_ptOverM",lead_p4.Pt()/mass);
    l.FillTree("pho1_sceta",pho1_sceta);
    l.FillTree("pho1_scphi",pho1_scphi);
    l.FillTree("pho1_e5x5",pho1_e5x5);
    l.FillTree("pho1_sigmaE",pho1_sigmaE);
    l.FillTree("pho1_scEnergy",pho1_scEnergy);
    l.FillTree("pho1_scEraw",pho1_scEraw);
    l.FillTree("pho1_scEpresh",pho1_scEpresh);
    l.FillTree("pho1_isconv",l.pho_isconv[lead]);

    l.FillTree("pho2_energy",sublead_p4.E());
    l.FillTree("pho2_energy_regr",l.pho_regr_energy[sublead]);
    l.FillTree("pho2_energy_noregr",pho2_energy_noregr);
    l.FillTree("pho2_pt",sublead_p4.Pt());
    l.FillTree("pho2_eta",sublead_p4.Eta());
    l.FillTree("pho2_phi",sublead_p4.Phi());
    l.FillTree("pho2_r9",sublead_r9);
    l.FillTree("pho2_phoidMva",phoid_mvaout_sublead);
    l.FillTree("pho2_ptOverM",sublead_p4.Pt()/mass);
    l.FillTree("pho2_sceta",pho2_sceta);
    l.FillTree("pho2_scphi",pho2_scphi);
    l.FillTree("pho2_e5x5",pho2_e5x5);
    l.FillTree("pho2_sigmaE",pho2_sigmaE);
    l.FillTree("pho2_scEnergy",pho2_scEnergy);
    l.FillTree("pho2_scEraw",pho2_scEraw);
    l.FillTree("pho2_scEpresh",pho2_scEpresh);
    l.FillTree("pho2_isconv",l.pho_isconv[sublead]);

    l.FillTree("sigmaMOverM",sigmaMrv/mass);
    l.FillTree("sigmaMOverM_wrongVtx",sigmaMwv/mass);
    l.FillTree("vtxProb",vtxProb);
    l.FillTree("cosDeltaPhi",cos_dphi);

    l.FillTree("dipho_mass",mass);
    l.FillTree("dipho_pt",ptHiggs);
    l.FillTree("dipho_mvaout",diphobdt_output);

    l.FillTree("nvtx",l.vtx_std_n);
    l.FillTree("vtx_x",vtx->x());
    l.FillTree("vtx_y",vtx->y());
    l.FillTree("vtx_z",vtx->z());

    l.FillTree("el1_eta",el1_eta);
    l.FillTree("el1_phi",el1_phi);
    l.FillTree("el2_eta",el2_eta);
    l.FillTree("el2_phi",el2_phi);

}

// ----------------------------------------------------------------------------------------------------
void MassFactorizedMvaAnalysis::fillZeeControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
						    const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
						    float phoid_mvaout_lead, float phoid_mvaout_sublead, 
						    float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
						    int diphoton_id, int category, int selectioncategory, float evweight, LoopAll & l )
{

    int lead = l.dipho_leadind[diphoton_id];
    int sublead = l.dipho_subleadind[diphoton_id];

    float mass = Higgs.M();
    float ptHiggs = Higgs.Pt();

    float cos_dphi = TMath::Cos(lead_p4.Phi()-sublead_p4.Phi());
    float pho1_sigmaE = energyCorrectedError[lead];
    float pho2_sigmaE = energyCorrectedError[sublead];

    l.FillHist("mass",0, mass, evweight);
    l.FillHist("mass",category+1, mass, evweight);
    l.FillHist("mass_basecat",selectioncategory, mass, evweight);
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
	l.FillHist("bdtout",selectioncategory+1,diphobdt_output,evweight);
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
	l.FillHist("eta",0, Higgs.Eta(), evweight);
	l.FillHist("eta",category+1, Higgs.Eta(), evweight);
	l.FillHist("rho",0, l.rho, evweight);

	l.FillHist2D("rhoVsNvtx",0, l.vtx_std_n, l.rho, evweight);

	if (fabs(lead_p4.Eta()) < 1.4442) {
	    l.FillHist("pho1_sigmaEOverE_EB",0,pho1_sigmaE/lead_p4.E(), evweight);
	    l.FillHist("pho1_sigmaEOverE_EB_up",0,(pho1_sigmaE*1.1)/lead_p4.E(), evweight);
	    l.FillHist("pho1_sigmaEOverE_EB_down",0,(pho1_sigmaE*0.9)/lead_p4.E(), evweight);
	} else {
	    l.FillHist("pho1_sigmaEOverE_EE",0,pho1_sigmaE/lead_p4.E(), evweight);
	    l.FillHist("pho1_sigmaEOverE_EE_up",0,(pho1_sigmaE*1.1)/lead_p4.E(), evweight);
	    l.FillHist("pho1_sigmaEOverE_EE_down",0,(pho1_sigmaE*0.9)/lead_p4.E(), evweight);
	}

	if (fabs(sublead_p4.Eta()) < 1.4442) {
	    l.FillHist("pho2_sigmaEOverE_EB",0,pho2_sigmaE/sublead_p4.E(), evweight);
	    l.FillHist("pho2_sigmaEOverE_EB_up",0,(pho2_sigmaE*1.1)/sublead_p4.E(), evweight);
	    l.FillHist("pho2_sigmaEOverE_EB_down",0,(pho2_sigmaE*0.9)/sublead_p4.E(), evweight);
	} else {
	    l.FillHist("pho2_sigmaEOverE_EE",0,pho2_sigmaE/sublead_p4.E(), evweight);
	    l.FillHist("pho2_sigmaEOverE_EE_up",0,(pho2_sigmaE*1.1)/sublead_p4.E(), evweight);
	    l.FillHist("pho2_sigmaEOverE_EE_down",0,(pho2_sigmaE*0.9)/sublead_p4.E(), evweight);
	}

	if (fabs(lead_p4.Eta()) < 1.4442) {
	    l.FillHist("pho1_phoidMva_EB",0,phoid_mvaout_lead, evweight);
	    if (l.vtx_std_n>15) {
		l.FillHist("pho1_phoidMva_EB_nvtxgt15",0,phoid_mvaout_lead, evweight);
	    } else {
		l.FillHist("pho1_phoidMva_EB_nvtxlt15",0,phoid_mvaout_lead, evweight);
	    }
	    l.FillHist("pho1_phoidMva_EB_up",0,(phoid_mvaout_lead+0.01), evweight);
	    l.FillHist("pho1_phoidMva_EB_down",0,(phoid_mvaout_lead-0.01), evweight);
	} else {
	    l.FillHist("pho1_phoidMva_EE",0,phoid_mvaout_lead, evweight);
	    if (l.vtx_std_n>15) {
		l.FillHist("pho1_phoidMva_EE_nvtxgt15",0,phoid_mvaout_lead, evweight);
	    } else {
		l.FillHist("pho1_phoidMva_EE_nvtxlt15",0,phoid_mvaout_lead, evweight);
	    }
	    l.FillHist("pho1_phoidMva_EE_up",0,(phoid_mvaout_lead+0.01), evweight);
	    l.FillHist("pho1_phoidMva_EE_down",0,(phoid_mvaout_lead-0.01), evweight);
	}

	if (fabs(sublead_p4.Eta()) < 1.4442) {
	    l.FillHist("pho2_phoidMva_EB",0,phoid_mvaout_sublead, evweight);
	    if (l.vtx_std_n>15) {
		l.FillHist("pho2_phoidMva_EB_nvtxgt15",0,phoid_mvaout_lead, evweight);
	    } else {
		l.FillHist("pho2_phoidMva_EB_nvtxlt15",0,phoid_mvaout_lead, evweight);
	    }
	    l.FillHist("pho2_phoidMva_EB_up",0,(phoid_mvaout_sublead+0.01), evweight);
	    l.FillHist("pho2_phoidMva_EB_down",0,(phoid_mvaout_sublead-0.01), evweight);
	} else {
	    l.FillHist("pho2_phoidMva_EE",0,phoid_mvaout_sublead, evweight);
	    if (l.vtx_std_n>15) {
		l.FillHist("pho2_phoidMva_EE_nvtxgt15",0,phoid_mvaout_lead, evweight);
	    } else {
		l.FillHist("pho2_phoidMva_EE_nvtxlt15",0,phoid_mvaout_lead, evweight);
	    }
	    l.FillHist("pho2_phoidMva_EE_up",0,(phoid_mvaout_sublead+0.01), evweight);
	    l.FillHist("pho2_phoidMva_EE_down",0,(phoid_mvaout_sublead-0.01), evweight);
	}

	if( mass>=83. && mass<=96. ) {
	    l.FillHist("nvtx_83to96",0, l.vtx_std_n, evweight);
	    l.FillHist("nvtx_83to96",category+1, l.vtx_std_n, evweight);
	}

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


	for (int i=0; i<2; i++) {

	    int iPhoton = (i==0) ? lead : sublead;
	    float ptpho = (i==0) ? lead_p4.Pt() : sublead_p4.Pt();
	    float r9 = (i==0) ? lead_r9 : sublead_r9;

	    float pfchargedisobad03=0.;
	    for(int ivtx=0; ivtx<l.vtx_std_n; ivtx++) {
		pfchargedisobad03=(*(l.pho_pfiso_mycharged03))[iPhoton][ivtx]>pfchargedisobad03?(*(l.pho_pfiso_mycharged03))[iPhoton][ivtx]:pfchargedisobad03;
	    }
	    float pfchargedisogood03 = (*(l.pho_pfiso_mycharged03))[iPhoton][l.dipho_vtxind[diphoton_id]];
	    float pfphotoniso03      = l.pho_pfiso_myphoton03[iPhoton];
  
	    float sieie        = l.pho_sieie[iPhoton];
	    float sieip        = l.pho_sieip[iPhoton];
	    float etawidth     = l.pho_etawidth[iPhoton];
	    float phiwidth     = l.sc_sphi[l.pho_scind[iPhoton]];
	    float s4ratio  = l.pho_s4ratio[iPhoton];
	    float ESEffSigmaRR = l.pho_ESEffSigmaRR[iPhoton];

	    if (l.pho_isEB[iPhoton]) {
		l.FillHist("pfchargedisogood03_EB",0,pfchargedisogood03, evweight);
		l.FillHist("pfchargedisobad03_EB",0,pfchargedisogood03, evweight);
		l.FillHist("pfphotoniso03_EB",0,pfphotoniso03, evweight);
		l.FillHist("pfchargedisogood03_tail_EB",0,pfchargedisogood03, evweight);
		l.FillHist("pfchargedisobad03_tail_EB",0,pfchargedisogood03, evweight);
		l.FillHist("pfphotoniso03_tail_EB",0,pfphotoniso03, evweight);
		l.FillHist("pfchargedisogood03_rel_EB",0,pfchargedisogood03/ptpho, evweight);
		l.FillHist("pfchargedisobad03_rel_EB",0,pfchargedisogood03/ptpho, evweight);
		l.FillHist("pfphotoniso03_rel_EB",0,pfphotoniso03/ptpho, evweight);
		l.FillHist("sieie_EB",0,sieie, evweight);
		l.FillHist("sieip_EB",0,sieip, evweight);
		l.FillHist("etawidth_EB",0,etawidth, evweight);
		l.FillHist("phiwidth_EB",0,phiwidth, evweight);
		l.FillHist("s4ratio_EB",0,s4ratio, evweight);
		l.FillHist("r9_EB",0,lead_r9, evweight);
		l.FillHist("rho_EB",0,l.rho, evweight);
	    } else {
		l.FillHist("pfchargedisogood03_EE",0,pfchargedisogood03, evweight);
		l.FillHist("pfchargedisobad03_EE",0,pfchargedisogood03, evweight);
		l.FillHist("pfphotoniso03_EE",0,pfphotoniso03, evweight);
		l.FillHist("pfchargedisogood03_tail_EE",0,pfchargedisogood03, evweight);
		l.FillHist("pfchargedisobad03_tail_EE",0,pfchargedisogood03, evweight);
		l.FillHist("pfphotoniso03_tail_EE",0,pfphotoniso03, evweight);
		l.FillHist("pfchargedisogood03_rel_EE",0,pfchargedisogood03/ptpho, evweight);
		l.FillHist("pfchargedisobad03_rel_EE",0,pfchargedisogood03/ptpho, evweight);
		l.FillHist("pfphotoniso03_rel_EE",0,pfphotoniso03/ptpho, evweight);
		l.FillHist("sieie_EE",0,sieie, evweight);
		l.FillHist("sieip_EE",0,sieip, evweight);
		l.FillHist("etawidth_EE",0,etawidth, evweight);
		l.FillHist("phiwidth_EE",0,phiwidth, evweight);
		l.FillHist("s4ratio_EE",0,s4ratio, evweight);
		l.FillHist("ESEffSigmaRR_EE",0,ESEffSigmaRR, evweight);
		l.FillHist("r9_EE",0,r9, evweight);
		l.FillHist("rho_EE",0,l.rho, evweight);
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
	    int cat = categoryFromBoundaries( bdtCategoryBoundaries, bdtout );
	    if( VBFevent && cat > -1 ) cat = bdtCategoryBoundaries.size();
	    return cat;
    } else std::cerr << "No BDT Philosophy known - " << bdtTrainingPhilosophy << std::endl;
}

void MassFactorizedMvaAnalysis::fillTrainTree(LoopAll &l, Int_t leadingPho, Int_t subleadingPho, Int_t vtx, float vtxProb, TLorentzVector &leadP4, TLorentzVector &subleadP4, float sigmaMrv, float sigmaMwv, float sigmaMeonly, const char* type, float photonID_1,float photonID_2){
    
    Float_t mva = 99.;
    TLorentzVector Higgs = leadP4+subleadP4;
    float leadPt    = leadP4.Pt();
    float subleadPt = subleadP4.Pt();
    float mass     = Higgs.M();
    float diphopt   = Higgs.Pt();

    l.FillTree("dmom",sigmaMrv/mass);
    l.FillTree("dmom_wrong_vtx",sigmaMwv/mass);
    l.FillTree("vtxprob",vtxProb);
    l.FillTree("ptom1",leadPt/mass);
    l.FillTree("ptom2",subleadPt/mass);
    l.FillTree("eta1",leadP4.Eta());
    l.FillTree("eta2",subleadP4.Eta());
    l.FillTree("dphi",TMath::Cos(leadP4.Phi()-subleadP4.Phi()));
    l.FillTree("ph1mva",photonID_1);
    l.FillTree("ph2mva",photonID_2);
}

void MassFactorizedMvaAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    if ( doEresolSmear ) eResolSmearer->resetRandom();
}

	
// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
