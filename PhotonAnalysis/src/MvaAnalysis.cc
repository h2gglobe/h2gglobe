#include "../interface/MvaAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::MvaAnalysis()  : 
    name_("MvaAnalysis"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    nMasses  = 9;
    doTraining = false;
    trainingMassHypothesis = 123.0;
    
}

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::~MvaAnalysis() 
{
    //Default
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Term(LoopAll& l) 
{
    // All Final Fits + rebinning are now done in Macros/FullMvaToolkit rather than in Term
    eventListText.close();
    PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
        cout << "InitRealMvaAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;

    eventListText.open(Form("%s",l.outputTextFileName.c_str()));

    if (doTraining) {
        nMasses  = 2;
        names[0]="_121.0";
        BDTnames[0]="_121";
        masses[0] = 121.;

        names[1]="_123.0";
        BDTnames[1]="_123";
        masses[1] = 123.;
    }
    else {
        names[0]="_105.0";
        BDTnames[0]="_105";
        masses[0] = 105.;

        names[1]="_110.0";
        BDTnames[1]="_110";
        masses[1] = 110.;

        names[2]="_115.0";
        BDTnames[2]="_115";
        masses[2] = 115.;

        names[3]="_120.0";
        BDTnames[3]="_120";
        masses[3] = 120.;

        names[4]="_125.0";
        BDTnames[4]="_125";
        masses[4] = 125.;

        names[5]="_130.0";
        BDTnames[5]="_130";
        masses[5] = 130.;

        names[6]="_135.0";
        BDTnames[6]="_135";
        masses[6] = 135.;

        names[7]="_140.0";
        BDTnames[7]="_140";
        masses[7] = 140.;

        names[8]="_150.0";
        BDTnames[8]="_150";
        masses[8] = 150.;
        /*
          names[9]="_121.0";
          BDTnames[9]="_121";
          masses[9] = 121.;

          names[10]="_123.0";
          BDTnames[10]="_123";
          masses[10] = 123.;
        */
    }

    // Make sure that we wont try to use more sidebands than available
    assert(numberOfSidebandsForAlgos+numberOfSidebandGaps <= numberOfSidebands);

    std::string outputfilename = (std::string) l.histFileName;
    //
    // These parameters are set in the configuration file
    std::cout
        << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << "MvaAnalysis " << "\n"
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

    cout << "using bdt philosophy:       " << bdtTrainingPhilosophy << endl;
    if (doTraining) cout << "Training ON:    saving trees" << endl;
    else cout << "Training OFF:     running analysis" << endl;
    if (!doTraining) cout << "Reading MVA of type:       " << MVAtype << endl; 
    // call the base class initializer
    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
        l.histoContainer[ind].setScale(1.);
    }

    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

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
        eResolSmearer->doEnergy(false);
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

    /*    if(doKFactorSmear2D) {
    // kFactor efficiency
    std::cerr << __LINE__ << std::endl; 
    kFactorSmearer2D = new KFactorSmearer2D( kfacHist );
    kFactorSmearer2D->name("kFactor2D");
    kFactorSmearer2D->init();
    genLevelSmearers_.push_back(kFactorSmearer2D);
    }
    */

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed

    // FIXME move these params to config file
    l.rooContainer->SetNCategories(1+(int)includeVBF);

    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;
    l.rooContainer->SaveRooDataHists(false);
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
    /*    if(doKFactorSmear2D && doKFactorSyst2D) {
          systGenLevelSmearers_.push_back(kFactorSmearer2D);
          std::vector<std::string> sys(1,kFactorSmearer2D->name());
          std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
          l.rooContainer->MakeSystematicStudy(sys,sys_t);
          }
    */

    // ----------------------------------------------------
    // ----------------------------------------------------
    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.045,1.00);
    // ----------------------------------------------------
    l.rooContainer->AddObservable("CMS_hgg_mass",massMin,massMax);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);


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

    if (doTraining){ // do nothing since Trees are setup inside TreeContainers
    }

    else{

        // override defaults if Zee Validation
        if (l.runZeeValidation){
            mHMaximum=90.0;
            mHMinimum=90.0;
        }

        l.rooContainer->AddObservable("BDT" ,-1.,1.);

        //Set up TMVA reader (only two variables)
        tmvaReader_= new TMVA::Reader();

        tmvaReader_->AddVariable("bdtoutput",&_bdtoutput);
        tmvaReader_->AddVariable("deltaMOverM", &_deltaMOverM);

        //Invariant Mass Spectra
        l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass",nDataBins);
        l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass" ,nDataBins);       

        int nBDTbins = 5000;    // The initial number of bins (ie before any binning algorithm is performed)


        // Usual datasets
        for (double mass=mHMinimum; mass<=mHMaximum; mass+=mHStep){

            //Gradient Boost
            l.rooContainer->CreateDataSet("BDT",Form("data_BDT_grad_%3.1f",mass)     ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT",Form("bkg_BDT_grad_%3.1f",mass)      ,nBDTbins);

            for (int sideband_i=1;sideband_i<=numberOfSidebands;sideband_i++){
                // Always create all of the sidebands, even if we skip some of them
                // later on for the sums

                l.rooContainer->CreateDataSet("BDT",Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass)  ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass) ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass) ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass),nBDTbins);

            }

        }

        // loop signal mass points signal datasets
        for (double sig=mHMinimum;sig<=mHMaximum;sig+=5.0){  // We are ignoring mass 105 for now
            if (sig==145.0) continue;
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_ggh_%3.1f",sig)      ,nBDTbins); 
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_vbf_%3.1f",sig)      ,nBDTbins); 
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_wzh_%3.1f",sig)      ,nBDTbins); 
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_tth_%3.1f",sig)      ,nBDTbins); 


            // Make the signal Systematic Sets
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_ggh_%3.1f",sig),-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_vbf_%3.1f",sig),-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_wzh_%3.1f",sig),-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_tth_%3.1f",sig),-1);


        }
        //TMVA Reader
        tmvaReader_->BookMVA("BDT_grad_123", mvaWeightsFolder+"/TMVAClassification_BDTgradMIT.weights.xml");
    }

    FillSignalLabelMap();

    if(PADEBUG) 
        cout << "InitRealMvaAnalysis END"<<endl;
}


// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------

int MvaAnalysis::SignalType(int cur_type){
    int i0 = -1;
    if (doTraining) {
        if (cur_type == -53 || cur_type == -54 ||  cur_type == -55 || cur_type == -56){//121 
            i0 = 0;}
        else if (cur_type == -57 || cur_type == -58 ||  cur_type == -59 || cur_type == -60){//123 
            i0 = 1;}
    } else{
        if (cur_type == -13 || cur_type == -14 ||  cur_type == -15 || cur_type == -16){//105 
            i0 = 0;}
        else if (cur_type == -17 || cur_type == -18 ||  cur_type == -19 || cur_type == -20){//110 
            i0 = 1;}
        else if (cur_type == -21 || cur_type == -22 ||  cur_type == -23 || cur_type == -24){//115 
            i0 = 2;}
        else if (cur_type == -25 || cur_type == -26 ||  cur_type == -27 || cur_type == -28){//120 
            i0 = 3;}
        else if (cur_type == -37 || cur_type == -38 ||  cur_type == -39 || cur_type == -40){//125 
            i0 = 4;}
        else if (cur_type == -29 || cur_type == -30 ||  cur_type == -31 || cur_type == -32){//130 
            i0 = 5;}
        else if (cur_type == -41 || cur_type == -42 ||  cur_type == -43 || cur_type == -44){//135 
            i0 = 6;}
        else if (cur_type == -33 || cur_type == -34 ||  cur_type == -35 || cur_type == -36){//140 
            i0 = 7;}
        //    else if (cur_type == -45 || cur_type == -46 ||  cur_type == -47 || cur_type == -48){//145 
        //        i0 = 8;}
        else if (cur_type == -49 || cur_type == -50 ||  cur_type == -51 || cur_type == -52){//150 
            i0 = 8;}
    }
    return i0;
}
// ----------------------------------------------------------------------------------------------------

void MvaAnalysis::FillSignalLabelMap(){

    // Basically A Map of the ID (type) to the signal's name which can be filled Now:
    signalLabels[-57]="ggh_123.0";
    signalLabels[-58]="vbf_123.0";
    signalLabels[-60]="wzh_123.0";
    signalLabels[-59]="tth_123.0";
    signalLabels[-53]="ggh_121.0";
    signalLabels[-54]="vbf_121.0";
    signalLabels[-56]="wzh_121.0";
    signalLabels[-55]="tth_121.0";
    signalLabels[-65]="ggh_160.0";
    signalLabels[-66]="vbf_160.0";
    signalLabels[-68]="wzh_160.0";
    signalLabels[-67]="tth_160.0";
    signalLabels[-61]="ggh_155.0";
    signalLabels[-62]="vbf_155.0";
    signalLabels[-64]="wzh_155.0";
    signalLabels[-63]="tth_155.0";
    signalLabels[-49]="ggh_150.0";
    signalLabels[-50]="vbf_150.0";
    signalLabels[-52]="wzh_150.0";
    signalLabels[-51]="tth_150.0";
    signalLabels[-45]="ggh_145.0";
    signalLabels[-46]="vbf_145.0";
    signalLabels[-48]="wzh_145.0";
    signalLabels[-47]="tth_145.0";
    signalLabels[-33]="ggh_140.0";
    signalLabels[-34]="vbf_140.0";
    signalLabels[-36]="wzh_140.0";
    signalLabels[-35]="tth_140.0";
    signalLabels[-41]="ggh_135.0";
    signalLabels[-42]="vbf_135.0";
    signalLabels[-44]="wzh_135.0";
    signalLabels[-43]="tth_135.0";
    signalLabels[-29]="ggh_130.0";
    signalLabels[-30]="vbf_130.0";
    signalLabels[-32]="wzh_130.0";
    signalLabels[-31]="tth_130.0";
    signalLabels[-37]="ggh_125.0";
    signalLabels[-38]="vbf_125.0";
    signalLabels[-40]="wzh_125.0";
    signalLabels[-39]="tth_125.0";
    signalLabels[-25]="ggh_120.0";
    signalLabels[-26]="vbf_120.0";
    signalLabels[-28]="wzh_120.0";
    signalLabels[-27]="tth_120.0";
    signalLabels[-21]="ggh_115.0";
    signalLabels[-22]="vbf_115.0";
    signalLabels[-24]="wzh_115.0";
    signalLabels[-23]="tth_115.0";
    signalLabels[-17]="ggh_110.0";
    signalLabels[-18]="vbf_110.0";
    signalLabels[-19]="wzh_110.0";
    signalLabels[-20]="tth_110.0";
    signalLabels[-13]="ggh_105.0";
    signalLabels[-14]="vbf_105.0";
    signalLabels[-16]="wzh_105.0";
    signalLabels[-15]="tth_105.0";
    signalLabels[-69]="ggh_100.0";
    signalLabels[-70]="vbf_100.0";
    signalLabels[-72]="wzh_100.0";
    signalLabels[-71]="tth_100.0";
}
// ----------------------------------------------------------------------------------------------------

std::string MvaAnalysis::GetSignalLabel(int id){

    // For the lazy man, can return a memeber of the map rather than doing it yourself
    std::map<int,std::string>::iterator it = signalLabels.find(id);

    if (it!=signalLabels.end()){
        return it->second;

    } else { 

        std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
        return "NULL";
    }

}
// ----------------------------------------------------------------------------------------------------

void MvaAnalysis::SetBDTInputVariables(TLorentzVector *lead_p4, TLorentzVector *sublead_p4, double lead_r9, double sublead_r9, MassResolution *massResolutionCalculator, double vtx_mva, double mass_hypothesis, double bdtoutput, double evweight, int cat){

    TLorentzVector Higgs = *lead_p4 + *sublead_p4;   
    float mass    = Higgs.M();

    _log_H_pt =  log10( Higgs.Pt());
    _H_eta = fabs(Higgs.Eta());
    _d_phi = fabs(lead_p4->DeltaPhi(*sublead_p4));
    _cos_d_phi = TMath::Cos(lead_p4->Phi()-sublead_p4->Phi());        
    _max_eta = max(fabs(lead_p4->Eta()),fabs(sublead_p4->Eta()));
    _min_r9  = min(lead_r9,sublead_r9);
    _pho1_eta = lead_p4->Eta();
    _pho2_eta = sublead_p4->Eta();

    _mgg = mass;
    _pho1_phi = lead_p4->Phi();
    _pho1_pt = lead_p4->Pt();
    _pho1_r9 = lead_r9;

    _pho2_phi = sublead_p4->Phi();
    _pho2_pt = sublead_p4->Pt();
    _pho2_r9 = sublead_r9;

    _H_pt = Higgs.Pt();
    _Ht = lead_p4->Pt()+sublead_p4->Pt();

    _d_eta = lead_p4->Eta()-sublead_p4->Eta();
    _mod_d_eta = fabs(lead_p4->Eta()-sublead_p4->Eta());
    _cos_theta_star = fabs(lead_p4->E()-sublead_p4->E())/Higgs.P();


    _vtx_prob = 1.-0.49*(vtx_mva+1.0);

    _wt= evweight;

    _pho1_ptOverM = lead_p4->Pt()/mass_hypothesis;
    _pho2_ptOverM = sublead_p4->Pt()/mass_hypothesis;
    _deltaMOverM = (mass-mass_hypothesis)/mass_hypothesis;

    double massResolution = massResolutionCalculator->massResolutionCorrVtx();
    _deltaMOverSigmaM = (mass-mass_hypothesis)/massResolution;
    _sigmaMOverM = massResolution/mass;
    _sigmaMOverM_wrongVtx = massResolutionCalculator->massResolutionWrongVtx()/mass;
    _H_ptOverM    = Higgs.Pt()/mass_hypothesis;
    _cat        = cat;
    _bdtoutput = bdtoutput;
}
// ----------------------------------------------------------------------------------------------------

void MvaAnalysis::SetBDTInputTree(TTree *tree){
    mvaFile_->cd();
    TBranch *b_log_H_pt             = tree->Branch("log_H_pt", &_log_H_pt , "log_H_pt/F");;
    TBranch *b_H_ptOverM            = tree->Branch("H_ptOverM", &_H_ptOverM , "_H_ptOverM/F");;
    TBranch *b_H_eta                = tree->Branch("H_eta", &_H_eta , "H_eta/F");;
    TBranch *b_d_phi                = tree->Branch("d_phi", &_d_phi,"d_phi/F");
    TBranch *b_max_eta              = tree->Branch("max_eta", &_max_eta , "max_eta/F");;
    TBranch *b_min_r9               = tree->Branch("min_r9", &_min_r9 , "min_r9/F");;
    TBranch *b_pho1_eta             = tree->Branch("pho1_eta", &_pho1_eta , "pho1_eta/F");
    TBranch *b_pho2_eta             = tree->Branch("pho2_eta", &_pho2_eta , "pho2_eta/F");
    TBranch *b_pho1_ptOverM         = tree->Branch("pho1_ptOverM", &_pho1_ptOverM , "pho1_ptOverM/F");;
    TBranch *b_pho2_ptOverM         = tree->Branch("pho2_ptOverM", &_pho2_ptOverM , "pho2_ptOverM/F");;
    TBranch *b_deltaMOverM          = tree->Branch("deltaMOverM", &_deltaMOverM,"deltaMOverM/F");
    TBranch *b_deltaMOverSigmaM     = tree->Branch("deltaMOverSigmaM", &_deltaMOverSigmaM,"deltaMOverSigmaM/F");
    TBranch *b_sigmaMOverM          = tree->Branch("sigmaMOverM", &_sigmaMOverM,"sigmaMOverM/F");
    TBranch *b_sigmaMOverM_wrongVtx = tree->Branch("sigmaMOverM_wrongVtx", &_sigmaMOverM_wrongVtx,"sigmaMOverM_wrongVtx/F");
    TBranch *b_mgg                  = tree->Branch("mgg", &_mgg, "mgg/F");
    TBranch *b_pho1_phi             = tree->Branch("pho1_phi", &_pho1_phi , "pho1_phi/F");
    TBranch *b_pho1_pt              = tree->Branch("pho1_pt", &_pho1_pt , "pho1_pt/F");;
    TBranch *b_pho1_r9              = tree->Branch("pho1_r9", &_pho1_r9 , "_pho1_r9/F");;
    TBranch *b_pho2_phi             = tree->Branch("pho2_phi", &_pho2_phi , "pho2_phi/F");
    TBranch *b_pho2_pt              = tree->Branch("pho2_pt", &_pho2_pt , "pho2_pt/F");;
    TBranch *b_pho2_r9              = tree->Branch("pho2_r9", &_pho2_r9 , "_pho2_r9/F");;
    TBranch *b_H_pt                 = tree->Branch("H_pt", &_H_pt , "H_pt/F");;
    TBranch *b_Ht                   = tree->Branch("Ht", &_Ht , "Ht/F");;
    TBranch *b_d_eta                = tree->Branch("d_eta", &_d_eta,"d_eta/F");
    TBranch *b_mod_d_eta            = tree->Branch("mod_d_eta", &_mod_d_eta,"mod_d_eta/F");
    TBranch *b_cos_theta_star       = tree->Branch("cos_theta_star", &_cos_theta_star , "cos_theta_star/F");;
    TBranch *b_vtx_prob             = tree->Branch("vtx_prob", &_vtx_prob, "vtx_prob/F");
    TBranch *b_wt                   = tree->Branch("wt", &_wt, "wt/F");
    TBranch *b_cat                  = tree->Branch("cat", &_cat, "cat/I");
    TBranch *b_sideband             = tree->Branch("sideband", &_sideband, "sideband/I");
    TBranch *b_bdtoutput            = tree->Branch("bdtoutput", &_bdtoutput, "bdtoutput/F");
}
// ----------------------------------------------------------------------------------------------------


float MvaAnalysis::tmvaGetVal(double mass, double mass_hypothesis, float kinematic_bdt){

    _deltaMOverM = (mass-mass_hypothesis)/mass_hypothesis;
    _bdtoutput = kinematic_bdt;
    return tmvaReader_->EvaluateMVA( "BDT_grad_123" );

}
// ----------------------------------------------------------------------------------------------------
int MvaAnalysis::GetBDTBoundaryCategory(float bdtout, bool isEB, bool VBFevent){

    if (VBFevent) {
        if (bdtout >= 0.05) return 1;
        return -1;
    } else {
        if (bdtout >= 0.05) return 0;
        return -1;
    }
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::fillTMVATrees(LoopAll& l,float mass,float diphotonMVA,int category,float evweight, int cur_type){

        // The event is exclusively in one of these bands for a given mH
        float mH = trainingMassHypothesis;

        // Fill Signal Window
        float sideband_boundaries[2];
        sideband_boundaries[0] = mH*(1-sidebandWidth);
        sideband_boundaries[1] = mH*(1+sidebandWidth);
        if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut

            // Signal Region
            l.FillHist("diphotonMVA",diphotonMVA);
            l.FillHist("mass",mass);
            l.FillHist("category",category);
            l.FillHist("weight",evweight);
            l.FillHist("deltaMoverM",(float)((mass-mH)/mH));
            return;
        }

        if (cur_type>0){ // Dont use signal in sidebands

            // Loop over N lower sidebands
            for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                double hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth);
                double mass_hypothesis_low     = (mH*(1.-signalRegionWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                double sideband_boundaries_low = mass_hypothesis_low*(1.-sidebandWidth);
                double sideband_boundaries_high= mass_hypothesis_low*(1.+sidebandWidth);

                if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                    l.FillHist("diphotonMVA",diphotonMVA);
                    l.FillHist("mass",mass);
                    l.FillHist("category",category);
                    l.FillHist("weight",evweight);
                    l.FillHist("deltaMoverM",(float)((mass-mass_hypothesis_low)/mass_hypothesis_low));
                    return;
                }
            }
            // Loop over N higher sidebands
            for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                double hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth);
                double mass_hypothesis_high     = (mH*(1.+signalRegionWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                double sideband_boundaries_low = mass_hypothesis_high*(1.-sidebandWidth);
                double sideband_boundaries_high= mass_hypothesis_high*(1.+sidebandWidth);

                if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                    l.FillHist("diphotonMVA",diphotonMVA);
                    l.FillHist("mass",mass);
                    l.FillHist("category",category);
                    l.FillHist("weight",evweight);
                    l.FillHist("deltaMoverM",(float)((mass-mass_hypothesis_high)/mass_hypothesis_high));
                    return;
                }
            }
        }

}
// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA, 
				    int category, float weight, bool isCorrectVertex) 
{


  if (doTraining){
    fillTMVATrees(l,mass,diphotonMVA,category,weight,cur_type);
  } else {

    // --- Fill invariant mass spectrum -------
    if (cur_type==0){  // Data
                l.rooContainer->InputDataPoint("data_mass",category,mass);
    } else if (cur_type>0){ // Background MC
                l.rooContainer->InputBinnedDataPoint("bkg_mass",category,mass,weight);
    }


    if (cur_type<0){ // signal MC

            std::string currentTypeSignalLabel = GetSignalLabel(cur_type);
            float mass_hypothesis = masses[SignalType(cur_type)];
            float sideband_boundaries[2];
            sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
            sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);
            if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                   float bdt_grad =  tmvaGetVal(mass,mass_hypothesis,diphotonMVA);
                   l.rooContainer->InputBinnedDataPoint("sig_BDT_grad_"+currentTypeSignalLabel ,category,bdt_grad,weight);
            
            }

     } else {

        for (double mH=mHMinimum; mH<=mHMaximum; mH+=mHStep){
            // Fill Signal Window
            float sideband_boundaries[2];
            sideband_boundaries[0] = mH*(1-sidebandWidth);
            sideband_boundaries[1] = mH*(1+sidebandWidth);
            if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                float bdt_grad =  tmvaGetVal(mass,mH,diphotonMVA);
                if (cur_type==0) l.rooContainer->InputBinnedDataPoint(Form("data_BDT_grad_%3.1f",mH),category,bdt_grad,weight);
                else if (cur_type > 0) l.rooContainer->InputBinnedDataPoint(Form("bkg_BDT_grad_%3.1f",mH) ,category,bdt_grad,weight);
            }
            
            // Loop over N lower sidebands
            for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                double hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth);
                double mass_hypothesis_low     = (mH*(1.-signalRegionWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                double sideband_boundaries_low = mass_hypothesis_low*(1.-sidebandWidth);
                double sideband_boundaries_high= mass_hypothesis_low*(1.+sidebandWidth);

                if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                   float bdt_grad =  tmvaGetVal(mass,mass_hypothesis_low,diphotonMVA);
                   if (cur_type==0) l.rooContainer->InputBinnedDataPoint(Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mH) ,category,bdt_grad,weight);
                   else if (cur_type > 0) l.rooContainer->InputBinnedDataPoint(Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mH) ,category,bdt_grad,weight);
                   break;

                }
            }
            // Loop over N higher sidebands
            for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                double hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth);
                double mass_hypothesis_high     = (mH*(1.+signalRegionWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                double sideband_boundaries_low = mass_hypothesis_high*(1.-sidebandWidth);
                double sideband_boundaries_high= mass_hypothesis_high*(1.+sidebandWidth);

                if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                   float bdt_grad =  tmvaGetVal(mass,mass_hypothesis_high,diphotonMVA);
                   if (cur_type==0) l.rooContainer->InputBinnedDataPoint(Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mH) ,category,bdt_grad,weight);
                   else if (cur_type > 0) l.rooContainer->InputBinnedDataPoint(Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mH) ,category,bdt_grad,weight);
                   break;
                }
            }

        }
    }
  }
}
// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::AccumulateSyst(int cur_type, float mass, float diphotonMVA, 
				  int category, float weight,
				  std::vector<double> & mass_errors,
				  std::vector<double> & mva_errors,
				  std::vector<int>    & categories,
				  std::vector<double> & weights)
{
    float mass_hypothesis = masses[SignalType(cur_type)];
    // define the sidebands
    float sideband_boundaries[2];
    sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
    sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);
	// Signal Window cut
	if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){
			double syst_bdt_grad = tmvaGetVal(mass,mass_hypothesis,diphotonMVA);
			categories.push_back(category);
			mva_errors.push_back(syst_bdt_grad);
			weights.push_back(weight);
	} else {
			categories.push_back(-1);
			mva_errors.push_back(-100);
			weights.push_back(0);
	}
}
// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::FillRooContainerSyst(LoopAll& l, const std::string &name, int cur_type,
			  std::vector<double> & mass_errors, std::vector<double> & mva_errors,
			  std::vector<int>    & categories, std::vector<double> & weights) 
{
        // Get Signal Label for the current type
        if (cur_type<0){
            std::string currentTypeSignalLabel = GetSignalLabel(cur_type);
		    l.rooContainer->InputSystematicSet("sig_BDT_grad_"+currentTypeSignalLabel,name,categories,mva_errors,weights);
        }
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
