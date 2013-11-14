#ifndef __PHOTONANALYSIS__
#define __PHOTONANALYSIS__

#include "BaseAnalysis.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/PhotonInfo.h"

#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include "TriggerSelection.h"
#include "EnergySmearer.h"

#include "MassResolution.h"

#include "TMVA/Reader.h"
#include "PhotonFix.h"
#include <stdio.h>
// #include "HiggsToGammaGamma/interface/GBRForest.h"
//#include "../../../../HiggsToGammaGamma/interface/GBRForest.h"
//#include "HiggsAnalysis/HiggsToGammaGamma/interface/GBRForest.h"

#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
//#include "HiggsAnalysis/GBRLikelihoodEGTools/interface/EGEnergyCorrectorSemiParm.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"
#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForest.h"
#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForestD.h"

class JetHandler;

// ------------------------------------------------------------------------------------
class PhotonAnalysis : public BaseAnalysis
{
    friend class LoopAll;
 public:

    PhotonAnalysis();
    virtual ~PhotonAnalysis();

    virtual const std::string & name() const { return name_; };

    // LoopAll analysis interface implementation
    virtual void Init(LoopAll&);
    virtual void Term(LoopAll&);

    virtual void ReducedOutputTree(LoopAll &l, TTree *);
    virtual void GetBranches(TTree *, std::set<TBranch *>& );

    virtual void FillReductionVariables(LoopAll& l, int jentry);
    virtual bool SelectEventsReduction(LoopAll&, int);

    virtual bool SkimEvents(LoopAll&, int);
    virtual bool SelectEvents(LoopAll&, int);
    virtual bool Analysis(LoopAll&, Int_t);

    virtual void ResetAnalysis();

    float zero_;
    void GetRegressionCorrections(LoopAll&);
    void GetRegressionCorrectionsV5(LoopAll&); // 8 TeV
    void GetRegressionCorrectionsV8(LoopAll&); // 7 TeV (V6-Barrel / V7-Endcap)
    void GetSinglePhotonRegressionCorrectionV6(LoopAll&,int,double *,double *);
    void GetSinglePhotonRegressionCorrectionV7(LoopAll&,int,double *,double *);

    //  void GetRegressionCorrections(LoopAll&);
    // Public parameters to be read from config file
    VertexAlgoParameters vtxAlgoParams;
    std::vector<std::string> vtxVarNames;
    std::vector<string> tmvaPerVtxVariables;

    bool reRunVtx, rematchConversions;

    std::string tmvaPerVtxMethod;
    std::string tmvaPerVtxWeights;
    std::string tmvaPerEvtMethod;
    std::string tmvaPerEvtWeights;

    bool useDefaultVertex;
    float forcedRho;
    bool applyPtoverM;
    float massMin, massMax;
    bool doTriggerSelection;
    bool useRunDTriggersForZee;
    std::vector<TriggerSelection> triggerSelections;

    // Options
    float phoidMvaCut;
    bool run7TeV4Xanalysis;
    bool includeVBF;
    bool includeVHhad;
    bool includeVHhadBtag;
    bool includeTTHlep;
    bool includeTTHhad;
    bool includeVHlep;
    bool includeVHlepPlusMet;
    int nElectronCategories;
    int nMuonCategories;
    bool includeVHmet;  //met at analysis step
    int nMetCategories;
    bool moriond2013MetCorrection;

    bool reRunCiCForData;
    bool reComputeCiCPF;
    bool skimOnDiphoN;

    float leadEtCut;
    float leadEtVBFCut;
    float leadEtVHhadCut;
    float leadEtVHhadBtagCut;
    float leadEtTTHlepCut;
    float leadEtTTHhadCut;
    float leadEtVHlepCut;
    float leadEtVHmetCut;  //met at analysis step
    float subleadEtCut;
    float subleadEtVBFCut;
    float subleadEtVHhadCut;
    float subleadEtVHhadBtagCut;
    float subleadEtVHlepCut;
    float subleadEtVHmetCut;  //met at analysis step
    float subleadEtTTHhadCut;
    float subleadEtTTHlepCut;
    int nVBFEtaCategories;
    int nVHhadEtaCategories;
    int nVBFDijetJetCategories;

    bool emulateBeamspot;
    double emulatedBeamspotWidth;
    double beamspotWidth;
    double beamspotSigma; // for massResolution
    double targetsigma;
    double sourcesigma;
    bool reweighBeamspot;
    bool saveBSTrees_;
    bool rescaleDZforVtxMVA;
    bool pickRandomVtx;
    bool randomizeHiggsPt;
    
    bool reweighPt;
    string ptreweightype;
    string ptreweighfilename;
    TFile *ptreweighfile;
    TH1F *ptreweighHistSM;
    TH1F *ptreweighHistGG;
    TH1F *ptreweighHistQQ;

    bool saveDatacardTrees_;
    double datacardTreeMass;
    bool saveSpinTrees_;
    bool runJetsForSpin;
    bool saveVBFTrees_;

    //vhhadronic cuts                                                                                                                                                           
    float ptgg_btag_cut,costhetastar_btag_cut,ptjet_btag_cut;
    float ptgg_0tag_cut,costhetastar_0tag_cut,ptjet_0tag_cut;
    float ptjet_loosecut;

    float deltaRPholep_cut;
    bool doMinvCut;

    float ptJets_ttH_thresh;
    int njets_tthHad_thresh;
    bool doDrGsfTrackCut;

    float drSC_lep;
    float    drGsf_lep;
    int isLep_ele,isLep_mu;
    int eleIndex,muIndex;

    bool doApplyEleVeto;
    bool removeBtagtth;

    bool createCS;
    bool doLooseLep;

    float diphobdt_output_Cut_TTHlep;
    float diphobdt_output_Cut_TTHhad;
    float diphobdt_output_Cut_VHhadBtag;
    float diphobdt_output_Cut_VHhad;
    float diphobdt_output_Cut_VHLep;
    float diphobdt_output_Cut_VHMet;

    bool optimizeMVA;

    // Preselection indexes
    float presel_scet1, presel_scet2, presel_maxeta;
    float presel_ecaliso_eb, presel_ecaliso_ee, presel_sieie_eb, presel_sieie_ee, presel_hoe;

    bool doEcorrectionSmear, doEcorrectionSyst;

    EnergySmearer::energySmearingParameters eSmearDataPars;
    std::string scale_offset_file;
    float scale_offset_EBHighR9         ;
    float scale_offset_EBLowR9          ;
    float scale_offset_EBm4HighR9         ;
    float scale_offset_EBm4LowR9          ;
    float scale_offset_EEHighR9         ;
    float scale_offset_EELowR9          ;
    float scale_offset_error_EBHighR9   ;
    float scale_offset_error_EBLowR9    ;
    float scale_offset_error_EBm4HighR9   ;
    float scale_offset_error_EBm4LowR9    ;
    float scale_offset_error_EEHighR9   ;
    float scale_offset_error_EELowR9    ;

    bool splitEscaleSyst, doStocasticSmearingSyst;
    std::string scale_offset_error_file, scale_offset_corr_error_file;
    EnergySmearer::energySmearingParameters eScaleCorrPars;
    void setupEscaleSmearer();
    void setupEscaleSyst(LoopAll &l);

    bool splitEresolSyst;
    std::string corr_smearing_file, smearing_file;
    EnergySmearer::energySmearingParameters eResolCorrPars;
    void setupEresolSmearer();
    void setupEresolSyst(LoopAll &l);

    EnergySmearer::energySmearingParameters eSmearPars;
    float smearing_sigma_EBHighR9       ;
    float smearing_sigma_EBLowR9        ;
    float smearing_sigma_EBm4HighR9       ;
    float smearing_sigma_EBm4LowR9        ;
    float smearing_sigma_EEHighR9       ;
    float smearing_sigma_EELowR9        ;
    float smearing_sigma_error_EBHighR9 ;
    float smearing_sigma_error_EBLowR9  ;
    float smearing_sigma_error_EBm4HighR9 ;
    float smearing_sigma_error_EBm4LowR9  ;
    float smearing_sigma_error_EEHighR9 ;
    float smearing_sigma_error_EELowR9  ;
    std::string mass_resol_file;
    EnergySmearer::energySmearingParameters massResoPars;

    std::vector<int> pho_acc;
    std::vector<int> pho_presel;
    std::vector<int> pho_presel_lead;
    std::vector<float> pho_et;
    // Other options
    bool runStatAnalysis;
    TString puHist, puMap, puTarget;//name of pileup reweighting histogram
    std::vector<TString> puTargets; 

    enum BkgCategory{promptprompt,promptfake,fakefake};
    bool keepPP, keepPF, keepFF;
    bool selectprocess;
    int processtoselect;

    std::string energyCorrectionMethod;
    //std::string massResolutionFileName;

    bool mvaVertexSelection, addConversionToMva;

    // PhotonFix
    std::string photonFixDat;
    std::string regressionFile;
    int regressionVersion;

    int   nEtaCategories, nR9Categories, nPtCategories, nPtOverMCategories, nVtxCategories;
    float R9CatBoundary;
    bool  usePUjetveto;
    //std::string photonFixDat;
    //std::string regressionFile;

    // n-1 plot for ClassicCats
    float sublead_r9;
    float sublead_isoOverEt;
    float sublead_badisoOverEt;
    float sublead_trkisooet;
    float sublead_sieie;
    float sublead_drtotk;
    float sublead_hovere;
    float sublead_mgg;

    // n-1 plots for VH electron tag
    float myEl_leptonSig ;
    float myEl_elpt      ;
    float myEl_oEsuboP   ;
    float myEl_D0        ;
    float myEl_DZ        ;
    float myEl_mishit    ;
    float myEl_conv      ;
    float myEl_detain    ;
    float myEl_dphiin    ;
    float myEl_sieie     ;
    float myEl_hoe       ;
    float myEl_drlead    ;
    float myEl_drsub     ;
    float myEl_melead    ;
    float myEl_meleadveto10;
    float myEl_meleadveto15;
    float myEl_mesub     ;
    float myEl_mesubveto5 ;
    float myEl_mesubveto10;
    float myEl_elvetolead;
    float myEl_elvetosub ;
    float myEl_reliso    ;
    float myEl_iso       ;
    float myEl_mvaTrig   ;
    float myEl_mvaNonTrig;
    float myEl_dZ_ee     ;
    float myEl_mass_ee   ;
    float myEl_inwindow_ee;
    float myEl_ptlead    ;
    float myEl_ptsub     ;
    float myEl_ptleadom  ;
    float myEl_ptsubom   ;
    float myEl_ptgg      ;
    float myEl_phomaxeta ;
    float myEl_sumpt3    ;
    float myEl_sumpt4    ;
    float myEl_dRtklead  ;
    float myEl_dRtksub   ;
    float myEl_MVAlead   ;
    float myEl_MVAsub    ;
    float myEl_CiClead   ;
    float myEl_CiCsub    ;
    float myEl_mgg       ;
    float myEl_MET       ;
    float myEl_METphi    ;
    float myEl_diphomva  ;
    float myEl_presellead ;
    float myEl_matchellead;
    float myEl_preselsub  ;
    float myEl_matchelsub ;
    float myEl_category ;
    float myEl_ElePho   ;
    float myEl_passelcuts ;

    // Di-pho MVA
    bool doDiphoMvaUpFront;
    bool useGbrDiphotonMva;
    std::string gbrDiphotonFile;
    std::string bdtTrainingPhilosophy;
    std::string bdtTrainingType;
    std::string photonLevelMvaUCSD  ;
    std::string eventLevelMvaUCSD   ;                    
    std::string photonLevel2011IDMVA_EB;
    std::string photonLevel2011IDMVA_EE;
    std::string eventLevelMvaMIT    ;
    std::string photonLevel2012IDMVA_EB;
    std::string photonLevel2012IDMVA_EE;
    std::string photonLevel2013IDMVA_EB;
    std::string photonLevel2013IDMVA_EE;
    std::string photonLevel2013_7TeV_IDMVA_EE;
    std::string photonLevel2013_7TeV_IDMVA_EB;
    std::vector<float> bdtCategoryBoundaries;

    // n-1 plots for VH hadronic tag 2011
    float  myVHhadLeadJPt;
    float  myVHhadSubJPt;
    float  myVHhaddEta;
    float  myVHhadZep;
    float  myVHhaddPhi;
    float  myVHhad_Mjj;
    float  myVHhad_Mgg;

    float myVBFDIPHObdt;
    float myVBFcombined;
    
    // n-1 plots for VBF tag 2011
    float  myVBFLeadJPt;
    float  myVBFSubJPt;
    float  myVBFLeadJEta;
    float  myVBFSubJEta;
    float  myVBFdEta;
    float  myVBFZep;
    float  myVBFdPhi;
    float  myVBFdPhiTrunc;
    float  myVBF_Mjj;
    float  myVBF_Mgg;
    float  myVBFLeadPhoPtOverM;
    float  myVBFSubPhoPtOverM;
    float  myVBFDiPhoPtOverM;
    float  myVBF_deltaPhiGamGam;
    float  myVBF_etaJJ;
    float  myVBF_MVA;
    float  myVBF_MVA0;
    float  myVBF_MVA1;
    float  myVBF_MVA2;
    float  myVBF_Pz;
    float  myVBF_S;
    float  myVBF_K1;
    float  myVBF_K2;

    // VBF Spin studies
    float myVBFSpin_Discriminant;
    float myVBFSpin_DeltaPhiJJ;
    float myVBFSpin_absDeltaPhiJJ;
    float myVBFSpin_CosThetaJ1;
    float myVBFSpin_absCosThetaJ1;
    float myVBFSpin_CosThetaJ2;
    float myVBFSpin_absCosThetaJ2;
    // Small deflection
    float myVBFSpin_CosThetaS;
    float myVBFSpin_absCosThetaS;
    float myVBFSpin_DeltaPhiJJS;
    float myVBFSpin_absDeltaPhiJJS;
    // Large deflection
    float myVBFSpin_CosThetaL;
    float myVBFSpin_absCosThetaL;
    float myVBFSpin_DeltaPhiJJL;
    float myVBFSpin_absDeltaPhiJJL;

    bool useGbrVbfMva;
    std::string gbrVbfFile, gbrVbfDiPhoFile;
    
    bool bookDiPhoCutsInVbf;
    bool mvaVbfSelection, mvaVbfUseDiPhoPt, mvaVbfUsePhoPt;
    bool combinedmvaVbfSelection;
    bool mvaVbfSpin;
    bool multiclassVbfSelection, vbfVsDiphoVbfSelection;
    TString mvaVbfDiphoWeights, mvaVbfDiphoMethod;
    TString mvaVbfWeights, mvaVbfMethod;
    TString mvaVbfSpinWeights, mvaVbfSpinMethod;
    std::vector<float> mvaVbfCatBoundaries;
    std::vector<float> multiclassVbfCatBoundaries0;
    std::vector<float> multiclassVbfCatBoundaries1;
    std::vector<float> multiclassVbfCatBoundaries2;

    // Smearings / corrections and systematics
    bool  doMCSmearing, doSystematics;

    bool recomputeBetas, recorrectJets, rerunJetMva, recomputeJetWp;
    bool applyJer, applyJecUnc, emulateJetResponse;
    float jerShift, jecShift;
    float jetResponseLumiStep;
    std::string jetHandlerCfg;
    bool applyBtagSF;
    bool  shiftBtagEffUp_bc;
    bool shiftBtagEffDown_bc;
    bool  shiftBtagEffUp_l;
    bool shiftBtagEffDown_l;
    bool applyLeptonSF;

    // progress
    int lastRun;
    int lastEvent;
    int lastLumi;

    // genLevels for calculating pdf errors post ws production (until PdfWeightSmearer works)
    float generatorPt_;
    float generatorY_;

 protected:
    void PreselectPhotons(LoopAll& l, int jentry);

    void SetNullHiggs(LoopAll& l);
    bool FindHiggsObjects(LoopAll& l);
    Bool_t GenMatchedPhoton(LoopAll& l, int ipho);

    bool ClassicCatsNm1Plots(LoopAll& l, int diphoton_nm1_id, float* smeared_pho_energy, float eventweight, float myweight);

    RooFuncReader *gbrVbfReader_, *gbrVbfDiphoReader_;
    
    // Exclusive tags
    TMVA::Reader *tmvaVbfDiphoReader_;

    int  categoryFromBoundaries(std::vector<float> & v, float val);
    int  categoryFromBoundaries2D(std::vector<float> & v1, std::vector<float> & v2, std::vector<float> & v3, float val1, float val2, float val3);
    
    bool VBFTag2013(int & ijet1, int & ijet2, LoopAll& l, int& diphotonVBF_id, float* smeared_pho_energy=0, bool vetodipho=true, bool kinonly=true, bool mvaselection=true, float eventweight=1, float myweight=1);
    bool FillDijetVariables(int & ijet1, int & ijet2, LoopAll& l, int diphoton_id, float* smeared_pho_energy=0,bool* jetid_flag=0, bool getAngles=0);
    // ICHEP2012
    bool VBFTag2012(int & ijet1, int & ijet2, LoopAll& l, int diphoton_id,
		    float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1,bool * jetid_flags=0);
    TMVA::Reader *tmvaVbfReader_;
    // Moriond 2012
    bool VBFTag2011(LoopAll& l, int diphoton_id, float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1);
    // VHhadronic Tag
    bool VHhadronicTag2011(LoopAll& l, int& diphoton_id, float* smeared_pho_energy=0, bool *jetid_flags=0,bool mvaselection=false,bool vetodipho=false, bool kinonly=false);
    // VHhadronic Tag
    bool VHhadronicTag2012(LoopAll& l, int& diphoton_id, float* smeared_pho_energy=0, bool *jetid_flags=0,bool mvaselection=false,bool vetodipho=false, bool kinonly=false);
    // VH category w btag
    bool VHhadronicBtag2012(LoopAll& l, int& diphoton_id, float* smeared_pho_energy=0, bool *jetid_flags=0,bool mvaselection=false,bool vetodipho=false, bool kinonly=false);
    //TTH leptonic category
    bool TTHleptonicTag2012(LoopAll& l, int& diphoton_id, float* smeared_pho_energy=0, bool *jetid_flags=0,bool mvaselection=false,bool vetodipho=false, bool kinonly=false);
    //TTH hadronic category
    bool TTHhadronicTag2012(LoopAll& l, int& diphoton_id, float* smeared_pho_energy=0, bool *jetid_flags=0, bool mvaselection=false,bool vetodipho=false, bool kinonly=false);
    //only one tth category for 7 TeV
    bool TTHTag7TeV(LoopAll& l, int& diphoton_id, float* smeared_pho_energy=0, bool *jetid_flags=0,bool mvaselection=false,bool vetodipho=false, bool kinonly=false);

    //btag syst
    float BtagReweight(LoopAll& l, bool shiftBtagEffUp_bc, bool shiftBtagEffDown_bc, bool shiftBtagEffUp_l, bool shiftBtagEffDown_l,int WP);
    float BtagReweight2013(LoopAll& l, bool shiftBtagEffUp_bc, bool shiftBtagEffDown_bc, bool shiftBtagEffUp_l, bool shiftBtagEffDown_l,int WP);

    //electrons-muons SF
    float ElectronSFReweight(LoopAll& l);
    float MuonSFReweight(LoopAll& l);

    int nBMC,nCMC,nLMC;
    int nBtagB,nBtagC,nBtagL;
    void computeBtagEff(LoopAll& l);

    // VBF Spin studies
    TMVA::Reader *tmvaVbfSpinReader_;

    // Moriond 2012
    bool ElectronTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1);
    // ICHEP2012
    bool ElectronTag2012(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, ofstream& lep_sync, bool nm1=false, float eventweight=1, float myweight=1);
    // HCP 2012
    bool ElectronTag2012B(LoopAll& l, int& diphotonVHlep_id, int& el_ind, int& elVtx, int& el_cat, float* smeared_pho_energy, ofstream& lep_sync, bool mvaselection=true, float phoidMvaCut=-0.2, float eventweight=1.0, std::vector<float>  smeared_pho_weight=std::vector<float>(), bool fillHist=false, bool vetodipho=false, bool kinonly=false);
    bool ElectronStudies2012B(LoopAll& l, float* smeared_pho_energy, bool mvaselection, float phoidMvaCut, float eventweight=1, float myweight=1, int jentry=-1);
    bool ElectronTagStudies2012(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, bool nm1=true, float eventweight=1, float myweight=1, int jentry=-1);
    void ZWithFakeGammaCS(LoopAll& l, float* smeared_pho_energy);
    void ControlPlotsElectronTag2012B(LoopAll& l, TLorentzVector lead_p4, TLorentzVector sublead_p4, int el_ind, float bdtoutput, float evweight, std::string label);
    void ControlPlotsMetTag2012B(LoopAll& l, TLorentzVector lead_p4, TLorentzVector sublead_p4, float bdtoutput, float evweight, std::string label);

    // Moriond 2012
    bool MuonTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1);
    // ~ ICHEP2012
    bool MuonTag2012(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, ofstream& lep_sync, bool nm1=false, float eventweight=1, float myweight=1);
    // HCP2012
    bool MuonTag2012B(LoopAll& l, int& diphotonVHlep_id, int& mu_ind, int& muVtx, int& mu_cat, float* smeared_pho_energy, ofstream& lep_sync, bool mvaselection=true, float phoidMvaCut=-0.2, float eventweight=1.0, std::vector<float>  smeared_pho_weight=std::vector<float>(), bool fillHist=false, bool vetodipho=false, bool kinonly=false);
    void ControlPlotsMuonTag2012B(LoopAll& l, TLorentzVector lead_p4, TLorentzVector sublead_p4, int mu_ind, float bdtoutput, float evweight, std::string label);


    // ~ ICHEP2012
    bool METTag2012(LoopAll& l, int& diphotonVHmet_id , float* smeared_pho_energy);  //met at analysis step
    void MetCorrections2012(LoopAll& l);
    void MetCorrections2012_Simple(LoopAll& l, TLorentzVector lead_p4, TLorentzVector sublead_p4);

    // HCP2012
    bool METTag2012B(LoopAll& l, int& diphotonVHmet_id, int& met_cat, float* smeared_pho_energy, ofstream& met_sync, bool mvaselection=true, float phoidMvaCut=-0.2, bool useUncor=true);


    int GenMatch(LoopAll& l, TLorentzVector* recop4);
    bool PhotonMatchElectron(LoopAll& l, TLorentzVector* pho_p4);
    bool PhotonMatchElectron(LoopAll& l, TLorentzVector* pho_p4, int& el_match_ind);
    bool HLTPhotonPreselection(LoopAll& l, TLorentzVector* pho_p4, int phoind);

    ofstream met_sync;
    ofstream lep_sync;

    // Pile-up reweighing
    void loadPuMap(const char * fname, TDirectory * dir, TH1 * target=0);
    void loadPuWeights(int typid, TDirectory * dir, TH1 * target=0);
    void load2DPuWeights(int typid, TDirectory* dir, std::vector<TH1*> target);
    float getPuWeight(int npu, int sample_type, SampleContainer* container, bool warnMe, int run=0);
    TH1 * puTargetHist;
    std::vector<TH1*> puTargetHists;

    std::string name_;

    // Beamsport reweighting
    float BeamspotReweight(double vtxZ, double genZ);
    void saveBSTrees(LoopAll &l, float evweight, int category, TLorentzVector Higgs, TVector3 *chosenVtx, TVector3 *genVtx, float diphobdt_output=-100.);

    // Pt reweighting
    float PtReweight(double pt, int cur_type);

    // Track systematics
    float ComputeEventScaleError(LoopAll& l, int ipho1, int ipho2, float & scale1, float & scale1_err, float & scale2, float & scale2_err);
    float ComputeEventSmearError(LoopAll& l, int ipho1, int ipho2, float & smear1, float & smear1_err, float & smear2, float & smear2_err);
    pair<double,double> ComputeNewSigmaMs(LoopAll &l, int ipho1, int ipho2, int ivtx, float syst_shift);
    void saveDatCardTree(LoopAll& l, int cur_type, int category, int inc_cat, float evweight, int ipho1, int ipho2, int ivtx, TLorentzVector lead_p4, TLorentzVector sublead_p4, bool isCutBased, string proc, double sigmaMrv=0., double sigmaMwv=0., double sigmaMeonly=0., float vtxProb=0., string trainPhi="", string bdtType="", float lead_id_mva=0., float sublead_id_mva=0.);
    
    // Save spin trees
    void saveSpinTree(LoopAll &l, int category, float evweight, TLorentzVector Higgs, TLorentzVector lead_p4, TLorentzVector sublead_p4, int ipho1, int ipho2, int diphoton_id, float vtxProb, bool isCorrectVertex);
    void saveSpinTree(LoopAll& l, int category, float evweight, TLorentzVector Higgs, TLorentzVector lead_p4, TLorentzVector sublead_p4, int ipho1, int ipho2, float diphobdt, double sigmaMrv, double sigmaMwv, double lead_sigmaE, double lead_sigmaE_nosmear, double sublead_sigmaE, double sublead_sigmaE_nosmear, float vtxProb, float lead_id_mva, float sublead_id_mva);

    // Save VBF trees
    void saveVBFTree(LoopAll& l, int category, float evweight,float diphobdt);

    
    // Vertex analysis
    HggVertexAnalyzer vtxAna_;
    HggVertexFromConversions vtxConv_;
    void reVertex(LoopAll & l);

    // Jets
    JetHandler * jetHandler_;
    void postProcessJets(LoopAll & l, int vtx=-1);
    void switchJetIdVertex(LoopAll &l, int ivtx);

    std::map<int, vector<double> > weights;
    std::map<int, std::vector<vector<double> > > rd_weights;
    int trigCounter_;

    // MC smearing and correction machinery
    void applyGenLevelSmearings(double & genLevWeight, const TLorentzVector & gP4, int npu, int sample_type, BaseGenLevelSmearer * sys=0, float systshift=0.);

    void applySinglePhotonSmearings(std::vector<float> & smeared_pho_energy, std::vector<float> & smeared_pho_r9, std::vector<float> & smeared_pho_weight,
				    int cur_type, const LoopAll & l, const float * energyCorrected, const float * energyCorrectedError,
				    BaseSmearer * sys=0, float syst_shift=0.);

    void fillDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs,
                  float & lead_r9, float & sublead_r9, TVector3 *& vtx, const float * energy,
                  const LoopAll & l, int diphoton_id);
    void fillDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs,
                  float & lead_r9, float & sublead_r9, TVector3 *& vtx, const float * energy,
                  const LoopAll & l, int leadind, int subleadind, int myvtx);

    void applyDiPhotonSmearings(TLorentzVector & Higgs, TVector3 & vtx, int category, int cur_type, const TVector3 & truevtx,
				float & evweight, float & idmva1, float & idmva2,
				BaseDiPhotonSmearer * sys=0, float syst_shift=0.);

    std::pair<TLorentzVector, TLorentzVector> GetVBF_IntermediateBoson(TLorentzVector& Pho1, TLorentzVector& Pho2, TLorentzVector& Jet1, TLorentzVector& Jet2);
    Double_t GetPerpendicularAngle(TLorentzVector& ref, TLorentzVector& v1, TLorentzVector& v2);
    void VBFAngles(TLorentzVector& gamma1, TLorentzVector& gamma2, TLorentzVector& J1, TLorentzVector& J2);
    
    double getCosThetaCS(TLorentzVector, TLorentzVector,int);
    double getCosThetaHX(TLorentzVector, TLorentzVector,int);

    std::vector<BaseSmearer *> photonSmearers_;
    std::vector<BaseSmearer *> systPhotonSmearers_;
    std::vector<BaseDiPhotonSmearer *> diPhotonSmearers_;
    std::vector<BaseDiPhotonSmearer *> systDiPhotonSmearers_;
    std::vector<BaseGenLevelSmearer *> genLevelSmearers_;
    std::vector<BaseGenLevelSmearer *> systGenLevelSmearers_;

    // common smearers
    void addResolSmearer(EnergySmearer * theSmear);
    EnergySmearer *eScaleDataSmearer ; // corrections for energy scale data
    EnergySmearer *eScaleSmearer, *eScaleCorrSmearer;      // corrections for energy scale  MC
    std::vector<EnergySmearer *> eScaleSmearers_;
    EnergySmearer *eResolSmearer, *eResolCorrSmearer;
    std::vector<std::pair<EnergySmearer*,EnergySmearerExtrapolation *> > eResolSmearers_;
    EnergySmearer *eCorrSmearer;      // corrections for energy scale  MC
    std::vector<float> corrected_pho_energy;
    std::vector<PhotonReducedInfo> photonInfoCollection;



    Float_t *energyCorrected;
    Float_t *energyCorrectedError;

    TMVA::Reader *tmvaPerVtxReader_;
    TMVA::Reader *tmvaPerEvtReader_;

    MassResolution *massResolutionCalculator;

    void VHLepTag2013(LoopAll& l, int & diphotonVHlep_id, bool & VHlep1event, bool & VHlep2event, bool mvaselection, int & mu_ind, int & muVtx, int VHmuevent_cat, int & el_ind, int & elVtx, int VHelevent_cat, float* smeared_pho_energy, float phoidMvaCut, float eventweight, std::vector<float> smeared_pho_weight, bool isSyst, bool vetodipho = false, bool kinonly = false);
    int VHNumberOfJets(LoopAll& l, int diphotonVHlep_id, int vertex, bool VHelevent_prov, bool VHmuevent_prov, int el_ind, int mu_ind, float* smeared_pho_energy);


    // For semi-parametric Regression 
    
    std::vector<float> _vals;
    
    HybridGBRForest *_foresteb;
    HybridGBRForest *_forestee;

    HybridGBRForestD *_forestDeb;
    HybridGBRForestD *_forestDee;

    RooRealVar *_mean;
    RooRealVar *_tgt;
    RooRealVar *_sigma;
    RooRealVar *_n1;
    RooRealVar *_n2;  
    
    RooAbsReal *_meanlim;
    RooAbsReal *_sigmalim;
    RooAbsReal *_n1lim;
    RooAbsReal *_n2lim;        
    RooAbsReal *alpha1;
    RooAbsReal *alpha2;        
    
    RooAbsPdf *_pdf;
    
    RooArgList _args;
    

    //TFile *fgbr;
    //GBRForest *fReadereb;
    //GBRForest *fReaderebvariance;
    //GBRForest *fReaderee;
    //GBRForest *fReadereevariance;
    std::pair<int, int> SelectBtaggedAndHighestPtJets(LoopAll& l,int diphoton_id,const TLorentzVector& leadpho,const TLorentzVector& subleadpho, Bool_t * jetid_flags=0);
    float getDiphoBDTOutput(LoopAll &l,int diphoton_id,TLorentzVector lead_p4, TLorentzVector sublead_p4,std::string bdtTrainingPhilosophy="MIT");
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
