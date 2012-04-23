#ifndef __PHOTONANALYSIS__
#define __PHOTONANALYSIS__

#include "BaseAnalysis.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/PhotonInfo.h"

#include "TriggerSelection.h"
#include "EnergySmearer.h"

#include "MassResolution.h"

#include "TMVA/Reader.h"
#include "PhotonFix.h"
// #include "HiggsToGammaGamma/interface/GBRForest.h"
//#include "../../../../HiggsToGammaGamma/interface/GBRForest.h"
//#include "HiggsAnalysis/HiggsToGammaGamma/interface/GBRForest.h"

// ------------------------------------------------------------------------------------
class PhotonAnalysis : public BaseAnalysis 
{
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
    virtual void Analysis(LoopAll&, Int_t);

    virtual void ResetAnalysis();
    
    //void GetRegressionCorrections(LoopAll&);  
    //  void GetRegressionCorrections(LoopAll&);    
    // Public parameters to be read from config file
    float zero_;
    VertexAlgoParameters vtxAlgoParams;  
    std::vector<std::string> vtxVarNames;
    bool useDefaultVertex;
    float forcedRho;
    bool applyPtoverM;  
    bool doTriggerSelection; 
    std::vector<TriggerSelection> triggerSelections;
    
    // Options
    bool includeVBF;
    bool includeVHhad;
    bool includeVHlep;
    bool reRunCiCForData;

    float leadEtCut;
    float leadEtVBFCut;
    float leadEtVHhadCut;
    float leadEtVHlepCut;
    float subleadEtCut;
    float subleadEtVBFCut;
    float subleadEtVHhadCut;
    float subleadEtVHlepCut;
    int nVBFEtaCategories;
    int nVHhadEtaCategories;
    
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

    std::vector<int> pho_acc;
    std::vector<int> pho_presel;
    std::vector<int> pho_presel_lead;
    std::vector<float> pho_et;
    // Other options
    bool runStatAnalysis;
    TString puHist, puMap, puTarget;//name of pileup reweighting histogram

    enum BkgCategory{promptprompt,promptfake,fakefake};
    bool keepPP, keepPF, keepFF;

    std::string scale_offset_error_file, smearing_file;
    std::string energyCorrectionMethod;
    //std::string massResolutionFileName;

    std::string tmvaPerVtxMethod;                           
    std::string tmvaPerVtxWeights;                  
    std::string tmvaPerEvtMethod;                   
    std::string tmvaPerEvtWeights;                  
                                                        
    bool mvaVertexSelection, addConversionToMva;     

    // PhotonFix
    std::string photonFixDat;
    std::string regressionFile;
    
    int   nEtaCategories, nR9Categories, nPtCategories;
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
  
    // n-1 plots for VH hadronic tag 2011
    float  myVHhadLeadJPt;
    float  myVHhadSubJPt;
    float  myVHhaddEta;
    float  myVHhadZep;
    float  myVHhaddPhi;
    float  myVHhad_Mjj;
    float  myVHhad_Mgg;
    
    // n-1 plots for VBF tag 2011 
    float  myVBFLeadJPt;
    float  myVBFSubJPt;
    float  myVBFdEta;
    float  myVBFZep;
    float  myVBFdPhi;
    float  myVBF_Mjj;
    float  myVBF_Mgg;

    // Smearings / corrections and systematics
    bool  doMCSmearing, doSystematics;

 protected:
    void PreselectPhotons(LoopAll& l, int jentry);
    float GetSmearSigma(float eta, float r9, int epoch=0);
    
    void SetNullHiggs(LoopAll& l);
    bool FindHiggsObjects(LoopAll& l);
    Bool_t GenMatchedPhoton(LoopAll& l, int ipho);
    
    bool ClassicCatsNm1Plots(LoopAll& l, int diphoton_nm1_id, float* smeared_pho_energy, float eventweight, float myweight);
    
    // Exclusive tags
    bool VBFTag2011(LoopAll& l, int diphoton_id, float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1);
    bool VHhadronicTag2011(LoopAll& l, int diphoton_id, float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1);
    bool ElectronTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1);
    bool MuonTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy=0, bool nm1=false, float eventweight=1, float myweight=1);
    
    // Pile-up reweighing
    void loadPuMap(const char * fname, TDirectory * dir, TH1 * target=0);
    void loadPuWeights(int typid, TDirectory * dir, TH1 * target=0);
    float getPuWeight(int npu, int sample_type, bool warnMe); 
	
    std::string name_;
    
    // Vertex analysis
    HggVertexAnalyzer vtxAna_;
    HggVertexFromConversions vtxConv_;
    
    std::map<int, vector<double> > weights;
    int trigCounter_;

    // MC smearing and correction machinery
    void applyGenLevelSmearings(double & genLevWeight, const TLorentzVector & gP4, int npu, int sample_type, BaseGenLevelSmearer * sys=0, float systshift=0.);
    
    void applySinglePhotonSmearings(std::vector<float> & smeared_pho_energy, std::vector<float> & smeared_pho_r9, std::vector<float> & smeared_pho_weight,
				    int cur_type, const LoopAll & l, const float * energyCorrected, const float * energyCorrectedError,
				    BaseSmearer * sys=0, float syst_shift=0.);
    
    void fillDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs, float & lead_r9, float & sublead_r9, TVector3 *& vtx, 
		      const float * energy, const LoopAll & l,  int diphoton_id);
    
    void applyDiPhotonSmearings(TLorentzVector & Higgs, TVector3 & vtx, int category, int cur_type, const TVector3 & truevtx, 
				float & evweight, float & idmva1, float & idmva2,
				BaseDiPhotonSmearer * sys=0, float syst_shift=0.);
    
    std::vector<BaseSmearer *> photonSmearers_;
    std::vector<BaseSmearer *> systPhotonSmearers_;
    std::vector<BaseDiPhotonSmearer *> diPhotonSmearers_;
    std::vector<BaseDiPhotonSmearer *> systDiPhotonSmearers_;
    std::vector<BaseGenLevelSmearer *> genLevelSmearers_;
    std::vector<BaseGenLevelSmearer *> systGenLevelSmearers_;
    
    // common smearers
    EnergySmearer *eScaleDataSmearer ; // corrections for energy scale data
    EnergySmearer *eScaleSmearer;      // corrections for energy scale  MC
    EnergySmearer *eCorrSmearer;      // corrections for energy scale  MC
    std::vector<float> corrected_pho_energy;
    std::vector<PhotonReducedInfo> photonInfoCollection;
    
    Float_t *energyCorrected;
    Float_t *energyCorrectedError;

    vector<string> tmvaPerVtxVariables_;
    TMVA::Reader *tmvaPerVtxReader_;
    TMVA::Reader *tmvaPerEvtReader_;

    MassResolution *massResolutionCalculator;

    int DiphotonMVASelection(LoopAll &l, HggVertexAnalyzer & vtxAna, Float_t & diphoMVA,  
                             Float_t minLeadingMVA=-0.3, Float_t minSubleadingMVA=-0.3, Float_t leadPtMin=30, 
                             Float_t subleadPtMin=20, std::string type="UCSD", int ncategories=7,
                             bool applyPtoverM=true, float *pho_energy_array=0, bool split=false);
    int DiphotonMVAEventClass(LoopAll &l, float diphoMVA, int nCat, std::string type, int EBEB=1);
    
    //TFile *fgbr;
    //GBRForest *fReadereb;
    //GBRForest *fReaderebvariance;
    //GBRForest *fReaderee;
    //GBRForest *fReadereevariance;      

};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
