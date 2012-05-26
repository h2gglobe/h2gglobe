#ifndef __MVAANALYSIS__
#define __MVAANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include "InterferenceSmearer.h"
#include "MassFactorizedMvaAnalysis.h"
#include "MassResolution.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class MvaAnalysis : public MassFactorizedMvaAnalysis
{
 public:
    
    MvaAnalysis();
    virtual ~MvaAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    void GetBranches(TTree *, std::set<TBranch *>& );
    
    void FillSignalLabelMap();
    int GetBDTBoundaryCategory(float,bool,bool);
    std::string GetSignalLabel(int);

    int SignalType(int);
    void SetBDTInputVariables(TLorentzVector*, TLorentzVector*, double, double, MassResolution* ,double, double, double, double, int cat = 0);
    void SetTree(TTree *);
    void SetBDTInputTree(TTree *);
    
   
    double mHMaximum, mHMinimum, mHStep; 
    float systRange;
    int   nSystSteps;   
    int   nEtaCategories, nR9Categories, nPtCategories;
    float massMin, massMax;
    int nDataBins;  
    float signalRegionWidth;
    float sidebandWidth;
    float sidebandShift;
    float massSidebandMin;
    float massSidebandMax;
    int numberOfSidebandGaps;
    int numberOfSidebands;
    int numberOfSidebandsForAlgos;
    
    std::string bdtTrainingPhilosophy;
    std::string photonLevelMvaUCSD  ;
    std::string eventLevelMvaUCSD   ;                    
    std::string photonLevelMvaMIT_EB;
    std::string photonLevelMvaMIT_EE;
    std::string eventLevelMvaMIT    ;
    std::string kfacHist;
    std::string MVAtype;

    int nMasses;

    bool doTraining;
    bool makeTrees;
    bool splitSignalSample;
    bool splitBackgroundSample;
    double trainingMassHypothesis;
    //int nMassPt;
    std::string names[9];
    std::string BDTnames[9];
    double masses[9];
    std::vector<float> bkg_masses; //This is for Nicks training method

    std::string mvaWeightsFolder;

    // Bin edges pre-defined
    bool rederiveOptimizedBinEdges;
    std::vector<double> VbfBinEdges_110, GradBinEdges_110, AdaBinEdges_110;
    std::vector<double> VbfBinEdges_115, GradBinEdges_115, AdaBinEdges_115;
    std::vector<double> VbfBinEdges_120, GradBinEdges_120, AdaBinEdges_120;
    std::vector<double> VbfBinEdges_125, GradBinEdges_125, AdaBinEdges_125;
    std::vector<double> VbfBinEdges_130, GradBinEdges_130, AdaBinEdges_130;
    std::vector<double> VbfBinEdges_135, GradBinEdges_135, AdaBinEdges_135;
    std::vector<double> VbfBinEdges_140, GradBinEdges_140, AdaBinEdges_140;
    std::vector<double> VbfBinEdges_150, GradBinEdges_150, AdaBinEdges_150;

 protected:

    float tmvaGetVal(double,double,float);
    void fillTMVATrees(LoopAll &l,float,float,int,float,int);

    virtual void FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA, int category, float weight, bool isCorrectVertex);

    virtual void AccumulateSyst(int cur_type, float mass, float diphotonMVA, int category, float weight,
				std::vector<double> & mass_errors,
				std::vector<double> & mva_errors,
				std::vector<int>    & categories,
				std::vector<double> & weights);
    virtual void FillRooContainerSyst(LoopAll& l, const std::string & name,int cur_type,
				      std::vector<double> & mass_errors, std::vector<double> & mva_errors,
				      std::vector<int>    & categories, std::vector<double> & weights);
    
    std::string name_;
    float nevents, sumwei, sumaccept, sumsmear, sumev; 
    
    int nPhotonCategories_;
    // Vertex analysis
    HggVertexAnalyzer vtxAna_;
    HggVertexFromConversions vtxConv_;
    
    // RooStuff
    //RooContainer *rooContainer;

    ofstream eventListText;

    std::map<int,std::string> signalLabels;

    //MVA variables
    float _log_H_pt;     
    float _H_ptOverM;     
    float _H_eta;        
    float _d_phi;        
    float _cos_d_phi;        
    float _max_eta;      
    float _min_r9;       
    float _pho1_eta;    
    float _pho2_eta;     
    float _pho1_ptOverM; 
    float _pho2_ptOverM; 
    float _sigmaMOverM;
    float _sigmaMOverM_wrongVtx;
    float _deltaMOverSigmaM;
    float _deltaMOverM; 
    float _mgg;          
    float _pho1_phi;     
    float _pho1_pt;      
    float _pho1_r9;      
    float _pho2_phi;     
    float _pho2_pt;      
    float _pho2_r9;      
    float _H_pt;         
    float _Ht;           
    float _d_eta;        
    float _mod_d_eta;    
    float _cos_theta_star;
    float _vtx_prob;           
    float _wt;           
    float _bdtoutput;
    int _cat;           
    int _sideband;           

    // extras for trees
    int _vbf;
    int _cur_type;

    // To make trees for toys
    TFile *treeFile_;
    TTree *dataTree_;
    TTree *bkgTree_;
    TTree *sigTree_[9]; 

    //vector<double> weights;
    TFile *kfacFile;
    
    TMVA::Reader * tmvaReader_;

    // Used to create the training samples
    TTree * signalTrainTree_[2];
    TTree * signalTestTree_[2];
    TTree * backgroundTrainTree_2pt_[2];
    TTree * backgroundTestTree_2pt_[2];
    TTree * backgroundTrainTree_7pt_[2];
    TTree * backgroundTestTree_7pt_[2];
    TFile * mvaFile_;
    
};

#endif


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
