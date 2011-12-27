#ifndef __STATANALYSISEXCLUSIVE__
#define __STATANALYSISEXCLUSIVE__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class StatAnalysisExclusive : public PhotonAnalysis 
{
public:
	
	StatAnalysisExclusive();
	virtual ~StatAnalysisExclusive();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	void Init(LoopAll&);
	void Term(LoopAll&);
	
	void GetBranches(TTree *, std::set<TBranch *>& );
	
	virtual bool SelectEvents(LoopAll&, int);
	virtual void Analysis(LoopAll&, Int_t);
	
	// Options
  bool includeVBF;
  bool includeVHad;
  bool includeVHlep;
  bool reRunCiCForData;

  float leadEtCut;
  float leadEtVBFCut;
  float leadEtVHadCut;
  float leadEtVHlepCut;
  float subleadEtCut;
  float subleadEtVBFCut;
  float subleadEtVHadCut;
  float subleadEtVHlepCut;
  int nVBFEtaCategories;
	int nVHadEtaCategories;
  std::string efficiencyFile;
	
	// EnergySmearer::energySmearingParameters eSmearPars; // gone to PhotonAnalysis GF
	EfficiencySmearer::efficiencySmearingParameters effSmearPars;
	DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

	double GetDifferentialKfactor(double, int);

	void FillSignalLabelMap();
	std::string GetSignalLabel(int) ;

  // for N-1
  float  myVBFLeadJPt;
  float  myVBFSubJPt;
  float  myVBFdEta;
  float  myVBFZep;
  float  myVBFdPhi;
  float  myVBF_Mjj;

  float  myVBF_Mgg;

  float  myVHadLeadJPt;
  float  myVHadSubJPt;
  float  myVHaddEta;
  float  myVHadZep;
  float  myVHaddPhi;
  float  myVHad_Mjj;

  float  myVHad_Mgg;

  float  myAllLeadJPt;
  float  myAllSubJPt;
  float  myAllLeadJEta;
  float  myAllSubJEta;
  float  myAlldEta;
  float  myAllZep;
  float  myAlldPhi;
  float  myAll_Mjj;

  float  myAll_nvtx;

  float  myAll_Mgg;
  float  myAllPtHiggs;
  float  myAll_phoet1;
  float  myAll_phoet2;
  float  myAll_phoetom1;
  float  myAll_phoetom2;

  float  myInclusive_Mgg;
  float  myInclusivePtHiggs;

  std::pair<int, int> Select2HighestPtJets(LoopAll&, TLorentzVector& leadpho, TLorentzVector& subleadpho, float jtLMinPt, float jtTMinPt, TLorentzVector* jet3=0);
  int RescaleJetEnergy(LoopAll&);

	bool  doMCSmearing;
	bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
	bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, doKFactorSmear;
	float systRange;
	int   nSystSteps;   
	int   nEtaCategories, nR9Categories, nPtCategories;
	float massMin, massMax;
	int nDataBins;	
	//float smearing_sigma_EBHighR9       ;
	//float smearing_sigma_EBLowR9        ;
	//float smearing_sigma_EEHighR9       ;
	//float smearing_sigma_EELowR9        ;
	//float smearing_sigma_error_EBHighR9 ;
	//float smearing_sigma_error_EBLowR9  ;
	//float smearing_sigma_error_EEHighR9 ;
	//float smearing_sigma_error_EELowR9  ;
	
	std::string kfacHist;

	TH1D *thm110,*thm120,*thm130,*thm140;
	int nMasses;

protected:
	std::vector<BaseSmearer *> photonSmearers_;
	std::vector<BaseSmearer *> systPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> diPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> systDiPhotonSmearers_;
	std::vector<BaseGenLevelSmearer *> genLevelSmearers_;
	std::vector<BaseGenLevelSmearer *> systGenLevelSmearers_;
	
	EnergySmearer /* *eScaleSmearer,*/ *eResolSmearer ; // moved to PhotonAnalysis GF 
	EfficiencySmearer *idEffSmearer, *r9Smearer;
	DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer;
	KFactorSmearer * kFactorSmearer;
	
	std::string name_;
	std::map<int,std::string> signalLabels;
	float nevents, sumwei, sumaccept, sumsmear, sumev; 
	
	int nCategories_;
	int nPhotonCategories_;
	int diPhoCounter_;
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	// RooStuff
	RooContainer *rooContainer;

	ofstream eventListText;
	ofstream eventListTextVBF;
	//vector<double> weights;
	TFile *kfacFile;

	//JIM


#define MPDEBUGCOPY 0
#define MPDEBUG 0

	LoopAll* ll;



TFile * opfile;
TTree * optree;
// Declaration of leaf types for jim's OLD ntuple

Float_t   t_j1pt;
Float_t   t_j1eta;
Float_t   t_j1phi;

Float_t   t_j2pt;
Float_t   t_j2eta;
Float_t   t_j2phi;

Float_t   t_j3pt;
Float_t   t_j3eta;
Float_t   t_j3phi;

Float_t   t_jjmass;
Float_t   t_jjdeta;
Float_t   t_jjzep;
Float_t   t_jjdphi;


Int_t           t_event;
Int_t           t_run;
Int_t           t_lumis;

Int_t           t_itype;
Int_t           t_processid;
Float_t         t_w;
Float_t         t_wpu;
Float_t         t_mass;
Float_t         t_dmom0;
Float_t         t_dmom;
Float_t         t_deltaM;
Int_t           t_category;
Int_t         t_diphocat2r92eta;
Int_t         t_diphosubcat4;
Int_t           t_barrel;
Float_t         t_diphor9;
Float_t         t_diphoeta;
Float_t         t_costhetastar;
Float_t         t_diphopt;
Float_t         t_diphopz;
Float_t         t_deltar;
Float_t         t_etamax;
Float_t         t_etamin;
Float_t         t_deta;
Int_t           t_leadcat;
Float_t         t_leadr9;
Float_t         t_leadeta;
Float_t         t_leadsceta;
Float_t         t_leadpt;
Int_t           t_leadgenmatch;
Int_t           t_leadbarrel;
Float_t         t_leadsee;
Float_t         t_leadspp;
Float_t         t_leadpi0nn;
Int_t           t_subleadcat;
Float_t         t_subleadr9;
Float_t         t_subleadeta;
Float_t         t_subleadsceta;
Float_t         t_subleadpt;
Int_t           t_subleadgenmatch;
Int_t           t_subleadbarrel;
Float_t         t_subleadsee;
Float_t         t_subleadspp;
Float_t         t_subleadpi0nn;
Int_t           t_leadtrkiso;
Float_t         t_leadecaliso;
Float_t         t_leadhcaliso;
Float_t         t_leadhovere;
Int_t           t_subleadtrkiso;
Float_t         t_subleadecaliso;
Float_t         t_subleadhcaliso;
Float_t         t_subleadhovere;
Float_t         t_leadtrkdeltar0;
Float_t         t_leadtrkdeltar15;
Float_t         t_leadtrkdeltar20;
Float_t         t_leadtrkdeltar25;
Float_t         t_leadtrkecone20;
Float_t         t_leadtrkecone25;
Float_t         t_leadtrkecone30;
Float_t         t_leadtrkecone35;
Float_t         t_leadtrkecone40;
Float_t         t_leadtrkecone45;
Float_t         t_leadecalhitsJ_060330;
Float_t         t_leadecal_048_32_16_018;
Float_t         t_leadecal_048_32_14_015;
Float_t         t_leadecal_048_32_18_018;
Float_t         t_leadecal_075_36_20_015;
Float_t         t_leadecal_060_32_20_018;
Float_t         t_leadecal_075_36_12_015;
Float_t         t_subleadecal_048_32_16_018;
Float_t         t_subleadecal_048_32_14_015;
Float_t         t_subleadecal_048_32_18_018;
Float_t         t_subleadecal_075_36_20_015;
Float_t         t_subleadecal_060_32_20_018;
Float_t         t_subleadecal_075_36_12_015;
Float_t         t_leadtrkplusecal;
Float_t         t_subleadtrkplusecal;
Float_t         t_subleadtrkdeltar0;
Float_t         t_subleadtrkdeltar15;
Float_t         t_subleadtrkdeltar20;
Float_t         t_subleadtrkdeltar25;
Float_t         t_subleadtrkecone20;
Float_t         t_subleadtrkecone25;
Float_t         t_subleadtrkecone30;
Float_t         t_subleadtrkecone35;
Float_t         t_subleadtrkecone40;
Float_t         t_subleadtrkecone45;
Float_t         t_subleadecalhitsJ_060330;

//new for the ntuple
Float_t t_rho;
Float_t t_nvtx;
Float_t t_pvtx;
Float_t t_sigma_mz;
Float_t t_leadphoidmitmva;
Float_t t_subleadphoidmitmva;
Float_t t_diphomitmva;
Float_t t_bsZ;

Float_t t_leadcutindex;
Float_t t_subleadcutindex;
Float_t t_leadci6cindex;
Float_t t_subleadci6cindex;
Float_t t_leadci6cpfindex;
Float_t t_subleadci6cpfindex;
Float_t t_leadci6cpfmva;
Float_t t_leadci6cpfmvacat;
Float_t t_leadci6cpfmvaptom;
Float_t t_leadci6cpfmvaptom2;
Float_t t_subleadci6cpfmva;
Float_t t_subleadci6cpfmvacat;
Float_t t_subleadci6cpfmvaptom;
Float_t t_subleadci6cpfmvaptom2;
Float_t t_diphomva;
Float_t t_diphomva2;
Float_t t_leadci4cindex;
Float_t t_subleadci4cindex;

////lead vars
Int_t t_leadpixel;
Float_t t_leadsieie;
Float_t t_leadtrkhollowdr03;
Float_t t_leadtrkhollowdr04;
Float_t t_leadtrksoliddr03;
Float_t t_leadtrksoliddr04;
Float_t t_leadecaldr03;
Float_t t_leadecaldr04;
Float_t t_leadhcaldr03;
Float_t t_leadhcaldr04;

Float_t t_lead_beta_recvtx_out30;
Float_t t_lead_beta_recvtx_out40;
Float_t t_lead_beta_badvtx_out30;
Float_t t_lead_beta_badvtx_out40;
Float_t t_sublead_beta_recvtx_out30;
Float_t t_sublead_beta_recvtx_out40;
Float_t t_sublead_beta_badvtx_out30;
Float_t t_sublead_beta_badvtx_out40;

Float_t t_lead_ecal_in035out20;
Float_t t_lead_ecal_in035out25;
Float_t t_lead_ecal_in035out30;
Float_t t_lead_ecal_in035out35;
Float_t t_lead_ecal_in035out40;
Float_t t_lead_ecal_in035out45;
Float_t t_lead_ecal_in035out50;
Float_t t_sublead_ecal_in035out20;
Float_t t_sublead_ecal_in035out25;
Float_t t_sublead_ecal_in035out30;
Float_t t_sublead_ecal_in035out35;
Float_t t_sublead_ecal_in035out40;
Float_t t_sublead_ecal_in035out45;
Float_t t_sublead_ecal_in035out50;

Float_t t_lead_ecal_rhocorr12_in035out20;
Float_t t_lead_ecal_rhocorr12_in035out25;
Float_t t_lead_ecal_rhocorr12_in035out30;
Float_t t_lead_ecal_rhocorr12_in035out35;
Float_t t_lead_ecal_rhocorr12_in035out40;
Float_t t_lead_ecal_rhocorr12_in035out45;
Float_t t_lead_ecal_rhocorr12_in035out50;
Float_t t_sublead_ecal_rhocorr12_in035out20;
Float_t t_sublead_ecal_rhocorr12_in035out25;
Float_t t_sublead_ecal_rhocorr12_in035out30;
Float_t t_sublead_ecal_rhocorr12_in035out35;
Float_t t_sublead_ecal_rhocorr12_in035out40;
Float_t t_sublead_ecal_rhocorr12_in035out45;
Float_t t_sublead_ecal_rhocorr12_in035out50;

//sublead vars
Int_t t_subleadpixel;
Float_t t_subleadsieie;
Float_t t_subleadtrkhollowdr03;
Float_t t_subleadtrkhollowdr04;
Float_t t_subleadtrksoliddr03;
Float_t t_subleadtrksoliddr04;
Float_t t_subleadecaldr03;
Float_t t_subleadecaldr04;
Float_t t_subleadhcaldr03;
Float_t t_subleadhcaldr04;


// default reproduced
Float_t t_lead_tkiso_recvtx_030_004_0015_02_01;
Float_t t_sublead_tkiso_recvtx_030_004_0015_02_01;
Float_t t_lead_tkiso_recvtx_030_000_0015_02_01;
Float_t t_sublead_tkiso_recvtx_030_000_0015_02_01;
Float_t t_lead_tkiso_recvtx_040_004_0015_02_01;
Float_t t_sublead_tkiso_recvtx_040_004_0015_02_01;
Float_t t_lead_tkiso_recvtx_040_000_0015_02_01;
Float_t t_sublead_tkiso_recvtx_040_000_0015_02_01;

//outer cone
Float_t t_lead_tkiso_recvtx_026_004_0000_02_01;
Float_t t_sublead_tkiso_recvtx_026_004_0000_02_01;
Float_t t_lead_tkiso_recvtx_028_004_0000_02_01;
Float_t t_sublead_tkiso_recvtx_028_004_0000_02_01;
Float_t t_lead_tkiso_recvtx_032_004_0000_02_01;
Float_t t_sublead_tkiso_recvtx_032_004_0000_02_01;
Float_t t_lead_tkiso_recvtx_034_004_0000_02_01;
Float_t t_sublead_tkiso_recvtx_034_004_0000_02_01;
Float_t t_lead_tkiso_recvtx_036_004_0000_02_01;
Float_t t_sublead_tkiso_recvtx_036_004_0000_02_01;
Float_t t_lead_tkiso_recvtx_026_004_0000_10_01;
Float_t t_sublead_tkiso_recvtx_026_004_0000_10_01;
Float_t t_lead_tkiso_recvtx_028_004_0000_10_01;
Float_t t_sublead_tkiso_recvtx_028_004_0000_10_01;
Float_t t_lead_tkiso_recvtx_032_004_0000_10_01;
Float_t t_sublead_tkiso_recvtx_032_004_0000_10_01;
Float_t t_lead_tkiso_recvtx_034_004_0000_10_01;
Float_t t_sublead_tkiso_recvtx_034_004_0000_10_01;
Float_t t_lead_tkiso_recvtx_036_004_0000_10_01;
Float_t t_sublead_tkiso_recvtx_036_004_0000_10_01;

//inner cone
Float_t t_lead_tkiso_recvtx_030_003_0000_02_01;
Float_t t_sublead_tkiso_recvtx_030_003_0000_02_01;
Float_t t_lead_tkiso_recvtx_030_005_0000_02_01;
Float_t t_sublead_tkiso_recvtx_030_005_0000_02_01;
Float_t t_lead_tkiso_recvtx_030_006_0000_02_01;
Float_t t_sublead_tkiso_recvtx_030_006_0000_02_01;
Float_t t_lead_tkiso_recvtx_030_005_0000_10_01;
Float_t t_sublead_tkiso_recvtx_030_005_0000_10_01;
Float_t t_lead_tkiso_recvtx_030_006_0000_10_01;
Float_t t_sublead_tkiso_recvtx_030_006_0000_10_01;


Float_t t_lead_tkiso_recvtx_030_000_0000_10_01;
Float_t t_lead_tkiso_recvtx_030_001_0000_10_01;
Float_t t_lead_tkiso_recvtx_030_002_0000_10_01;
Float_t t_lead_tkiso_recvtx_030_003_0000_10_01;
Float_t t_lead_tkiso_recvtx_035_000_0000_10_01;
Float_t t_lead_tkiso_recvtx_035_001_0000_10_01;
Float_t t_lead_tkiso_recvtx_035_002_0000_10_01;
Float_t t_lead_tkiso_recvtx_035_003_0000_10_01;
Float_t t_lead_tkiso_recvtx_040_000_0000_10_01;
Float_t t_lead_tkiso_recvtx_040_001_0000_10_01;
Float_t t_lead_tkiso_recvtx_040_002_0000_10_01;
Float_t t_lead_tkiso_recvtx_040_003_0000_10_01;
Float_t t_sublead_tkiso_recvtx_030_000_0000_10_01;
Float_t t_sublead_tkiso_recvtx_030_001_0000_10_01;
Float_t t_sublead_tkiso_recvtx_030_002_0000_10_01;
Float_t t_sublead_tkiso_recvtx_030_003_0000_10_01;
Float_t t_sublead_tkiso_recvtx_035_000_0000_10_01;
Float_t t_sublead_tkiso_recvtx_035_001_0000_10_01;
Float_t t_sublead_tkiso_recvtx_035_002_0000_10_01;
Float_t t_sublead_tkiso_recvtx_035_003_0000_10_01;
Float_t t_sublead_tkiso_recvtx_040_000_0000_10_01;
Float_t t_sublead_tkiso_recvtx_040_001_0000_10_01;
Float_t t_sublead_tkiso_recvtx_040_002_0000_10_01;
Float_t t_sublead_tkiso_recvtx_040_003_0000_10_01;

Float_t t_lead_tkiso_genvtx_030_000_0000_10_01;
Float_t t_lead_tkiso_genvtx_030_001_0000_10_01;
Float_t t_lead_tkiso_genvtx_030_002_0000_10_01;
Float_t t_lead_tkiso_genvtx_030_003_0000_10_01;
Float_t t_lead_tkiso_genvtx_035_000_0000_10_01;
Float_t t_lead_tkiso_genvtx_035_001_0000_10_01;
Float_t t_lead_tkiso_genvtx_035_002_0000_10_01;
Float_t t_lead_tkiso_genvtx_035_003_0000_10_01;
Float_t t_lead_tkiso_genvtx_040_000_0000_10_01;
Float_t t_lead_tkiso_genvtx_040_001_0000_10_01;
Float_t t_lead_tkiso_genvtx_040_002_0000_10_01;
Float_t t_lead_tkiso_genvtx_040_003_0000_10_01;
Float_t t_sublead_tkiso_genvtx_030_000_0000_10_01;
Float_t t_sublead_tkiso_genvtx_030_001_0000_10_01;
Float_t t_sublead_tkiso_genvtx_030_002_0000_10_01;
Float_t t_sublead_tkiso_genvtx_030_003_0000_10_01;
Float_t t_sublead_tkiso_genvtx_035_000_0000_10_01;
Float_t t_sublead_tkiso_genvtx_035_001_0000_10_01;
Float_t t_sublead_tkiso_genvtx_035_002_0000_10_01;
Float_t t_sublead_tkiso_genvtx_035_003_0000_10_01;
Float_t t_sublead_tkiso_genvtx_040_000_0000_10_01;
Float_t t_sublead_tkiso_genvtx_040_001_0000_10_01;
Float_t t_sublead_tkiso_genvtx_040_002_0000_10_01;
Float_t t_sublead_tkiso_genvtx_040_003_0000_10_01;

Float_t t_lead_tkiso_badvtx_030_000_0000_10_01;
Float_t t_lead_tkiso_badvtx_030_001_0000_10_01;
Float_t t_lead_tkiso_badvtx_030_002_0000_10_01;
Float_t t_lead_tkiso_badvtx_030_003_0000_10_01;
Float_t t_lead_tkiso_badvtx_035_000_0000_10_01;
Float_t t_lead_tkiso_badvtx_035_001_0000_10_01;
Float_t t_lead_tkiso_badvtx_035_002_0000_10_01;
Float_t t_lead_tkiso_badvtx_035_003_0000_10_01;
Float_t t_lead_tkiso_badvtx_040_000_0000_10_01;
Float_t t_lead_tkiso_badvtx_040_001_0000_10_01;
Float_t t_lead_tkiso_badvtx_040_002_0000_10_01;
Float_t t_lead_tkiso_badvtx_040_003_0000_10_01;
Float_t t_sublead_tkiso_badvtx_030_000_0000_10_01;
Float_t t_sublead_tkiso_badvtx_030_001_0000_10_01;
Float_t t_sublead_tkiso_badvtx_030_002_0000_10_01;
Float_t t_sublead_tkiso_badvtx_030_003_0000_10_01;
Float_t t_sublead_tkiso_badvtx_035_000_0000_10_01;
Float_t t_sublead_tkiso_badvtx_035_001_0000_10_01;
Float_t t_sublead_tkiso_badvtx_035_002_0000_10_01;
Float_t t_sublead_tkiso_badvtx_035_003_0000_10_01;
Float_t t_sublead_tkiso_badvtx_040_000_0000_10_01;
Float_t t_sublead_tkiso_badvtx_040_001_0000_10_01;
Float_t t_sublead_tkiso_badvtx_040_002_0000_10_01;
Float_t t_sublead_tkiso_badvtx_040_003_0000_10_01;

Float_t t_lead_tkiso_badvtx_rhocorr314_030_000_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_030_001_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_030_002_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_030_003_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_035_000_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_035_001_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_035_002_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_035_003_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_040_000_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_040_001_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_040_002_0000_10_01;
Float_t t_lead_tkiso_badvtx_rhocorr314_040_003_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_030_000_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_030_001_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_030_002_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_030_003_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_035_000_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_035_001_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_035_002_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_035_003_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_040_000_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_040_001_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_040_002_0000_10_01;
Float_t t_sublead_tkiso_badvtx_rhocorr314_040_003_0000_10_01;

//eta strip
Float_t t_lead_tkiso_recvtx_030_004_0010_02_01;
Float_t t_sublead_tkiso_recvtx_030_004_0010_02_01;
Float_t t_lead_tkiso_recvtx_030_004_0020_02_01;
Float_t t_sublead_tkiso_recvtx_030_004_0020_02_01;
Float_t t_lead_tkiso_recvtx_030_004_0010_10_01;
Float_t t_sublead_tkiso_recvtx_030_004_0010_10_01;
Float_t t_lead_tkiso_recvtx_030_004_0020_10_01;
Float_t t_sublead_tkiso_recvtx_030_004_0020_10_01;

// deltaz
Float_t t_lead_tkiso_recvtx_030_004_0000_06_01;
Float_t t_sublead_tkiso_recvtx_030_004_0000_06_01;
Float_t t_lead_tkiso_recvtx_030_004_0000_10_01;
Float_t t_sublead_tkiso_recvtx_030_004_0000_10_01;
Float_t t_lead_tkiso_recvtx_030_004_0000_15_01;
Float_t t_sublead_tkiso_recvtx_030_004_0000_15_01;
Float_t t_lead_tkiso_recvtx_030_004_0000_99_01;
Float_t t_sublead_tkiso_recvtx_030_004_0000_99_01;
Float_t t_lead_tkiso_recvtx_040_004_0000_06_01;
Float_t t_sublead_tkiso_recvtx_040_004_0000_06_01;
Float_t t_lead_tkiso_recvtx_040_004_0000_10_01;
Float_t t_sublead_tkiso_recvtx_040_004_0000_10_01;
Float_t t_lead_tkiso_recvtx_040_004_0000_99_01;
Float_t t_sublead_tkiso_recvtx_040_004_0000_99_01;

Float_t t_lead_tkiso_badvtx_030_004_0000_02_01;
Float_t t_sublead_tkiso_badvtx_030_004_0000_02_01;
Float_t t_lead_tkiso_badvtx_030_004_0000_06_01;
Float_t t_sublead_tkiso_badvtx_030_004_0000_06_01;
Float_t t_lead_tkiso_badvtx_030_004_0000_10_01;
Float_t t_sublead_tkiso_badvtx_030_004_0000_10_01;
Float_t t_lead_tkiso_badvtx_030_004_0000_15_01;
Float_t t_sublead_tkiso_badvtx_030_004_0000_15_01;
Float_t t_lead_tkiso_badvtx_030_004_0000_99_01;
Float_t t_sublead_tkiso_badvtx_030_004_0000_99_01;
Float_t t_lead_tkiso_badvtx_040_004_0000_02_01;
Float_t t_sublead_tkiso_badvtx_040_004_0000_02_01;
Float_t t_lead_tkiso_badvtx_040_004_0000_06_01;
Float_t t_sublead_tkiso_badvtx_040_004_0000_06_01;
Float_t t_lead_tkiso_badvtx_040_004_0000_10_01;
Float_t t_sublead_tkiso_badvtx_040_004_0000_10_01;
Float_t t_lead_tkiso_badvtx_040_004_0000_99_01;
Float_t t_sublead_tkiso_badvtx_040_004_0000_99_01;

Float_t t_lead_tkiso_recvtx_030_000_0000_06_01;
Float_t t_sublead_tkiso_recvtx_030_000_0000_06_01;
Float_t t_lead_tkiso_recvtx_030_000_0000_99_01;
Float_t t_sublead_tkiso_recvtx_030_000_0000_99_01;
Float_t t_lead_tkiso_recvtx_040_000_0000_06_01;
Float_t t_sublead_tkiso_recvtx_040_000_0000_06_01;
Float_t t_lead_tkiso_recvtx_040_000_0000_99_01;
Float_t t_sublead_tkiso_recvtx_040_000_0000_99_01;

Float_t t_lead_tkiso_badvtx_030_000_0000_02_01;
Float_t t_sublead_tkiso_badvtx_030_000_0000_02_01;
Float_t t_lead_tkiso_badvtx_030_000_0000_06_01;
Float_t t_sublead_tkiso_badvtx_030_000_0000_06_01;
Float_t t_lead_tkiso_badvtx_030_000_0000_99_01;
Float_t t_sublead_tkiso_badvtx_030_000_0000_99_01;
Float_t t_lead_tkiso_badvtx_040_000_0000_02_01;
Float_t t_sublead_tkiso_badvtx_040_000_0000_02_01;
Float_t t_lead_tkiso_badvtx_040_000_0000_06_01;
Float_t t_sublead_tkiso_badvtx_040_000_0000_06_01;
Float_t t_lead_tkiso_badvtx_040_000_0000_99_01;
Float_t t_sublead_tkiso_badvtx_040_000_0000_99_01;

// dxy
Float_t t_lead_tkiso_recvtx_030_004_0000_10_02;
Float_t t_sublead_tkiso_recvtx_030_004_0000_10_02;

Float_t t_genvtxz;
Float_t t_badvtxz;
Float_t t_recvtxz;
Float_t t_genvtx_sumpt;
Float_t t_recvtx_sumpt;
Float_t t_badvtx_sumpt;

// pfiso
Float_t t_lead_pfiso_charged03;
Float_t t_lead_pfiso_photon03;
Float_t t_lead_pfiso_neutral03;
Float_t t_lead_pfiso_charged04;
Float_t t_lead_pfiso_photon04;
Float_t t_lead_pfiso_neutral04;
Float_t t_lead_pfiso_charged_badvtx_03;
Float_t t_lead_pfiso_photon_badvtx_03;
Float_t t_lead_pfiso_neutral_badvtx_03;
Float_t t_lead_pfiso_charged_badvtx_04;
Float_t t_lead_pfiso_photon_badvtx_04;
Float_t t_lead_pfiso_neutral_badvtx_04;
Float_t t_sublead_pfiso_charged03;
Float_t t_sublead_pfiso_photon03;
Float_t t_sublead_pfiso_neutral03;
Float_t t_sublead_pfiso_charged04;
Float_t t_sublead_pfiso_photon04;
Float_t t_sublead_pfiso_neutral04;
Float_t t_sublead_pfiso_charged_badvtx_03;
Float_t t_sublead_pfiso_photon_badvtx_03;
Float_t t_sublead_pfiso_neutral_badvtx_03;
Float_t t_sublead_pfiso_charged_badvtx_04;
Float_t t_sublead_pfiso_photon_badvtx_04;
Float_t t_sublead_pfiso_neutral_badvtx_04;

//delta-r to track
Float_t t_lead_drtotk_10_06;
Float_t t_sublead_drtotk_10_06;
Float_t t_lead_drtotk_15_06;
Float_t t_sublead_drtotk_15_06;
Float_t t_lead_drtotk_20_06;
Float_t t_sublead_drtotk_20_06;
Float_t t_lead_drtotk_25_06;
Float_t t_sublead_drtotk_25_06;
Float_t t_lead_drtotkworstvtx_10_06;
Float_t t_sublead_drtotkworstvtx_10_06;
Float_t t_lead_drtotkworstvtx_15_06;
Float_t t_sublead_drtotkworstvtx_15_06;
Float_t t_lead_drtotkworstvtx_20_06;
Float_t t_sublead_drtotkworstvtx_20_06;
Float_t t_lead_drtotkworstvtx_25_06;
Float_t t_sublead_drtotkworstvtx_25_06;

Float_t t_lead_drtotk_10_10;
Float_t t_sublead_drtotk_10_10;
Float_t t_lead_drtotk_15_10;
Float_t t_sublead_drtotk_15_10;
Float_t t_lead_drtotk_20_10;
Float_t t_sublead_drtotk_20_10;
Float_t t_lead_drtotk_25_10;
Float_t t_sublead_drtotk_25_10;
Float_t t_lead_drtotkworstvtx_10_10;
Float_t t_sublead_drtotkworstvtx_10_10;
Float_t t_lead_drtotkworstvtx_15_10;
Float_t t_sublead_drtotkworstvtx_15_10;
Float_t t_lead_drtotkworstvtx_20_10;
Float_t t_sublead_drtotkworstvtx_20_10;
Float_t t_lead_drtotkworstvtx_25_10;
Float_t t_sublead_drtotkworstvtx_25_10;

Float_t t_lead_drtotk_10_99;
Float_t t_sublead_drtotk_10_99;
Float_t t_lead_drtotk_15_99;
Float_t t_sublead_drtotk_15_99;
Float_t t_lead_drtotk_20_99;
Float_t t_sublead_drtotk_20_99;
Float_t t_lead_drtotk_25_99;
Float_t t_sublead_drtotk_25_99;
Float_t t_lead_drtotkworstvtx_10_99;
Float_t t_sublead_drtotkworstvtx_10_99;
Float_t t_lead_drtotkworstvtx_15_99;
Float_t t_sublead_drtotkworstvtx_15_99;
Float_t t_lead_drtotkworstvtx_20_99;
Float_t t_sublead_drtotkworstvtx_20_99;
Float_t t_lead_drtotkworstvtx_25_99;
Float_t t_sublead_drtotkworstvtx_25_99;

	//END JIM
	// Jim's Functions BEGIN
void setPhotonCuts4();
void setPhotonCuts6();
void setPhotonCuts6pf();
void setDiphoCuts();
int diphoCutLevel(int leadind, int subleadind, int vtxind);
float diphoWeight(int c4, int dipholev);
int photonCutLevel4(float r9, float eta, float et, float isosum, float isosumbad, float tkiso, float sieie, float hoe, float drtotk);
int photonCutLevel6(float r9, float eta, float et, float isosum, float isosumbad, float tkiso, float sieie, float hoe, float drtotk);
int photonCutLevel6pf(float r9, float eta, float et, float isosum, float isosumbad, float chargediso, float neutraliso, float sieie, float hoe, float drtotk);

int diphoSubCategory(int, int, int);
float diphoSubCategory(int, int);
float              cutdiphoptom[20][4];
float                   cutdmom[20][4];
float       cutfsubleadcutindex[20][4];
float          cutfleadcutindex[20][4];
float                 cutetamax[20][4];
float                 cutetamin[20][4];
float         cutsubleadptomass[20][4];
float            cutleadptomass[20][4];
float                cutsumptom[20][4];

float       cutsubleadisosumoet6c[20][6];
float    cutsubleadisosumoetbad6c[20][6];
float     cutsubleadtrkisooetom6c[20][6];
float           cutsubleadsieie6c[20][6];
float          cutsubleadhovere6c[20][6];
float              cutsubleadr96c[20][6];
float   cutsublead_drtotk_25_996c[20][6];
bool photonCutsSet6;

float       cutpfisosumoet6c[20][6];
float    cutpfisosumoetbad6c[20][6];
float        cutpfchisooet6c[20][6];
float    cutpfhcalisooetom6c[20][6];
float           cutpfsieie6c[20][6];
float          cutpfhovere6c[20][6];
float              cutpfr96c[20][6];
float   cutpf_drtotk_25_996c[20][6];
bool photonCutsSet6pf;

float       cutsubleadisosumoet[20][4];
float    cutsubleadisosumoetbad[20][4];
float     cutsubleadtrkisooetom[20][4];
float           cutsubleadsieie[20][4];
float          cutsubleadhovere[20][4];
float              cutsubleadr9[20][4];
float   cutsublead_drtotk_25_99[20][4];
bool photonCutsSet4;

                 // Jim's Functions END

enum phoCiCIDLevel { phoNOCUTS, phoLOOSE, phoMEDIUM, phoTIGHT, phoSUPERTIGHT, phoHYPERTIGHT1, phoHYPERTIGHT2, phoHYPERTIGHT3, phoHYPERTIGHT4, phoNCUTLEVELS };

float Ci6C_cut_lead_isosumoet[phoNCUTLEVELS][6];
float Ci6C_cut_lead_isosumoetbad[phoNCUTLEVELS][6];
float Ci6C_cut_lead_trkisooet[phoNCUTLEVELS][6];
float Ci6C_cut_lead_sieie[phoNCUTLEVELS][6];
float Ci6C_cut_lead_hovere[phoNCUTLEVELS][6];
float Ci6C_cut_lead_r9[phoNCUTLEVELS][6];
float Ci6C_cut_lead_drtotk_25_99[phoNCUTLEVELS][6];
float Ci6C_cut_lead_pixel[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_isosumoet[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_isosumoetbad[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_trkisooet[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_sieie[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_hovere[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_r9[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_drtotk_25_99[phoNCUTLEVELS][6];
float Ci6C_cut_sublead_pixel[phoNCUTLEVELS][6];

float Ci4C_cut_lead_isosumoet[phoNCUTLEVELS][4];
float Ci4C_cut_lead_isosumoetbad[phoNCUTLEVELS][4];
float Ci4C_cut_lead_trkisooet[phoNCUTLEVELS][4];
float Ci4C_cut_lead_sieie[phoNCUTLEVELS][4];
float Ci4C_cut_lead_hovere[phoNCUTLEVELS][4];
float Ci4C_cut_lead_r9[phoNCUTLEVELS][4];
float Ci4C_cut_lead_drtotk_25_99[phoNCUTLEVELS][4];
float Ci4C_cut_lead_pixel[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_isosumoet[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_isosumoetbad[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_trkisooet[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_sieie[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_hovere[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_r9[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_drtotk_25_99[phoNCUTLEVELS][4];
float Ci4C_cut_sublead_pixel[phoNCUTLEVELS][4];


Int_t allevtvtxcount[21];
Int_t goodevtvtxcount1[21];
Int_t goodevtvtxcount14[21];
Int_t goodevtvtxcount17[21];
Int_t goodevtvtxcount2[21];

//WAS IN GENERAL FUNCTIONS

TMVA::Reader* tmvaReader, *tmvaReader1, *tmvaReader2, *tmvaReader3, *tmvaReader_dipho, *tmvaReader_dipho2;

Float_t tmva_sieie, tmva_goodpf_iso, tmva_badpf_iso, tmva_drtotk, tmva_ptom;
Float_t tmva_hoe, tmva_tkisopf, tmva_r9, tmva_pt, tmva_eta, tmva_isLeading;
Int_t tmva_cat;
Float_t tmva_etamax, tmva_subleadptomass, tmva_diphoptom, tmva_sumptom;
 Float_t tmva_subleadmva, tmva_leadmva, tmva_leadeta, tmva_subleadeta, tmva_leadr9, tmva_subleadr9, tmva_dmom;

void SetBDT();
Float_t BDT(Int_t jentry, Int_t iPhoton, Int_t vtx);
Float_t BDT_categorized(Int_t jentry, Int_t iPhoton, Int_t vtx, float isLead);
Float_t BDT_ptom(Int_t jentry, Int_t iPhoton, Int_t vtx, float mass);
Float_t BDT_ptom2(Int_t jentry, Int_t iPhoton, Int_t vtx, float mass);
 Float_t BDT_dipho(Int_t jentry, Int_t ilindex, Int_t islindex, float lmva, float slmva, float diphopt, float mass, float dmom);
Float_t BDT_dipho2(Int_t jentry, Int_t ilindex, Int_t islindex, float lmva, float slmva, float diphopt, float mass);

 void HggBookOptree();
 void SetOutputNtupleVariables(int jentry, int itype, int leadind, int subleadind, int vtxind, float mass, TLorentzVector *leadp4, TLorentzVector *subleadp4, float evweight, float pileupWeight, TLorentzVector * jet1=0, TLorentzVector * jet2=0, TLorentzVector * jet3=0);


 //float GetShiftValue(int run, float thiseta, float r9, int epoch=0);

float GetSmearSigma(float thiseta, float r9, int epoch=0);

std::vector<Float_t> arrayToVector(size_t length, const Float_t *data);

//r9cat+nr9*etaCat
std::vector<float> smear_nov14;
std::vector<float> smearErr_nov14;

// runCat+nRunCat*r9cat +nRunCat*nr9Cat*etaCat
//std::vector<float> shift_nov14;
//std::vector<float> shiftErr_nov14;

 Int_t GenIndexHgg(TLorentzVector *p_p4);

 Int_t itype_jim(int itype);

 
std::pair<int,int> DiphotonCiCSelectionIndicesJim( phoCiCIDLevel LEADCUTLEVEL = phoLOOSE, phoCiCIDLevel SUBLEADCUTLEVEL = phoLOOSE, Float_t leadPtMin = 30, Float_t subleadPtMin = 20, int ncat=6, int vtxind=-1, bool applyPtoverM=false);
// for a photon index, applies all levels of cuts and returns the index to the highest cut level passed (can do lead and sublead - same for now)
int   PhotonCiCSelectionLevelJim( int photon_index, int ncat=6, int vtxind=-1, int doSublead=1, int diphoind=-1, int print=0);
int   PhotonCiCpfSelectionLevelJim( int photon_index, int ncat=6, int vtxind=-1, int doSublead=1, int diphoind=-1, int print=0);

};

#endif
