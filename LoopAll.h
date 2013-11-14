#ifndef LoopAll_h
#define LoopAll_h

#include "CommonParameters.h"

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMacro.h"
#include "TClonesArray.h"

#include <sstream>
#include <set>
#include <string>
#include <list>
#include <map>

#include <THStack.h>
#include <TLegend.h>
#include <TAxis.h>

#include "TMVA/Reader.h"
#include "TStopwatch.h"

class BaseAnalysis;

#include "HistoContainer.h"
#include "CounterContainer.h"
#include "SampleContainer.h"
#include "TreeContainer.h"
#include "Cut.h"
#include "branchdef/Limits.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/HggVertexFromConversions.h"
#include "VertexAnalysis/interface/VertexAlgoParameters.h"
#include "VertexAnalysis/interface/PhotonInfo.h"
#include "VertexAnalysis/interface/VertexAlgoParameters.h"
#include "Macros/Normalization_8TeV.h"
#include "RooFuncReader.h"

#define BRANCH_DICT(NAME) branchDict[# NAME] = branch_info_t(& b ## _ ## NAME, & LoopAll::SetBranchAddress ## _ ## NAME, & LoopAll::Branch ## _ ## NAME )

#define DEBUG 1

class LoopAll {
 public :
  TTree          *fChain;

#include "branchdef/branchdef.h"
#include "branchdef/treedef.h"
#include "branchdef/setbranchaddress.h"
#include "branchdef/treebranch.h"

  std::vector<HistoContainer> histoContainer;
  std::vector<CounterContainer> counterContainer;
  std::vector<SampleContainer> sampleContainer;
  std::vector<Cut> cutContainer;
  //std::vector<TreeContainer> treeContainer;	 
  std::map<std::string, std::vector<TreeContainer> > treeContainer;	 
  
  RooContainer *rooContainer;
  int sqrtS;
  Normalization_8TeV * normalizer();
  
  LoopAll(TTree *tree=0);
  virtual ~LoopAll();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void   Init(Int_t typerunpass, TTree *tree);
  virtual void   InitReal(Int_t typerunpass);
  virtual void   TermReal(Int_t typerunpass);
  //virtual void   Loop(Double_t a);
  virtual Bool_t Notify();
  virtual void   Show(Long64_t entry = -1);
  virtual void   InitHistos();

  void LoopAndFillHistos(TString treename="event");
  void MergeContainers();
  //void WriteHist();  
  //void WriteFits();  
  //void WriteCounters();  

  /** @param type is the type of the analysis:
             0 for looper.py (filling histograms)
             1 for reduce.py (reduction step)
      @param n is the name of the output file 
  */
  void SetTypeRun(int type, const char* n);

  //void SetOutputNames(const char* n, const char* n2="");
  void StoreProcessedLumis(TTree * tree);
  void AddFile(std::string,int);
  void ReadInput(int t=0);
  
  void SetSubJob(bool);

  /** adds a new entry to sampleContainer and returns a reference
      to this entry. */
  SampleContainer & DefineSamples(const char *filesshortnam,
				  int type, int histtoindfromfiles, int histoplotit,
				  int nred, long long ntot, float intlumi,
				  float lumi, float xsec, float kfactor,
				  float scale, bool ignoreEvWeight=false,
				  int forceVersion=0, 
				  bool addnevents=false, TString pileup="");

  void Term(); 

  void checkDuty(int n, float thr, float start) { checkBench=n; benchThr=thr; benchStart=start; }
  
  int checkBench;
  TStopwatch stopWatch;
  float benchThr, benchStart;
  bool is_subjob;

  std::vector<TMacro*> configFiles;
  std::vector<std::string> files;
  std::vector<int> itype;
  //int lumireal[MAXFILES];
  int nfiles;
  float intlumi;
  float intlumi_;

  // Zee Validaton Flag
  bool runZeeValidation;
  bool applyEcalIsoPresel;
  bool makeDummyTrees;
  Float_t * pho_r9_cic;
  std::string cicVersion;
  bool usePFCiC;
  
  std::vector<TTree*> Trees;
  std::vector<TTree*> LumiTrees;
  std::vector<TFile*> Files;
  std::vector<TTree*> TreesPar;


  TTree * outputTree;
  TTree * outputTreeLumi;
  TTree * outputTreePar;
  TTree * plotvartree;
  TTree * inputfiletree;
  
  std::vector<TH1*> globalHistos;
  void AddGlobalHisto(TH1 * x ) { x->SetDirectory(0); globalHistos.push_back(x); }
  TH1D  * pileup;


  TFile * outputFile;
  TString outputFileName;
  TString histFileName;
  std::string outputTextFileName;
  Int_t makeOutputTree;

  char inputFilesName[1024];
  char countersName[1024];

  enum runtypes { kFill=0, kReduce=1, kFillReduce=2 };
  Int_t typerun;
  Int_t typeread;

  std::map<int,int> type2HistVal;

  Int_t        current; //current file
  Int_t        current_sample_index; //current file
  // global parameters
  Int_t tot_events, sel_events, type, version, reductions;
  bool createCS_;


  float diPhotonBDTOutput;

  std::vector<std::string>  globalCountersNames; 
  std::vector<int> globalCounters; 
  std::vector<int> fileGlobalCounters;
  
  Int_t        outputParTot_Events, outputParSel_Events, outputParType, outputParVersion, outputParReductions, outputParRed_Events[20];
  std::vector<std::string>* outputParParameters;
  std::string* outputParJobMaker;

  Int_t currentindexfiles;
  
  std::vector<Float_t> counters;
  std::vector<Float_t> countersred;
   
  TFile *hfile;
  Int_t outputEvents;

  void Loop(Int_t);
  void WriteHist();
  void WriteFits();
  void WriteCounters();
  
  int FillAndReduce(int);
  
  //void BookHistos();
  void AddCut(char*,int,int,int,float*,float*);
  void InitCounters();

  void AddCounter(int, const char*, const char*, const char*, const char*);

  int ApplyCut(int, float, int);
  int ApplyCut(std::string, float, int);
  void FillCutPlots(int icat, int cutset, std::string postfix, float histweight, float countweight);
  
  void myPrintCounters();
  void myPrintCountersNew();

  virtual void   BookHisto(int, int, int, int, int, int,
                           float, float, float, float, const char*,
			   const char* xaxis="", const char* yaxis="");

  
  // Cut down (flat) trees for MVA Training 
  void InitTrees(std::string);
  void BookTreeBranch(std::string name, int type, std::string dirName="");
  template <class T> void BookExternalTreeBranch(const char * name, T* addr, std::string dirName) {
	  for(unsigned int ind=0; ind<treeContainer[dirName].size(); ind++) {
		  treeContainer[dirName][ind].AddExternalBranch<T>(name,addr);
	  }
  }

  template <class T> void BookExternalTreeBranch(const char * name, T* addr, const char * type, std::string dirName) {
	  for(unsigned int ind=0; ind<treeContainer[dirName].size(); ind++) {
		  treeContainer[dirName][ind].AddExternalBranch<T>(name,addr,type);
	  }
  }

  template <class T> void BookExternalTreeBranch(const char * name, T* addr, int bufsize, int splitlevel, std::string dirName) {
	  for(unsigned int ind=0; ind<treeContainer[dirName].size(); ind++) {
		  treeContainer[dirName][ind].AddExternalBranch<T>(name,addr,bufsize,splitlevel);
	  }
  }
    

  void FillTreeContainer(std::string dir="");

  void FillTree(std::string name, float x, std::string dirName="");
  void FillTree(std::string name, double x, std::string dirName="");
  void FillTree(std::string name, int x, std::string dirName="");
  void FillTree(std::string name, unsigned int x, std::string dirName="");
  void FillTree(std::string name, std::string x, std::string dirName="");
  void FillTree(std::string name, bool x, std::string dirName="");
 
  void WritePI();
  void AddCut2(char*, int, int, int, float*, float*, int, int, int, float, float, char*, char*);

  int ApplyCut(int icut, int icat);
  int ApplyCut(TString cutname, int icat);

  float GetCutValue(TString cutname, int icat, int highcut=0);

  //int ApplyCut(int icut, float var, int icat);
  //int ApplyCut(TString cutname, float var, int icat);

  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCut(int icut, int * passcategory); //returns the number of categories
  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCut(TString cutname, int * passcategory); //returns the number of categories
  
  int ApplyCutsFill(int icat, int cutset, int & ncutsapplied, int & ncutspassed,  int & ncutsfailed, float histweight=1.0, float countweight=1.0);
  int ApplyCuts(int icat, int cutset, int & ncutsapplied, int & ncutspassed,  int & ncutsfailed);
  int ApplyCutsFill(int icat, int cutset, float histweight=1.0, float countweight=1.0);
  int ApplyCuts(int icat, int cutset);
  
  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCuts(int cutset, int * passcategory, int * ncutsapplied, int * ncutspassed,  int * ncutsfailed); //returns the number of categories
  //DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
  int ApplyCuts(int cutset, int * passcategory); //returns the number of categories
  
  int SetCutVariables(int i, float * variables);
  int SetCutVariables(TString cutname, float * variables);  

  void FillHist(std::string, float);
  void FillHist2D(std::string, float, float);

  void FillHist(std::string, int, float, float wt = 1.0); 
  void FillHist2D(std::string, int, float, float, float wt = 1.0);

  void FillCounter(std::string name, float weight=1., int cat=0);

  BaseAnalysis* AddAnalysis(BaseAnalysis*); 
  template<class T> T * GetAnalysis( const std::string & name ) {
	  std::vector<BaseAnalysis*>::iterator it=find( analyses.begin(), analyses.end(), name);
	  if( it != analyses.end() ) { return dynamic_cast<T*>( *it ); }
	  return 0;
  }
  
  /// void SkimBranch(const std::string & name)   { skimBranchNames.insert(name);  };
  void InputBranch(const std::string & name, int typ)  { inputBranchNames.insert(std::pair<std::string,int>(name,typ)); };
  void OutputBranch(const std::string & name) { if( find(outputBranchNames.begin(), outputBranchNames.end(), name)==outputBranchNames.end() ) { outputBranchNames.push_back(name); } };

  void GetBranches(std::map<std::string,int> & names, std::set<TBranch *> & branches);
  void SetBranchAddresses(std::map<std::string,int> & names);
  void Branches(std::list<std::string> & names);
  void GetEntry(std::set<TBranch *> & branches, int jentry);

  bool CheckLumiSelection( int run, int lumi );
  bool CheckEventList( int run, int lumi, int event );

  void StoreConfigFile(std::string);

#ifndef __CINT__
  typedef void (LoopAll::*branch_io_t) (TTree *);
  struct branch_info_t {
	  branch_info_t(TBranch ** b=0, branch_io_t r=0, branch_io_t w=0 ) :
		  branch(b), read(r), write(w) 
		  {};
	  TBranch ** branch;
	  branch_io_t read, write;
  };
  typedef std::map<std::string, branch_info_t> dict_t;
  dict_t branchDict;
#endif

  //// std::set<std::string> skimBranchNames;
  //// std::set<TBranch *> skimBranches; 
  std::map<std::string,int> inputBranchNames;
  std::set<TBranch *> inputBranches; 
  std::list<std::string> outputBranchNames;
  
  /** list of the analyses to be performed */
  std::vector<BaseAnalysis*> analyses;

  float pfTkIsoWithVertex(int phoindex, int vtxInd, float dRmax, float dRvetoBarrel, float dRvetoEndcap, float ptMin, float dzMax, float dxyMax, int pfToUse=1);
  float pfEcalIso(int phoindex, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStripBarrel, float etaStripEndcap, 
		  float thrBarrel, float thrEndcaps, int pfToUse=4);

  RooFuncReader *funcReader_dipho_MIT;
  TMVA::Reader *tmvaReaderID_UCSD, * tmvaReader_dipho_UCSD;
  TMVA::Reader *tmvaReaderID_MIT_Barrel, *tmvaReaderID_MIT_Endcap;
  TMVA::Reader *tmvaReader_dipho_MIT;
  TMVA::Reader *tmvaReaderID_Single_Barrel, *tmvaReaderID_Single_Endcap;
  TMVA::Reader *tmvaReaderID_2013_Barrel, *tmvaReaderID_2013_Endcap;
  TMVA::Reader *tmvaReaderID_2013_7TeV_MIT_Barrel, *tmvaReaderID_2013_7TeV_MIT_Endcap;

  Float_t photonIDMVA(Int_t, Int_t, TLorentzVector &, const char*);

  Float_t tmva_photonid_pfchargedisogood03;
  Float_t tmva_photonid_pfchargedisobad03;
  Float_t tmva_photonid_pfphotoniso03;
  Float_t tmva_photonid_pfneutraliso03;
  Float_t tmva_photonid_sieie;
  Float_t tmva_photonid_sieip;
  Float_t tmva_photonid_etawidth;
  Float_t tmva_photonid_phiwidth;
  Float_t tmva_photonid_r9;
  Float_t tmva_photonid_scrawe;
  Float_t tmva_photonid_s4ratio;
  Float_t tmva_photonid_lambdaratio;
  Float_t tmva_photonid_sceta;
  Float_t tmva_photonid_eventrho;
  Float_t tmva_photonid_ESEffSigmaRR;

  Float_t tmva_id_ucsd_sieie;
  Float_t tmva_id_ucsd_goodpf_iso;
  Float_t tmva_id_ucsd_badpf_iso;
  Float_t tmva_id_ucsd_drtotk;
  Float_t tmva_id_ucsd_hoe;
  Float_t tmva_id_ucsd_tkisopf;
  Float_t tmva_id_ucsd_r9;
  Float_t tmva_id_ucsd_ptom;
  Float_t tmva_id_ucsd_eta;
  Float_t tmva_id_ucsd_isLeading;
  
  Float_t tmva_dipho_UCSD_subleadptomass;
  Float_t tmva_dipho_UCSD_diphoptom;
  Float_t tmva_dipho_UCSD_sumptom;
  Float_t tmva_dipho_UCSD_subleadmva;
  Float_t tmva_dipho_UCSD_leadmva;
  Float_t tmva_dipho_UCSD_leadeta;
  Float_t tmva_dipho_UCSD_subleadeta;
  Float_t tmva_dipho_UCSD_leadr9;
  Float_t tmva_dipho_UCSD_subleadr9;
  Float_t tmva_dipho_UCSD_dmom;
  Float_t tmva_dipho_UCSD_diphocat2r92eta;

  Float_t tmva_id_mit_hoe;
  Float_t tmva_id_mit_sieie;
  Float_t tmva_id_mit_tiso1;
  Float_t tmva_id_mit_tiso3;
  Float_t tmva_id_mit_tiso2;
  Float_t tmva_id_mit_r9;
  Float_t tmva_id_mit_ecal;
  Float_t tmva_id_mit_hcal;
  Float_t tmva_id_mit_e5x5;
  Float_t tmva_id_mit_etawidth;
  Float_t tmva_id_mit_phiwidth;
  Float_t tmva_id_mit_sieip;
  Float_t tmva_id_mit_sipip;
  Float_t tmva_id_mit_nvtx;
  Float_t tmva_id_mit_sceta;
  Float_t tmva_id_mit_preshower;
 
  std::vector<Float_t> tmva_dipho_MIT_buf;
  std::map<int,std::vector<Float_t> > tmva_dipho_MIT_cache;
  Float_t *tmva_dipho_MIT_dmom;
  Float_t *tmva_dipho_MIT_dmom_wrong_vtx;
  Float_t *tmva_dipho_MIT_vtxprob;
  Float_t *tmva_dipho_MIT_ptom1;
  Float_t *tmva_dipho_MIT_ptom2;
  Float_t *tmva_dipho_MIT_eta1;
  Float_t *tmva_dipho_MIT_eta2;
  Float_t *tmva_dipho_MIT_dphi;
  Float_t *tmva_dipho_MIT_ph1mva;
  Float_t *tmva_dipho_MIT_ph2mva;

void GlobeCtIsol(int, TLorentzVector*, float, float, float, Int_t&, Float_t&, Float_t&, Float_t&, Float_t&);
int GlobeMatchIsl(TLorentzVector*, Float_t&);


enum eIDLevel {UltraLoose, VeryLoose, Loose, Medium, Tight, SuperTight, HyperTight1, HyperTight2, HyperTight3, HyperTight4, Robust};
int ElectronClassification(int);
std::pair<bool, bool> ElectronId(int, eIDLevel); 
void eIDInfo(Int_t, Int_t&, Int_t&,Int_t eIDMaxLevel=10);
Float_t sipCalculator(int);

// Match the Photon with the merged collection of conversions
PhotonInfo fillPhotonInfos(int p1, int useAllConvs=2, float * energy=0);
int matchPhotonToConversion(int lpho, int useAllConvs=2);
double phiNorm (float &phi);
double etaTransformation(  float EtaParticle , float Zvertex);

vector<double> generate_flat10_weights(TH1D* data_npu_estimated);


//----------------------------------------
// Vertex analysis
//----------------------------------------
void vertexAnalysis(HggVertexAnalyzer & vtxAna,  PhotonInfo pho1, PhotonInfo pho2);
//std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, int p1, int p2, std::vector<std::string> & vtxVarNames);

/** @return the indices of the vertices ranked by some algorithm */
std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, PhotonInfo & pho1, PhotonInfo & pho2,
				 std::vector<std::string> & vtxVarNames, 					  
				 bool useMva=false, TMVA::Reader * reader=0, std::string tmvaMethod="");

//----------------------------------------
bool FindMCLeptons(int index, int& mc1, int& mc2, int& pho, int leptonType=11);
bool FindMCHiggsPhotons(int& higgsind, int& mc1, int& mc2, int& i1, int& i2  );
bool FindMCVBF(int higgsind, int& vbfq1, int& vbfq2 );
bool FindMCVH(int higgsind, int& vh, int& vh1, int& vh2 );
TLorentzVector GetHiggs() 
{
	TLorentzVector  gP4(0,0,0,0);
	assert(gh_higgs_p4 != 0 || gp_p4 != 0);
	if( gh_higgs_p4 != 0 && gh_higgs_p4->At(0) != 0 ) {
	    gP4 = *((TLorentzVector*)gh_higgs_p4->At(0));
	} else {
	    for (int gi=0;gi<gp_n;gi++){
		if (gp_pdgid[gi]==25){
		    gP4 = *((TLorentzVector*)gp_p4->At(gi));
		    break;
		}
	    }
	}
	return gP4; 
};

/** @return the photon four momentum calculated with respect to the given interaction vertex (??) */
TLorentzVector get_pho_p4(int ipho, int ivtx, const float *pho_energy_array=0) const ;
TLorentzVector get_pho_p4(int ipho, TVector3 * vtx, const float * energy=0) const ;
void set_pho_p4(int ipho, int ivtx, float *pho_energy_array=0);
double get_pho_zposfromconv(TVector3 convvtx, TVector3 superclustervtx, TVector3 beamSpot);
// end vertex analysis 

void FillCICPFInputs();
void FillCICInputs();
void FillCIC();

// CiC SELECTION CODE BEGIN - SSIMON

// defines photon CiC ID cuts for all cut levels
enum phoCiCIDLevel { phoNOCUTS=0, phoLOOSE, phoMEDIUM, phoTIGHT, phoSUPERTIGHT, phoHYPERTIGHT1, phoHYPERTIGHT2, phoHYPERTIGHT3, phoHYPERTIGHT4, phoHYPERTIGHT5, phoHYPERTIGHT6, phoHYPERTIGHT7, phoHYPERTIGHT8, phoHYPERTIGHT9, phoNCUTLEVELS };
enum phoCiCCuts { phoISOSUMOET=0,  phoISOSUMOETBAD,   phoTRKISOOETOM,   phoSIEIE,   phoHOVERE,   phoR9,   phoDRTOTK_25_99,   phoPIXEL, phoNCUTS };
enum phoCiC6Categories { phoCiC6EBhighR9=0, phoCiC6EBmidR9, phoCiC6EBlowR9, phoCiC6EEhighR9, phoCiC6EEmidR9, phoCiC6EElowR9, phoCiC6NCATEGORIES };
enum phoCiC4Categories { phoCiC4EBhighR9=0, phoCiC4EBlowR9, phoCiC4EEhighR9, phoCiC4EElowR9, phoCiC4NCATEGORIES };

/** @param cutlevel           (input) is the required level of the cuts (e.g. phoSUPERTIGHT) 
    @param cic6_cuts_lead     (output) will be filled with the CIC6 cut values for the LEADING photon
    @param cic6_cuts_sublead, (output) will be filled with the CIC6 cut values for the SUBLEADING photon
    @param cic4_cuts_lead     (output) will be filled with the CIC4 cut values for the LEADING photon
    @param cic4_cuts_sublead  (output) will be filled with the CIC4 cut values for the SUBLEADING photon
*/
 void SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cic6_cuts_lead, float * cic6_cuts_sublead, float * cic4_cuts_lead, float * cic4_cuts_sublead, float*, float*);

//----------------------------------------------------------------------
// cut values for cuts in categories selection
//
// first index is the level ('tightness') of the selection,
// second index is the category index.
//
// See LoopAll::SetPhotonCutsInCategories(..) in GeneralFunctions_cc.h
//   for the actual values.
// initialized in PhotonAnalysis::Init(..)
//
// The (highest) level of selection which a photon passes is
// calculated in LoopAll::PhotonCiCSelectionLevel(..)
//----------------------------------------------------------------------
// 6 categories
//   leading photon
float cic6_cut_lead_isosumoet[phoNCUTLEVELS][6];
float cic6_cut_lead_isosumoetbad[phoNCUTLEVELS][6];
float cic6_cut_lead_trkisooet[phoNCUTLEVELS][6];
float cic6_cut_lead_sieie[phoNCUTLEVELS][6];
float cic6_cut_lead_hovere[phoNCUTLEVELS][6];
float cic6_cut_lead_r9[phoNCUTLEVELS][6];
float cic6_cut_lead_drtotk_25_99[phoNCUTLEVELS][6];
float cic6_cut_lead_pixel[phoNCUTLEVELS][6];

//   subleading photon
float cic6_cut_sublead_isosumoet[phoNCUTLEVELS][6];
float cic6_cut_sublead_isosumoetbad[phoNCUTLEVELS][6];
float cic6_cut_sublead_trkisooet[phoNCUTLEVELS][6];
float cic6_cut_sublead_sieie[phoNCUTLEVELS][6];
float cic6_cut_sublead_hovere[phoNCUTLEVELS][6];
float cic6_cut_sublead_r9[phoNCUTLEVELS][6];
float cic6_cut_sublead_drtotk_25_99[phoNCUTLEVELS][6];
float cic6_cut_sublead_pixel[phoNCUTLEVELS][6];

// 4 categories
//   leading photon
float cic4_cut_lead_isosumoet[phoNCUTLEVELS][4];
float cic4_cut_lead_isosumoetbad[phoNCUTLEVELS][4];
float cic4_cut_lead_trkisooet[phoNCUTLEVELS][4];
float cic4_cut_lead_sieie[phoNCUTLEVELS][4];
float cic4_cut_lead_hovere[phoNCUTLEVELS][4];
float cic4_cut_lead_r9[phoNCUTLEVELS][4];
float cic4_cut_lead_drtotk_25_99[phoNCUTLEVELS][4];
float cic4_cut_lead_pixel[phoNCUTLEVELS][4];

//   subleading photon
float cic4_cut_sublead_isosumoet[phoNCUTLEVELS][4];
float cic4_cut_sublead_isosumoetbad[phoNCUTLEVELS][4];
float cic4_cut_sublead_trkisooet[phoNCUTLEVELS][4];
float cic4_cut_sublead_sieie[phoNCUTLEVELS][4];
float cic4_cut_sublead_hovere[phoNCUTLEVELS][4];
float cic4_cut_sublead_r9[phoNCUTLEVELS][4];
float cic4_cut_sublead_drtotk_25_99[phoNCUTLEVELS][4];
float cic4_cut_sublead_pixel[phoNCUTLEVELS][4];

float pfisoOffset;
//   leading photon
float cic4pf_cut_lead_isosumoet[phoNCUTLEVELS][4];
float cic4pf_cut_lead_isosumoetbad[phoNCUTLEVELS][4];
float cic4pf_cut_lead_trkisooet[phoNCUTLEVELS][4];
float cic4pf_cut_lead_sieie[phoNCUTLEVELS][4];
float cic4pf_cut_lead_hovere[phoNCUTLEVELS][4];
float cic4pf_cut_lead_r9[phoNCUTLEVELS][4];
float cic4pf_cut_lead_drtotk_25_99[phoNCUTLEVELS][4];
float cic4pf_cut_lead_pixel[phoNCUTLEVELS][4];

//   subleading photon
float cic4pf_cut_sublead_isosumoet[phoNCUTLEVELS][4];
float cic4pf_cut_sublead_isosumoetbad[phoNCUTLEVELS][4];
float cic4pf_cut_sublead_trkisooet[phoNCUTLEVELS][4];
float cic4pf_cut_sublead_sieie[phoNCUTLEVELS][4];
float cic4pf_cut_sublead_hovere[phoNCUTLEVELS][4];
float cic4pf_cut_sublead_r9[phoNCUTLEVELS][4];
float cic4pf_cut_sublead_drtotk_25_99[phoNCUTLEVELS][4];
float cic4pf_cut_sublead_pixel[phoNCUTLEVELS][4];

//----------------------------------------------------------------------

/** loops through photons and returns indices to two photons passing desired selection .
    if more than one diphoton passes, returns pair with highest lead photon pt, or if lead is same, with highest sublead pt. */
int DiphotonCiCSelection( phoCiCIDLevel LEADCUTLEVEL = phoLOOSE, 
                          phoCiCIDLevel SUBLEADCUTLEVEL = phoLOOSE, 
                          Float_t leadPtMin = 30, 
                          Float_t subleadPtMin = 20, 
                          int ncategories=6, 
                          bool applyPtoverM=false, 
                          float *pho_energy_array=0, 
                          bool split=false, int fixedvtx=-1, std::vector<bool> veto_indices=std::vector<bool>(false),
                          std::vector<int> cutsbycat=std::vector<int>(0));


int DiphotonMITPreSelection(const char * type, Float_t leadPtMin, Float_t subleadPtMin, Float_t phoidMvaCut, bool applyPtoverM, float *pho_energy_array=0, bool vetodipho=false, bool kinonly=false, float dipho_BDT_Cut=-100,int fixedvtx=-1, bool split=false, std::vector<bool> veto_indices=std::vector<bool>(false));
float DiphotonMITPreSelectionPerDipho(const char * type, int idipho, Float_t leadPtMin, Float_t subleadPtMin, Float_t phoidMvaCut, bool applyPtoverM, float *pho_energy_array=0, int fixedvtx=-1, bool split=false, bool kinonly=false, std::vector<bool> veto_indices=std::vector<bool>(false));
int DiphotonMITPreSelection2011(Float_t leadPtMin, Float_t subleadPtMin, Float_t phoidMvaCut, bool applyPtoverM, float *pho_energy_array=0, bool kinonly=false);

/** for a photon index, applies all levels of cuts and returns the
    index to the highest cut level passed (can do lead and sublead -
    same for now). This should be compared to one of the phoCiCIDLevel
    enum constants (e.g. phoSUPERTIGHT).

    @param ph_passcut will contain flags which cuts the given photon
           passes. The first index is the level of the selection
           looked at, the second index is the index of the cut. 
           Will be automatically resized.

    @param doSublead zero, if the cuts for the leading photon
           should be applied, non-zero if the cuts for the 
           subleading photon should be applied.

 */
int   PhotonCiCPFSelectionLevel( int photon_index, int vertex_index, std::vector<std::vector<bool> > & ph_passcut, int ncategories=6, int doSublead=1, float *pho_energy_array=0);

int   PhotonCiCSelectionLevel( int photon_index, int vertex_index, std::vector<std::vector<bool> > & ph_passcut, int ncategories=6, int doSublead=1, float *pho_energy_array=0);

bool   PhotonMITPreSelection( int photon_index, int vertex_index,float *pho_energy_array=0);
bool   PhotonMITPreSelection2011( int photon_index, int vertex_index,float *pho_energy_array=0);

// Functions to calculate variables used in CiC selection
Float_t DeltaRToTrack(Int_t photon_ind=-1, Int_t vtxind=-1, Float_t PtMin=1., Float_t dzmax=0.2, Float_t dxymax=0.1, int maxlosthits=0);
Float_t IsoEcalHitsSumEtNumCrystal( TVector3 *calopos, Float_t innerConeDR, Float_t outerConeDR, Float_t stripEtaHalfWidth, Float_t stripHalfLength=99.);
std::pair<Int_t, Float_t> WorstSumTrackPtInCone(int ipho, Int_t returnVtxIndex=0, Float_t PtMin=0, Float_t OuterConeRadius=0.3, Float_t InnerConeRadius=0.04, Float_t EtaStripHalfWidth=0.015, Float_t dzmax=0.2, Float_t dxymax=0.1, bool Zee_validation=false);
Float_t SumTrackPtInCone(TLorentzVector *photon_p4, Int_t vtxind, Float_t PtMin=0, Float_t OuterConeRadius=0.3, Float_t InnerConeRadius=0.04, Float_t EtaStripHalfWidth=0.015, Float_t dzmax=0.2, Float_t dxymax=0.1, bool Zee_validation=false, Int_t pho_ind=-1);


//----------------------------------------------------------------------
/** photon category functions (r9 and eta) */
int PhotonCategory(int photonindex, int n_r9cat=3, int n_etacat=2) { 


  // example: n_r9cat = 2 and n_etacat = 2
  //  -> return value 0 high R9, barrel
  //                  1 low R9,  barrel
  //                  2 high R9, endcap
  //                  3 low R9,  endcap
  return PhotonR9Category(photonindex,n_r9cat) + n_r9cat * PhotonEtaCategory(photonindex,n_etacat);
}

//----------------------------------------------------------------------

/** @return the photon R9 category number */
Int_t PhotonR9Category(int photonindex, int n_r9cat=3, float r9boundary=0.94) { 
  if(photonindex < 0) return -1;
  if(n_r9cat<2)return 0;
  int r9cat=0;
  float r9 = pho_r9_cic[photonindex];
  if(n_r9cat==3) {
    r9cat = (Int_t)(r9<0.94) + (r9<0.9);// 0, 1, or 2 (high r9 --> low r9)
  } else if(n_r9cat==2) {
    r9cat = (Int_t)(r9<r9boundary);// 0, 1(high r9 --> low r9)
  }
  return r9cat;
}

//----------------------------------------------------------------------

/** @return the photon eta category number */
int PhotonEtaCategory(int photonindex, int n_etacat=4) {
  if(photonindex < 0) return -1;
  if(n_etacat<2)return 0;
  int etacat;
  // Float_t eta = fabs(((TVector3*)pho_calopos->At(photonindex))->Eta());
  Float_t eta = fabs(((TVector3*)sc_xyz->At(pho_scind[photonindex]))->Eta());
  if(n_etacat==4) {
    etacat = (Int_t)(eta>0.9) + (Int_t)(eta>1.479) + (Int_t)(eta>2.1);   // 0, 1, 2, or 3 (barrel --> endcap)
  } else if(n_etacat==2) {
    etacat = (Int_t)(eta>1.479);   // 0, 1 (barrel --> endcap)
  }
  return  etacat;
}

//----------------------------------------------------------------------
//diphoton category functions ( r9, eta, and diphoton pt)
int DiphotonCategory(Int_t leadind, Int_t subleadind, float pTh, float pToMh,  int n_etacat=4,int n_r9cat=3, float r9boundary=0.94, int n_pThcat=0, int n_pToMhcat=0, int nVtxCategories=0, int nvtx=0, float vtxMva=-1.) {
  Int_t r9cat  =  TMath::Max(PhotonR9Category(leadind,n_r9cat,r9boundary),PhotonR9Category(subleadind,n_r9cat,r9boundary));
  Int_t etacat =  TMath::Max(PhotonEtaCategory(leadind,n_etacat),PhotonEtaCategory(subleadind,n_etacat));
  Int_t pThcat =  DiphotonPtCategory(pTh,n_pThcat);
  Int_t pToMhcat =  DiphotonPtOverMCategory(pToMh,n_pToMhcat);
  Int_t vtxCat =  DiphotonVtxCategory(nVtxCategories,nvtx);
  return  (r9cat + n_r9cat*etacat + (n_r9cat*n_etacat)*pThcat) + (n_r9cat*n_etacat*(n_pThcat>0?n_pThcat:1))*pToMhcat + (n_r9cat*n_etacat*(n_pThcat>0?n_pThcat:1)*(n_pToMhcat>0?n_pToMhcat:1))*vtxCat;  // (n_r9cat*c_etacat*n_pThcat) categories
}

//----------------------------------------------------------------------
int DijetSubCategory(float mjj, float leadPt, float subledPt, float ncat)
{
	if( ncat == 1 ) { return 0; }
	return ( mjj < 500. || subledPt < 30. );
}


//----------------------------------------------------------------------

int DiphotonVtxCategory(int nVtxCategories, int nvtx)
{
	int cat=0;
	if(nVtxCategories==3) {
		cat = (Int_t)(nvtx > 18) + (Int_t)(nvtx > 13);
	} else if (nVtxCategories==2) {
		cat = (Int_t)(nvtx > 15);
	}
	/// cout << "DiphotonVtxCategory " << cat << " " << vtxMva << " " << nVtxCategories << endl;
	return  cat;
}

int DiphotonVtxCategory(float vtxMva, int nVtxCategories)
{
	int cat=0;
	if(nVtxCategories==2) {
		cat = (Int_t)(vtxMva > -0.8);
	} else if (nVtxCategories>0) {
		cat = (Int_t)(vtxMva > -0.8) + (Int_t)(vtxMva > -0.55);
	}
	/// cout << "DiphotonVtxCategory " << cat << " " << vtxMva << " " << nVtxCategories << endl;
	return  cat;
}

int DiphotonPtCategory(double pTh, int n_pThcat=0) {
  if(n_pThcat<2)return 0;
  int pThcat=0;
  if(n_pThcat == 2) {
    pThcat = (Int_t)(pTh < 40.);
  } else if (n_pThcat == 3) {
    pThcat = (Int_t)((pTh < 70.) + (pTh < 40.));
  }
  return pThcat;
}

int DiphotonPtOverMCategory(double pToMh, int n_pToMhcat=0) {
  if(n_pToMhcat<2) return 0;
  int pToMhcat=0;
  if(n_pToMhcat == 2) {
    pToMhcat = (Int_t)(pToMh < 40./125.);
  } else if (n_pToMhcat == 3) {
    pToMhcat = (Int_t)((pToMh < 70./125.) + (pToMh < 40./125.));
  }
  return pToMhcat;
}
// CiC SELECTION CODE END - SSIMON

// Functions movec from Tools.h
double DeltaPhi(double,double);

///
/// Missing and addition branches
///

// Higgs and company 
Int_t gh_gen2reco1;
Int_t gh_gen2reco2;
Int_t gh_vbfq1_pdgid;
Int_t gh_vbfq2_pdgid;
Int_t gh_vh_pdgid;
Int_t gh_vh1_pdgid;
Int_t gh_vh2_pdgid;
TClonesArray *gh_higgs_p4;
TClonesArray *gh_pho1_p4;
TClonesArray *gh_pho2_p4;
TClonesArray *gh_vbfq1_p4;
TClonesArray *gh_vbfq2_p4;
TClonesArray *gh_vh1_p4;
TClonesArray *gh_vh2_p4;
//TClonesArray *METcorrected;  //met at analysis step

Float_t rho;
Int_t gv_n; 
TClonesArray * gv_pos; 
Int_t pu_n; 
std::vector<float> * pu_zpos;
std::vector<float> * pu_sumpt_lowpt;
std::vector<float> * pu_sumpt_highpt;
std::vector<int> * pu_ntrks_lowpt;
std::vector<int> * pu_ntrks_highpt;

bool pho_idmva_cached;
float pho_idmva[MAX_PHOTONS][MAX_VERTICES];

#define MAX_DIPHOTONS 50
Int_t dipho_n;
Int_t dipho_leadind[MAX_DIPHOTONS];
Int_t dipho_subleadind[MAX_DIPHOTONS];
Int_t dipho_vtxind[MAX_DIPHOTONS];
Float_t dipho_sumpt[MAX_DIPHOTONS];
Bool_t dipho_sel[MAX_DIPHOTONS];
Float_t dipho_BDT[MAX_DIPHOTONS];
Bool_t pho_genmatched[MAX_PHOTONS];
Float_t pho_regr_energy_otf[MAX_PHOTONS];
Float_t pho_regr_energyerr_otf[MAX_PHOTONS];
Float_t pho_ESEffSigmaRR[MAX_PHOTONS];
Float_t pho_s4ratio[MAX_PHOTONS];
//// Float_t dipho_leadet[MAX_DIPHOTONS];
//// Float_t dipho_subleadet[MAX_DIPHOTONS];
//// Float_t dipho_leadeta[MAX_DIPHOTONS];
//// Float_t dipho_subleadeta[MAX_DIPHOTONS];
//// Int_t dipho_leadci6cindex[MAX_DIPHOTONS];
//// Int_t dipho_subleadci6cindex[MAX_DIPHOTONS];
//// Int_t dipho_leadci4cindex[MAX_DIPHOTONS];
//// Int_t dipho_subleadci4cindex[MAX_DIPHOTONS];
//// Float_t dipho_mass[MAX_DIPHOTONS];
//// Float_t dipho_pt[MAX_DIPHOTONS];
//// Float_t dipho_eta[MAX_DIPHOTONS];
//// Float_t dipho_phi[MAX_DIPHOTONS];
//// Float_t dipho_cts[MAX_DIPHOTONS];


//correctMETinRED
Float_t shiftMET_pt;
Float_t shiftMET_phi;
Float_t smearMET_pt;
Float_t smearMET_phi;
Float_t shiftscaleMET_pt;
Float_t shiftscaleMET_phi;
Float_t shiftsmearMET_pt;
Float_t shiftsmearMET_phi;
Float_t correctedpfMET;
Float_t correctedpfMET_phi;

Float_t shiftMET_eta;
Float_t shiftMET_e;
Float_t shiftscaleMET_eta;
Float_t shiftscaleMET_e;



TBranch *b_gh_gen2reco1;
TBranch *b_gh_gen2reco2;
TBranch *b_gh_vbfq1_pdgid;
TBranch *b_gh_vbfq2_pdgid;
TBranch *b_gh_vh_pdgid;
TBranch *b_gh_vh1_pdgid;
TBranch *b_gh_vh2_pdgid;
//TBranch *b_METcorrected;  //met at analysis step
TBranch *b_gh_higgs_p4;
TBranch *b_gh_pho1_p4;
TBranch *b_gh_pho2_p4;
TBranch *b_gh_vbfq1_p4;
TBranch *b_gh_vbfq2_p4;
TBranch *b_gh_vh1_p4;
TBranch *b_gh_vh2_p4;
TBranch * b_dipho_n;
TBranch * b_dipho_leadind;
TBranch * b_dipho_subleadind;
TBranch * b_dipho_vtxind;
TBranch * b_dipho_sumpt;
TBranch * b_pho_genmatched;
TBranch * b_pho_regr_energy_otf;
TBranch * b_pho_regr_energyerr_otf;
//// TBranch * b_dipho_leadet;
//// TBranch * b_dipho_subleadet;
//// TBranch * b_dipho_leadeta;
//// TBranch * b_dipho_subleadeta;
//// TBranch * b_dipho_leadci6cindex;
//// TBranch * b_dipho_subleadci6cindex;
//// TBranch * b_dipho_leadci4cindex;
//// TBranch * b_dipho_subleadci4cindex;
//// TBranch * b_dipho_mass;
//// TBranch * b_dipho_pt;
//// TBranch * b_dipho_eta;
//// TBranch * b_dipho_phi;
//// TBranch * b_dipho_cts;


//correctMETinRED
TBranch * b_shiftMET_pt;
TBranch * b_shiftMET_phi;
TBranch * b_smearMET_pt;
TBranch * b_smearMET_phi;
TBranch * b_shiftscaleMET_pt;
TBranch * b_shiftscaleMET_phi;
TBranch * b_shiftsmearMET_pt;
TBranch * b_shiftsmearMET_phi;
TBranch * b_correctedpfMET;
TBranch * b_correctedpfMET_phi;

TBranch * b_shiftMET_eta;
TBranch * b_shiftMET_e;
TBranch * b_shiftscaleMET_eta;
TBranch * b_shiftscaleMET_e;

//----------------------------------------
// photon and diphoton vertex selection
//----------------------------------------
/** calculated in PhotonAnalysis.cc: vertex selected
    for the diphoton pair with highest sum of Pt */
int vtx_std_sel;

std::vector<int> *  dipho_vtx_std_sel;

/** calculated e.g. in PhotonAnalysis.cc */
std::vector<std::vector<int> > * vtx_std_ranked_list;
std::vector<float> * vtx_std_evt_mva;
// std::vector<int> * vtx_std_ranked_list;

//----------------------------------------

// CiC inputs
std::vector<std::vector<float> >* pho_mitmva;
std::vector<std::vector<float> >* pho_tkiso_recvtx_030_002_0000_10_01;
Float_t pho_tkiso_badvtx_040_002_0000_10_01[MAX_PHOTONS];
Float_t pho_pfiso_charged_badvtx_04[MAX_PHOTONS];
Int_t pho_pfiso_charged_badvtx_id[MAX_PHOTONS];
Int_t pho_tkiso_badvtx_id[MAX_PHOTONS];
std::vector<std::vector<float> >* pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01;
Float_t pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01[MAX_PHOTONS];
Int_t pho_ZeeVal_tkiso_badvtx_id[MAX_PHOTONS];
Float_t pho_drtotk_25_99[MAX_PHOTONS];

bool runCiC;

/** cut levels of CIC photon identification (stored in the tree).
    The first index is the index of the photon object (0..pho_n-1),
    the second index is the index of the vertex (0..vtx_std_n-1) 
    with respect to which the id is calculated. 

    The contents are the return values of the function
    PhotonCiCSelectionLevel(..)
*/

Int_t mu_glo_hasgsftrack[MAX_MUONS];
  
std::vector<std::vector<Short_t> >* pho_cic6cutlevel_lead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic6passcuts_lead;
std::vector<std::vector<Short_t> >* pho_cic6cutlevel_sublead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic6passcuts_sublead;
std::vector<std::vector<Short_t> >* pho_cic4cutlevel_lead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic4passcuts_lead;
std::vector<std::vector<Short_t> >* pho_cic4cutlevel_sublead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic4passcuts_sublead;

std::vector<std::vector<Short_t> >* pho_cic4pfcutlevel_lead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic4pfpasscuts_lead;
std::vector<std::vector<Short_t> >* pho_cic4pfcutlevel_sublead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic4pfpasscuts_sublead;

std::vector<std::vector<Short_t> >* pho_cutlevel_lead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_passcuts_lead;
std::vector<std::vector<Short_t> >* pho_cutlevel_sublead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_passcuts_sublead;


// Indices of conversions matching the photons
std::vector<int> * pho_matchingConv;
TBranch *b_pho_matchingConv;

TBranch *b_gv_n;
TBranch *b_gv_pos;
TBranch *b_pu_n;
TBranch *b_pu_zpos;
TBranch *b_pu_sumpt_lowpt;
TBranch *b_pu_sumpt_highpt;
TBranch *b_pu_ntrks_lowpt;
TBranch *b_pu_ntrks_highpt;

TBranch *b_vtx_std_sel;
TBranch *b_vtx_std_ranked_list;
TBranch *b_vtx_std_evt_mva;

TBranch * b_pho_mitmva;
TBranch * b_pho_pfiso_charged_badvtx_04;
TBranch * b_pho_pfiso_charged_badvtx_id;
TBranch * b_pho_tkiso_recvtx_030_002_0000_10_01;
TBranch * b_pho_tkiso_badvtx_040_002_0000_10_01;
TBranch * b_pho_tkiso_badvtx_id;
TBranch * b_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01;
TBranch * b_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01;
TBranch * b_pho_ZeeVal_tkiso_badvtx_id;
TBranch * b_pho_drtotk_25_99;

TBranch * b_mu_glo_hasgsftrack;
  
TBranch * b_pho_cic6cutlevel_lead;
TBranch * b_pho_cic6passcuts_lead;
TBranch * b_pho_cic6cutlevel_sublead;
TBranch * b_pho_cic6passcuts_sublead;
TBranch * b_pho_cic4cutlevel_lead;
TBranch * b_pho_cic4passcuts_lead;
TBranch * b_pho_cic4cutlevel_sublead;
TBranch * b_pho_cic4passcuts_sublead;

TBranch * b_pho_cic4pfcutlevel_lead;
TBranch * b_pho_cic4pfpasscuts_lead;
TBranch * b_pho_cic4pfcutlevel_sublead;
TBranch * b_pho_cic4pfpasscuts_sublead;

TBranch * b_pho_cutlevel_lead;
TBranch * b_pho_passcuts_lead;
TBranch * b_pho_cutlevel_sublead;
TBranch * b_pho_passcuts_sublead;

Float_t mc_et[MAX_ELECTRONS];
Float_t mc_eta[MAX_ELECTRONS];
Float_t mc_phi[MAX_ELECTRONS];
Float_t fsr_et[MAX_ELECTRONS];
Float_t fsr_eta[MAX_ELECTRONS];
Float_t fsr_phi[MAX_ELECTRONS];
TLorentzVector* higgs;
TBranch* b_mc_et,*b_mc_eta, *b_mc_phi, *b_higgs, *b_fsr_et, *b_fsr_eta, *b_fsr_phi;

void DefineUserBranches();

void Branch_mc_et(TTree* tree) { tree->Branch("mc_et", &mc_et, "mc_et[100]/F"); };
void Branch_mc_eta(TTree* tree) { tree->Branch("mc_eta", &mc_eta, "mc_eta[100]/F"); };
void Branch_mc_phi(TTree* tree) { tree->Branch("mc_phi", &mc_phi, "mc_phi[100]/F"); };
void Branch_fsr_et(TTree* tree) { tree->Branch("fsr_et", &fsr_et, "fsr_et[100]/F"); };
void Branch_fsr_eta(TTree* tree) { tree->Branch("fsr_eta", &fsr_eta, "fsr_eta[100]/F"); };
void Branch_fsr_phi(TTree* tree) { tree->Branch("fsr_phi", &fsr_phi, "fsr_phi[100]/F"); };
void Branch_higgs(TTree* tree) { tree->Branch("higgs", "TLorentzVector", &higgs, 32000, 0); };

void SetBranchAddress_mc_et(TTree * tree) { tree->SetBranchAddress("mc_et", &mc_et, &b_mc_et); }; 
void SetBranchAddress_mc_eta(TTree * tree) { tree->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta); }; 
void SetBranchAddress_mc_phi(TTree * tree) { tree->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi); }; 
void SetBranchAddress_fsr_et(TTree * tree) { tree->SetBranchAddress("fsr_et", &fsr_et, &b_fsr_et); }; 
void SetBranchAddress_fsr_eta(TTree * tree) { tree->SetBranchAddress("fsr_eta", &fsr_eta, &b_fsr_eta); }; 
void SetBranchAddress_fsr_phi(TTree * tree) { tree->SetBranchAddress("fsr_phi", &fsr_phi, &b_fsr_phi); }; 
void SetBranchAddress_higgs(TTree * tree) { tree->SetBranchAddress("higgs", &higgs, &b_higgs); }; 

// I/O
//
void Branch_gh_gen2reco1(TTree * tree) { tree->Branch("gh_gen2reco1",&gh_gen2reco1, "gh_gen2reco1/I");  }; 
void Branch_gh_gen2reco2(TTree * tree) { tree->Branch("gh_gen2reco2",&gh_gen2reco2, "gh_gen2reco2/I");  }; 
void Branch_gh_vbfq1_pdgid(TTree * tree) { tree->Branch("gh_vbfq1_pdgid",&gh_vbfq1_pdgid, "gh_vbfq1_pdgid/I");  }; 
void Branch_gh_vbfq2_pdgid(TTree * tree) { tree->Branch("gh_vbfq2_pdgid",&gh_vbfq2_pdgid, "gh_vbfq2_pdgid/I");  }; 
void Branch_gh_vh_pdgid(TTree * tree) { tree->Branch("gh_vh_pdgid",&gh_vh_pdgid, "gh_vh_pdgid/I");  }; 
void Branch_gh_vh1_pdgid(TTree * tree) { tree->Branch("gh_vh1_pdgid",&gh_vh1_pdgid, "gh_vh1_pdgid/I");  }; 
void Branch_gh_vh2_pdgid(TTree * tree) { tree->Branch("gh_vh2_pdgid",&gh_vh2_pdgid, "gh_vh2_pdgid/I");  }; 
//void Branch_METcorrected(TTree * tree) { tree->Branch("METcorrected", "TClonesArray",&METcorrected, 32000, 0); }; //met at analysis step
void Branch_gh_higgs_p4(TTree * tree) { tree->Branch("gh_higgs_p4", "TClonesArray",&gh_higgs_p4, 32000, 0); }; 
void Branch_gh_pho1_p4(TTree * tree) { tree->Branch("gh_pho1_p4", "TClonesArray",&gh_pho1_p4, 32000, 0); }; 
void Branch_gh_pho2_p4(TTree * tree) { tree->Branch("gh_pho2_p4", "TClonesArray",&gh_pho2_p4, 32000, 0); }; 
void Branch_gh_vbfq1_p4(TTree * tree) { tree->Branch("gh_vbfq1_p4", "TClonesArray",&gh_vbfq1_p4, 32000, 0); }; 
void Branch_gh_vbfq2_p4(TTree * tree) { tree->Branch("gh_vbfq2_p4", "TClonesArray",&gh_vbfq2_p4, 32000, 0); }; 
void Branch_gh_vh1_p4(TTree * tree) { tree->Branch("gh_vh1_p4", "TClonesArray",&gh_vh1_p4, 32000, 0); }; 
void Branch_gh_vh2_p4(TTree * tree) { tree->Branch("gh_vh2_p4", "TClonesArray",&gh_vh2_p4, 32000, 0); }; 
// vertex branches
void Branch_vtx_std_sel(TTree * tree) { tree->Branch("vtx_std_sel", &vtx_std_sel, "vtx_std_sel/I"); }; 
void Branch_vtx_std_evt_mva(TTree * tree) { tree->Branch("vtx_std_evt_mva", "std::vector<float>", &vtx_std_evt_mva); }; 
void Branch_vtx_std_ranked_list(TTree * tree) { tree->Branch("vtx_std_ranked_list", "std::vector<std::vector<int> >", &vtx_std_ranked_list); }; 
void Branch_pho_matchingConv(TTree * tree) { tree->Branch("pho_matchingConv", "std::vector<int>", &pho_matchingConv); }; 


void SetBranchAddress_gh_gen2reco1(TTree * tree) { tree->SetBranchAddress("gh_gen2reco1", &gh_gen2reco1, &b_gh_gen2reco1); }; 
void SetBranchAddress_gh_gen2reco2(TTree * tree) { tree->SetBranchAddress("gh_gen2reco2", &gh_gen2reco2, &b_gh_gen2reco2); }; 
void SetBranchAddress_gh_vbfq1_pdgid(TTree * tree) { tree->SetBranchAddress("gh_vbfq1_pdgid", &gh_vbfq1_pdgid, &b_gh_vbfq1_pdgid); }; 
void SetBranchAddress_gh_vbfq2_pdgid(TTree * tree) { tree->SetBranchAddress("gh_vbfq2_pdgid", &gh_vbfq2_pdgid, &b_gh_vbfq2_pdgid); }; 
void SetBranchAddress_gh_vh_pdgid(TTree * tree) { tree->SetBranchAddress("gh_vh_pdgid", &gh_vh_pdgid, &b_gh_vh_pdgid); }; 
void SetBranchAddress_gh_vh1_pdgid(TTree * tree) { tree->SetBranchAddress("gh_vh1_pdgid", &gh_vh1_pdgid, &b_gh_vh1_pdgid); }; 
void SetBranchAddress_gh_vh2_pdgid(TTree * tree) { tree->SetBranchAddress("gh_vh2_pdgid", &gh_vh2_pdgid, &b_gh_vh2_pdgid); }; 
//void SetBranchAddress_METcorrected(TTree * tree) { tree->SetBranchAddress("METcorrected", &METcorrected, &b_METcorrected); }; //met at analysis step
void SetBranchAddress_gh_higgs_p4(TTree * tree) { tree->SetBranchAddress("gh_higgs_p4", &gh_higgs_p4, &b_gh_higgs_p4); }; 
void SetBranchAddress_gh_pho1_p4(TTree * tree) { tree->SetBranchAddress("gh_pho1_p4", &gh_pho1_p4, &b_gh_pho1_p4); }; 
void SetBranchAddress_gh_pho2_p4(TTree * tree) { tree->SetBranchAddress("gh_pho2_p4", &gh_pho2_p4, &b_gh_pho2_p4); }; 
void SetBranchAddress_gh_vbfq1_p4(TTree * tree) { tree->SetBranchAddress("gh_vbfq1_p4", &gh_vbfq1_p4, &b_gh_vbfq1_p4); }; 
void SetBranchAddress_gh_vbfq2_p4(TTree * tree) { tree->SetBranchAddress("gh_vbfq2_p4", &gh_vbfq2_p4, &b_gh_vbfq2_p4); }; 
void SetBranchAddress_gh_vh1_p4(TTree * tree) { tree->SetBranchAddress("gh_vh1_p4", &gh_vh1_p4, &b_gh_vh1_p4); }; 
void SetBranchAddress_gh_vh2_p4(TTree * tree) { tree->SetBranchAddress("gh_vh2_p4", &gh_vh2_p4, &b_gh_vh2_p4); }; 
void SetBranchAddress_vtx_std_sel(TTree * tree) { tree->SetBranchAddress("vtx_std_sel", &vtx_std_sel, &b_vtx_std_sel); }; 
void SetBranchAddress_vtx_std_evt_mva(TTree * tree) { tree->SetBranchAddress("vtx_std_evt_mva", &vtx_std_evt_mva, &b_vtx_std_evt_mva); };
void SetBranchAddress_vtx_std_ranked_list(TTree * tree) { tree->SetBranchAddress("vtx_std_ranked_list", &vtx_std_ranked_list, &b_vtx_std_ranked_list); };
void SetBranchAddress_pho_matchingConv(TTree * tree) { tree->SetBranchAddress("pho_matchingConv", &pho_matchingConv, &b_pho_matchingConv); }; 

void Branch_dipho_n(TTree * tree) { tree->Branch("dipho_n", &dipho_n, "dipho_n/I"); };
void Branch_dipho_leadind(TTree * tree) { tree->Branch("dipho_leadind", &dipho_leadind, "dipho_leadind[dipho_n]/I"); };
void Branch_dipho_subleadind(TTree * tree) { tree->Branch("dipho_subleadind", &dipho_subleadind, "dipho_subleadind[dipho_n]/I"); };
void Branch_dipho_vtxind(TTree * tree) { tree->Branch("dipho_vtxind", &dipho_vtxind, "dipho_vtxind[dipho_n]/I"); };
void Branch_dipho_sumpt(TTree * tree) { tree->Branch("dipho_sumpt", &dipho_sumpt, "dipho_sumpt[dipho_n]/F"); };
void Branch_pho_genmatched(TTree *tree){tree->Branch("pho_genmatched", &pho_genmatched, "pho_genmatched[pho_n]/O");};
void Branch_pho_regr_energy_otf(TTree *tree){tree->Branch("pho_regr_energy_otf", &pho_regr_energy_otf, "pho_regr_energy_otf[pho_n]/F");};
void Branch_pho_regr_energyerr_otf(TTree *tree){tree->Branch("pho_regr_energyerr_otf", &pho_regr_energyerr_otf, "pho_regr_energyerr_otf[pho_n]/F");};
//// void Branch_dipho_leadet(TTree * tree) { tree->Branch("dipho_leadet", &dipho_leadet, "dipho_leadet[dipho_n]/F"); };
//// void Branch_dipho_subleadet(TTree * tree) { tree->Branch("dipho_subleadet", &dipho_subleadet, "dipho_subleadet[dipho_n]/F"); };
//// void Branch_dipho_leadeta(TTree * tree) { tree->Branch("dipho_leadeta", &dipho_leadeta, "dipho_leadeta[dipho_n]/F"); };
//// void Branch_dipho_subleadeta(TTree * tree) { tree->Branch("dipho_subleadeta", &dipho_subleadeta, "dipho_subleadeta[dipho_n]/F"); };
//// void Branch_dipho_leadci6cindex(TTree * tree) { tree->Branch("dipho_leadci6cindex", &dipho_leadci6cindex, "dipho_leadci6cindex[dipho_n]/I"); };
//// void Branch_dipho_subleadci6cindex(TTree * tree) { tree->Branch("dipho_subleadci6cindex", &dipho_subleadci6cindex, "dipho_subleadci6cindex[dipho_n]/I"); };
//// void Branch_dipho_leadci4cindex(TTree * tree) { tree->Branch("dipho_leadci4cindex", &dipho_leadci4cindex, "dipho_leadci4cindex[dipho_n]/I"); };
//// void Branch_dipho_subleadci4cindex(TTree * tree) { tree->Branch("dipho_subleadci4cindex", &dipho_subleadci4cindex, "dipho_subleadci4cindex[dipho_n]/I"); };
//// void Branch_dipho_mass(TTree * tree) { tree->Branch("dipho_mass", &dipho_mass, "dipho_mass[dipho_n]/F"); };
//// void Branch_dipho_pt(TTree * tree) { tree->Branch("dipho_pt", &dipho_pt, "dipho_pt[dipho_n]/F"); };
//// void Branch_dipho_eta(TTree * tree) { tree->Branch("dipho_eta", &dipho_eta, "dipho_eta[dipho_n]/F"); };
//// void Branch_dipho_phi(TTree * tree) { tree->Branch("dipho_phi", &dipho_phi, "dipho_phi[dipho_n]/F"); };
//// void Branch_dipho_cts(TTree * tree) { tree->Branch("dipho_cts", &dipho_cts, "dipho_cts[dipho_n]/F"); };

void SetBranchAddress_dipho_n(TTree * tree) { tree->SetBranchAddress("dipho_n", &dipho_n, &b_dipho_n); };
void SetBranchAddress_dipho_leadind(TTree * tree) { tree->SetBranchAddress("dipho_leadind", &dipho_leadind, &b_dipho_leadind); };
void SetBranchAddress_dipho_subleadind(TTree * tree) { tree->SetBranchAddress("dipho_subleadind", &dipho_subleadind, &b_dipho_subleadind); };
void SetBranchAddress_dipho_vtxind(TTree * tree) { tree->SetBranchAddress("dipho_vtxind", &dipho_vtxind, &b_dipho_vtxind); };
void SetBranchAddress_dipho_sumpt(TTree * tree) { tree->SetBranchAddress("dipho_sumpt", &dipho_sumpt, &b_dipho_sumpt); };
void SetBranchAddress_pho_genmatched(TTree * tree) { tree->SetBranchAddress("pho_genmatched", &pho_genmatched, &b_pho_genmatched); };
void SetBranchAddress_pho_regr_energy_otf(TTree * tree) { tree->SetBranchAddress("pho_regr_energy_otf", &pho_regr_energy_otf, &b_pho_regr_energy_otf); };
void SetBranchAddress_pho_regr_energyerr_otf(TTree * tree) { tree->SetBranchAddress("pho_regr_energyerr_otf", &pho_regr_energyerr_otf, &b_pho_regr_energyerr_otf); };
//// void SetBranchAddress_dipho_leadet(TTree * tree) { tree->SetBranchAddress("dipho_leadet", &dipho_leadet, &b_dipho_leadet); };
//// void SetBranchAddress_dipho_subleadet(TTree * tree) { tree->SetBranchAddress("dipho_subleadet", &dipho_subleadet, &b_dipho_subleadet); };
//// void SetBranchAddress_dipho_leadeta(TTree * tree) { tree->SetBranchAddress("dipho_leadeta", &dipho_leadeta, &b_dipho_leadeta); };
//// void SetBranchAddress_dipho_subleadeta(TTree * tree) { tree->SetBranchAddress("dipho_subleadeta", &dipho_subleadeta, &b_dipho_subleadeta); };
//// void SetBranchAddress_dipho_leadci6cindex(TTree * tree) { tree->SetBranchAddress("dipho_leadci6cindex", &dipho_leadci6cindex, &b_dipho_leadci6cindex); };
//// void SetBranchAddress_dipho_subleadci6cindex(TTree * tree) { tree->SetBranchAddress("dipho_subleadci6cindex", &dipho_subleadci6cindex, &b_dipho_subleadci6cindex); };
//// void SetBranchAddress_dipho_leadci4cindex(TTree * tree) { tree->SetBranchAddress("dipho_leadci4cindex", &dipho_leadci4cindex, &b_dipho_leadci4cindex); };
//// void SetBranchAddress_dipho_subleadci4cindex(TTree * tree) { tree->SetBranchAddress("dipho_subleadci4cindex", &dipho_subleadci4cindex, &b_dipho_subleadci4cindex); };
//// void SetBranchAddress_dipho_mass(TTree * tree) { tree->SetBranchAddress("dipho_mass", &dipho_mass, &b_dipho_mass); };
//// void SetBranchAddress_dipho_pt(TTree * tree) { tree->SetBranchAddress("dipho_pt", &dipho_pt, &b_dipho_pt); };
//// void SetBranchAddress_dipho_eta(TTree * tree) { tree->SetBranchAddress("dipho_eta", &dipho_eta, &b_dipho_eta); };
//// void SetBranchAddress_dipho_phi(TTree * tree) { tree->SetBranchAddress("dipho_phi", &dipho_phi, &b_dipho_phi); };
//// void SetBranchAddress_dipho_cts(TTree * tree) { tree->SetBranchAddress("dipho_cts", &dipho_cts, &b_dipho_cts); };


//correctMETinRED
void Branch_shiftMET_pt(TTree * tree) { tree->Branch("shiftMET_pt", &shiftMET_pt, "shiftMET_pt[dipho_n]/F"); };
void Branch_shiftMET_phi(TTree * tree) { tree->Branch("shiftMET_phi", &shiftMET_phi, "shiftMET_phi[dipho_n]/F"); };
void Branch_smearMET_pt(TTree * tree) { tree->Branch("smearMET_pt", &smearMET_pt, "smearMET_pt[dipho_n]/F"); };
void Branch_smearMET_phi(TTree * tree) { tree->Branch("smearMET_phi", &smearMET_phi, "smearMET_phi[dipho_n]/F"); };
void Branch_shiftscaleMET_pt(TTree * tree) { tree->Branch("shiftscaleMET_pt", &shiftscaleMET_pt, "shiftscaleMET_pt[dipho_n]/F"); };
void Branch_shiftscaleMET_phi(TTree * tree) { tree->Branch("shiftscaleMET_phi", &shiftscaleMET_phi, "shiftscaleMET_phi[dipho_n]/F"); };
void Branch_shiftsmearMET_pt(TTree * tree) { tree->Branch("shiftsmearMET_pt", &shiftsmearMET_pt, "shiftsmearMET_pt[dipho_n]/F"); };
void Branch_shiftsmearMET_phi(TTree * tree) { tree->Branch("shiftsmearMET_phi", &shiftsmearMET_phi, "shiftsmearMET_phi[dipho_n]/F"); };
void Branch_correctedpfMET(TTree * tree) { tree->Branch("correctedpfMET", &correctedpfMET, "correctedpfMET[dipho_n]/F"); };
void Branch_correctedpfMET_phi(TTree * tree) { tree->Branch("correctedpfMET_phi", &correctedpfMET_phi, "correctedpfMET_phi[dipho_n]/F"); };

void SetBranchAddress_shiftMET_pt(TTree * tree) { tree->SetBranchAddress("shiftMET_pt", &shiftMET_pt, &b_shiftMET_pt); };
void SetBranchAddress_shiftMET_phi(TTree * tree) { tree->SetBranchAddress("shiftMET_phi", &shiftMET_phi, &b_shiftMET_phi); };
void SetBranchAddress_smearMET_pt(TTree * tree) { tree->SetBranchAddress("smearMET_pt", &smearMET_pt, &b_smearMET_pt); };
void SetBranchAddress_smearMET_phi(TTree * tree) { tree->SetBranchAddress("smearMET_phi", &smearMET_phi, &b_smearMET_phi); };
void SetBranchAddress_shiftscaleMET_pt(TTree * tree) { tree->SetBranchAddress("shiftscaleMET_pt", &shiftscaleMET_pt, &b_shiftscaleMET_pt); };
void SetBranchAddress_shiftscaleMET_phi(TTree * tree) { tree->SetBranchAddress("shiftscaleMET_phi", &shiftscaleMET_phi, &b_shiftscaleMET_phi); };
void SetBranchAddress_shiftsmearMET_pt(TTree * tree) { tree->SetBranchAddress("shiftsmearMET_pt", &shiftsmearMET_pt, &b_shiftsmearMET_pt); };
void SetBranchAddress_shiftsmearMET_phi(TTree * tree) { tree->SetBranchAddress("shiftsmearMET_phi", &shiftsmearMET_phi, &b_shiftsmearMET_phi); };
void SetBranchAddress_correctedpfMET(TTree * tree) { tree->SetBranchAddress("correctedpfMET", &correctedpfMET, &b_correctedpfMET); };
void SetBranchAddress_correctedpfMET_phi(TTree * tree) { tree->SetBranchAddress("correctedpfMET_phi", &correctedpfMET_phi, &b_correctedpfMET_phi); };

void Branch_shiftMET_eta(TTree * tree) { tree->Branch("shiftMET_eta", &shiftMET_eta, "shiftMET_eta[dipho_n]/F"); };
void Branch_shiftMET_e(TTree * tree) { tree->Branch("shiftMET_e", &shiftMET_e, "shiftMET_e[dipho_n]/F"); };
void Branch_shiftscaleMET_eta(TTree * tree) { tree->Branch("shiftscaleMET_eta", &shiftscaleMET_eta, "shiftscaleMET_eta[dipho_n]/F"); };
void Branch_shiftscaleMET_e(TTree * tree) { tree->Branch("shiftscaleMET_e", &shiftscaleMET_e, "shiftscaleMET_e[dipho_n]/F"); };

void SetBranchAddress_shiftMET_eta(TTree * tree) { tree->SetBranchAddress("shiftMET_eta", &shiftMET_eta, &b_shiftMET_eta); };
void SetBranchAddress_shiftMET_e(TTree * tree) { tree->SetBranchAddress("shiftMET_e", &shiftMET_e, &b_shiftMET_e); };
void SetBranchAddress_shiftscaleMET_eta(TTree * tree) { tree->SetBranchAddress("shiftscaleMET_eta", &shiftscaleMET_eta, &b_shiftscaleMET_eta); };
void SetBranchAddress_shiftscaleMET_e(TTree * tree) { tree->SetBranchAddress("shiftscaleMET_e", &shiftscaleMET_e, &b_shiftscaleMET_e); };

// ID branches
void Branch_pho_mitmva(TTree * tree) { tree->Branch("pho_mitmva", "std::vector<std::vector<float> >", &pho_mitmva); }; 

void Branch_pho_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->Branch("pho_tkiso_recvtx_030_002_0000_10_01", "std::vector<std::vector<float> >", &pho_tkiso_recvtx_030_002_0000_10_01); }; 
void Branch_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->Branch("pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01", "std::vector<std::vector<float> >", &pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01); }; 
void Branch_pho_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->Branch("pho_tkiso_badvtx_040_002_0000_10_01", &pho_tkiso_badvtx_040_002_0000_10_01, "pho_tkiso_badvtx_040_002_0000_10_01[pho_n]/F" ); };
void Branch_pho_pfiso_charged_badvtx_04(TTree * tree) { tree->Branch("pho_pfiso_charged_badvtx_04", &pho_pfiso_charged_badvtx_04, "pho_pfiso_charged_badvtx_04[pho_n]/F" ); };
void Branch_pho_pfiso_charged_badvtx_id(TTree * tree) { tree->Branch("pho_pfiso_charged_badvtx_id", &pho_pfiso_charged_badvtx_id, "pho_pfiso_charged_badvtx_id[pho_n]/I" ); };

void Branch_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->Branch("pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01", &pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01, "pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01[pho_n]/F" ); };

void Branch_pho_tkiso_badvtx_id(TTree * tree) { tree->Branch("pho_tkiso_badvtx_id", &pho_tkiso_badvtx_id, "pho_tkiso_badvtx_id[pho_n]/I" ); };
void Branch_pho_ZeeVal_tkiso_badvtx_id(TTree * tree) { tree->Branch("pho_ZeeVal_tkiso_badvtx_id", &pho_ZeeVal_tkiso_badvtx_id, "pho_ZeeVal_tkiso_badvtx_id[pho_n]/I" ); };
void Branch_pho_drtotk_25_99(TTree * tree) { tree->Branch("pho_drtotk_25_99", &pho_drtotk_25_99, "pho_drtotk_25_99[pho_n]/F" ); };

void Branch_mu_glo_hasgsftrack(TTree * tree) { tree->Branch("mu_glo_hasgsftrack", &mu_glo_hasgsftrack, "mu_glo_hasgsftrack[mu_glo_n]/I" ); };
void SetBranchAddress_mu_glo_hasgsftrack(TTree * tree) { tree->SetBranchAddress("mu_glo_hasgsftrack", &mu_glo_hasgsftrack, &b_mu_glo_hasgsftrack); };

void SetBranchAddress_pho_mitmva(TTree * tree) { tree->SetBranchAddress("pho_mitmva", &pho_mitmva, &b_pho_mitmva); }; 

void SetBranchAddress_pho_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_tkiso_recvtx_030_002_0000_10_01", &pho_tkiso_recvtx_030_002_0000_10_01, &b_pho_tkiso_recvtx_030_002_0000_10_01); }; 
void SetBranchAddress_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01", &pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01, &b_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01); }; 
void SetBranchAddress_pho_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_tkiso_badvtx_040_002_0000_10_01", &pho_tkiso_badvtx_040_002_0000_10_01, &b_pho_tkiso_badvtx_040_002_0000_10_01); };

void SetBranchAddress_pho_pfiso_charged_badvtx_04(TTree * tree) { tree->SetBranchAddress("pho_pfiso_charged_badvtx_04", &pho_pfiso_charged_badvtx_04, &b_pho_pfiso_charged_badvtx_04); };
void SetBranchAddress_pho_pfiso_charged_badvtx_id(TTree * tree) { tree->SetBranchAddress("pho_pfiso_charged_badvtx_id", &pho_pfiso_charged_badvtx_id, &b_pho_pfiso_charged_badvtx_id); };

void SetBranchAddress_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01", &pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01, &b_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01); };
void SetBranchAddress_pho_tkiso_badvtx_id(TTree * tree) { tree->SetBranchAddress("pho_tkiso_badvtx_id", &pho_tkiso_badvtx_id, &b_pho_tkiso_badvtx_id); };
void SetBranchAddress_pho_ZeeVal_tkiso_badvtx_id(TTree * tree) { tree->SetBranchAddress("pho_ZeeVal_tkiso_badvtx_id", &pho_ZeeVal_tkiso_badvtx_id, &b_pho_ZeeVal_tkiso_badvtx_id); };
void SetBranchAddress_pho_drtotk_25_99(TTree * tree) { tree->SetBranchAddress("pho_drtotk_25_99", &pho_drtotk_25_99, &b_pho_drtotk_25_99); };

// These are missing in branchdef
void Branch_gv_pos(TTree * tree) { tree->Branch("gv_pos", "TClonesArray",&gv_pos, 32000, 0); }; 
void Branch_gv_n(TTree * tree) { tree->Branch("gv_n",&gv_n, "gv_n/I"); }; 
void Branch_pu_n(TTree * tree) { tree->Branch("pu_n", &pu_n, "pu_n/I"); };
void Branch_pu_zpos		(TTree * tree) { tree->Branch("pu_zpos", "std::vector<float>", &pu_zpos		); }; 
void Branch_pu_sumpt_lowpt	(TTree * tree) { tree->Branch("pu_sumpt_lowpt", "std::vector<float>", &pu_sumpt_lowpt	); }; 
void Branch_pu_sumpt_highpt	(TTree * tree) { tree->Branch("pu_sumpt_highpt", "std::vector<float>", &pu_sumpt_highpt	); }; 
void Branch_pu_ntrks_lowpt	(TTree * tree) { tree->Branch("pu_ntrks_lowpt", "std::vector<int>", &pu_ntrks_lowpt	); }; 
void Branch_pu_ntrks_highpt  (TTree * tree) { tree->Branch("pu_ntrks_highpt" ,  "std::vector<int>", &pu_ntrks_highpt ); }; 

void Branch_pho_cic6cutlevel_lead(TTree * tree) { tree->Branch("pho_cic6cutlevel_lead", "std::vector<std::vector<Short_t> >", &pho_cic6cutlevel_lead); };
void Branch_pho_cic6passcuts_lead(TTree * tree) { tree->Branch("pho_cic6passcuts_lead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic6passcuts_lead); };
void Branch_pho_cic6cutlevel_sublead(TTree * tree) { tree->Branch("pho_cic6cutlevel_sublead", "std::vector<std::vector<Short_t> >", &pho_cic6cutlevel_sublead); };
void Branch_pho_cic6passcuts_sublead(TTree * tree) { tree->Branch("pho_cic6passcuts_sublead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic6passcuts_sublead); };
void Branch_pho_cic4cutlevel_lead(TTree * tree) { tree->Branch("pho_cic4cutlevel_lead", "std::vector<std::vector<Short_t> >", &pho_cic4cutlevel_lead); };
void Branch_pho_cic4passcuts_lead(TTree * tree) { tree->Branch("pho_cic4passcuts_lead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic4passcuts_lead); };
void Branch_pho_cic4cutlevel_sublead(TTree * tree) { tree->Branch("pho_cic4cutlevel_sublead", "std::vector<std::vector<Short_t> >", &pho_cic4cutlevel_sublead); };
void Branch_pho_cic4passcuts_sublead(TTree * tree) { tree->Branch("pho_cic4passcuts_sublead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic4passcuts_sublead); };


void Branch_pho_cic4pfcutlevel_lead(TTree * tree) { tree->Branch("pho_cic4pfcutlevel_lead", "std::vector<std::vector<Short_t> >", &pho_cic4pfcutlevel_lead); };
void Branch_pho_cic4pfpasscuts_lead(TTree * tree) { tree->Branch("pho_cic4pfpasscuts_lead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic4pfpasscuts_lead); };
void Branch_pho_cic4pfcutlevel_sublead(TTree * tree) { tree->Branch("pho_cic4pfcutlevel_sublead", "std::vector<std::vector<Short_t> >", &pho_cic4pfcutlevel_sublead); };
void Branch_pho_cic4pfpasscuts_sublead(TTree * tree) { tree->Branch("pho_cic4pfpasscuts_sublead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic4pfpasscuts_sublead); };


void Branch_pho_cutlevel_lead(TTree * tree) { tree->Branch("pho_cutlevel_lead", "std::vector<std::vector<Short_t> >", &pho_cutlevel_lead); };
void Branch_pho_passcuts_lead(TTree * tree) { tree->Branch("pho_passcuts_lead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_passcuts_lead); };
void Branch_pho_cutlevel_sublead(TTree * tree) { tree->Branch("pho_cutlevel_sublead", "std::vector<std::vector<Short_t> >", &pho_cutlevel_sublead); };
void Branch_pho_passcuts_sublead(TTree * tree) { tree->Branch("pho_passcuts_sublead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_passcuts_sublead); };


void SetBranchAddress_gv_n(TTree * tree) { tree->SetBranchAddress("gv_n", &gv_n, &b_gv_n); }; 
void SetBranchAddress_gv_pos(TTree * tree) { tree->SetBranchAddress("gv_pos", &gv_pos, &b_gv_pos); }; 
void SetBranchAddress_pu_n(TTree * tree) { tree->SetBranchAddress("pu_n", &pu_n, &b_pu_n); };
void SetBranchAddress_pu_zpos(TTree * tree) { tree->SetBranchAddress("pu_zpos", &pu_zpos, &b_pu_zpos); };
void SetBranchAddress_pu_sumpt_lowpt(TTree * tree) { tree->SetBranchAddress("pu_sumpt_lowpt", &pu_sumpt_lowpt, &b_pu_sumpt_lowpt); };
void SetBranchAddress_pu_sumpt_highpt(TTree * tree) { tree->SetBranchAddress("pu_sumpt_highpt", &pu_sumpt_highpt, &b_pu_sumpt_highpt); };
void SetBranchAddress_pu_ntrks_lowpt(TTree * tree) { tree->SetBranchAddress("pu_ntrks_lowpt", &pu_ntrks_lowpt, &b_pu_ntrks_lowpt); };
void SetBranchAddress_pu_ntrks_highpt(TTree * tree) { tree->SetBranchAddress("pu_ntrks_highpt", &pu_ntrks_highpt, &b_pu_ntrks_highpt); };

void SetBranchAddress_pho_cic6cutlevel_lead(TTree * tree) { tree->SetBranchAddress("pho_cic6cutlevel_lead", &pho_cic6cutlevel_lead, &b_pho_cic6cutlevel_lead ); };
void SetBranchAddress_pho_cic6passcuts_lead(TTree * tree) { tree->SetBranchAddress("pho_cic6passcuts_lead", &pho_cic6passcuts_lead, &b_pho_cic6passcuts_lead ); };
void SetBranchAddress_pho_cic6cutlevel_sublead(TTree * tree) { tree->SetBranchAddress("pho_cic6cutlevel_sublead", &pho_cic6cutlevel_sublead, &b_pho_cic6cutlevel_sublead ); };
void SetBranchAddress_pho_cic6passcuts_sublead(TTree * tree) { tree->SetBranchAddress("pho_cic6passcuts_sublead", &pho_cic6passcuts_sublead, &b_pho_cic6passcuts_sublead ); };
void SetBranchAddress_pho_cic4cutlevel_lead(TTree * tree) { tree->SetBranchAddress("pho_cic4cutlevel_lead", &pho_cic4cutlevel_lead, &b_pho_cic4cutlevel_lead ); };
void SetBranchAddress_pho_cic4passcuts_lead(TTree * tree) { tree->SetBranchAddress("pho_cic4passcuts_lead", &pho_cic4passcuts_lead, &b_pho_cic4passcuts_lead ); };
void SetBranchAddress_pho_cic4cutlevel_sublead(TTree * tree) { tree->SetBranchAddress("pho_cic4cutlevel_sublead", &pho_cic4cutlevel_sublead, &b_pho_cic4cutlevel_sublead ); };
void SetBranchAddress_pho_cic4passcuts_sublead(TTree * tree) { tree->SetBranchAddress("pho_cic4passcuts_sublead", &pho_cic4passcuts_sublead, &b_pho_cic4passcuts_sublead ); };


void SetBranchAddress_pho_cic4pfcutlevel_lead(TTree * tree) { tree->SetBranchAddress("pho_cic4pfcutlevel_lead", &pho_cic4pfcutlevel_lead, &b_pho_cic4pfcutlevel_lead ); };
void SetBranchAddress_pho_cic4pfpasscuts_lead(TTree * tree) { tree->SetBranchAddress("pho_cic4pfpasscuts_lead", &pho_cic4pfpasscuts_lead, &b_pho_cic4pfpasscuts_lead ); };
void SetBranchAddress_pho_cic4pfcutlevel_sublead(TTree * tree) { tree->SetBranchAddress("pho_cic4pfcutlevel_sublead", &pho_cic4pfcutlevel_sublead, &b_pho_cic4pfcutlevel_sublead ); };
void SetBranchAddress_pho_cic4pfpasscuts_sublead(TTree * tree) { tree->SetBranchAddress("pho_cic4pfpasscuts_sublead", &pho_cic4pfpasscuts_sublead, &b_pho_cic4pfpasscuts_sublead ); };


void SetBranchAddress_pho_cutlevel_lead(TTree * tree) { tree->SetBranchAddress("pho_cutlevel_lead", &pho_cutlevel_lead, &b_pho_cutlevel_lead ); };
void SetBranchAddress_pho_passcuts_lead(TTree * tree) { tree->SetBranchAddress("pho_passcuts_lead", &pho_passcuts_lead, &b_pho_passcuts_lead ); };
void SetBranchAddress_pho_cutlevel_sublead(TTree * tree) { tree->SetBranchAddress("pho_cutlevel_sublead", &pho_cutlevel_sublead, &b_pho_cutlevel_sublead ); };
void SetBranchAddress_pho_passcuts_sublead(TTree * tree) { tree->SetBranchAddress("pho_passcuts_sublead", &pho_passcuts_sublead, &b_pho_passcuts_sublead ); };

void doJetMatching(TClonesArray & reco, TClonesArray & gen, Bool_t * match_flag, Bool_t * match_vbf_flag, Bool_t * match_b_flag, Bool_t * match_c_flag,
		   Bool_t * match_l_flag, Float_t * match_pt, Float_t * match_dr, Float_t maxDr=0.4 );		 

std::pair<int, int> Select2HighestPtJets(TLorentzVector& leadpho, TLorentzVector& subleadpho, Bool_t * jetid_flags=0);
int RescaleJetEnergy(bool force=false);

//Moriond2012
int MuonSelection(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind);
//ICHEP2012
int MuonSelection2012(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind);
bool MuonPhotonCuts2012(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector* thismu);
//HCP2012
int MuonSelection2012B(float muptcut=20.);
bool MuonPhotonCuts2012B(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector* thismu,float deltaRcut=1.0);
bool MuonLooseID2012(int indmu);
bool MuonTightID2012(int indmu, int vtxind=-1);
bool MuonIsolation2012(int indmu, float mupt, bool isTight=false);
int FindMuonVertex(int mu_ind);

//Moriond2012
int ElectronSelection(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind);
bool ElectronPhotonCuts(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector& ele);
//ICHEP2012
int ElectronSelection2012(TLorentzVector& pho1, TLorentzVector& pho2, int vtxind,  bool phodepend=true);
bool ElectronLooseEGammaID(int electronindex, int vertexindex=-1);
bool ElectronTightEGammaID(int electronindex, int vertexindex=-1);
//HCP2012
int ElectronSelectionMVA2012(float elptcut=20.);
bool ElectronMVACuts(int el_ind, int vtx_ind=-1);
bool ElectronPhotonCuts2012B(TLorentzVector& pho1, TLorentzVector& pho2, TLorentzVector& ele, bool includeVHlepPlusMet=false,float deltaRcut=1.0);
int FindElectronVertex(int el_ind);
void PhotonsToVeto(TLorentzVector* veto_p4, float drtoveto, std::vector<bool>& vetos, bool drtotkveto,float drgsftoveto=1.0);

int ElectronSelectionMVA2012_nocutOnMVA(float elptcut);
bool ElectronMVACuts_nocutOnMVA(int el_ind, int vtx_ind);
//HCP2012
TLorentzVector METCorrection2012B(TLorentzVector lead_p4, TLorentzVector sublead_p4, bool moriond2013MetCorrection);
//bool METAnalysis2012B(float MET);
bool METAnalysis2012B(TLorentzVector lead_p4, TLorentzVector sublead_p4, bool useUncor, bool doMETCleaning=true, bool moriond2013MetCorrection=false);
bool METCleaning2012B(TLorentzVector& lead_p4, TLorentzVector& sublead_p4, TLorentzVector& myMet);

//~ICHEP2012
//correct met
double ErrEt( double Et, double Eta);
TLorentzVector shiftMet(TLorentzVector *uncormet, bool isMC);
TLorentzVector correctMet_Simple(TLorentzVector& pho_lead, TLorentzVector& pho_sublead, TLorentzVector *uncormet, bool smearing, bool scale);


void SetAllMVA();
void FillMuonGsfTracks();
Float_t diphotonMVA(Int_t diphoton_id, Int_t leadingPho, Int_t subleadingPho, Int_t vtx, float vtxProb, TLorentzVector &leadP4, TLorentzVector &subleadP4, float sigmaMrv, float sigmaMwv, float sigmaMeonly, const char* idType, const char* bdtType, float photonID_1,float photonID_2);
float getDmOverDz(Int_t, Int_t, Float_t*);
Float_t deltaMassVtx(Int_t, Int_t, Float_t);

int IPhi(double phi){
  if (phi < -999) return 0;
  do {
	phi+=2*TMath::Pi();
  } while (phi<0);
  return int(TMath::Ceil(phi/(2 *TMath::Pi()/(18*4*5))));
};
int IEta (double eta){
  return int(TMath::Ceil((eta/1.479)*17*5) + (eta<0 ? -1 : 0));
};

void getIetaIPhi(int phoid, int & ieta, int & iphi ) const ;
bool CheckSphericalPhoton(int ieta, int iphi) const;
bool CheckSphericalPhoton(int phoind) const;

void VHNewLeptonCategorization(bool & VHlep1event, bool & VHlep2event, int diphotonVHlep_id, int vertex, bool VHelevent_prov, bool VHmuevent_prov, int el_ind, int mu_ind, float* smeared_pho_energy, float METcut, bool moriond2013MetCorrection);
void VHTwoMuonsEvents(bool & VHlep1event, bool & VHlep2event, int & diphotonVHlep_id, int & muVtx, float* smeared_pho_energy, float leadEtVHlepCut, float subleadEtVHlepCut, bool applyPtoverM, bool mvaselection, float diphobdt_output_Cut_VHLep, float phoidMvaCut, bool vetodipho, bool kinonly, const char * type);
void VHTwoElectronsEvents(bool & VHlep1event, bool & VHlep2event, int & diphotonVHlep_id, int & elVtx, float* smeared_pho_energy, float leadEtVHlepCut, float subleadEtVHlepCut, bool applyPtoverM, bool mvaselection, float diphobdt_output_Cut_VHLep, float phoidMvaCut, bool vetodipho, bool kinonly, const char * type);
 
private:
  Float_t photonIDMVA2012(Int_t, Int_t, TLorentzVector &, const char*);
  Float_t photonIDMVA2013(Int_t, Int_t, TLorentzVector &, const char*);
  Float_t photonIDMVA2013_7TeV(Int_t, Int_t, TLorentzVector &, const char*);
  Float_t photonIDMVA2011(Int_t, Int_t, TLorentzVector &, const char*);

  Normalization_8TeV *signalNormalizer;

#ifdef NewFeatures
#include "Marco/plotInteractive_h.h"
#endif

};


#endif
