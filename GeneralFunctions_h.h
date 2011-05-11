
void GlobeCtIsol(int, TLorentzVector*, float, float, float, Int_t&, Float_t&, Float_t&, Float_t&, Float_t&);
int GlobeMatchIsl(TLorentzVector*, Float_t&);


enum eIDLevel {UltraLoose, VeryLoose, Loose, Medium, Tight, SuperTight, HyperTight1, HyperTight2, HyperTight3, HyperTight4, Robust};
int ElectronClassification(int);
std::pair<bool, bool> ElectronId(int, eIDLevel); 
void eIDInfo(Int_t, Int_t&, Int_t&,Int_t eIDMaxLevel=10);
Float_t sipCalculator(int);

// Vertex analysis
void vertexAnalysis(HggVertexAnalyzer & vtxAna, int pho1, int pho2);
std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, int p1, int p2, std::vector<std::string> & vtxVarNames);

TLorentzVector get_pho_p4(int ipho, int ivtx);
void set_pho_p4(int ipho, int ivtx);

// end vertex analysis 

void FillCICInputs();
void FillCIC();

// CiC SELECTION CODE BEGIN - SSIMON

// defines photon CiC ID cuts for all cut levels
enum phoCiCIDLevel { phoLOOSE=0, phoMEDIUM, phoTIGHT, phoSUPERTIGHT, phoHYPERTIGHT1, phoHYPERTIGHT2, phoHYPERTIGHT3, phoHYPERTIGHT4, phoNCUTLEVELS };
enum phoCiCCuts { phoISOSUMOET=0,  phoISOSUMOETBAD,   phoTRKISOOETOM,   phoSIEIE,   phoHOVERE,   phoR9,   phoDRTOTK_25_99,   phoPIXEL, phoNCUTS };
enum phoCiCCategories { phoEBhighR9=0, phoEBmidR9, phoEBlowR9, phoEEhighR9, phoEEmidR9, phoEElowR9, phoNCATEGORIES };
void SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cuts_lead, float * cuts_sublead);
float cut_lead_isosumoet[phoNCUTLEVELS][6];
float cut_lead_isosumoetbad[phoNCUTLEVELS][6];
float cut_lead_trkisooet[phoNCUTLEVELS][6];
float cut_lead_sieie[phoNCUTLEVELS][6];
float cut_lead_hovere[phoNCUTLEVELS][6];
float cut_lead_r9[phoNCUTLEVELS][6];
float cut_lead_drtotk_25_99[phoNCUTLEVELS][6];
float cut_lead_pixel[phoNCUTLEVELS][6];
float cut_sublead_isosumoet[phoNCUTLEVELS][6];
float cut_sublead_isosumoetbad[phoNCUTLEVELS][6];
float cut_sublead_trkisooet[phoNCUTLEVELS][6];
float cut_sublead_sieie[phoNCUTLEVELS][6];
float cut_sublead_hovere[phoNCUTLEVELS][6];
float cut_sublead_r9[phoNCUTLEVELS][6];
float cut_sublead_drtotk_25_99[phoNCUTLEVELS][6];
float cut_sublead_pixel[phoNCUTLEVELS][6];

// loops through photons and returns indices to two photons passing desired selection 
// if more than one diphoton passes, returns pair with highest lead photon pt, or if lead is same, with highest sublead pt.
std::pair<int,int> DiphotonCiCSelection( phoCiCIDLevel LEADCUTLEVEL = phoLOOSE, phoCiCIDLevel SUBLEADCUTLEVEL = phoLOOSE, Float_t leadPtMin = 30, Float_t subleadPtMin = 20, bool applyPtoverM=false);

// for a photon index, applies all levels of cuts and returns the index to the highest cut level passed (can do lead and sublead - same for now)
int   PhotonCiCSelectionLevel( int photon_index, std::vector<std::vector<bool> > & ph_passcut, int doSublead=1);


// Functions to calculate variables used in CiC selection
Float_t DeltaRToTrack(Int_t photon_ind=-1, Int_t vtxind=-1, Float_t PtMin=1., Float_t dzmax=0.2, Float_t dxymax=0.1, int maxlosthits=0);
Float_t IsoEcalHitsSumEtNumCrystal( TVector3 *calopos, Float_t innerConeDR, Float_t outerConeDR, Float_t stripEtaHalfWidth, Float_t stripHalfLength=99.);
Float_t WorstSumTrackPtInCone(TLorentzVector *photon_p4, Int_t returnVtxIndex=0, Float_t PtMin=0, Float_t OuterConeRadius=0.3, Float_t InnerConeRadius=0.04, Float_t EtaStripHalfWidth=0.015, Float_t dzmax=0.2, Float_t dxymax=0.1);
Float_t SumTrackPtInCone(TLorentzVector *photon_p4, Int_t vtxind, Float_t PtMin=0, Float_t OuterConeRadius=0.3, Float_t InnerConeRadius=0.04, Float_t EtaStripHalfWidth=0.015, Float_t dzmax=0.2, Float_t dxymax=0.1);

//photon category functions (r9 and eta)
int PhotonCategory(int photonindex, int n_r9cat=3, int n_etacat=2) { 
  return PhotonR9Category(photonindex,n_r9cat) + n_r9cat*PhotonEtaCategory(photonindex,n_etacat);
}
Int_t PhotonR9Category(int photonindex, int n_r9cat=3) { 
  if(photonindex < 0) return -1;
  if(n_r9cat<2)return 0;
  int r9cat=0;
  float r9 = pho_r9[photonindex];
  if(n_r9cat==3) {
    r9cat = (Int_t)(r9<0.95) + (r9<0.9);// 0, 1, or 2 (high r9 --> low r9)
  } else if(n_r9cat==2) {
    r9cat = (Int_t)(r9<0.93);// 0, 1(high r9 --> low r9)
  }
  return r9cat;
}
int PhotonEtaCategory(int photonindex, int n_etacat=4) {
  if(photonindex < 0) return -1;
  if(n_etacat<2)return 0;
  int etacat;
  Float_t eta = fabs(((TVector3*)pho_calopos->At(photonindex))->Eta());
  if(n_etacat==4) {
    etacat = (Int_t)((eta>0.9) + (eta>1.479) + (eta>2.1));   // 0, 1, 2, or 3 (barrel --> endcap)
  } else if(n_etacat==2) {
    etacat = (Int_t)(eta>1.479);   // 0, 1 (barrel --> endcap)
  }
  return  etacat;
}
//diphoton category functions ( r9, eta, and diphoton pt)
int DiphotonCategory(Int_t leadind, Int_t subleadind, int n_r9cat=3, int n_etacat=4, int n_pThcat=0) {
  Int_t r9cat  =  TMath::Max(PhotonR9Category(leadind,n_r9cat),PhotonR9Category(subleadind,n_r9cat));
  Int_t etacat =  TMath::Max(PhotonEtaCategory(leadind,n_etacat),PhotonEtaCategory(subleadind,n_etacat));
  Int_t pThcat =  DiphotonPtCategory(leadind,subleadind,n_pThcat);
  return  (r9cat + n_r9cat*etacat + (n_r9cat*n_etacat)*pThcat);  // (n_r9cat*c_etacat*n_pThcat) categories
}
int DiphotonPtCategory(Int_t leadind, Int_t subleadind, int n_pThcat=0) {
  if(leadind < 0 || subleadind < 0) return -1;
  if(n_pThcat<2)return 0;
  int pThcat=0;
  TLorentzVector * leadp4 = (TLorentzVector*)pho_p4->At(leadind);
  TLorentzVector * subleadp4 = (TLorentzVector*)pho_p4->At(subleadind);
  double pTh = (*leadp4 + *subleadp4).Pt();
  if(n_pThcat == 2) {
    pThcat = (Int_t)(pTh < 40.);
  } else if (n_pThcat == 3) {
    pThcat = (Int_t)((pTh < 50.) + (pTh < 25.));
  }
  return pThcat;
}
// CiC SELECTION CODE END - SSIMON


///
/// Missing and addition branches
///
Float_t rho;
Int_t gv_n; 
TClonesArray * gv_pos; 
Int_t pu_n; 
std::vector<float> * pu_zpos;
std::vector<float> * pu_sumpt_lowpt;
std::vector<float> * pu_sumpt_highpt;
std::vector<int> * pu_ntrks_lowpt;
std::vector<int> * pu_ntrks_highpt;

// Vertex choiche
std::vector<int> * vtx_std_ranked_list;
int vtx_std_sel;

// CiC inputs
std::vector<std::vector<float> >* pho_tkiso_recvtx_030_002_0000_10_01;
Float_t pho_tkiso_badvtx_040_002_0000_10_01[MAX_PHOTONS];
Float_t pho_drtotk_25_99[MAX_PHOTONS];

std::vector<Short_t>* pho_cutlevel_lead;
std::vector<std::vector<UInt_t> >* pho_passcuts_lead;
std::vector<Short_t>* pho_cutlevel_sublead;
std::vector<std::vector<UInt_t> >* pho_passcuts_sublead;

TBranch *b_rho;
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

TBranch * b_pho_tkiso_recvtx_030_002_0000_10_01;
TBranch * b_pho_tkiso_badvtx_040_002_0000_10_01;
TBranch * b_pho_drtotk_25_99;

TBranch * b_pho_cutlevel_lead;
TBranch * b_pho_passcuts_lead;
TBranch * b_pho_cutlevel_sublead;
TBranch * b_pho_passcuts_sublead;

void DefineUserBranches();

// I/O
//
// vertex branches
void Branch_vtx_std_sel(TTree * tree) { tree->Branch("vtx_std_sel", &vtx_std_sel, "vtx_std_sel/I"); }; 
void Branch_vtx_std_ranked_list(TTree * tree) { tree->Branch("vtx_std_ranked_list", "std::vector<int>", &vtx_std_ranked_list); }; 

void SetBranchAddress_vtx_std_sel(TTree * tree) { tree->SetBranchAddress("vtx_std_sel", &vtx_std_sel, &b_vtx_std_sel); }; 
void SetBranchAddress_vtx_std_ranked_list(TTree * tree) { tree->SetBranchAddress("std_ranked_list", &vtx_std_ranked_list, &b_vtx_std_ranked_list); }; 

// ID branches
void Branch_pho_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->Branch("pho_tkiso_recvtx_030_002_0000_10_01", "std::vector<std::vector<float> >", &pho_tkiso_recvtx_030_002_0000_10_01); }; 
void Branch_pho_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->Branch("pho_tkiso_badvtx_040_002_0000_10_01", &pho_tkiso_badvtx_040_002_0000_10_01, "pho_tkiso_badvtx_040_002_0000_10_01[pho_n]/F" ); };
void Branch_pho_drtotk_25_99(TTree * tree) { tree->Branch("pho_drtotk_25_99", &pho_drtotk_25_99, "pho_drtotk_25_99[pho_n]/F" ); };

void SetBranchAddress_pho_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_tkiso_recvtx_030_002_0000_10_01", &pho_tkiso_recvtx_030_002_0000_10_01, &b_pho_tkiso_recvtx_030_002_0000_10_01); }; 
void SetBranchAddress_pho_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_tkiso_badvtx_040_002_0000_10_01", &pho_tkiso_badvtx_040_002_0000_10_01, &b_pho_tkiso_badvtx_040_002_0000_10_01); };
void SetBranchAddress_pho_drtotk_25_99(TTree * tree) { tree->SetBranchAddress("pho_drtotk_25_99", &pho_drtotk_25_99, &b_pho_drtotk_25_99); };

// These are missing in branchdef
void Branch_rho(TTree * tree) { tree->Branch("rho", &rho, "rho/F"); }; 
void Branch_gv_pos(TTree * tree) { tree->Branch("gv_pos", "TClonesArray",&gv_pos, 32000, 0); }; 
void Branch_gv_n(TTree * tree) { tree->Branch("gv_n",&gv_n, "gv_n/I"); }; 
void Branch_pu_n(TTree * tree) { tree->Branch("pu_n", &pu_n, "pu_n/I"); };
void Branch_pu_zpos		(TTree * tree) { tree->Branch("pu_zpos", "std::vector<float>", &pu_zpos		); }; 
void Branch_pu_sumpt_lowpt	(TTree * tree) { tree->Branch("pu_sumpt_lowpt", "std::vector<float>", &pu_sumpt_lowpt	); }; 
void Branch_pu_sumpt_highpt	(TTree * tree) { tree->Branch("pu_sumpt_highpt", "std::vector<float>", &pu_sumpt_highpt	); }; 
void Branch_pu_ntrks_lowpt	(TTree * tree) { tree->Branch("pu_ntrks_lowpt", "std::vector<int>", &pu_ntrks_lowpt	); }; 
void Branch_pu_ntrks_highpt  (TTree * tree) { tree->Branch("pu_ntrks_highpt" ,  "std::vector<int>", &pu_ntrks_highpt ); }; 

void Branch_pho_cutlevel_lead(TTree * tree) { tree->Branch("pho_cutlevel_lead", "std::vector<Short_t>", &pho_cutlevel_lead); };
void Branch_pho_passcuts_lead(TTree * tree) { tree->Branch("pho_passcuts_lead", "std::vector<std::vector<UInt_t> >", &pho_passcuts_lead); };
void Branch_pho_cutlevel_sublead(TTree * tree) { tree->Branch("pho_cutlevel_sublead", "std::vector<Short_t>", &pho_cutlevel_sublead); };
void Branch_pho_passcuts_sublead(TTree * tree) { tree->Branch("pho_passcuts_sublead", "std::vector<std::vector<UInt_t> >", &pho_passcuts_sublead); };

void SetBranchAddress_rho(TTree * tree) { tree->SetBranchAddress("rho", &rho, &b_rho); }; 
void SetBranchAddress_gv_n(TTree * tree) { tree->SetBranchAddress("gv_n", &gv_n, &b_gv_n); }; 
void SetBranchAddress_gv_pos(TTree * tree) { tree->SetBranchAddress("gv_pos", &gv_pos, &b_gv_pos); }; 
void SetBranchAddress_pu_n(TTree * tree) { tree->SetBranchAddress("pu_n", &pu_n, &b_pu_n); };
void SetBranchAddress_pu_zpos(TTree * tree) { tree->SetBranchAddress("pu_zpos", &pu_zpos, &b_pu_zpos); };
void SetBranchAddress_pu_sumpt_lowpt(TTree * tree) { tree->SetBranchAddress("pu_sumpt_lowpt", &pu_sumpt_lowpt, &b_pu_sumpt_lowpt); };
void SetBranchAddress_pu_sumpt_highpt(TTree * tree) { tree->SetBranchAddress("pu_sumpt_highpt", &pu_sumpt_highpt, &b_pu_sumpt_highpt); };
void SetBranchAddress_pu_ntrks_lowpt(TTree * tree) { tree->SetBranchAddress("pu_ntrks_lowpt", &pu_ntrks_lowpt, &b_pu_ntrks_lowpt); };
void SetBranchAddress_pu_ntrks_highpt(TTree * tree) { tree->SetBranchAddress("pu_ntrks_highpt", &pu_ntrks_highpt, &b_pu_ntrks_highpt); };

void SetBranchAddress_pho_cutlevel_lead(TTree * tree) { tree->SetBranchAddress("pho_cutlevel_lead", &pho_cutlevel_lead, &b_pho_cutlevel_lead ); };
void SetBranchAddress_pho_passcuts_lead(TTree * tree) { tree->SetBranchAddress("pho_passcuts_lead", &pho_passcuts_lead, &b_pho_passcuts_lead ); };
void SetBranchAddress_pho_cutlevel_sublead(TTree * tree) { tree->SetBranchAddress("pho_cutlevel_sublead", &pho_cutlevel_sublead, &b_pho_cutlevel_sublead ); };
void SetBranchAddress_pho_passcuts_sublead(TTree * tree) { tree->SetBranchAddress("pho_passcuts_sublead", &pho_passcuts_sublead, &b_pho_passcuts_sublead ); };

