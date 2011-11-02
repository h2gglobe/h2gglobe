
void GlobeCtIsol(int, TLorentzVector*, float, float, float, Int_t&, Float_t&, Float_t&, Float_t&, Float_t&);
int GlobeMatchIsl(TLorentzVector*, Float_t&);


enum eIDLevel {UltraLoose, VeryLoose, Loose, Medium, Tight, SuperTight, HyperTight1, HyperTight2, HyperTight3, HyperTight4, Robust};
int ElectronClassification(int);
std::pair<bool, bool> ElectronId(int, eIDLevel); 
void eIDInfo(Int_t, Int_t&, Int_t&,Int_t eIDMaxLevel=10);
Float_t sipCalculator(int);

// Match the Photon with the merged collection of conversions
PhotonInfo fillPhotonInfos(int p1, bool useAllConvs);
int matchPhotonToConversion(int); 
double phiNorm (float &phi);
double etaTransformation(  float EtaParticle , float Zvertex);

vector<double> generate_flat10_weights(TH1D* data_npu_estimated);


// Vertex analysis
void vertexAnalysis(HggVertexAnalyzer & vtxAna,  PhotonInfo pho1, PhotonInfo pho2);
//std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, int p1, int p2, std::vector<std::string> & vtxVarNames);
std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, PhotonInfo & pho1, PhotonInfo & pho2,
				 std::vector<std::string> & vtxVarNames, 					  
				 bool useMva=false, TMVA::Reader * reader=0, std::string tmvaMethod="");

TLorentzVector get_pho_p4(int ipho, int ivtx, float *pho_energy_array=0);
TLorentzVector get_pho_p4(int ipho, TVector3 * vtx, float * energy=0);
void set_pho_p4(int ipho, int ivtx, float *pho_energy_array=0);

// end vertex analysis 

void FillCICInputs();
void FillCIC();

// CiC SELECTION CODE BEGIN - SSIMON

// defines photon CiC ID cuts for all cut levels
enum phoCiCIDLevel { phoNOCUTS=0, phoLOOSE, phoMEDIUM, phoTIGHT, phoSUPERTIGHT, phoHYPERTIGHT1, phoHYPERTIGHT2, phoHYPERTIGHT3, phoHYPERTIGHT4, phoNCUTLEVELS };
enum phoCiCCuts { phoISOSUMOET=0,  phoISOSUMOETBAD,   phoTRKISOOETOM,   phoSIEIE,   phoHOVERE,   phoR9,   phoDRTOTK_25_99,   phoPIXEL, phoNCUTS };
enum phoCiC6Categories { phoCiC6EBhighR9=0, phoCiC6EBmidR9, phoCiC6EBlowR9, phoCiC6EEhighR9, phoCiC6EEmidR9, phoCiC6EElowR9, phoCiC6NCATEGORIES };
enum phoCiC4Categories { phoCiC4EBhighR9=0, phoCiC4EBlowR9, phoCiC4EEhighR9, phoCiC4EElowR9, phoCiC4NCATEGORIES };
void SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cic6_cuts_lead, float * cic6_cuts_sublead, float * cic4_cuts_lead, float * cic4_cuts_sublead);
float cic6_cut_lead_isosumoet[phoNCUTLEVELS][6];
float cic6_cut_lead_isosumoetbad[phoNCUTLEVELS][6];
float cic6_cut_lead_trkisooet[phoNCUTLEVELS][6];
float cic6_cut_lead_sieie[phoNCUTLEVELS][6];
float cic6_cut_lead_hovere[phoNCUTLEVELS][6];
float cic6_cut_lead_r9[phoNCUTLEVELS][6];
float cic6_cut_lead_drtotk_25_99[phoNCUTLEVELS][6];
float cic6_cut_lead_pixel[phoNCUTLEVELS][6];
float cic6_cut_sublead_isosumoet[phoNCUTLEVELS][6];
float cic6_cut_sublead_isosumoetbad[phoNCUTLEVELS][6];
float cic6_cut_sublead_trkisooet[phoNCUTLEVELS][6];
float cic6_cut_sublead_sieie[phoNCUTLEVELS][6];
float cic6_cut_sublead_hovere[phoNCUTLEVELS][6];
float cic6_cut_sublead_r9[phoNCUTLEVELS][6];
float cic6_cut_sublead_drtotk_25_99[phoNCUTLEVELS][6];
float cic6_cut_sublead_pixel[phoNCUTLEVELS][6];

float cic4_cut_lead_isosumoet[phoNCUTLEVELS][4];
float cic4_cut_lead_isosumoetbad[phoNCUTLEVELS][4];
float cic4_cut_lead_trkisooet[phoNCUTLEVELS][4];
float cic4_cut_lead_sieie[phoNCUTLEVELS][4];
float cic4_cut_lead_hovere[phoNCUTLEVELS][4];
float cic4_cut_lead_r9[phoNCUTLEVELS][4];
float cic4_cut_lead_drtotk_25_99[phoNCUTLEVELS][4];
float cic4_cut_lead_pixel[phoNCUTLEVELS][4];
float cic4_cut_sublead_isosumoet[phoNCUTLEVELS][4];
float cic4_cut_sublead_isosumoetbad[phoNCUTLEVELS][4];
float cic4_cut_sublead_trkisooet[phoNCUTLEVELS][4];
float cic4_cut_sublead_sieie[phoNCUTLEVELS][4];
float cic4_cut_sublead_hovere[phoNCUTLEVELS][4];
float cic4_cut_sublead_r9[phoNCUTLEVELS][4];
float cic4_cut_sublead_drtotk_25_99[phoNCUTLEVELS][4];
float cic4_cut_sublead_pixel[phoNCUTLEVELS][4];

// loops through photons and returns indices to two photons passing desired selection 
// if more than one diphoton passes, returns pair with highest lead photon pt, or if lead is same, with highest sublead pt.
int DiphotonCiCSelection( phoCiCIDLevel LEADCUTLEVEL = phoLOOSE, phoCiCIDLevel SUBLEADCUTLEVEL = phoLOOSE, Float_t leadPtMin = 30, Float_t subleadPtMin = 20, int ncategories=6, bool applyPtoverM=false, float *pho_energy_array=0);

// for a photon index, applies all levels of cuts and returns the index to the highest cut level passed (can do lead and sublead - same for now)
int   PhotonCiCSelectionLevel( int photon_index, int vertex_index, std::vector<std::vector<bool> > & ph_passcut, int ncategories=6, int doSublead=1, float *pho_energy_array=0);


// Functions to calculate variables used in CiC selection
Float_t DeltaRToTrack(Int_t photon_ind=-1, Int_t vtxind=-1, Float_t PtMin=1., Float_t dzmax=0.2, Float_t dxymax=0.1, int maxlosthits=0);
Float_t IsoEcalHitsSumEtNumCrystal( TVector3 *calopos, Float_t innerConeDR, Float_t outerConeDR, Float_t stripEtaHalfWidth, Float_t stripHalfLength=99.);
std::pair<Int_t, Float_t> WorstSumTrackPtInCone(int ipho, Int_t returnVtxIndex=0, Float_t PtMin=0, Float_t OuterConeRadius=0.3, Float_t InnerConeRadius=0.04, Float_t EtaStripHalfWidth=0.015, Float_t dzmax=0.2, Float_t dxymax=0.1);
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
    r9cat = (Int_t)(r9<0.94) + (r9<0.9);// 0, 1, or 2 (high r9 --> low r9)
  } else if(n_r9cat==2) {
    r9cat = (Int_t)(r9<0.94);// 0, 1(high r9 --> low r9)
  }
  return r9cat;
}
int PhotonEtaCategory(int photonindex, int n_etacat=4) {
  if(photonindex < 0) return -1;
  if(n_etacat<2)return 0;
  int etacat;
  // Float_t eta = fabs(((TVector3*)pho_calopos->At(photonindex))->Eta());
  Float_t eta = fabs(((TVector3*)sc_xyz->At(pho_scind[photonindex]))->Eta());
  if(n_etacat==4) {
    etacat = (Int_t)((eta>0.9) + (eta>1.479) + (eta>2.1));   // 0, 1, 2, or 3 (barrel --> endcap)
  } else if(n_etacat==2) {
    etacat = (Int_t)(eta>1.479);   // 0, 1 (barrel --> endcap)
  }
  return  etacat;
}
//diphoton category functions ( r9, eta, and diphoton pt)
int DiphotonCategory(Int_t leadind, Int_t subleadind, float pTh, int n_r9cat=3, int n_etacat=4, int n_pThcat=0, int nVtxCategories=0, float vtxMva=-1.) {
  Int_t r9cat  =  TMath::Max(PhotonR9Category(leadind,n_r9cat),PhotonR9Category(subleadind,n_r9cat));
  Int_t etacat =  TMath::Max(PhotonEtaCategory(leadind,n_etacat),PhotonEtaCategory(subleadind,n_etacat));
  Int_t pThcat =  DiphotonPtCategory(pTh,n_pThcat);
  Int_t vtxCat =  DiphotonVtxCategory(vtxMva,nVtxCategories);
  return  (r9cat + n_r9cat*etacat + (n_r9cat*n_etacat)*pThcat) + (n_r9cat*n_etacat*(n_pThcat>0?n_pThcat:1))*vtxCat;  // (n_r9cat*c_etacat*n_pThcat) categories
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
    pThcat = (Int_t)((pTh < 50.) + (pTh < 25.));
  }
  return pThcat;
}
// CiC SELECTION CODE END - SSIMON

// Functions movec from Tools.h
double DeltaPhi(double,double);

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

#define MAX_DIPHOTONS 50
Int_t dipho_n;
Int_t dipho_leadind[MAX_DIPHOTONS];
Int_t dipho_subleadind[MAX_DIPHOTONS];
Int_t dipho_vtxind[MAX_DIPHOTONS];
Float_t dipho_sumpt[MAX_DIPHOTONS];
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

TBranch * b_dipho_n;
TBranch * b_dipho_leadind;
TBranch * b_dipho_subleadind;
TBranch * b_dipho_vtxind;
TBranch * b_dipho_sumpt;
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

// Vertex choiche
int vtx_std_sel;
std::vector<int> *  dipho_vtx_std_sel;
std::vector<std::vector<int> > * vtx_std_ranked_list;
std::vector<float> * vtx_std_evt_mva;
// std::vector<int> * vtx_std_ranked_list;

// CiC inputs
std::vector<std::vector<float> >* pho_tkiso_recvtx_030_002_0000_10_01;
Float_t pho_tkiso_badvtx_040_002_0000_10_01[MAX_PHOTONS];
Int_t pho_tkiso_badvtx_id[MAX_PHOTONS];
Float_t pho_drtotk_25_99[MAX_PHOTONS];

bool runCiC;
std::vector<std::vector<Short_t> >* pho_cic6cutlevel_lead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic6passcuts_lead;
std::vector<std::vector<Short_t> >* pho_cic6cutlevel_sublead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic6passcuts_sublead;
std::vector<std::vector<Short_t> >* pho_cic4cutlevel_lead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic4passcuts_lead;
std::vector<std::vector<Short_t> >* pho_cic4cutlevel_sublead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_cic4passcuts_sublead;

std::vector<std::vector<Short_t> >* pho_cutlevel_lead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_passcuts_lead;
std::vector<std::vector<Short_t> >* pho_cutlevel_sublead;
std::vector<std::vector<std::vector<UInt_t> > >* pho_passcuts_sublead;


// Indices of conversions matching the photons
std::vector<int> * pho_matchingConv;
TBranch *b_pho_matchingConv;

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
TBranch *b_vtx_std_evt_mva;

TBranch * b_pho_tkiso_recvtx_030_002_0000_10_01;
TBranch * b_pho_tkiso_badvtx_040_002_0000_10_01;
TBranch * b_pho_tkiso_badvtx_id;
TBranch * b_pho_drtotk_25_99;

TBranch * b_pho_cic6cutlevel_lead;
TBranch * b_pho_cic6passcuts_lead;
TBranch * b_pho_cic6cutlevel_sublead;
TBranch * b_pho_cic6passcuts_sublead;
TBranch * b_pho_cic4cutlevel_lead;
TBranch * b_pho_cic4passcuts_lead;
TBranch * b_pho_cic4cutlevel_sublead;
TBranch * b_pho_cic4passcuts_sublead;

TBranch * b_pho_cutlevel_lead;
TBranch * b_pho_passcuts_lead;
TBranch * b_pho_cutlevel_sublead;
TBranch * b_pho_passcuts_sublead;

void DefineUserBranches();

// I/O
//
// vertex branches
void Branch_vtx_std_sel(TTree * tree) { tree->Branch("vtx_std_sel", &vtx_std_sel, "vtx_std_sel/I"); }; 
void Branch_vtx_std_evt_mva(TTree * tree) { tree->Branch("vtx_std_evt_mva", "std::vector<float>", &vtx_std_evt_mva); }; 
void Branch_vtx_std_ranked_list(TTree * tree) { tree->Branch("vtx_std_ranked_list", "std::vector<std::vector<int> >", &vtx_std_ranked_list); }; 
void Branch_pho_matchingConv(TTree * tree) { tree->Branch("pho_matchingConv", "std::vector<int>", &pho_matchingConv); }; 


void SetBranchAddress_vtx_std_sel(TTree * tree) { tree->SetBranchAddress("vtx_std_sel", &vtx_std_sel, &b_vtx_std_sel); }; 
void SetBranchAddress_vtx_std_evt_mva(TTree * tree) { tree->SetBranchAddress("vtx_std_evt_mva", &vtx_std_evt_mva, &b_vtx_std_evt_mva); };
void SetBranchAddress_vtx_std_ranked_list(TTree * tree) { tree->SetBranchAddress("vtx_std_ranked_list", &vtx_std_ranked_list, &b_vtx_std_ranked_list); };
void SetBranchAddress_pho_matchingConv(TTree * tree) { tree->SetBranchAddress("pho_matchingConv", &pho_matchingConv, &b_pho_matchingConv); }; 

void Branch_dipho_n(TTree * tree) { tree->Branch("dipho_n", &dipho_n, "dipho_n/I"); };
void Branch_dipho_leadind(TTree * tree) { tree->Branch("dipho_leadind", &dipho_leadind, "dipho_leadind[dipho_n]/I"); };
void Branch_dipho_subleadind(TTree * tree) { tree->Branch("dipho_subleadind", &dipho_subleadind, "dipho_subleadind[dipho_n]/I"); };
void Branch_dipho_vtxind(TTree * tree) { tree->Branch("dipho_vtxind", &dipho_vtxind, "dipho_vtxind[dipho_n]/I"); };
void Branch_dipho_sumpt(TTree * tree) { tree->Branch("dipho_sumpt", &dipho_sumpt, "dipho_sumpt[dipho_n]/F"); };
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



// ID branches
void Branch_pho_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->Branch("pho_tkiso_recvtx_030_002_0000_10_01", "std::vector<std::vector<float> >", &pho_tkiso_recvtx_030_002_0000_10_01); }; 
void Branch_pho_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->Branch("pho_tkiso_badvtx_040_002_0000_10_01", &pho_tkiso_badvtx_040_002_0000_10_01, "pho_tkiso_badvtx_040_002_0000_10_01[pho_n]/F" ); };
void Branch_pho_tkiso_badvtx_id(TTree * tree) { tree->Branch("pho_tkiso_badvtx_id", &pho_tkiso_badvtx_id, "pho_tkiso_badvtx_id[pho_n]/I" ); };
void Branch_pho_drtotk_25_99(TTree * tree) { tree->Branch("pho_drtotk_25_99", &pho_drtotk_25_99, "pho_drtotk_25_99[pho_n]/F" ); };

void SetBranchAddress_pho_tkiso_recvtx_030_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_tkiso_recvtx_030_002_0000_10_01", &pho_tkiso_recvtx_030_002_0000_10_01, &b_pho_tkiso_recvtx_030_002_0000_10_01); }; 
void SetBranchAddress_pho_tkiso_badvtx_040_002_0000_10_01(TTree * tree) { tree->SetBranchAddress("pho_tkiso_badvtx_040_002_0000_10_01", &pho_tkiso_badvtx_040_002_0000_10_01, &b_pho_tkiso_badvtx_040_002_0000_10_01); };
void SetBranchAddress_pho_tkiso_badvtx_id(TTree * tree) { tree->SetBranchAddress("pho_tkiso_badvtx_id", &pho_tkiso_badvtx_id, &b_pho_tkiso_badvtx_id); };
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

void Branch_pho_cic6cutlevel_lead(TTree * tree) { tree->Branch("pho_cic6cutlevel_lead", "std::vector<std::vector<Short_t> >", &pho_cic6cutlevel_lead); };
void Branch_pho_cic6passcuts_lead(TTree * tree) { tree->Branch("pho_cic6passcuts_lead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic6passcuts_lead); };
void Branch_pho_cic6cutlevel_sublead(TTree * tree) { tree->Branch("pho_cic6cutlevel_sublead", "std::vector<std::vector<Short_t> >", &pho_cic6cutlevel_sublead); };
void Branch_pho_cic6passcuts_sublead(TTree * tree) { tree->Branch("pho_cic6passcuts_sublead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic6passcuts_sublead); };
void Branch_pho_cic4cutlevel_lead(TTree * tree) { tree->Branch("pho_cic4cutlevel_lead", "std::vector<std::vector<Short_t> >", &pho_cic4cutlevel_lead); };
void Branch_pho_cic4passcuts_lead(TTree * tree) { tree->Branch("pho_cic4passcuts_lead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic4passcuts_lead); };
void Branch_pho_cic4cutlevel_sublead(TTree * tree) { tree->Branch("pho_cic4cutlevel_sublead", "std::vector<std::vector<Short_t> >", &pho_cic4cutlevel_sublead); };
void Branch_pho_cic4passcuts_sublead(TTree * tree) { tree->Branch("pho_cic4passcuts_sublead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_cic4passcuts_sublead); };

void Branch_pho_cutlevel_lead(TTree * tree) { tree->Branch("pho_cutlevel_lead", "std::vector<std::vector<Short_t> >", &pho_cutlevel_lead); };
void Branch_pho_passcuts_lead(TTree * tree) { tree->Branch("pho_passcuts_lead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_passcuts_lead); };
void Branch_pho_cutlevel_sublead(TTree * tree) { tree->Branch("pho_cutlevel_sublead", "std::vector<std::vector<Short_t> >", &pho_cutlevel_sublead); };
void Branch_pho_passcuts_sublead(TTree * tree) { tree->Branch("pho_passcuts_sublead", "std::vector<std::vector<std::vector<UInt_t> > >", &pho_passcuts_sublead); };


void SetBranchAddress_rho(TTree * tree) { tree->SetBranchAddress("rho", &rho, &b_rho); }; 
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

void SetBranchAddress_pho_cutlevel_lead(TTree * tree) { tree->SetBranchAddress("pho_cutlevel_lead", &pho_cutlevel_lead, &b_pho_cutlevel_lead ); };
void SetBranchAddress_pho_passcuts_lead(TTree * tree) { tree->SetBranchAddress("pho_passcuts_lead", &pho_passcuts_lead, &b_pho_passcuts_lead ); };
void SetBranchAddress_pho_cutlevel_sublead(TTree * tree) { tree->SetBranchAddress("pho_cutlevel_sublead", &pho_cutlevel_sublead, &b_pho_cutlevel_sublead ); };
void SetBranchAddress_pho_passcuts_sublead(TTree * tree) { tree->SetBranchAddress("pho_passcuts_sublead", &pho_passcuts_sublead, &b_pho_passcuts_sublead ); };
