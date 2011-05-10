
void GlobeCtIsol(int, TLorentzVector*, float, float, float, Int_t&, Float_t&, Float_t&, Float_t&, Float_t&);
int GlobeMatchIsl(TLorentzVector*, Float_t&);


enum eIDLevel {UltraLoose, VeryLoose, Loose, Medium, Tight, SuperTight, HyperTight1, HyperTight2, HyperTight3, HyperTight4, Robust};
int ElectronClassification(int);
std::pair<bool, bool> ElectronId(int, eIDLevel); 
void eIDInfo(Int_t, Int_t&, Int_t&,Int_t eIDMaxLevel=10);
Float_t sipCalculator(int);

void vertexAnalysis(HggVertexAnalyzer & vtxAna, int pho1, int pho2);
//std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, int p1, int p2, std::vector<std::string> & vtxVarNames);
std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, int p1, int p2, std::vector<std::string> & vtxVarNames);

TLorentzVector get_pho_p4(int ipho, int ivtx);

std::pair<Float_t,Int_t> TrackIsoHgg(TLorentzVector *photon_p4, Int_t vtxind, Float_t PtMin=0, 
				     Float_t OuterConeRadius=0.3, Float_t InnerConeRadius=0.04, Float_t EtaStripHalfWidth=0.015, 
				     Float_t dzmax=0.2, Float_t dxymax=0.1);

void fillTrackIsolation(float tkIso_ptmin,float tkIso_outerCone,float tkIso_innerCone,float tkIso_etaStripHalfW,float tkIso_dzmax,float tkIso_dxymax);
