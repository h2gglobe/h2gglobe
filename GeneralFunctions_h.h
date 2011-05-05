
void GlobeCtIsol(int, TLorentzVector*, float, float, float, Int_t&, Float_t&, Float_t&, Float_t&, Float_t&);
int GlobeMatchIsl(TLorentzVector*, Float_t&);


enum eIDLevel {UltraLoose, VeryLoose, Loose, Medium, Tight, SuperTight, HyperTight1, HyperTight2, HyperTight3, HyperTight4, Robust};
int ElectronClassification(int);
std::pair<bool, bool> ElectronId(int, eIDLevel); 
void eIDInfo(Int_t, Int_t&, Int_t&,Int_t eIDMaxLevel=10);
Float_t sipCalculator(int);

void vertexAnalysis(HggVertexAnalyzer & vtxAna, int pho1, int pho2);
std::vector<int> vertexSelection(HggVertexAnalyzer & vtxAna, int p1, int p2, std::vector<std::string> & vtxVarNames);

TLorentzVector get_pho_p4(int ipho, int ivtx);
