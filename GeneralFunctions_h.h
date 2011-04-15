
void GlobeCtIsol(int, TLorentzVector*, float, float, float, Int_t&, Float_t&, Float_t&, Float_t&, Float_t&);
int GlobeMatchIsl(TLorentzVector*, Float_t&);


enum eIDLevel {UltraLoose, VeryLoose, Loose, Medium, Tight, SuperTight, HyperTight1, HyperTight2, HyperTight3, HyperTight4, Robust};
int ElectronClassification(int);
std::pair<bool, bool> ElectronId(int, eIDLevel); 
void eIDInfo(Int_t, Int_t&, Int_t&,Int_t eIDMaxLevel=10);
Float_t sipCalculator(int);

void vertexAnalysis(int pho1, int pho2);

