
#ifndef __TapAnalysis__
#define __TapAnalysis__

#include "PhotonAnalysis/interface/StatAnalysis.h"

#include "TMVA/Reader.h"

// ------------------------------------------------------------------------------------
class TapAnalysis : public StatAnalysis {
 public:
  TapAnalysis();
  virtual ~TapAnalysis();
  
  virtual const std::string & name() const { return name_; };
  
  // LoopAll analysis interface implementation
  virtual void Init(LoopAll&);
  virtual void Term(LoopAll&);
  
  virtual void ReducedOutputTree(LoopAll &l, TTree *);
  virtual void GetBranches(TTree *, std::set<TBranch *>& );

  virtual void ResetAnalysis();
  
  virtual void FillReductionVariables(LoopAll& l, int jentry);   
  virtual bool SelectEventsReduction(LoopAll&, int);
  
  virtual bool SkimEvents(LoopAll&, int);
  virtual bool SelectEvents(LoopAll&, int);
  virtual bool Analysis(LoopAll&, Int_t);

  void FillFlatTree(LoopAll&, Int_t, Int_t, Int_t, Int_t, Int_t,
		    Float_t, Float_t, Float_t, Float_t);

  bool ElectronId(LoopAll&, Int_t, Int_t, std::string, Float_t);
  Float_t PhotonId(LoopAll&, Int_t, Int_t, std::string, Float_t);
  bool checkEventHLT(LoopAll&, std::vector<std::string>);
  std::vector<std::pair<int, int> > TPPairs(LoopAll&, std::vector<int> tags, std::vector<int> probes, int type, int chargePairing);
  Int_t ChooseVertex(LoopAll&, int, bool);
  int PhotonIDCategory(LoopAll&, int, int);
  bool CiCPhotonIDPF(LoopAll& l, int nCategories, int photonindex, int chosenVtx, int IDlevel);
  TLorentzVector get_pho_p4(LoopAll& l, Int_t ipho, int ivtx);

  float GetR9Weight(LoopAll& l, Int_t ipho);
  int GetPrescaleWeight(int run, int lumi);
  void antiRescaleClusterVariables(LoopAll& l);

  Float_t forElectrons;
  Float_t cutPFMET, cutETTag, cutETProbe;
  std::string selectionTypeTag, selectionTypeProbe, selectionTypeToMeasure;
  Float_t cutSelectionTag, cutSelectionProbe, cutSelectionToMeasure;
  Float_t minZMass, maxZMass;
  Float_t applyPreselection;
  bool hltPrescaleWeight;

 protected:
  std::string name_;
  TMVA::Reader *tmvaReaderID_Single_Barrel, *tmvaReaderID_Single_Endcap;
  std::vector<std::string> hltPaths;
  std::vector<std::string> hltPathsDE;

  std::vector<float> r9Weight;
  std::vector<int> runList, lumiList, scaleList;
};

#endif
