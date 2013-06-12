#ifndef __PhotonPlusJetVertexAnalysis__
#define __PhotonPlusJetVertexAnalysis__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "PhotonAnalysis/interface/PhotonAnalysis.h"
#include "PhotonAnalysis/interface/StatAnalysis.h"
#include "PhotonJetAnalysis/interface/PhotonJetAnalysis.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include "PdfWeightSmearer.h"
#include "InterferenceSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------

class PhotonPlusJetVertexAnalysis : public StatAnalysis 
{
 public:
    
    PhotonPlusJetVertexAnalysis();
    virtual ~PhotonPlusJetVertexAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    bool SkimEvents(LoopAll& l, int jentry);
    void FillReductionVariables(LoopAll& l, int jentry);   
    void ReducedOutputTree(LoopAll &l, TTree * outputTree);
    bool SelectEventsReduction(LoopAll&, int);

    bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
		      bool & isCorrectVertex, float &kinematic_bdtout,
		      bool isSyst=false, 
		      float syst_shift=0., bool skipSelection=false,
		      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 

    void fillControlPlots( LoopAll & l, TLorentzVector photon, TLorentzVector jet, HggVertexAnalyzer vtxAna_, int phoind, float evweight);

 private:   

};

#endif
