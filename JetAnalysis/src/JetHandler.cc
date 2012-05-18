#include "JetAnalysis/interface/JetHandler.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "LoopAll.h"

// ---------------------------------------------------------------------------------------------------------------
JetHandler::JetHandler(const std::string & cfg, const std::string & pset, LoopAll & l):
    l_(l)
{
    const edm::ParameterSet & myPset = edm::readPSetsFrom(cfg)->getParameter<edm::ParameterSet>(pset);
    
    cutbased = new PileupJetIdAlgo(myPset.getParameter<edm::ParameterSet>("cutbased"));
    simple   = new PileupJetIdAlgo(myPset.getParameter<edm::ParameterSet>("simple"));
    full     = new PileupJetIdAlgo(myPset.getParameter<edm::ParameterSet>("full"));
}

// ---------------------------------------------------------------------------------------------------------------
JetHandler::~JetHandler()
{
    delete cutbased;
    delete simple;
    delete full;
}

// ---------------------------------------------------------------------------------------------------------------
void compute_beta(int ijet, int ivx)
{
    //////////// TVector3 * tkpos= (TVector3 *) tk_vtx_pos->At(itk);
    //////////// /// double deltaz = fabs(vtxpos->Z() - tkpos->Z()); 
    //////////// double deltaz = fabs( (tkpos->Z()-vtxpos->Z()) - ( (tkpos->X()-vtxpos->X())*tkp4->Px() + (tkpos->Y()-vtxpos->Y())*tkp4->Py() )/tkp4->Pt() * tkp4->Pz()/tkp4->Pt() );

}

// ---------------------------------------------------------------------------------------------------------------
void compute_mva(int ijet, int ivtx)
{
    
}

// ---------------------------------------------------------------------------------------------------------------
void compute_wp(int ijet, int ivtx)
{
    
}

// ---------------------------------------------------------------------------------------------------------------
void recompute_jec(int ijet, int ivtx)
{
    
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
