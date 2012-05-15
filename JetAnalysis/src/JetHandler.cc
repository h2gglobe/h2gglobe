#include "JetAnalysis/interface/JetHandler.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

JetHandler::JetHandler(const std::string & cfg, const std::string & pset)
{
    const edm::ParameterSet & myPset = edm::readPSetsFrom(cfg)->getParameter<edm::ParameterSet>(pset);
    
    jetIdAlgo = new PileupJetIdAlgo(myPset.getParameter<edm::ParameterSet>("jetIdAlgo"));
}

JetHandler::~JetHandler()
{
    delete jetIdAlgo;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
