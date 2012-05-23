#include "JetAnalysis/interface/JetHandler.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "LoopAll.h"

// ---------------------------------------------------------------------------------------------------------------
JetHandler::JetHandler(const std::string & cfg, LoopAll & l):
    l_(l)
{
    boost::shared_ptr<edm::ParameterSet> myPset = edm::readConfig(cfg);
    
    cutbased = new PileupJetIdAlgo(myPset->getParameter<edm::ParameterSet>("cutbased"));
    simple   = new PileupJetIdAlgo(myPset->getParameter<edm::ParameterSet>("simple"));
    full     = new PileupJetIdAlgo(myPset->getParameter<edm::ParameterSet>("full"));
    
    
    jecCorData_ = 0;
    jecCorMc_   = 0;
    std::vector<std::string> dataJEC = myPset->getParameter<std::vector<std::string> >("dataJEC");
    std::vector<std::string> mcJEC   = myPset->getParameter<std::vector<std::string> >("mcJEC");
    
    for(std::vector<std::string>::iterator corr = dataJEC.begin(); corr!=dataJEC.end(); ++corr) {
	jetCorParsData_.push_back( JetCorrectorParameters(*corr) );
    }

    for(std::vector<std::string>::iterator corr = mcJEC.begin(); corr!=mcJEC.end(); ++corr) {
	jetCorParsMc_.push_back( JetCorrectorParameters(*corr) );
    }
    
    if( ! jetCorParsData_.empty() ) {
	jecCorData_ = new FactorizedJetCorrector(jetCorParsData_);
    }

    if( ! jetCorParsMc_.empty() ) {
	jecCorMc_ = new FactorizedJetCorrector(jetCorParsMc_);
    }
    
}

// ---------------------------------------------------------------------------------------------------------------
JetHandler::~JetHandler()
{
    delete cutbased;
    delete simple;
    delete full;
}

// ---------------------------------------------------------------------------------------------------------------
float dz(TLorentzVector * tkp4, TVector3 * tkpos, TVector3 * vtxpos)
{
    return fabs( (tkpos->Z()-vtxpos->Z()) - 
		 ( (tkpos->X()-vtxpos->X())*tkp4->Px() + (tkpos->Y()-vtxpos->Y())*tkp4->Py() )/tkp4->Pt() * tkp4->Pz()/tkp4->Pt() );
}

// ---------------------------------------------------------------------------------------------------------------
void JetHandler::compute_betas(int ijet, int vtx)
{
    float sumTkPt = 0.;
    float & beta = (*l_.jet_algoPF1_beta_ext)[ijet][vtx];
    float & betaStar = (*l_.jet_algoPF1_betaStar_ext)[ijet][vtx];
    float & betaStarClassic = (*l_.jet_algoPF1_betaStarClassic_ext)[ijet][vtx];
    beta = 0., betaStar = 0., betaStarClassic = 0.;

    const std::vector<unsigned short> & vtx_tracks = (*l_.vtx_std_tkind)[vtx];
    TVector3 * vtxpos = (TVector3*)l_.vtx_std_xyz->At(vtx);
    
    std::cout << l_.jet_algoPF1_tkind << " " << l_.jet_algoPF1_tkind->size() << std::endl;
    std::cout << l_.jet_algoPF2_tkind << " " << l_.jet_algoPF2_tkind->size() << std::endl;
    std::cout << l_.jet_algoPF3_tkind << " " << l_.jet_algoPF3_tkind->size() << std::endl;
    const std::vector<unsigned short> & jet_tracks = (*l_.jet_algoPF1_tkind)[ijet];
    for(std::vector<unsigned short>::const_iterator itrack=jet_tracks.begin(); itrack!=jet_tracks.end(); ++itrack ) {
	bool inVtx0 = find( vtx_tracks.begin(),  vtx_tracks.end(), *itrack ) != vtx_tracks.end();
	bool inAnyOther = false;
	
	TVector3 * tkpos= (TVector3 *) l_.tk_vtx_pos->At(*itrack);
	TLorentzVector * tkp4= (TLorentzVector *) l_.tk_p4->At(*itrack);
	float tkpt = tkp4->Pt();
	
	float dZ0 = dz(tkp4, tkpos, vtxpos);
	float dZ = dZ0;
	
	for(int ivtx=0; ivtx<l_.vtx_std_n; ++ivtx) {
	    bool isVtx0 = (ivtx == vtx);
	    if( ! isVtx0 && ! inAnyOther ) {
		const std::vector<unsigned short> & ivtx_tracks = (*l_.vtx_std_tkind)[ivtx];
		inAnyOther = find( ivtx_tracks.begin(),  ivtx_tracks.end(), *itrack ) != ivtx_tracks.end();
	    }
	    TVector3 * ivtxpos = (TVector3*)l_.vtx_std_xyz->At(ivtx);
	    dZ = std::min(dZ, dz(tkp4, tkpos, ivtxpos));
	}
	
	sumTkPt += tkpt;
	if( ! inVtx0 && inAnyOther ) {
	    betaStarClassic += tkpt; 
	}
	if( dZ0 < 0.2 ) {
	    beta += tkpt;
	} else if( dZ < 0.2 ) {
	    betaStar += tkpt;
	}
    }
    
    if( sumTkPt != 0. ) {
	beta     /= sumTkPt;
	betaStar /= sumTkPt;
	betaStarClassic /= sumTkPt;
    } else {
	assert( beta == 0. && betaStar == 0. && betaStarClassic == 0. );
    }
    
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
