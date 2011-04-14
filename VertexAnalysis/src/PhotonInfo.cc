#include "../interface/PhotonInfo.h"

#include <algorithm>
#include <numeric>
#include <iterator>

#include <assert.h>
#include "TVectorD.h"


using namespace std;


PhotonInfo::PhotonInfo(const TVector3 & caloPosition, float energy) :
	caloPosition_(caloPosition), energy_(energy)
{
}

PhotonInfo::PhotonInfo(const TVector3 & caloPosition, 
		       const TVector3 &  bs, 
		       const TVector3 &  convVtx, 
		       float energy,
		       int   iDet,
		       int   nTracks,
		       bool  convVtxValid,
		       float convVtxChi2Prob,
		       float EoP
 ):
  caloPosition_(caloPosition), 
  beamSpot_(bs),
  conversionVertex_(convVtx),
  energy_(energy),
  iDet_(iDet),
  nTracks_(nTracks),
  convVtxValid_(convVtxValid),
  convVtxChi2Prob_(convVtxChi2Prob),
  EoP_(EoP)
{
}
	
// -------------------------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector PhotonInfo::p4(float vtxx, float vtxy, float vtxz) const
{
	TVector3 vPos(vtxx,vtxy,vtxz);
	TVector3 direction = caloPosition() - vPos;
	TVector3 p = direction.Unit() * energy();
	TLorentzVector p4(p.x(),p.y(),p.z(),energy());
	return p4;
}

bool  PhotonInfo::isAConversion() {
  bool isAConversion=false;
  if (  nTracks_ == 2  &&  convVtxValid_ &&   convVtxChi2Prob_ > 0.0005 )  isAConversion=true;

  return isAConversion;

}
