#include "../interface/PhotonInfo.h"

#include <algorithm>
#include <numeric>
#include <iterator>

#include <assert.h>
#include "TVectorD.h"


using namespace std;


PhotonInfo::PhotonInfo(int id, const TVector3 & caloPosition, float energy) :
	id_(id), caloPosition_(caloPosition), energy_(energy), nTracks_(0), isFake_(false)
{
}

PhotonInfo::PhotonInfo(int id,
		       const TVector3 & caloPosition, 
		       const TVector3 &  bs, 
		       const TVector3 &  convVtx,
               const TVector3 & refittedMomentum,
		       float energy,
		       int   iDet,
		       int   nTracks,
		       bool  convVtxValid,
		       float convVtxChi2Prob,
		       float EoP
 ):
	id_(id),
  caloPosition_(caloPosition), 
  beamSpot_(bs),
  conversionVertex_(convVtx),
  refittedMomentum_(refittedMomentum),
  energy_(energy),
  iDet_(iDet),
  nTracks_(nTracks),
  convVtxValid_(convVtxValid),
  convVtxChi2Prob_(convVtxChi2Prob),
	EoP_(EoP),
	isFake_(false)
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

bool  PhotonInfo::isAConversion() const {
  bool isAConversion=false;
  if (  nTracks_ == 2  &&  convVtxValid_ &&   convVtxChi2Prob_ > 0.000001 )  isAConversion=true;

  return isAConversion;

}
