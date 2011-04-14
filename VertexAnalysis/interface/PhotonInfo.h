#ifndef hgg_PhotonInfo_h
#define hgg_PhotonInfo_h

#include <vector>
#include <map>
#include <string>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMatrixDSym.h"

class PhotonInfo
{
public:
  PhotonInfo(const TVector3 & caloPosition, float energy);
  PhotonInfo(const TVector3 & caloPosition, 
	 const TVector3 &  bs, 
	 const TVector3 &  convVtx, 
	 float energy);
	
  
  const TVector3 & beamSpot() const { return beamSpot_; }
  const TVector3 & conversionVertex() const { return conversionVertex_; }
  const TVector3 & caloPosition() const { return caloPosition_; }
  float      energy()       const { return energy_; }
  TLorentzVector p4(float vtxx, float vtxy, float vtxz) const;
	
protected:
  TVector3 caloPosition_;
  TVector3 beamSpot_;
  TVector3 conversionVertex_;
  float energy_;
};


#endif


