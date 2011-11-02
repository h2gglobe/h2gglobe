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

  PhotonInfo(int id, const TVector3 & caloPosition, float energy);
  PhotonInfo(int id, const TVector3 & caloPosition, 
	     const TVector3 &  bs, 
	     const TVector3 &  convVtx, 
	     float energy,
	     int iDet,
	     int nTracks,
	     bool convVtxValid,
	     float convVtxChi2Prob,
	     float EoP
	     );
	
	int id() const { return id_; }  
	const TVector3 & beamSpot() const { return beamSpot_; }
	const TVector3 & conversionVertex() const { return conversionVertex_; }
	const TVector3 & caloPosition() const { return caloPosition_; }
	const int iDet() const {return iDet_;}
	const float      energy()       const { return energy_; }
	TLorentzVector p4(float vtxx, float vtxy, float vtxz) const;
	bool isAConversion() const;
	bool isFake() const { return isFake_; };
	void isFake(bool x=true) { isFake_ = x; };
	
protected:
	int id_;
	TVector3 caloPosition_;
	TVector3 beamSpot_;
	TVector3 conversionVertex_;
	float energy_;
	int iDet_;
	int nTracks_;
	bool convVtxValid_;
	float convVtxChi2Prob_;
	float EoP_;
	bool isFake_;
};


#endif


