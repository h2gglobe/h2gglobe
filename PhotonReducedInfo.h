#ifndef PhotonReducedInfo_h
#define PhotonReducedInfo_h

#include "TVector3.h"
#include "TLorentzVector.h"

class PhotonReducedInfo
{
public:

  
  PhotonReducedInfo(const TVector3 & caloPosition, float energy, int iDet, float r9);

  const TVector3 & caloPosition() const { return caloPosition_; }
  int iDet() const {return iDet_;}
  float energy()       const { return energy_; }
  float r9()       const { return r9_; }
  TLorentzVector p4(float vtxx, float vtxy, float vtxz) const;

  void setEnergy(float energy) {energy_=energy; };
  void setCaloPosition(const TVector3 & caloPosition) { caloPosition_=caloPosition; };
  void setR9(float r9) {r9_=r9; };
  void setDet(int det) { iDet_=det; };

protected:

  TVector3 caloPosition_;
  float energy_;
  int iDet_;
  float r9_;
};


#endif


