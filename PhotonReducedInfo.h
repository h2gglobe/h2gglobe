#ifndef PhotonReducedInfo_h
#define PhotonReducedInfo_h

#include "TVector3.h"
#include "TLorentzVector.h"

class PhotonReducedInfo
{
public:
  
  PhotonReducedInfo();
  PhotonReducedInfo(const TVector3 & caloPosition, float energy, float corrEnergy, int iDet, float r9, bool passId, float corrEnergyErr=0);

  
  PhotonReducedInfo(const PhotonReducedInfo &obj){copy_(obj);}
  PhotonReducedInfo & operator= (const PhotonReducedInfo &obj);

  const TVector3 & caloPosition() const { return caloPosition_; }
  int iDet() const {return iDet_;}
  float energy()       const { return energy_; }
  float corrEnergy()       const { return corrEnergy_; }
  float corrEnergyErr()       const { return corrEnergyErr_; }
  float r9()       const { return r9_; }
  bool passId() const { return passId_; }
  bool isSphericalPhoton() const { return sphericalPhoton_; }
  
  TLorentzVector p4(float vtxx, float vtxy, float vtxz) const;

  void setEnergy(float energy) {energy_=energy; };
  void setCorrEnergy(float corrEnergy) {corrEnergy_=corrEnergy; };
  void setCorrEnergyErr(float corrEnergyErr) {corrEnergyErr_=corrEnergyErr; };
  void setCaloPosition(const TVector3 & caloPosition) { caloPosition_=caloPosition; };
  void setR9(float r9) {r9_=r9; };
  void setDet(int det) { iDet_=det; };
  void setSphericalPhoton(bool issph){sphericalPhoton_= issph;};

protected:

  TVector3 caloPosition_;
  float energy_;
  float corrEnergy_;
  float corrEnergyErr_;
  int iDet_;
  float r9_;
  bool passId_;
  bool sphericalPhoton_;
  void copy_(const PhotonReducedInfo &);
};


#endif


