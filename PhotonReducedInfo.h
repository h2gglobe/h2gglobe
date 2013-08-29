#ifndef PhotonReducedInfo_h
#define PhotonReducedInfo_h

#include "TVector3.h"
#include "TLorentzVector.h"
#include <iostream>

#include "BaseSmearer.h"

class PhotonReducedInfo
{
public:
	
	// PhotonReducedInfo();
  PhotonReducedInfo(const TVector3 & caloPosition, float energy, float corrEnergy, int iDet, float r9, bool passId, float corrEnergyErr=0);

  
  // PhotonReducedInfo(const PhotonReducedInfo &obj){copy_(obj);}
  // PhotonReducedInfo & operator= (const PhotonReducedInfo &obj);

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

  unsigned int nSmearingSeeds() { return smearingSeeds_.size(); }
  int smearingSeed(int ised=0) { return smearingSeeds_[ised];  };
  void addSmearingSeed(int seed) { return smearingSeeds_.push_back(seed);  };
  void dump();
  void reset();
  
  void cacheVal(int id, const BaseSmearer * smearer, float val) { 
	  cache_.resize(BaseSmearer::nRegisteredSmerers(),std::make_pair((const BaseSmearer *)0,0.)); 
	  cache_[id] = std::make_pair(smearer,val); 
  };
  bool hasCachedVal(int id) { return cache_.size() > id && cache_[id].first != 0; };
  const std::pair<const BaseSmearer *, float> & cachedVal(int id) { return cache_[id]; };
  void dumpCache() {
	  std::cout << "Photon " << this;
	  for(std::vector<std::pair<const BaseSmearer *, float> >::iterator it=cache_.begin(); it!=cache_.end(); ++it) {
		  std::cout << "\n " << (it->first == 0 ? " - " : it->first->name() ) << " " << it->second;
	  }
	  std::cout << std::endl;
  };

protected:

  TVector3 caloPosition_, rawCaloPosition_;
  float energy_;
  float corrEnergy_;
  float corrEnergyErr_;
  int iDet_;
  float r9_;
  bool passId_;
  bool sphericalPhoton_;
  float rawEnergy_, rawCorrEnergy_, rawR9_, rawCorrEnergyErr_;

  std::vector<int> smearingSeeds_;
  std::vector<std::pair<const BaseSmearer *, float> > cache_;
  /// void copy_(const PhotonReducedInfo &);
};


#endif


