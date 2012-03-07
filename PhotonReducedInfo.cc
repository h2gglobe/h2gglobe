#include "PhotonReducedInfo.h"


using namespace std;


PhotonReducedInfo::PhotonReducedInfo(const TVector3 & caloPosition, float energy, float corrEnergy, int iDet, float r9, bool passId, float corrEnergyErr) :
  caloPosition_(caloPosition), energy_(energy), corrEnergy_(corrEnergy), iDet_(iDet), r9_(r9), passId_(passId), corrEnergyErr_(corrEnergyErr),sphericalPhoton_(false)
{
}

PhotonReducedInfo& PhotonReducedInfo::operator=(const PhotonReducedInfo &obj){
	copy_(obj);
	return *this;
}

void PhotonReducedInfo::copy_(const PhotonReducedInfo &obj){
	
	caloPosition_=obj.caloPosition_;
	energy_=obj.energy_;
	corrEnergy_=obj.corrEnergy_;
	corrEnergyErr_=obj.corrEnergyErr_;
	iDet_=obj.iDet_;
	r9_=obj.r9_;
	passId_=obj.passId_;
	sphericalPhoton_=obj.sphericalPhoton_;

}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector PhotonReducedInfo::p4(float vtxx, float vtxy, float vtxz) const
{
	TVector3 vPos(vtxx,vtxy,vtxz);
	TVector3 direction = caloPosition() - vPos;
	TVector3 p = direction.Unit() * energy();
	TLorentzVector p4(p.x(),p.y(),p.z(),energy());
	return p4;
}


