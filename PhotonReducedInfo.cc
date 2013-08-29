#include "PhotonReducedInfo.h"


using namespace std;


PhotonReducedInfo::PhotonReducedInfo(const TVector3 & caloPosition, float energy, float corrEnergy, int iDet, float r9, bool passId, float corrEnergyErr) :
	caloPosition_(caloPosition), rawCaloPosition_(caloPosition), energy_(energy), corrEnergy_(corrEnergy), iDet_(iDet), r9_(r9), passId_(passId), corrEnergyErr_(corrEnergyErr),sphericalPhoton_(false),
	rawEnergy_(energy), rawCorrEnergy_(corrEnergy), rawR9_(r9),rawCorrEnergyErr_(corrEnergyErr)
{
}

//// PhotonReducedInfo& PhotonReducedInfo::operator=(const PhotonReducedInfo &obj){
//// 	copy_(obj);
//// 	return *this;
//// }
//// 
//// void PhotonReducedInfo::copy_(const PhotonReducedInfo &obj){
//// 	
//// 	caloPosition_=obj.caloPosition_;
//// 	energy_=obj.energy_;
//// 	corrEnergy_=obj.corrEnergy_;
//// 	corrEnergyErr_=obj.corrEnergyErr_;
//// 	iDet_=obj.iDet_;
//// 	r9_=obj.r9_;
//// 	passId_=obj.passId_;
//// 	sphericalPhoton_=obj.sphericalPhoton_;
//// 	smearingSeeds_ = obj.smearingSeeds_;
//// 	
//// 	rawCaloPosition_ =obj.rawCaloPosition_ ;
//// 	rawEnergy_       =obj.rawEnergy_       ;
//// 	rawCorrEnergy_   =obj.rawCorrEnergy_   ;
//// 	rawR9_           =obj.rawR9_           ;
//// 	rawCorrEnergyErr_=obj.rawCorrEnergyErr_;
//// 
//// }

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector PhotonReducedInfo::p4(float vtxx, float vtxy, float vtxz) const
{
	TVector3 vPos(vtxx,vtxy,vtxz);
	TVector3 direction = caloPosition() - vPos;
	TVector3 p = direction.Unit() * energy();
	TLorentzVector p4(p.x(),p.y(),p.z(),energy());
	return p4;
}

void PhotonReducedInfo::dump(){
	std::cout << " --- Dumping Reduced Photon Info --- " << std::endl;
	std::cout << "SC Position - " <<std::endl;
	caloPosition_.Dump();
	std::cout << "Current Energy - " << energy_<<std::endl;
	std::cout << "Corrected Energy - " << corrEnergy_ <<std::endl;
	std::cout << "Resolution - " << corrEnergyErr_<<std::endl;
	std::cout << "Calo IDET - " <<iDet_<<std::endl;
	std::cout << "r9 - " << r9_<<std::endl;
	std::cout << "Pass ID - " << passId_<<std::endl;
	std::cout << "Special Photon - " << sphericalPhoton_<<std::endl;
	std::cout << " ----------------------------------- " << std::endl;

}

void PhotonReducedInfo::reset()
{
	caloPosition_  = rawCaloPosition_ ;
	energy_        = rawEnergy_       ;
	corrEnergy_    = rawCorrEnergy_   ;
	r9_            = rawR9_           ;
	corrEnergyErr_ = rawCorrEnergyErr_;
}
