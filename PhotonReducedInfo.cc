#include "PhotonReducedInfo.h"


using namespace std;


PhotonReducedInfo::PhotonReducedInfo(const TVector3 & caloPosition, float energy,int iDet, float r9, bool passId) :
	caloPosition_(caloPosition), energy_(energy), iDet_(iDet), r9_(r9), passId_(passId)
{
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


