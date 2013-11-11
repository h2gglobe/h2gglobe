#include "PhotonIdMVASmearer.h"
#include "PhotonReducedInfo.h"
#include "LoopAll.h"

PhotonIdMVASmearer::PhotonIdMVASmearer(LoopAll * lo, float syst_size, const std::string &type) :
	l_(lo), syst_size_(syst_size), type_(type)
{
}

PhotonIdMVASmearer::~PhotonIdMVASmearer()
{
	
}

bool PhotonIdMVASmearer::smearPhoton(PhotonReducedInfo &info, float & weight, int run, float syst_shift) const
{
	int ipho = info.iPho();
	/// if( syst_shift != 0. ) {
	/// 	std::cout << "PhotonIdMVASmearer::smearPhoton " << ipho << " " << syst_shift << std::endl;
	/// }
	
	std::vector<int> vtxs;
	for(int idipho=0; idipho<l_->dipho_n; ++idipho) {
		if(l_->dipho_leadind[idipho] == ipho || l_->dipho_subleadind[idipho] == ipho ) {
			vtxs.push_back(l_->dipho_vtxind[idipho]);
		}
	}
	
	l_->pho_idmva_cached = false;
	for(int ivtx=0; ivtx<l_->vtx_std_n; ++ivtx) {
		if( find(vtxs.begin(),vtxs.end(),ivtx) != vtxs.end() ) {
			TVector3 * vtx = (TVector3*)l_->vtx_std_xyz->At(ivtx);
			TLorentzVector p4 = info.p4(vtx->X(),vtx->Y(),vtx->Z());
			l_->pho_idmva[ipho][ivtx] = l_->photonIDMVA(ipho, ivtx, p4, type_.c_str()) + syst_shift*syst_size_;
			/// std::cout << "PhotonIdMVASmearer::smearPhoton " << ipho << " " << ivtx <<  " "  << l_->pho_idmva[ipho][ivtx] << std::endl;
		} else {
			l_->pho_idmva[ipho][ivtx] = -2;
		}
	}
	l_->pho_idmva_cached = true;
	
	return true;
}
