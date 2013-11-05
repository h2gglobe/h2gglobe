#ifndef __PHOTONIDMVASMEARER__
#define __PHOTONIDMVASMEARER__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include <utility>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

class LoopAll;
class PhotonReducedInfo;
class TRandom3;


// ------------------------------------------------------------------------------------
class PhotonIdMVASmearer : public BaseSmearer
{
public:

	PhotonIdMVASmearer(LoopAll * lo, float sys_shift, const std::string & type);
	virtual ~PhotonIdMVASmearer();
  
	virtual const std::string & name() const { return name_; };
  
	virtual bool smearPhoton(PhotonReducedInfo &, float & weight, int run, float syst_shift) const;
	
	void name(const std::string & x) { name_ = x; };

protected:
	LoopAll * l_;
	float syst_size_;
	std::string name_, type_;
	
};


#endif
