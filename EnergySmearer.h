#ifndef __ENERGYSMEARER__
#define __ENERGYSMEARER__

#include "BaseSmearer.h"
#include <string>
#include <map>

class PhotonReducedInfo;
class TRandom3;

// ------------------------------------------------------------------------------------
class EnergySmearer : public BaseSmearer
{
public:

  struct energySmearingParameters
  {
    int n_categories;
    std::string categoryType;
    std::string parameterSetName;

    typedef std::map<std::string,float> parameterMap;
    typedef std::map<std::string,float>::iterator parameterMapIt;
    typedef std::map<std::string,float>::const_iterator parameterMapConstIt;
 
    // Scale offset and smearing error should be espressed as a relative value
    // Example: scale_offset["EB"]=1.002 , smearing_sigma["EB"]=0.01  

    std::map<std::string,float> scale_offset;
    std::map<std::string,float> scale_offset_error;

    std::map<std::string,float> smearing_sigma;
    std::map<std::string,float> smearing_sigma_error;
  };

  EnergySmearer(const energySmearingParameters& par) ;
  virtual ~EnergySmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearPhoton(PhotonReducedInfo &) const;
  
  energySmearingParameters  myParameters_;
  
 protected:
  std::string photonCategory(PhotonReducedInfo &) const;
  
  std::string name_;
  TRandom3* rgen_;
};

#endif
