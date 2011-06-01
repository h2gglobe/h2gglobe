#ifndef __ENERGYSMEARER__
#define __ENERGYSMEARER__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include "TFile.h"
#include "TGraphAsymmErrors.h"

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

          std::string efficiency_file;
  };
  
  EnergySmearer(const energySmearingParameters& par) ;
  virtual ~EnergySmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearPhoton(PhotonReducedInfo &, float & weight, float syst_shift) const;
  
  void name(const std::string & x) { name_ = x; };
  
  void scaleOrSmear(bool x) { scaleOrSmear_=x; }; 

  void doEnergy(bool x) { doEnergy_=x; }; 
   
  void doEfficiencies(bool x) { doEfficiencies_=x; }; 
  
  void setEffName(std::string x) { effName_ =x; };

  bool initEfficiency();

  energySmearingParameters  myParameters_;
  
 protected:
  bool doEnergy_, scaleOrSmear_, doEfficiencies_;

  std::string photonCategory(PhotonReducedInfo &) const;

  double getWeight(double pt, std::string theCategory, float syst_shift) const;
  
  std::string   name_;
  TRandom3     *rgen_;
  std::string   effName_;
  TFile        *theEfficiencyFile_; 
  std::map<std::string,TGraphAsymmErrors*> smearing_eff_graph_;
};

#endif
