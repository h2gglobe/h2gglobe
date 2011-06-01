#ifndef __EFFICIENCYSMEARER__
#define __EFFICIENCYSMEARER__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include "TFile.h"
#include "TGraphAsymmErrors.h"

class PhotonReducedInfo;
class TRandom3;

// ------------------------------------------------------------------------------------
class EfficiencySmearer : public BaseSmearer
{
public:

  struct efficiencySmearingParameters
  {
	  int n_categories;
	  std::string categoryType;
	  std::string parameterSetName;
	  
	  typedef std::map<std::string,float> parameterMap;
	  typedef std::map<std::string,float>::iterator parameterMapIt;
	  typedef std::map<std::string,float>::const_iterator parameterMapConstIt;
	  
          std::string efficiency_file;
  };
  
  EfficiencySmearer(const efficiencySmearingParameters& par) ;
  virtual ~EfficiencySmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearPhoton( PhotonReducedInfo & pho, float & weight, float syst_shift=0. ) const;
  
  void name(const std::string & x) { name_ = x; };

  void setEffName(std::string x) { effName_ =x; };
  void doPhoId(bool x) { doPhoId_ = true; };

  bool init();

  efficiencySmearingParameters  myParameters_;
  
 protected:
  std::string photonCategory(PhotonReducedInfo &) const;

  double getWeight(double pt, std::string theCategory, float syst_shift) const;
  
  std::string   name_;
  TRandom3     *rgen_;
  std::string   effName_;
  bool doPhoId_;
  TFile        *theEfficiencyFile_; 
  std::map<std::string,TGraphAsymmErrors*> smearing_eff_graph_;
};

#endif
