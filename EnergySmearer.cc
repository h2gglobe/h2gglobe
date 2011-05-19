#include "EnergySmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

EnergySmearer::EnergySmearer(const energySmearingParameters& par) : myParameters_(par)
{
  rgen_ = new TRandom3(0);
  name_="EnergySmearer_"+ par.categoryType + "_" + par.parameterSetName;
  //Checking consistency of input parameters
  assert( myParameters_.n_categories == myParameters_.scale_offset.size() );
  assert( myParameters_.n_categories == myParameters_.scale_offset_error.size() );
  assert( myParameters_.n_categories == myParameters_.smearing_sigma.size() );
  assert( myParameters_.n_categories == myParameters_.smearing_sigma_error.size() );
  assert( ( myParameters_.categoryType == "EBEE" && myParameters_.n_categories == 2) ||
	  ( myParameters_.categoryType == "2CatR9_EBEE" && myParameters_.n_categories == 4) );
}

EnergySmearer::~EnergySmearer()
{
  delete rgen_;
}

std::string EnergySmearer::photonCategory(PhotonReducedInfo & aPho) const
{
  std::string myCategory="";
  if (myParameters_.categoryType=="2CatR9_EBEE")
    {
      if (aPho.iDet()==1)
	myCategory+="EB";
      else
	myCategory+="EE";
      
      if (aPho.r9()>=0.94)
	myCategory+="HighR9";
      else
	myCategory+="LowR9";
    }
  else if (myParameters_.categoryType=="EBEE")
    {
      if (aPho.iDet()==1)
	myCategory+="EB";
      else
	myCategory+="EE";
    } 
  else
    {
      std::cout << "Unknown categorization. No category name is returned" << std::endl;
    }
  return myCategory;
}

bool EnergySmearer::smearPhoton(PhotonReducedInfo & aPho) const
{
  std::string category=photonCategory(aPho);
  
  if (category == "")
    {
      std::cout << "No category has been found associated with this photon. Giving Up" << std::endl;
      return false;
    }

  energySmearingParameters::parameterMapConstIt it=myParameters_.scale_offset.find(category);

  if ( it == myParameters_.scale_offset.end())
    {
      std::cout << "Category was not found in the configuration. Giving Up" << std::endl;
      return false;
    }
  
  float scale_offset=it->second;
  //  float scale_offset_error=myParameters_.scale_offset_error[category];
  
  float smearing_sigma=myParameters_.smearing_sigma.find(category)->second;
  //  float smearing_sigma_error=myParameters_.smearing_sigma_error[category];

  //now smearing energy with correct value
  float newEnergy= aPho.energy() * rgen_->Gaus(scale_offset,smearing_sigma);
  aPho.setEnergy(newEnergy);

  return true;
}
