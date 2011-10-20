#include "EnergySmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

EnergySmearer::EnergySmearer(const energySmearingParameters& par) : myParameters_(par), scaleOrSmear_(true), doCorrections_(false)
{
  rgen_ = new TRandom3(0);
  name_="EnergySmearer_"+ par.categoryType + "_" + par.parameterSetName;
  //Checking consistency of input parameters
  std::cerr << myParameters_.categoryType << " " <<  myParameters_.n_categories << std::endl;
  assert( myParameters_.n_categories == myParameters_.smearing_sigma.size() );
  assert( myParameters_.n_categories == myParameters_.smearing_sigma_error.size() );
  assert( ( myParameters_.categoryType == "EBEE" && myParameters_.n_categories == 2 ) ||
	  ( myParameters_.categoryType == "2CatR9_EBEE" && myParameters_.n_categories == 4 ) ||
	  ( myParameters_.categoryType == "2CatR9_EBEE_ByRun" && myParameters_.n_categories == 4 ) ||
	  ( myParameters_.categoryType == "2CatR9_EBEBm4EE" && myParameters_.n_categories == 6 ) ||
	  ( myParameters_.categoryType == "2CatR9_EBEBm4EE_ByRun" && myParameters_.n_categories == 6 )
	  );
  if( myParameters_.byRun ) {
    for(energySmearingParameters::eScaleVectorIt it=myParameters_.scale_offset_byrun.begin(); it!=myParameters_.scale_offset_byrun.end();
	++it ) {
      assert( myParameters_.n_categories == it->scale_offset.size() );
      assert( myParameters_.n_categories == it->scale_offset_error.size() );
    }
  } else {
    assert( myParameters_.n_categories == myParameters_.scale_offset.size() );
    assert( myParameters_.n_categories == myParameters_.scale_offset_error.size() );
  }
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
  else if (myParameters_.categoryType=="2CatR9_EBEBm4EE")
    {
      if (aPho.iDet()==1 && fabs(aPho.caloPosition().PseudoRapidity())      < 1.)
	myCategory+="EB";
      else if (aPho.iDet()==1 && fabs(aPho.caloPosition().PseudoRapidity()) > 1.)
	myCategory+="EBm4";
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


float EnergySmearer::getScaleOffset(int run, const std::string & category) const
{
  const std::map<std::string, float> * scale_offset =  &(myParameters_.scale_offset);
  
  if( myParameters_.byRun ) {
    scale_offset = &(find(myParameters_.scale_offset_byrun.begin(),myParameters_.scale_offset_byrun.end(),run)->scale_offset) ;
  }
  energySmearingParameters::parameterMapConstIt it=scale_offset->find(category);
  
  if ( it == scale_offset->end())
    {
      std::cout << "Category was not found in the configuration. Giving Up" << std::endl;
      return false;
    }
  
  return 1. + it->second;
  
}

bool EnergySmearer::smearPhoton(PhotonReducedInfo & aPho, float & weight, int run, float syst_shift) const
{
  std::string category=photonCategory(aPho);
    
  if (category == "")
    {
      std::cout << "No category has been found associated with this photon. G<iving Up" << std::endl;
      return false;
    }

  float scale_offset   = getScaleOffset(run, category);
  float smearing_sigma = myParameters_.smearing_sigma.find(category)->second;

  /////////////////////// smearing or re-scaling photon energy ///////////////////////////////////////////
  
  float newEnergy=aPho.energy();
  /////////////////////// apply MC-based photon energy corrections ///////////////////////////////////////////
  if (  doCorrections_ ) {
    // corrEnergy is the corrected photon energy
    newEnergy = aPho.corrEnergy() + syst_shift * myParameters_.corrRelErr * (aPho.corrEnergy() - aPho.energy());
  }

  if( scaleOrSmear_ ) {
	  scale_offset   += syst_shift * myParameters_.scale_offset_error.find(category)->second;
	  newEnergy *=  scale_offset;
  } else {
	  smearing_sigma += syst_shift * myParameters_.smearing_sigma_error.find(category)->second;
	  newEnergy *=  rgen_->Gaus(1.,smearing_sigma);
  }
  
  //std::cout << "doCorrections: " << doCorrections_ << " ene: " <<  aPho.energy() 
  //<< " corr ene: " <<  aPho.corrEnergy() << " newEne: " << newEnergy << " syst_shift: " 
  //<< syst_shift << std::endl; 

  assert( newEnergy != 0. );
  aPho.setEnergy(newEnergy);
  
  /////////////////////// changing weigh of photon according to efficiencies ///////////////////////////////////////////
  //////////////////////  if you're doing corrections, don't touch the weights ////////////////////////////////////////
  if(doEfficiencies_ && (!doCorrections_) ) {
    if( !smearing_eff_graph_.empty()  ){
      weight = getWeight( ( aPho.energy() / cosh(aPho.caloPosition().PseudoRapidity()) ) ,category, syst_shift);
    }
  }

  return true;
}



bool EnergySmearer::initEfficiency() 
{

  // if map is not empty, yuo're initilized and happy..
  if( !smearing_eff_graph_.empty() ){
     std:cout << "initialization of efficiencies already done; proceed with usage. " << std::endl;
    return true;
  }

  //otherwise, get smearing functions from file and set up map
  std::cout << "\n>>>initializing one efficiency for photon re-weighting; " <<  std::endl;
  
  // do basic sanity checks first
  if(doEnergy_){
    std::cout << "*** Initializing reweighting for efficiencies; energy smearing active TOO, do you want them both? " << std::endl;  return false; }
  if(!doEfficiencies_){
    std::cout << "*** Initializing reweighting for efficiencies - BUT doEfficiencies_ is set to false; doing nothing. " << std::endl;  return false; }
  if( effName_.empty()){
    std::cout << "you're initialzinfg reweighting for efficiency but effName_ is empty ; doing nothing. " << std::endl;  return false; }
  if( myParameters_.efficiency_file.empty()){
    std::cout << "you're initialzinfg reweighting for efficiency: " << effName_  << " but input file with TGraphErrors is not specified; doing nothing. " << std::endl;  return false; }
  
  theEfficiencyFile_ = new TFile(myParameters_.efficiency_file.c_str());

  // initialize formulas for the four categories; 
  std::string effTmpName; std::string photonCat; TGraphAsymmErrors *graphTmp, *graphClone;

  photonCat =  std::string("EBHighR9");
  effTmpName = effName_+std::string("_")+photonCat; graphTmp = (TGraphAsymmErrors*) theEfficiencyFile_->Get(effTmpName.c_str());   // smearing_eff_graph_[photonCat]=graphTmp;  
  graphClone=(TGraphAsymmErrors*)graphTmp->Clone();    smearing_eff_graph_[photonCat]=graphClone;    

  photonCat =  std::string("EBLowR9");
  effTmpName = effName_+std::string("_")+photonCat; graphTmp = (TGraphAsymmErrors*) theEfficiencyFile_->Get(effTmpName.c_str());   // smearing_eff_graph_[photonCat]=graphTmp;  
  graphClone=(TGraphAsymmErrors*)graphTmp->Clone();    smearing_eff_graph_[photonCat]=graphClone;    

  photonCat =  std::string("EEHighR9");
  effTmpName = effName_+std::string("_")+photonCat; graphTmp = (TGraphAsymmErrors*) theEfficiencyFile_->Get(effTmpName.c_str());   // smearing_eff_graph_[photonCat]=graphTmp;  
  graphClone=(TGraphAsymmErrors*)graphTmp->Clone();    smearing_eff_graph_[photonCat]=graphClone;    

  photonCat =  std::string("EELowR9");
  effTmpName = effName_+std::string("_")+photonCat; graphTmp = (TGraphAsymmErrors*) theEfficiencyFile_->Get(effTmpName.c_str());   // smearing_eff_graph_[photonCat]=graphTmp;  
  graphClone=(TGraphAsymmErrors*)graphTmp->Clone();    smearing_eff_graph_[photonCat]=graphClone;    

  delete  theEfficiencyFile_;
  return true;

}


double EnergySmearer::getWeight(double pt, std::string theCategory, float syst_shift) const
{
  std::map<std::string,TGraphAsymmErrors*>::const_iterator theIter = smearing_eff_graph_.find(theCategory);
  if( theIter != smearing_eff_graph_.end()  ) {

    // determine the pair of bins between which  you interpolate
    int numPoints = ( theIter->second )->GetN();
    double x, y;
    int myBin = -1;
    for (int bin=0; bin<numPoints; bin++ ){
      ( theIter->second )->GetPoint(bin, x, y);
      if(pt > x) {
	myBin = bin; }
      else break;
    }
    int binLow, binHigh; 
    if(myBin == -1)                      {binHigh = 0; binLow=0;}
    else if (myBin == (numPoints-1))     {binHigh = numPoints-1; binLow=numPoints-1;}
    else {binLow=myBin; binHigh=myBin+1;}


    // get hold of efficiency ratio and error at either points
    // low-high refer to the points ; up-down refers to the errors 
    double xLow, yLow;    double xHigh, yHigh;
    ( theIter->second )->GetPoint(binLow, xLow, yLow);
    ( theIter->second )->GetPoint(binHigh, xHigh, yHigh);

    double errLowYup    = ( theIter->second )->GetErrorYhigh(binLow);
    double errLowYdown  = ( theIter->second )->GetErrorYlow(binLow);
    double errHighYup   = ( theIter->second )->GetErrorYhigh(binHigh);
    double errHighYdown = ( theIter->second )->GetErrorYlow(binHigh);

    double theErrorLow, theErrorHigh;
    if(syst_shift>0) {theErrorLow = errLowYup;   theErrorHigh = errHighYup;}
    else             {theErrorLow = errLowYdown; theErrorHigh = errHighYdown;}
    
    double theWeight, theError;
    theWeight = yLow + (yHigh-yLow) / (xHigh-xLow) * (pt-xLow);
    theError  = theErrorLow + (theErrorHigh-theErrorLow) / (xHigh-xLow) * (pt-xLow);
     
    return  ( theWeight + (theError*syst_shift));
  }
  else {
    std::cout << "category asked: " << theCategory << " was not found - which is a problem. Returning weight 1. " << std::endl;
    return 1.;
  }
  
}
