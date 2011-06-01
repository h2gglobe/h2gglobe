#include "EfficiencySmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

EfficiencySmearer::EfficiencySmearer(const efficiencySmearingParameters& par) : myParameters_(par), doPhoId_(false)
{
  rgen_ = new TRandom3(0);
  name_="EfficiencySmearer_"+ par.categoryType + "_" + par.parameterSetName;
  
  assert( ( myParameters_.categoryType == "2CatR9_EBEE" && myParameters_.n_categories == 4) );
}

EfficiencySmearer::~EfficiencySmearer()
{
  delete rgen_;
}

std::string EfficiencySmearer::photonCategory(PhotonReducedInfo & aPho) const
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

bool EfficiencySmearer::smearPhoton(PhotonReducedInfo & aPho, float & weight, float syst_shift) const
{
  std::string category=photonCategory(aPho);
    
  if (category == "")
    {
      std::cout << "No category has been found associated with this photon. G<iving Up" << std::endl;
      return false;
    }

  /////////////////////// changing weigh of photon according to efficiencies ///////////////////////////////////////////
  assert( !smearing_eff_graph_.empty() );
  if( ! doPhoId_ || aPho.passId() ) {
	  weight = getWeight( ( aPho.energy() / cosh(aPho.caloPosition().PseudoRapidity()) ) ,category, syst_shift);
  }
  
  return true;
}



bool EfficiencySmearer::init() 
{

  // if map is not empty, yuo're initilized and happy..
  if( !smearing_eff_graph_.empty() ){
     std:cout << "initialization of efficiencies already done; proceed with usage. " << std::endl;
    return true;
  }

  //otherwise, get smearing functions from file and set up map
  std::cout << "\n>>>initializing one efficiency for photon re-weighting; " <<  std::endl;
  
  // do basic sanity checks first
  if( effName_.empty()){
    std::cout << "you're initialzinfg reweighting for efficiency but effName_ is empty" << std::endl;  assert(false); }
  if( myParameters_.efficiency_file.empty()){
    std::cout << "you're initialzinfg reweighting for efficiency: " << effName_  << " but input file with TGraphErrors is not specified; doing nothing. " << std::endl;  assert(false); }
  
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


double EfficiencySmearer::getWeight(double pt, std::string theCategory, float syst_shift) const
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
