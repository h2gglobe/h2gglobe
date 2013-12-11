#include "DiPhoEfficiencySmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

DiPhoEfficiencySmearer::DiPhoEfficiencySmearer(const diPhoEfficiencySmearingParameters& par) : myParameters_(par)
{
  rgen_ = new TRandom3(0);
  name_="DiPhoEfficiencySmearer_"+ par.categoryType + "_" + par.parameterSetName;
  passFailWeights_=false;
  doVtxEff_       =false;
  doMvaIdEff_	  =false;
  idMVASystSize_  =par.idMVASystSize;
}

DiPhoEfficiencySmearer::~DiPhoEfficiencySmearer()
{
  delete rgen_;
}

bool DiPhoEfficiencySmearer::smearDiPhoton( TLorentzVector & p4, TVector3 & selVtx, float & weight, const int & category, 
					    const int & genMassPoint, const TVector3 & trueVtx, float & idMVA1, float & idMVA2, float syst_shift) const
{
  if (doMvaIdEff_){
      idMVA1 += syst_shift*idMVASystSize_; 
      idMVA2 += syst_shift*idMVASystSize_;
  } else {
    std::string cat=Form("cat%d", category);
  
    if (cat == "")
      {
       std::cout << effName_ <<" No category has been found associated with this diphoton. Giving Up" << std::endl;
       return false;
      }
  
    
    /////////////////////// changing weigh of photon according to efficiencies ///////////////////////////////////////////
    assert( ! smearing_eff_graph_.empty() );
  
    if( doVtxEff_ ) {
      if( (selVtx - trueVtx).Mag() < 1. ) {
        cat += "_pass"; 
      }
      else {
        cat += "_fail";
        syst_shift *=-1;      // shift for pass and fail need to be in opposite directions
      }
    }
  
    weight = getWeight( p4.Pt(), cat, syst_shift );
    weight = weight<0 ? 0 : weight;
  }
  
  return true;
}



bool DiPhoEfficiencySmearer::init() 
{

  // if map is not empty, yuo're initilized and happy..
  if( !smearing_eff_graph_.empty() ){
  std:cout << "initialization of DI-photon efficiency smearer " << effName_ << " already done; proceed with usage. " << std::endl;
    return true;
  }
  if( doVtxEff_ ) { passFailWeights_ = true; }
  else            { passFailWeights_ = false;}

  //otherwise, get smearing functions from file and set up map
  std::cout << "\n>>>initializing one efficiency for DI-photon re-weighting: " << effName_ <<  std::endl;
  
  // do basic sanity checks first
  if( effName_.empty()){
    std::cout << effName_ << " - you're initializing reweighting for efficiency but effName_ is empty" << std::endl;  assert(false); }
  if( myParameters_.efficiency_file.empty()){
    std::cout <<  effName_ <<  "- you're initializing reweighting for efficiency: " << effName_  << " but input file with TGraphErrors is not specified; doing nothing. " << std::endl;  assert(false); }
  
  theDiPhoEfficiencyFile_ = TFile::Open(myParameters_.efficiency_file.c_str());

  // initialize formulas for the di-photon categories; 
  for( int ii=0; ii<myParameters_.n_categories; ++ii ) {
    std::string cat = Form("cat%d", ii);
    //std::cout << "GF for efficiency: " << effName_ << "  passFailWeights_: " <<  passFailWeights_ << "  doVtxEff_: " <<  doVtxEff_ << std::endl;
    if( passFailWeights_ ) {
      // std::cout << "GF cloning: " << ( (effName_+"_"+cat+"_pass").c_str() ) << std::endl;
      smearing_eff_graph_[cat+"_pass"]=(TGraphAsymmErrors*) theDiPhoEfficiencyFile_->Get((effName_+"_"+cat+"_pass").c_str())->Clone();
      smearing_eff_graph_[cat+"_fail"]=(TGraphAsymmErrors*) theDiPhoEfficiencyFile_->Get((effName_+"_"+cat+"_fail").c_str())->Clone();
      std::cerr << "DiPhoEfficiencySmearerc " << cat+"_pass" << std::endl;
    } else {
      //std::cout << "GF cloning: " << ( (effName_+"_"+cat).c_str() ) << std::endl;
      smearing_eff_graph_[cat]=(TGraphAsymmErrors*) theDiPhoEfficiencyFile_->Get((effName_+"_"+cat).c_str())->Clone();
    }
  }

  theDiPhoEfficiencyFile_->Close();
  return true;

}


double DiPhoEfficiencySmearer::getWeight(double pt, std::string theCategory, float syst_shift) const
{
    static int nwarnings = 10;
  std::map<std::string,TGraphAsymmErrors*>::const_iterator theIter = smearing_eff_graph_.find(theCategory);
  if( theIter != smearing_eff_graph_.end()  ) {

    // determine the pair of bins between which  you interpolate
    int numPoints = ( theIter->second )->GetN();
    double x, y, xPrevious;
    int myBin = -1;
    xPrevious =-1e9;
    for (int bin=0; bin<numPoints; bin++ ){
      ( theIter->second )->GetPoint(bin, x, y);
      assert( xPrevious < x ) ;     // points in TGraphAsymmErrors must be in increasing order
      xPrevious=x;
      if(pt > x) {
	myBin = bin; }
      else break;
    }
    int binLow, binHigh; bool atBoundary(false);
    if      (myBin == -1)               {binHigh = 0; binLow=0; atBoundary=true;}
    else if (myBin == (numPoints-1))    {binHigh = numPoints-1; binLow=numPoints-1; atBoundary=true;}
    else                                {binLow=myBin; binHigh=myBin+1;}

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
    //           if you're NOT at the boundaris of TGraphAsymmErrors, linearly interpolate values and errors
    if(!atBoundary) {
      theWeight = yLow + (yHigh-yLow) / (xHigh-xLow) * (pt-xLow);
      theError  = theErrorLow + (theErrorHigh-theErrorLow) / (xHigh-xLow) * (pt-xLow); 
    } else { 
	//  if instead you ARE at the boundaris of TGraphAsymmErrors, collapse on first or last of its points  
	if(myBin == (numPoints-1)) {
	  theWeight = yHigh; 	  theError  = theErrorHigh;	}
	else if (myBin == -1) 	  {
	  theWeight = yLow; 	  theError  = theErrorLow;  	}
	else  {   std::cout <<  effName_ << " ** you claim to be at boundaries of TGraphAsymmErrors but your not! This is a problem " << std::endl;}
    }
    float ret = theWeight + theError*syst_shift; 
    return ret;
  }
  else {
      std::cout <<  effName_ << "- category asked: " << theCategory << " was not found - which is a problem. Returning weight 1. " << std::endl;
      return 1.;  
  }
  
}
