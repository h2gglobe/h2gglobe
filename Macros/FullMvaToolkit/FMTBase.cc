#include <iomanip>
#include <fstream>
#include "TMath.h"
#include "boost/lexical_cast.hpp"

#include "FMTBase.h"

using namespace std;

FMTBase::FMTBase(int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, bool includeVBF, bool includeLEP, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int, vector<double> > > AllBinEdges, bool verbose):
	
	mHMinimum_(mHMinimum),
	mHMaximum_(mHMaximum),
	mHStep_(mHStep),
	massMin_(massMin),
	massMax_(massMax),
	nDataBins_(nDataBins),
	
	signalRegionWidth_(signalRegionWidth),
	sidebandWidth_(sidebandWidth),
	numberOfSidebands_(numberOfSidebands),
	numberOfSidebandsForAlgos_(numberOfSidebandsForAlgos),
	numberOfSidebandGaps_(numberOfSidebandGaps),
	massSidebandMin_(massSidebandMin),
	massSidebandMax_(massSidebandMax),

	includeVBF_(includeVBF),
	includeLEP_(includeLEP),

	systematics_(systematics),

	rederiveOptimizedBinEdges_(rederiveOptimizedBinEdges),

  verbose_(verbose)
{
	BinEdges_ = AllBinEdges[0] ;
	if (includeVBF) VBFBinEdges_ = AllBinEdges[1];
	if (includeLEP) LEPBinEdges_ = AllBinEdges[2];

}

void FMTBase::checkMCMass(int mass){
	vector<int> theMasses = getMCMasses();
	bool valid=false;
	for (vector<int>::iterator m=theMasses.begin(); m!=theMasses.end(); m++){
		if (mass==*m) valid=true;
	}
	if (!valid){
		cerr << "Mass " << mass << " is not a valid rebinning mass. You must pick a MC mass which can only be [";
		for (int i=0; i<getNumMCMasses()-1; i++){
			cerr << (getMCMasses())[i] << ",";
		}
		cerr << (getMCMasses())[getNumMCMasses()-1] << "]. Bailing out." << endl;
		exit(1);
	}
}

// FMTBase::getters
int FMTBase::getmHMinimum(){
	return mHMinimum_;
}
int FMTBase::getmHMaximum(){
	return mHMaximum_;
}
double FMTBase::getmHStep(){
	return mHStep_;
}
double FMTBase::getmassMin(){
	return massMin_;
}
double FMTBase::getmassMax(){
	return massMax_;
}
int FMTBase::getnDataBins(){
	return nDataBins_;
}

double FMTBase::getsignalRegionWidth(){
	return signalRegionWidth_;
}
double FMTBase::getsidebandWidth(){
	return sidebandWidth_;
}
int FMTBase::getnumberOfSidebands(){
	return numberOfSidebands_;
}
int FMTBase::getnumberOfSidebandsForAlgos(){
	return numberOfSidebandsForAlgos_;
}
int FMTBase::getnumberOfSidebandGaps(){
	return numberOfSidebandGaps_;
}
double FMTBase::getmassSidebandMin(){
	return massSidebandMin_;
}
double FMTBase::getmassSidebandMax(){
	return massSidebandMax_;
}

pair<int,int> FMTBase::getNsidebandsUandD(double mass){
    int nLower=numberOfSidebandsForAlgos_;
    int nHigher=numberOfSidebandsForAlgos_;
    bool dropLow=false;
    bool dropHigh=false;
    
    // loop lower sidebands
    for (int sidebandNumb=numberOfSidebandGaps_+1; sidebandNumb<=numberOfSidebandsForAlgos_+numberOfSidebandGaps_; sidebandNumb++){
      double hypothesisModifier = (1.-sidebandWidth_)/(1.+sidebandWidth_);
      double sidebandCenter = mass*(1.-signalRegionWidth_)/(1.+sidebandWidth_)*TMath::Power(hypothesisModifier,sidebandNumb-1);
      double sidebandHigh = sidebandCenter*(1.+sidebandWidth_);
      double sidebandLow = sidebandCenter*(1.-sidebandWidth_);
      if (sidebandLow<massSidebandMin_){
        nLower--;
        nHigher++;
        dropLow=true;
      }
    }

    // loop higher sidebands
    for (int sidebandNumb=numberOfSidebandGaps_+1; sidebandNumb<=numberOfSidebandsForAlgos_+numberOfSidebandGaps_; sidebandNumb++){
      double hypothesisModifier = (1.+sidebandWidth_)/(1.-sidebandWidth_);
      double sidebandCenter = mass*(1.+signalRegionWidth_)/(1.-sidebandWidth_)*TMath::Power(hypothesisModifier,sidebandNumb-1);
      double sidebandHigh = sidebandCenter*(1.+sidebandWidth_);
      double sidebandLow = sidebandCenter*(1.-sidebandWidth_);
      if (sidebandHigh>massSidebandMax_){
        nHigher--;
        nLower++;
        dropHigh=true;
      }
    }
    if (dropLow && dropHigh) cerr << "There must be a problem as both low and high sidebands are being dropped" << endl;
    if (nLower+nHigher!=2*numberOfSidebandsForAlgos_) cerr << "There must be a problem as the total number of sidebands should equal " << 2*numberOfSidebandsForAlgos_ << " not " << nLower+nHigher << endl;

    cout << "Mass " << mass << " has " << nLower << " lower sidebands and " << nHigher << " higher sidebands" << endl;
   
    pair<int,int> result(nLower,nHigher);
    return result;
}

vector<double> FMTBase::getLowerSidebandEdges(double mass){
  
	vector<double> lowerEdges;
  int nLower = (getNsidebandsUandD(mass)).first;

  for (int sidebandNumb=numberOfSidebandGaps_+1; sidebandNumb<=nLower+numberOfSidebandGaps_; sidebandNumb++){
    double hypothesisModifier = (1.-sidebandWidth_)/(1.+sidebandWidth_);
    double sidebandCenter = mass*(1.-signalRegionWidth_)/(1.+sidebandWidth_)*TMath::Power(hypothesisModifier,sidebandNumb-1);
    
    if (sidebandNumb==numberOfSidebandGaps_+1) lowerEdges.push_back(sidebandCenter*(1.+sidebandWidth_));
    lowerEdges.push_back(sidebandCenter*(1.-sidebandWidth_));
  }
  return lowerEdges;
}

vector<double> FMTBase::getUpperSidebandEdges(double mass){
  
  vector<double> upperEdges;
  int nHigher = (getNsidebandsUandD(mass)).second;

  for (int sidebandNumb=numberOfSidebandGaps_+1; sidebandNumb<=nHigher+numberOfSidebandGaps_; sidebandNumb++){
    double hypothesisModifier = (1.+sidebandWidth_)/(1.-sidebandWidth_);
    double sidebandCenter = mass*(1.+signalRegionWidth_)/(1.-sidebandWidth_)*TMath::Power(hypothesisModifier,sidebandNumb-1);

    if (sidebandNumb==numberOfSidebandGaps_+1) upperEdges.push_back(sidebandCenter*(1.-sidebandWidth_));
    upperEdges.push_back(sidebandCenter*(1.+sidebandWidth_));
  }
  return upperEdges;
}

bool FMTBase::getincludeVBF(){
	return includeVBF_;
}
bool FMTBase::getincludeLEP(){
	return includeLEP_;
}
const int FMTBase::getNcats(){
	int ncats=1;
	if (includeVBF_) ncats++;
	if (includeLEP_) ncats++;
	return ncats;
}

vector<string> FMTBase::getsystematics(){
	return systematics_;
}
int FMTBase::getNsystematics(){
	return systematics_.size();
}

bool FMTBase::getrederiveOptimizedBinEdges(){
	return rederiveOptimizedBinEdges_;
}

vector<map<int, vector<double> > > FMTBase::getAllBinEdges(){
	vector<map<int, vector<double> > > theEdges;
	theEdges.push_back(getBinEdges());
	theEdges.push_back(getVBFBinEdges());
	theEdges.push_back(getLEPBinEdges());
	return theEdges;
}

map<int, vector<double> > FMTBase::getBinEdges(){
	if (BinEdges_.size() != getNumMCMasses()){
		cerr << "WARNING - the number of bin edges (" << BinEdges_.size() << ") does not match the number of MC masses (" << getNumMCMasses() << ")" << endl;
	}
	return BinEdges_;
}

map<int, vector<double> > FMTBase::getVBFBinEdges(){
	if (VBFBinEdges_.size() != getNumMCMasses()){
		cerr << "WARNING - the number of bin edges (" << VBFBinEdges_.size() << ") does not match the number of MC masses (" << getNumMCMasses() << ")" << endl;
	}
	return VBFBinEdges_;
}

map<int, vector<double> > FMTBase::getLEPBinEdges(){
	if (LEPBinEdges_.size() != getNumMCMasses()){
		cerr << "WARNING - the number of bin edges (" << LEPBinEdges_.size() << ") does not match the number of MC masses (" << getNumMCMasses() << ")" << endl;
	}
	return LEPBinEdges_;
}

vector<vector<double> > FMTBase::getAllBinEdges(int mass){
	vector<vector<double> > theEdges;
	theEdges.push_back(getBinEdges(mass));
	theEdges.push_back(getVBFBinEdges(mass));
	theEdges.push_back(getLEPBinEdges(mass));
	return theEdges;
}

vector<double> FMTBase::getBinEdges(int mass){
	if (BinEdges_.find(mass)==BinEdges_.end()){
		cerr << "WARNING -- FMTBase::getBinEdges() -- bin edges not found for mass " << mass << ". Bailing out" << endl;
		exit(1);
	}
	else{ 
		vector<double> theEdges = (BinEdges_.find(mass))->second;
		return theEdges;
	}
}

vector<double> FMTBase::getVBFBinEdges(int mass){
	if (VBFBinEdges_.find(mass)==VBFBinEdges_.end()){
		cerr << "WARNING -- FMTBase::getVBFBinEdges() -- bin edges not found for mass " << mass << ". Bailing out" << endl;
		exit(1);
	}
	else{ 
		vector<double> theEdges = (VBFBinEdges_.find(mass))->second;
		return theEdges;
	}
}

vector<double> FMTBase::getLEPBinEdges(int mass){
	if (LEPBinEdges_.find(mass)==LEPBinEdges_.end()){
		cerr << "WARNING -- FMTBase::getLEPBinEdges() -- bin edges not found for mass " << mass << ". Bailing out" << endl;
		exit(1);
	}
	else{ 
		vector<double> theEdges = (LEPBinEdges_.find(mass))->second;
		return theEdges;
	}
}

vector<double> FMTBase::getMHMasses(int mass){
	vector<double> theMasses;
	checkMCMass(mass);
	double m_low=0;
	double m_high=0;
	if (mass==mHMinimum_){
		m_low=0.;
		m_high=2.5;
	}
	else if (mass==mHMaximum_){
		m_low=-2.5;
		m_high=mHStep_;
		if (mass==150) m_low=-5.0;
	}
	else if (mass==140){
		m_low=-2.5;
		m_high=5.0;
	}
	else {
		m_low=-2.5;
		m_high=2.5;
	}
	for (double m=m_low; m<m_high-(mHStep_/2.); m+=mHStep_){
		theMasses.push_back(mass+m);
	}
	return theMasses;
}

vector<int> FMTBase::getUandDMCMasses(int mass){
  vector<int> theMasses;
  if (mass==140){
    if (135>=mHMinimum_) theMasses.push_back(135);
    theMasses.push_back(140);
    if (150<=mHMaximum_) theMasses.push_back(150);
  }
  else {
    if (mass-5>=mHMinimum_) theMasses.push_back(mass-5);
    theMasses.push_back(mass);
    if (mass+5<=mHMaximum_) theMasses.push_back(mass+5);
  }
  return theMasses;
}

pair<int,int> FMTBase::getInterpMasses(double mass){
  
  double bestdiff=1000.;
  double nbestdiff=1000.;
  int nearest;
  int nnearest;
  vector<int> mcMasses = getMCMasses();
  for (vector<int>::iterator mcM=mcMasses.begin(); mcM!=mcMasses.end(); mcM++){
    double diff = TMath::Abs((mass+0.001)-boost::lexical_cast<double>(*mcM));
    if (diff<bestdiff){
      bestdiff=diff;
      nearest=*mcM;
    }
  }
  for (vector<int>::iterator mcM=mcMasses.begin(); mcM!=mcMasses.end(); mcM++){
    double diff = TMath::Abs((mass+0.001)-boost::lexical_cast<double>(*mcM));
    if (diff>bestdiff && diff<nbestdiff){
      nbestdiff=diff;
      nnearest=*mcM;
    }
  }
  pair<int,int> nearestPair(nearest,nnearest);
  return nearestPair;
}

vector<int> FMTBase::getMCMasses(){
	vector<int> theMasses;
	for (int mH=mHMinimum_; mH<=mHMaximum_; mH+=5){
		if (mH==145) continue;
		theMasses.push_back(mH);
	}
	return theMasses;
}
int FMTBase::getNumMCMasses(){
	vector<int> theMasses = getMCMasses();
	return theMasses.size();
}

// FMTBase::setters
void FMTBase::setmHMinimum(int mHMinimum){
	mHMinimum_=mHMinimum;
}
void FMTBase::setmHMaximum(int mHMaximum){
	mHMaximum_=mHMaximum;
}
void FMTBase::setmHStep(double mHStep){
	mHStep_=mHStep;
}
void FMTBase::setmassMin(double massMin){
	massMin_=massMin;
}
void FMTBase::setmassMax(double massMax){
	massMax_=massMax;
}
void FMTBase::setnDataBins(int nDataBins){
	nDataBins_=nDataBins;
}

void FMTBase::setsignalRegionWidth(double signalRegionWidth){
	signalRegionWidth_=signalRegionWidth;
}
void FMTBase::setsidebandWidth(double sidebandWidth){
	sidebandWidth_=sidebandWidth;
}
void FMTBase::setnumberOfSidebands(int numberOfSidebands){
	numberOfSidebands_=numberOfSidebands;
}
void FMTBase::setnumberOfSidebandsForAlgos(int numberOfSidebandsForAlgos){
	numberOfSidebandsForAlgos_=numberOfSidebandsForAlgos;
}
void FMTBase::setnumberOfSidebandGaps(int numberOfSidebandGaps){
	numberOfSidebandGaps_=numberOfSidebandGaps;
}
void FMTBase::setmassSidebandMin(double massSidebandMin){
	massSidebandMin_=massSidebandMin;
}
void FMTBase::setmassSidebandMax(double massSidebandMax){
	massSidebandMax_=massSidebandMax;
}

void FMTBase::setincludeVBF(bool includeVBF){
	includeVBF_=includeVBF;
}
void FMTBase::setincludeLEP(bool includeLEP){
	includeLEP_=includeLEP;
}

void FMTBase::setsystematics(vector<string> systematics){
	systematics_=systematics;
}
void FMTBase::setsystematic(string systematic){
	if (find(systematics_.begin(), systematics_.end(), systematic)!=systematics_.end()){
		// already in
		return;
	}
	else {
		systematics_.push_back(systematic);
	}
}

void FMTBase::setrederiveOptimizedBinEdges(bool rederiveOptimizedBinEdges){
	rederiveOptimizedBinEdges_=rederiveOptimizedBinEdges;
}

void FMTBase::setAllBinEdges(vector<map<int,vector<double> > > theEdges){
	setBinEdges(*(theEdges.begin()));
	setVBFBinEdges(*(theEdges.begin()+1));
	setLEPBinEdges(*(theEdges.begin()+2));
}

void FMTBase::setBinEdges(int mass, vector<double> BinEdges){
	if (BinEdges_.find(mass)==BinEdges_.end()){
		BinEdges_.insert(pair<int,vector<double> >(mass,BinEdges));
	}
	else {
		BinEdges_[mass]=BinEdges;
	}
}
void FMTBase::setVBFBinEdges(int mass, vector<double> VBFBinEdges){
	if (VBFBinEdges_.find(mass)==VBFBinEdges_.end()){
		VBFBinEdges_.insert(pair<int,vector<double> >(mass,VBFBinEdges));
	}
	else {
		VBFBinEdges_[mass]=VBFBinEdges;
	}
}
void FMTBase::setLEPBinEdges(int mass, vector<double> LEPBinEdges){
	if (LEPBinEdges_.find(mass)==LEPBinEdges_.end()){
		LEPBinEdges_.insert(pair<int,vector<double> >(mass,LEPBinEdges));
	}
	else {
		LEPBinEdges_[mass]=LEPBinEdges;
	}
}
void FMTBase::setBinEdges(map<int,vector<double> > BinEdges){
	BinEdges_=BinEdges;
}
void FMTBase::setVBFBinEdges(map<int,vector<double> > VBFBinEdges){
	VBFBinEdges_=VBFBinEdges;
}
void FMTBase::setLEPBinEdges(map<int,vector<double> > LEPBinEdges){
	LEPBinEdges_=LEPBinEdges;
}

void FMTBase::printRunOptions(string filename){

	ofstream out;
	bool outFileReq;
	if (filename=="0") outFileReq=false;
	else {
		outFileReq=true;
		out.open(filename.c_str());
	}
	ostream &outFile = (outFileReq ? out : cout);
	
	out << "Running with following options:" << endl;
 	out << "\tmHMinimum                 " << mHMinimum_  << endl;               
 	out << "\tmHMaximum                 " << mHMaximum_  << endl;               
 	out << "\tmHStep                    " << mHStep_  << endl;               
 	out << "\tmassMin                   " << massMin_  << endl;               
 	out << "\tmassMax                   " << massMax_  << endl;               
 	out << "\tnDataBins                 " << nDataBins_  << endl;               
 	out << "\tsignalRegionWidth         " << signalRegionWidth_  << endl;               
 	out << "\tsidebandWidth             " << sidebandWidth_  << endl;               
 	out << "\tnumberOfSidebands         " << numberOfSidebands_  << endl;               
 	out << "\tnumberOfSidebandsForAlgos " << numberOfSidebandsForAlgos_  << endl;               
 	out << "\tnumberOfSidebandGaps      " << numberOfSidebandGaps_  << endl;               
 	out << "\tmassSidebandMin           " << massSidebandMin_  << endl;               
 	out << "\tmassSidebandMax           " << massSidebandMax_  << endl;               
 	out << "\tincludeVBF                " << includeVBF_  << endl;               
 	out << "\tincludeLEP                " << includeLEP_  << endl;               
  out << "\tsystematics               [";
  if (systematics_.size()>0) {
    for (int i=0; i<systematics_.size()-1; i++) out << systematics_[i] << ",";
    out << systematics_[systematics_.size()-1] << "]" << endl;
  }
	else out << "]" << endl;
  out << "\trederiveOptimizedBinEdges " << rederiveOptimizedBinEdges_ << endl;
  out << "\tBinEdges                  " << endl;
  for (int m=110; m<=150; m+=5){
    if (m==145) continue;
    out << "\t     BinEdges_" << m << "           [";
    if (BinEdges_[m].size()>0) {
      for (int i=0; i<BinEdges_[m].size()-1; i++) out << BinEdges_[m][i] << ",";
      out << BinEdges_[m][BinEdges_[m].size()-1];
    }
		out << "]" << endl;
  }
  out << "\tVBFBinEdges                  " << endl;
  for (int m=110; m<=150; m+=5){
    if (m==145) continue;
    out << "\t     VBFBinEdges_" << m << "        [";
    if (VBFBinEdges_[m].size()>0) {
      for (int i=0; i<VBFBinEdges_[m].size()-1; i++) out << VBFBinEdges_[m][i] << ",";
      out << VBFBinEdges_[m][VBFBinEdges_[m].size()-1];
    }
		out << "]" << endl;
  }
  out << "\tLEPBinEdges                  " << endl;
  for (int m=110; m<=150; m+=5){
    if (m==145) continue;
    out << "\t     LEPBinEdges_" << m << "        [";
    if (LEPBinEdges_[m].size()>0) {
      for (int i=0; i<LEPBinEdges_[m].size()-1; i++) out << LEPBinEdges_[m][i] << ",";
      out << LEPBinEdges_[m][LEPBinEdges_[m].size()-1];
    }
		out << "]" << endl;
  }
	out.close();
}

void FMTBase::checkHisto(TH1F *h){
  cout << "Histogram - " << left << setw(60) << h->GetName() << " e: " << h->GetEntries() << " b: " << h->GetNbinsX() << " i: " << h->Integral() << endl;
}

void FMTBase::write(TFile *tFile, TObject *obj){
  tFile->cd();
  obj->Write(obj->GetName(),TObject::kWriteDelete);
}

