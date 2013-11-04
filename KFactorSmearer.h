#ifndef __KFACTORSMEARER__
#define __KFACTORSMEARER__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include "TFile.h"
#include "TGraphAsymmErrors.h"

class Normalization_8TeV;

// ------------------------------------------------------------------------------------
class KFactorSmearer : public BaseGenLevelSmearer
{
public:

  std::string efficiency_file;
  
  //   downId and upId set according to prescription of Jun 3, 2011  
  KFactorSmearer(const std::string &efficiency_file, Normalization_8TeV * norm, unsigned int downId =1, unsigned int upId=6 ) ;
  virtual ~KFactorSmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift=0. ) const ;

  void name(const std::string & x) { name_ = x; };

  void setKFName(const std::string & x) { KFName_ =x; };

  bool init();
  
 protected:
  
  double getWeight( const TLorentzVector & p4, const int nPu, const int & genMassPoint, float syst_shift=0.) const;
  
  std::string   name_;
  std::string   KFName_;
  Normalization_8TeV * norm_;
  TFile        *theKFactorFile_; 
  std::map< int,std::vector<TH1*> > kFactorSmearers_;
  void   readMassPoint(int mass, int uId, int dId );
  double getKFactor(int genMassPoint, int id, double gPT ) const;

  unsigned int downId, upId; 
};

#endif
