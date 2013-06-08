#ifndef __PTSPINSMEARER__
#define __PTSPINSMEARER__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include "TFile.h"
#include "TGraphAsymmErrors.h"

// ------------------------------------------------------------------------------------
class PtSpinSmearer : public BaseGenLevelSmearer
{
public:

  std::string efficiency_file;
  
  //   downId and upId set according to prescription of Jun 3, 2011  
  PtSpinSmearer(const std::string &efficiency_file) ;
  virtual ~PtSpinSmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift=0. ) const ;

  void name(const std::string & x) { name_ = x; };

  void setPTSpinName(const std::string & x) { PTSpinName_ =x; };

  bool init();
  
 protected:
  
  double getWeight( const TLorentzVector & p4, const int nPu, const int & genMassPoint, const std::string & genType, float syst_shift=0.) const;
  
  std::string   name_;
  std::string   PTSpinName_;
  TFile        *thePtSpinFile_; 
  std::map< int,std::vector<TH1*> > ptSpinSmearersSM_;
  std::map< int,std::vector<TH1*> > ptSpinSmearersGGGRAV_;
  std::map< int,std::vector<TH1*> > ptSpinSmearersQQGRAV_;
  void   readMassPoint(int mass);
  double getPtSpin(int genMassPoint, std::string genType, int id, double gPT ) const;
};

#endif
