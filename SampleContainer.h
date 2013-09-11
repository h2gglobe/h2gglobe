#ifndef SAMPLECONTAINER
#define SAMPLECONTAINER

#include <string>
#include <map>
#include <vector>

class SampleContainer {
	
  static float defaultextw;
  
 public:	
  SampleContainer(const float * extw=0);
  ~SampleContainer();
  
  void computeWeight(float);

  /** adds a lumi section range to 'goodLumis' (typically
      taken from a 'json' file containing the list of 
      certified luminosity sections).

      @param run the run to be added 
      @param lumi1 the first lumi section to be added
      @param lumi2 the last lumi secetion to be added */
  void addGoodLumi(int run, int lumi1, int lumi2 );

  void addEventToList(int run, int lumi, int event );
  
  bool isdata() const { return itype == 0; };
  float weight() const { return ( (extweight!=0 && *extweight > 0 && ! isdata()) ? (*extweight)*intweight : intweight); };
  
  int itype;
  int ind;
  int histoplotit;
  std::string filesshortnam;
  long long int ntot;
  int nred;
  float lumi; 
  float xsec;
  float kfactor; 
  float scale;
  int forceVersion;
  float lumireal;
  bool hasLumiSelection, hasEventList;
  std::map<int, std::vector<std::pair<int,int> > > goodLumis;
  std::map<int, std::vector<std::pair<int,int> > > eventList;
  
  std::string pileup;
  
 private:
  const float * extweight;
  float intweight;


};

#endif
