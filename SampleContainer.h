#ifndef SAMPLECONTAINER
#define SAMPLECONTAINER

#include <string>
#include <map>
#include <vector>

class SampleContainer {

 public:
  SampleContainer();
  ~SampleContainer();
  
  void computeWeight(float);
  void addGoodLumi(int run, int lumi1, int lumi2 );
  void addEventToList(int run, int lumi, int event );

     
  float weight;
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
  float lumireal;
  bool hasLumiSelection, hasEventList;
  std::map<int, std::vector<std::pair<int,int> > > goodLumis;
  std::map<int, std::vector<std::pair<int,int> > > eventList;

  
 private:

};

#endif
