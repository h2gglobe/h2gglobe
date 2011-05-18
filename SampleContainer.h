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
  bool hasLumiSelection;
  std::map<int, std::vector<std::pair<int,int> > > goodLumis;

  
 private:

};

#endif
