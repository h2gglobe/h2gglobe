#ifndef FMTTree_h
#define FMTTree_h

#include "FMTFit.h"
#include "FMTBase.h"

using namespace std;

class FMTTree : public FMTBase {
	
	public:
    FMTTree(double,bool,int, int, double, double, double, int, double, double, int, int, int, double, double, int, bool, int, bool, int, vector<string>, bool, vector<map<int,vector<double> > >, bool verbose=false);
    ~FMTTree(){};

   
   private:
    float mass_;
    float bdtoutput_;
    float weight_;
    int category_;
    vector<pair<float,float> > massSyst_;
    vector<pair<float,float> > bdtouputSyst_;
    vector<pair<float,float> > weightSyst_;
    vector<pair<int,int> > categorySyst_;

};

#endif
