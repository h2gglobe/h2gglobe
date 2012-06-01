#ifndef FMTBase_h
#define FMTBase_h
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "TH1F.h"
#include "TFile.h"
#include "TObject.h"

using namespace std;

class FMTBase {

	public:
		FMTBase(){};
		FMTBase(int, int, double, double, double, int, double, double, int, int, int, double, double, bool, int, bool, int, vector<string>, bool, vector<map<int,vector<double> > >, bool verbose=false); 
		~FMTBase(){};

		void checkMCMass(int);
		template <class T>
		void printVec(vector<T>);
    template <class T>
    string returnVecAsString(vector<T>);

		// getters
		int getmHMinimum();
		int getmHMaximum();
		double getmHStep();
		double getmassMin();
		double getmassMax();
		int getnDataBins();
		
		double getsignalRegionWidth();
		double getsidebandWidth();
		int getnumberOfSidebands();
		int getnumberOfSidebandsForAlgos();
		int getnumberOfSidebandGaps();
		double getmassSidebandMin();
		double getmassSidebandMax();

		pair<int,int> getNsidebandsUandD(double);
		vector<double> getLowerSidebandEdges(double);
		vector<double> getUpperSidebandEdges(double);
    vector<string> getProdTypes();
	
		bool getincludeVBF();
    int getnVBFCategories();
		bool getincludeLEP();
    int getnLEPCategories();
		const int getNcats();

		vector<string> getsystematics();
		int getNsystematics();

		bool getrederiveOptimizedBinEdges();
		vector<map<int,vector<double> > > getAllBinEdges();
		map<int, vector<double> > getBinEdges();
		map<int, vector<double> > getVBFBinEdges();
		map<int, vector<double> > getLEPBinEdges();
		vector<vector<double> > getAllBinEdges(int);
		vector<double> getBinEdges(int);
		vector<double> getVBFBinEdges(int);
		vector<double> getLEPBinEdges(int);

		vector<double> getMHMasses(int);
    vector<int> getUandDMCMasses(int);
    pair<int,int> getInterpMasses(double);
		vector<int> getMCMasses();
    vector<double> getAllMH();
		int getNumMCMasses();
    int getNumMHMasses();

		// setters
		void setmHMinimum(int);
		void setmHMaximum(int);
		void setmHStep(double);
		void setmassMin(double);
		void setmassMax(double);
		void setnDataBins(int);
		
		void setsignalRegionWidth(double);
		void setsidebandWidth(double);
		void setnumberOfSidebands(int);
		void setnumberOfSidebandsForAlgos(int);
		void setnumberOfSidebandGaps(int);
		void setmassSidebandMin(double);
		void setmassSidebandMax(double);
	
		void setincludeVBF(bool);
    void setnVBFCategories(int);
		void setincludeLEP(bool);
    void setnLEPCategories(int);

		void setsystematics(vector<string>);
		void setsystematic(string);

		void setrederiveOptimizedBinEdges(bool);
		void setAllBinEdges(vector<map<int,vector<double> > >);
		void setBinEdges(int, vector<double>);
		void setVBFBinEdges(int, vector<double>);
		void setLEPBinEdges(int, vector<double>);
		void setBinEdges(map<int,vector<double> >);
		void setVBFBinEdges(map<int,vector<double> >);
		void setLEPBinEdges(map<int,vector<double> >);

    bool isIncCat(int);
    bool isVBFCat(int);
    bool isLEPCat(int);
		void printRunOptions(string filename="0");
    void dumpDatFile(string filename="mvaanalysis.dat");
		void checkHisto(TH1F*);

		void write(TFile*,TObject*);

  protected:
    bool verbose_;

	private:
		int mHMinimum_;
		int mHMaximum_;
		double mHStep_;
		double massMin_;
		double massMax_;
		int nDataBins_;
		
		double signalRegionWidth_;
		double sidebandWidth_;
		int numberOfSidebands_;
		int numberOfSidebandsForAlgos_;
		int numberOfSidebandGaps_;
		double massSidebandMin_;
		double massSidebandMax_;
	
		bool includeVBF_;
    int nVBFCategories_;
		bool includeLEP_;
    int nLEPCategories_;

		vector<string> systematics_;

		bool rederiveOptimizedBinEdges_;
		map<int, vector<double> > BinEdges_;
		map<int, vector<double> > VBFBinEdges_;
		map<int, vector<double> > LEPBinEdges_;

};

template <class T>
string FMTBase::returnVecAsString(vector<T> vec){
  string result="";
  for (typename vector<T>::iterator it=vec.begin(); it!=vec.end()-1; it++) result+=Form("%1.4f,",*it);
  result+=Form("%1.4f",*(vec.end()-1));
  return result;
}

template <class T>
void FMTBase::printVec(vector<T> vec){
	for (typename vector<T>::iterator it=vec.begin(); it!=vec.end()-1; it++) cout << *it << ",";
	cout << *(vec.end()-1);
}

#endif
