#ifndef FMTSetup_h
#define FMTSetup_h

#include "FMTBase.h"
#include "FMTFit.h"
#include "FMTRebin.h"

class FMTSetup : public FMTBase {
	
	public:
		FMTSetup(string="0");
		~FMTSetup();

		void OptionParser(int argc, char *argv[]);
		void ReadRunConfig();

		void CheckRunOptions();
		void checkAllHistos();

		template <class T> 
		vector<T> getVecFromString(string);

		vector<double> getBinEdgesFromString(string);

    double getLumiFromWorkspace();

		template <class T>
		T getOptFromConfig(string);
		template <class T>
		T getStepSize();
		
		void organiseVectors(vector<int>&, vector<double>&);

		void printPassedOptions();

		void runRebinning();
		void runFitting();
		void createCorrBkgModel();
		void interpolateBDT();
    void makePlots();
		void writeDataCards();
		void publishToWeb();
		void runCombine();

		void cleanUp();

	private:
		
		string filename_;
		string outfilename_;

		string fitString_;
		string rebinString_;
		vector<double> fitMasses_;
		vector<int> rebinMasses_;
		bool all_;
		bool fit_;
		bool rebin_;
    bool catByHand_;
		bool skipRebin_;
		bool justRebin_;
		bool binEdges_;
    bool dumpDatFile_;
		bool bkgModel_;
		bool interp_;
		bool datacards_;
		bool diagnose_;
		bool blinding_;
		bool web_;
		bool runCombine_;
		bool checkHistos_;
    bool noPlot_;
		bool runSB_;

		int tempmHMin_;
		int tempmHMax_;
		double tempmHStep_;

		string webDir_;
		string datFil_;
    string dumpDatFil_;

    double userLumi_;

		FMTRebin *rebinner;
		bool cleaned;

};

template <class T>
T FMTSetup::getStepSize(){
	return getmHStep();
}

template <>
inline int FMTSetup::getStepSize(){
	return 5;
}

template <class T>
vector<T> FMTSetup::getVecFromString(string name){
  vector<T> result;
  
  if (typeid(T)!=typeid(float) && typeid(T)!=typeid(double) && typeid(T)!=typeid(int)){
    cout << typeid(T).name() << " is not a valid type. Bailing out " << endl;
    exit(1);
  }
  T diff = getStepSize<T>();

  while (name.size()>0){
    if (name.find(",")==string::npos && name.find("-")==string::npos) {
      result.push_back(boost::lexical_cast<T>(name)); //niether
      name = "";
    }
    else if (name.find(",")==string::npos && name.find("-")!=string::npos){ // only "-"
      T lower = boost::lexical_cast<T>(name.substr(0,name.find("-")));
      T upper = boost::lexical_cast<T>(name.substr(name.find("-")+1,string::npos));
      for (T it=lower; it<=upper; it+=diff) result.push_back(it);
      name="";
    }
    else { // mixture of "," and "-"
      // first split string at ","
      string sub = name.substr(0,name.find(","));
      if (sub.find("-")==string::npos) { // not found "-" within
        result.push_back(boost::lexical_cast<T>(sub));
        name = name.substr(name.find(",")+1,string::npos);
      }
      else {
        T lower = boost::lexical_cast<T>(sub.substr(0,sub.find("-")));
        T upper = boost::lexical_cast<T>(sub.substr(sub.find("-")+1,string::npos));
        for (T it=lower; it<=upper; it+=diff) result.push_back(it);
        name = name.substr(name.find(",")+1,string::npos);
      }
    }
  }

  return result;
}

template <class T>
T FMTSetup::getOptFromConfig(string name){
  string res = name.substr(name.find("=")+1,string::npos);
  if (res.find("=")!=string::npos){
    cout << "String: " << name << " contains more than one = sign. Bailing out." << endl;
    exit(1);
  }
  return boost::lexical_cast<T>(res);
}

#endif

