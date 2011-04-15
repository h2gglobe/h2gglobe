#ifndef hgg_VertexAnalyzer_h
#define hgg_VertexAnalyzer_h

#include <vector>
#include <map>
#include <string>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMatrixDSym.h"

namespace TMVA { class Reader; }

class VertexInfoAdapter;
class PhotonInfo;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
struct VertexAlgoParameters {
	bool rescaleTkPtByError;
	float trackCountThr;
	bool highPurityOnly;
	float maxD0Signif, maxDzSignif;
	/// float maxD0, maxDz;
	/// float minPt;
};


// -------------------------------------------------------------------------------------------------------------------------------------------------------------
/**
 * 
 * \class HggVertexAnalyzer
 * Implements vertex analysis for Higgs to gamma gamma. 
 * The class provides vertexes classification through (T)MVA classifiers and ranks combination.
 *
 * Vertexes and tracks informations are read through the VertexInfoAdapter interface and the PhotonInfo class.
 * An implementation of the VertexInfoAdapter interface that can be used with (almost) any NTuple is provided 
 * by the TupleVertexInfo class.
 * 
 * usage example:
 * <code>
 * class MyAnalzer {
 * 
 *     // algo parameters
 *     VertexAlgoParameters vtxAlgoParams_;
 *     vector<string> rankVariables_, tmvaVariables_;
 *
 *     TMVA::Reader tmvaReader_;
 *     string tmvaMethod_;
 *
 *     // branches buffers
 *     inf nvtx_;
 *     float  vtxx_[100], vtxy_[100], vtxz_[100];
 *     
 *     int ntracks_;
 *     float tkpx_[100], tkpy_[100], tkpz_[100], tkPtErr_[100], tkVtxId_[100], tkWeight_[100];
 *     
 *     float phocalox_[30], phocaloy_[30], phocaloz_[30], phoen_[30];
 *     
 *     void init() {
 *           vtxAlgoParams_.rescaleTkPtByError = true;}
 *           // variables order matters to resolve ties
 *           rankVariables_.push_back("ptbal"), rankVariables.push_back("ptasym"), rankVariables.push_back("logsumpt2");
 *                    
 *           mvaMethod_ = "BDTMethod"; 
 *           // tmva has its own order
 *           tmvaVariables_.push_back("ptbal"), tmvaVariables.push_back("ptasym"), tmvaVariables.push_back("logsumpt2");
 *           tmvaReader_ = new TMVA::Reader( "!Color:!Silent" );
 *           HggVertexAnalyzer::bookVariables( *tmvaReader_, tmvaVariables_ );
 *           tmvaReader_->BookMVA( tmvaMethod_, "<path to weights>" ); 
 *     }
 *     
 *     void analyze() {
 *     
 *           TupleVertexInfo vinfo(nvtx_, vtxx_, vtxy_, vtxz_, 
 *     			   ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_);
 *           
 *           PhotonInfo pho1(TVector3(phocalox_[p1],phocaloy_[p1],phocaloz_[p1]),phoen_[p1]), 
 *                           pho2(TVector3(phocalox_[p2],phocaloy_[p2],phocaloz_[p2]),phoen_[p2]);
 *     
 *          
 *           HggVertexAnalyzer vAna(vtxAlgoParams_,nvtx);
 *     
 *           vAna.analyze(vinfo,pho1,pho2);
 *     
 *           vecto<int> rankprod = vAna.rankprod(rankVariables_);
 *           cout << "\n\nRanks product" << endl;
 *           cout << "best vertex " << rankprod[0] << endl;
 *           for(int ii=0; ii<nvtx_; ++ii) {
 *                 int vtxrank = find(rankprod.begin(), rankprod.end(), ii) - rankprod.begin();
 *                 cout << "vertx " << ii << " rank " << vtxrank << " " << vAna.ptbal(ii) << " " << vAna.ptasym(ii) << " " << vAna.logsumpt2(ii) << endl;
 *     
 *           }
 * 
 *           vecto<int> ranktmva = vAna.rankprod(*tmvaReader_,tmvaMethod_);
 *           cout << "\n\n" << tmvaMethod_ << endl;
 *           cout << "best vertex " << rankprod[0] << endl;
 *           for(int ii=0; ii<nvtx_; ++ii) {
 *                 int vtxrank = find(ranktmva.begin(), ranktmva.end(), ii) - ranktmva.begin();
 *                 cout << "vertx " << ii << " rank " << vtxrank << " " << vAna.ptbal(ii) << " " << vAna.ptasym(ii) << " " << vAna.logsumpt2(ii) << endl;
 *           }
 *     }
 * 
 * };
 * </code>
 */
class HggVertexAnalyzer
{
public:
	typedef VertexAlgoParameters AlgoParameters;

	HggVertexAnalyzer(AlgoParameters ap, int nvtx=40);


// CINT doesn't like function pointers
#ifndef __CINT__ 
	typedef double (HggVertexAnalyzer::*getter_t) (int) const;
	typedef std::map<std::string,std::pair<getter_t,bool> > dict_t;
	
	static dict_t dictionary_;
	static void fillDictionary();

	// TMVA interface
	static std::vector<getter_t> varmeths_;
#endif
	static std::vector<float> vars_;
	
	static const double spherPwr_;
	
	static void bookVariables(TMVA::Reader & reader, const std::vector<std::string> & vars);
	static void bookSpectators(TMVA::Reader & reader, const std::vector<std::string> & vars);
	void fillVariables(int iv);
	
	// Rank vertexes
	std::vector<int> rank(std::string method);
#ifndef __CINT__	
	std::vector<int> rank(getter_t method, bool sign);
#endif
	std::vector<int> rank(TMVA::Reader &reader, const std::string & method);
	std::vector<int> ranksum(const std::vector<std::string> & vars);
	std::vector<int> rankBorda(const std::vector<std::string> & vars);
	std::vector<int> rankprod(const std::vector<std::string> & vars);
	std::vector<int> rankreciprocal(const std::vector<std::string> & vars);
	std::vector<int> rankPairwise(const std::vector<std::string> & vars);

	void analyze(const VertexInfoAdapter &, const PhotonInfo & pho1, const PhotonInfo & pho2);

	void preselection(const std::vector<int> &ps) { preselection_ = ps; }

	// getters
	double mva(int i)    const { return 	mva_[i]; };	

	double diphopt(int i)    const { return diPhoton_[i].Pt(); };	
	double nch(int i)    const { return 	nch_[i]; };	
	double ptmax(int i)  const { return 	ptmax_[i]; };	
	double sumpt(int i)  const { return 	sumpt_[i]; };	
	double ptvtx(int i)  const { return 	vtxPt_[i].Mod(); };	
	double acosA(int i)  const { return  	acosA_[i]; }; 	
	double ptasym(int i) const { return 	ptasym_[i]; };	
	double ptbal(int i)  const { return 	ptbal_[i]; };	
	
	double nchthr(int i) const { return 	nchthr_[i]; };	
	double ptmax3(int i) const { return 	ptmax3_[i]; };	
	double thrust(int i) const { return 	thrust_[i]; };	

	double sumweight(int i) const { return sumweight_[i]; };	
	double logsumpt2(int i) const { return log(sumpt2_[i]); };	
	double ptratio(int i)   const { return ptratio_[i]; };	
	double pzasym(int i)    const { return pzasym_[i]; };	
	
	double spher(int i) const { return 	spher_[i]; };	
	double tspher(int i) const { return 	tspher_[i]; };
	double aplan(int i) const { return 	aplan_[i]; };	
	double threejetC(int i) const { return 	threejetC_[i]; };
	double fourjetD(int i) const { return 	fourjetD_[i]; };
	double sumpr(int i) const { return 	sumpr_[i]; };	
	
	double sumawy(int i) const { return 	sumawy_[i]; };	
	double sumtrv(int i) const { return 	sumtrv_[i]; };	
	double sumtwd(int i) const { return 	sumtwd_[i]; };	
	double awytwdasym(int i) const { return awytwdasym_[i]; };
		
private:
	AlgoParameters params_;
	int nvtx_;
	
	std::vector<int> preselection_;

	std::vector<double> mva_;
	
	// buffers
	std::vector<TLorentzVector> diPhoton_;
	
	std::vector<double> ptbal_;
	std::vector<double> thrust_;
	std::vector<double> sumpt_;
	std::vector<double> sumpt2_;
	std::vector<double> sumawy_;
	std::vector<double> sumtwd_;
	std::vector<double> sumtrv_;
	std::vector<double> sumweight_;
	std::vector<double> ptmax_;
	std::vector<double> nchthr_;
	std::vector<double> nch_;
	std::vector<std::vector<double> > tksPt_;
	std::vector<TMatrixDSym> sphers_;
	std::vector<double> sumpr_;
	std::vector<double> spher_;
	std::vector<double> tspher_;
	std::vector<double> aplan_;
	std::vector<double> threejetC_;
	std::vector<double> fourjetD_;

	std::vector<TVector3> vtxP_;
	std::vector<TVector2> vtxPt_;

	std::vector<TVector2> diPhotonPt_;
	std::vector<double> diPhotonPz_;
	
	std::vector<double> acosA_;
	std::vector<double> ptasym_;

	std::vector<double> ptmax3_;

	std::vector<double> ptratio_;
	std::vector<double> pzasym_;

	std::vector<double> awytwdasym_;

};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
class VertexInfoAdapter 
{
public:

	virtual int nvtx() const = 0;

	virtual float vtxx(int) const = 0;
	virtual float vtxy(int) const = 0;
	virtual float vtxz(int) const = 0;

	virtual int ntracks() const = 0;
	
	virtual float tkpx(int) const = 0;
	virtual float tkpy(int) const = 0;
	virtual float tkpz(int) const = 0;
	
	virtual float tkPtErr(int) const = 0;
	virtual int   tkVtxId(int) const = 0;

	virtual float tkWeight(int) const = 0;
	
	virtual float tkd0(int) const = 0;
	virtual float tkd0Err(int) const = 0;

	virtual float tkdz(int) const = 0;
	virtual float tkdzErr(int) const = 0;

	virtual bool tkIsHighPurity(int) const = 0;


	virtual ~VertexInfoAdapter();
};


// -------------------------------------------------------------------------------------------------------------------------------------------------------------
class TupleVertexInfo : public VertexInfoAdapter
{
	TupleVertexInfo(int nvtx, float * vtxx, float * vtxy, float * vtxz, 
			   int ntracks, float * tkpx, float * tkpy, float * tkpz,
			   float * tkPtErr, int * tkVtxId, float * tkWeight,
			   float * tkd0, float * tkd0Err, float * tkdz, float * tkdzErr,
			   bool * tkIsHighPurity
		);
	
	virtual int nvtx() const    { return nvtx_; };
	virtual int ntracks() const { return ntracks_; };
	
	virtual float tkpx(int ii) const { return tkpx_ != 0 ? tkpx_[ii] : 0.; };
	virtual float tkpy(int ii) const { return tkpx_ != 0 ? tkpy_[ii] : 0.; };
	virtual float tkpz(int ii) const { return tkpx_ != 0 ? tkpz_[ii] : 0.; };
	
	virtual float tkPtErr(int ii) const { return tkPtErr_  != 0 ? tkPtErr_[ii] : 999.; };
	virtual int   tkVtxId(int ii) const { return tkVtxId_  != 0 ? tkVtxId_[ii] : 999.; };

	virtual float tkWeight(int ii) const { return tkWeight_ != 0 ? tkWeight_[ii] : 0.; };
	
	virtual float vtxx(int ii) const { return vtxx_ != 0 ? vtxx_[ii] : 0.; };
	virtual float vtxy(int ii) const { return vtxy_ != 0 ? vtxy_[ii] : 0.; };
	virtual float vtxz(int ii) const { return vtxz_ != 0 ? vtxz_[ii] : 0.; };

	virtual float tkd0(int ii) const { return tkd0_ != 0 ? tkd0_[ii] : 0.; };
	virtual float tkd0Err(int ii) const { return tkd0Err_ != 0 ? tkd0Err_[ii] : 0.; };

	virtual float tkdz(int ii) const { return tkdz_ != 0 ? tkdz_[ii] : 0.; };
	virtual float tkdzErr(int ii) const { return tkdzErr_ != 0 ? tkdzErr_[ii] : 0.; };

	virtual bool tkIsHighPurity(int ii) const { return tkIsHighPurity_ != 0 ? tkIsHighPurity_[ii] : 0.; };

	//// virtual float tkWeight(int ii, int jj) const { return -1.;}  // FIXME
	//// virtual float tkd0(int ii, int jj) const { return 0.; } // FIXME
	//// virtual float tkd0Err(int ii, int jj) const { return 0.; };  // FIXME
	//// 
	//// virtual float tkdz(int ii, int jj) const { return 0.; };  // FIXME
	//// virtual float tkdzErr(int ii, int jj) const { return 0.; };  // FIXME


	virtual ~TupleVertexInfo();
	
private:
	int nvtx_;
	float * vtxx_;
	float * vtxy_;
	float * vtxz_;

	int ntracks_;
	float * tkpx_;
	float * tkpy_;
	float * tkpz_;
	float * tkPtErr_;
	int * tkVtxId_;	
	float * tkWeight_;
	float * tkd0_;
	float * tkd0Err_;
	float * tkdz_;
	float * tkdzErr_;

	bool * tkIsHighPurity_;
};

#endif

// Local Variables:
// mode: c++       
// mode: sensitive 
// c-basic-offset: 8
// End:             
