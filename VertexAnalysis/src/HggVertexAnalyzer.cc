#include "../interface/HggVertexAnalyzer.h"
#include "../interface/PhotonInfo.h"

#include <algorithm>
#include <numeric>
#include <iterator>

#include <assert.h>

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"

#include "TMVA/Reader.h"

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
const double HggVertexAnalyzer::spherPwr_(1.5);

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
HggVertexAnalyzer::dict_t HggVertexAnalyzer::dictionary_;
vector<float> HggVertexAnalyzer::vars_;
vector<HggVertexAnalyzer::getter_t> HggVertexAnalyzer::varmeths_;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void HggVertexAnalyzer::fillDictionary()
{
	//// dictionary_["mva"] =   make_pair(&HggVertexAnalyzer::mva,false);
	
	dictionary_["diphopt"]   = make_pair(&HggVertexAnalyzer::diphopt,false);
	dictionary_["nch"]   = make_pair(&HggVertexAnalyzer::nch,false);
	dictionary_["ptmax"] = make_pair(&HggVertexAnalyzer::ptmax,false);
	dictionary_["sumpt"] = make_pair(&HggVertexAnalyzer::sumpt,false);
	dictionary_["ptvtx"] = make_pair(&HggVertexAnalyzer::ptvtx,false);
	dictionary_["acosA"] = make_pair(&HggVertexAnalyzer::acosA,false);
	dictionary_["ptasym"] = make_pair(&HggVertexAnalyzer::ptasym,false);
	dictionary_["ptbal"] = make_pair(&HggVertexAnalyzer::ptbal,false);
	
	dictionary_["nchthr"] = make_pair(&HggVertexAnalyzer::nchthr,false);
	dictionary_["ptmax3"] = make_pair(&HggVertexAnalyzer::ptmax3,false);
	dictionary_["thrust"] = make_pair(&HggVertexAnalyzer::thrust,false);
	
	dictionary_["sumweight"] = make_pair(&HggVertexAnalyzer::sumweight,false);
	dictionary_["logsumpt2"] = make_pair(&HggVertexAnalyzer::logsumpt2,false);
	dictionary_["log(sumPt2)"] = make_pair(&HggVertexAnalyzer::logsumpt2,false);
	dictionary_["ptratio"] =   make_pair(&HggVertexAnalyzer::ptratio,false);
	dictionary_["pzasym"] =    make_pair(&HggVertexAnalyzer::pzasym,false);
	
	dictionary_["spher"] = make_pair(&HggVertexAnalyzer::spher,false);
	dictionary_["aplan"] = make_pair(&HggVertexAnalyzer::aplan,false);
	dictionary_["sumpr"] = make_pair(&HggVertexAnalyzer::sumpr,false);
	
	dictionary_["sumawy"] = make_pair(&HggVertexAnalyzer::sumawy,false);
	dictionary_["sumtrv"] = make_pair(&HggVertexAnalyzer::sumtrv,false);
	dictionary_["sumtwd"] = make_pair(&HggVertexAnalyzer::sumtwd,false);
	dictionary_["awytwdasym"] = make_pair( &HggVertexAnalyzer::awytwdasym,false);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
HggVertexAnalyzer::HggVertexAnalyzer(AlgoParameters ap, int nvtx) :
	params_(ap),
	nvtx_(nvtx),
	
	mva_(nvtx),

	diPhoton_(nvtx),

	ptbal_(nvtx,0.),
	thrust_(nvtx,0.),
	sumpt_(nvtx,0.),
	sumpt2_(nvtx,0.),
	sumawy_(nvtx,0.),
	sumtwd_(nvtx,0.),
	sumtrv_(nvtx,0.),
	sumweight_(nvtx,0.),
	ptmax_(nvtx,0.),
	nchthr_(nvtx,0.),
	nch_(nvtx,0.),
	tksPt_(nvtx, vector<double>(1)),
	sphers_(nvtx,TMatrixDSym(3)),
	sumpr_(nvtx,0.),
	spher_(nvtx,0.),
	tspher_(nvtx,0.),
	aplan_(nvtx,0.),
	threejetC_(nvtx,0.),
	fourjetD_(nvtx,0.),
	
	vtxP_(nvtx,0.),
	vtxPt_(nvtx),
	
	diPhotonPt_(nvtx),
	diPhotonPz_(nvtx),

	acosA_(nvtx),
	ptasym_(nvtx),
	
	ptmax3_(nvtx),
	
	ptratio_(nvtx),
	pzasym_(nvtx),
	
	awytwdasym_(nvtx)
{
	
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void HggVertexAnalyzer::bookVariables(TMVA::Reader & reader,const std::vector<std::string> & order)
{
	if( dictionary_.empty() ) { fillDictionary(); }
	vars_.resize(order.size(),0.);
	varmeths_.resize(order.size(),0);
	for(size_t ivar=0; ivar<order.size(); ++ivar) {
		reader.AddVariable( order[ivar], &vars_[ivar] );
		varmeths_[ivar] = dictionary_[order[ivar]].first;
	}
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void HggVertexAnalyzer::bookSpectators(TMVA::Reader & reader,const std::vector<std::string> & order)
{
	if( dictionary_.empty() ) { fillDictionary(); }
	for(size_t ivar=0; ivar<order.size(); ++ivar) {
		vars_.push_back(0.); 
		reader.AddSpectator( order[ivar], &*vars_.rbegin() );
		varmeths_.push_back( dictionary_[order[ivar]].first );
	}
}


// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void HggVertexAnalyzer::fillVariables(int iv)
{
	for(size_t ivar=0; ivar<varmeths_.size(); ++ivar) {
		vars_[ivar] = ((*this).*(varmeths_[ivar]))(iv);
	}
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
class RankHelper
{
public:
	RankHelper( HggVertexAnalyzer & vAna, vector<pair< HggVertexAnalyzer::getter_t, bool> > m) :
		vAna_(vAna), method_(m) {};
	RankHelper( HggVertexAnalyzer & vAna, HggVertexAnalyzer::getter_t m, bool sign) :
		vAna_(vAna), method_(vector<pair<HggVertexAnalyzer::getter_t,bool> >(1,make_pair(m,sign))) {};

	bool operator() (int lh,int rh) {

		if (lh == rh){
			return true;
		}

		vector<pair< HggVertexAnalyzer::getter_t, bool> >::iterator imeth;
		for(imeth = method_.begin(); imeth != method_.end(); ++imeth){
			double lhv = (vAna_.*imeth->first)(lh);
			double rhv = (vAna_.*imeth->first)(rh);
			if( lhv != rhv ){
				return imeth->second ? lhv < rhv : lhv > rhv;
			}
		}
		
		return lh < rh;
	}
	
private:
	HggVertexAnalyzer & vAna_;
	vector<pair< HggVertexAnalyzer::getter_t,bool> > method_;
};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> HggVertexAnalyzer::rank(HggVertexAnalyzer::getter_t method, bool sign)
{
	std::vector<int> vtxs = preselection_;
	RankHelper helper(*this, method, sign );
	sort(vtxs.begin(),vtxs.end(),helper);

	return vtxs;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> HggVertexAnalyzer::rank(std::string method)
{
	dict_t::value_type::second_type m = dictionary_[method];
	return rank(m.first,m.second);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> HggVertexAnalyzer::rank(TMVA::Reader &reader, const std::string & method)
{
	std::vector<int> vtxs = preselection_;
	for(int ii=0; ii<nvtx_; ++ii) { 
		fillVariables(ii);
		mva_[ii] = reader.EvaluateMVA(method);
	}
	RankHelper helper(*this,&HggVertexAnalyzer::mva,false);
	sort(vtxs.begin(),vtxs.end(),helper);

	return vtxs;	
}


// -------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> HggVertexAnalyzer::ranksum(const vector<string> & vars)
{
	std::vector<int> vtxs = preselection_;

	for(int ii=0; ii<nvtx_; ++ii) {
		mva_[ii] = 0.;
	}
	// fill the rank sum
	std::vector<int> vrank(nvtx_);
	std::vector<pair<HggVertexAnalyzer::getter_t, bool> > meths(1,make_pair(&HggVertexAnalyzer::mva,true));
	for( vector<string>::const_iterator ivar=vars.begin(); ivar!=vars.end(); ++ivar) {
		meths.push_back(dictionary_[*ivar]);
		vrank = rank( *ivar );
		for(size_t ii=0; ii<vtxs.size(); ++ii) {
			int ivert = vtxs[ii];
			mva_[ivert] += (double)(find( vrank.begin(), vrank.end(), ivert) - vrank.begin());
		}
	}
	RankHelper helper(*this,meths);
	sort(vtxs.begin(),vtxs.end(),helper);

	return vtxs;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> HggVertexAnalyzer::rankprod(const vector<string> & vars)
{
	std::vector<int> vtxs = preselection_;

	for(int ii=0; ii<nvtx_; ++ii) {
		mva_[ii] = 1.;
	}
	// fill the rank sum
	std::vector<int> vrank(nvtx_);
	std::vector<pair<HggVertexAnalyzer::getter_t, bool> > meths(1,make_pair(&HggVertexAnalyzer::mva,true));
	for( vector<string>::const_iterator ivar=vars.begin(); ivar!=vars.end(); ++ivar) {
		meths.push_back(dictionary_[*ivar]);
		vrank = rank( *ivar );
		for(size_t ii=0; ii<vtxs.size(); ++ii) {
			int ivert = vtxs[ii];
			int rank = find( vrank.begin(), vrank.end(), ivert) - vrank.begin(); 
			mva_[ivert] *= 1. + (double)(rank);
			/// if( mva_[ivert] < 0. ) {
			/// 	cerr << "HggVertexAn  : negative rankprod " << ivert << " " << rank << " " << vrank.size() << " " << mva_[ivert] << endl;
			/// }
		}
	}
	for(int ii=0; ii<nvtx_; ++ii) {
		mva_[ii] = pow( mva_[ii], 1./(double)vars.size() );
	}
	RankHelper helper(*this,meths);
	sort(vtxs.begin(),vtxs.end(),helper);

	return vtxs;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> HggVertexAnalyzer::rankreciprocal(const vector<string> & vars)
{
	std::vector<int> vtxs = preselection_;

	for(int ii=0; ii<nvtx_; ++ii) {
		mva_[ii] = 0.;
	}
	// fill the rank sum
	std::vector<int> vrank(nvtx_);
	std::vector<pair<HggVertexAnalyzer::getter_t, bool> > meths(1,make_pair(&HggVertexAnalyzer::mva,false));
	for( vector<string>::const_iterator ivar=vars.begin(); ivar!=vars.end(); ++ivar) {
		meths.push_back(dictionary_[*ivar]);
		vrank = rank( *ivar );
		for(size_t ii=0; ii<vtxs.size(); ++ii) {
			int ivert = vtxs[ii];
			mva_[ivert] += pow((double)(1 + (find( vrank.begin(), vrank.end(), ivert) - vrank.begin())),-2);
		}
	}
	RankHelper helper(*this,meths);
	sort(vtxs.begin(),vtxs.end(),helper);

	return vtxs;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> HggVertexAnalyzer::rankBorda(const vector<string> & vars)
{
	std::vector<int> vtxs = preselection_;

	for(int ii=0; ii<nvtx_; ++ii) {
		mva_[ii] = 0.;
	}
	// fill the rank sum
	std::vector<int> vrank(nvtx_);
	std::vector<pair<HggVertexAnalyzer::getter_t, bool> > meths(1,make_pair(&HggVertexAnalyzer::mva,false));
	for( vector<string>::const_iterator ivar=vars.begin(); ivar!=vars.end(); ++ivar) {
		meths.push_back(dictionary_[*ivar]);
		vrank = rank( *ivar );
		for(size_t ii=0; ii<vtxs.size(); ++ii) {
			int ivert = vtxs[ii];
			mva_[ivert] += (double)(vtxs.size() - (find( vrank.begin(), vrank.end(), ivert) - vrank.begin()));
		}
	}
	RankHelper helper(*this,meths);
	sort(vtxs.begin(),vtxs.end(),helper);

	return vtxs;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
vector<int> HggVertexAnalyzer::rankPairwise(const vector<string> & vars)
{
	vector<int> vtxs = preselection_;

	for(int ii=0; ii<nvtx_; ++ii) {
		mva_[ii] = 0.;
	}
	// fill the rank sum
	vector<int> vrank(nvtx_);
	vector<vector<int> > comp(nvtx_,vector<int>(nvtx_,0));
	vector<pair<HggVertexAnalyzer::getter_t, bool> > meths(1,make_pair(&HggVertexAnalyzer::mva,false));
	for( vector<string>::const_iterator ivar=vars.begin(); ivar!=vars.end(); ++ivar) {
		meths.push_back(dictionary_[*ivar]);
		vrank = rank( *ivar );
		for(size_t i=0; i<vtxs.size(); ++i) {
			int v1 = vtxs[i];
			for(size_t j=0; j<i; ++j) {
				int v2 = vtxs[j];
				int rankv1 = (int)(find( vrank.begin(), vrank.end(), v1) - vrank.begin());
				int rankv2 = (int)(find( vrank.begin(), vrank.end(), v2) - vrank.begin());
				if(rankv1 < rankv2){
					comp[v1][v2]++;
				}else{
					comp[v2][v1]++;
				}
			}
		}
	}

	for(size_t i=0; i<vtxs.size(); ++i) {
		int v1 = vtxs[i];
		for(size_t j=0; j<i; ++j) {
			int v2 = vtxs[j];
			if( comp[v1][v2] > comp[v2][v1] ){
				mva_[v1] += 1.0;
			}else if ( comp[v1][v2] < comp[v2][v1] ) {
				mva_[v2] += 1.0;
			}else{
				mva_[v1] += 0.5;
				mva_[v2] += 0.5;
			}
		}
	}

	RankHelper helper(*this,meths);
	sort(vtxs.begin(),vtxs.end(),helper);

	return vtxs;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void HggVertexAnalyzer::analyze(const VertexInfoAdapter & e, const PhotonInfo & p1, const PhotonInfo & p2)
{
	int const nvtx = e.nvtx();
	nvtx_ = nvtx;
	
	// initilise
	mva_.resize(nvtx,0.);
	ptbal_.resize(nvtx,0.);
	thrust_.resize(nvtx,0.);
	sumpt_.resize(nvtx,0.);
	sumpt2_.resize(nvtx,0.);
	sumawy_.resize(nvtx,0.);
	sumtwd_.resize(nvtx,0.);
	sumtrv_.resize(nvtx,0.);
	sumweight_.resize(nvtx,0.);
	ptmax_.resize(nvtx,0.);
	nchthr_.resize(nvtx,0.);
	nch_.resize(nvtx,0.);
	vtxP_.resize(nvtx,0.);
	tksPt_.resize(nvtx, vector<double>(1));
	sphers_.resize(nvtx,TMatrixDSym(3));
	sumpr_.resize(nvtx,0.);
	spher_.resize(nvtx,0.);
	tspher_.resize(nvtx,0.);
	aplan_.resize(nvtx,0.);
	threejetC_.resize(nvtx,0.);
	fourjetD_.resize(nvtx,0.);
	
	diPhotonPt_.resize(nvtx);
	vtxPt_.resize(nvtx);
	diPhotonPz_.resize(nvtx);
	
	acosA_.resize(nvtx);
	ptasym_.resize(nvtx);
	
	ptmax3_.resize(nvtx);
	thrust_.resize(nvtx);
	
	ptratio_.resize(nvtx);
	pzasym_.resize(nvtx);
	
	awytwdasym_.resize(nvtx);
	
	diPhoton_.resize(nvtx);
	for(int i=0; i<nvtx; ++i) {
		diPhoton_[i] = 
			p1.p4(e.vtxx(i),e.vtxy(i),e.vtxz(i)) +
			p2.p4(e.vtxx(i),e.vtxy(i),e.vtxz(i));
	}

	//calculating loop over tracks
	for(int i=0; i<e.ntracks(); ++i) {
		
		if( ( params_.highPurityOnly && !e.tkIsHighPurity(i)  )
		    || fabs(e.tkd0(i)/e.tkd0Err(i)) > params_.maxD0Signif 
		    || fabs(e.tkd0(i)/e.tkd0Err(i)) > params_.maxDzSignif ) { 
			continue; 
		}
		
		const TVector3 tkPVec(e.tkpx(i),e.tkpy(i),e.tkpz(i));

		int vid = e.tkVtxId(i);
		assert(vid >= 0 && vid < nvtx);

		TVector2 tkPtVec = tkPVec.XYvector();
		double tkPt = tkPtVec.Mod();
		/// cout << "tkPt "<< tkPt << " tkPtErr " << i << " " << e.tkPtErr(i) << endl;
		const double modpt = tkPt > e.tkPtErr(i) ? tkPt - e.tkPtErr(i)  : 0.;
		if( modpt == 0. ) { continue; }
		
		// correct track pt a la POG
		if( params_.rescaleTkPtByError ) {
			const double ptcorr = modpt/tkPt;
			tkPtVec *= ptcorr;
			tkPt = modpt;
		}
		
		if(tkPt > params_.trackCountThr) nchthr_[vid] += 1;
		nch_[vid] += 1;
		ptbal_[vid] -= tkPtVec * diPhoton_[vid].Vect().XYvector().Unit();
		sumpt_[vid] += tkPt;
		sumpt2_[vid] += tkPtVec.Mod2();
		double cosTk = tkPVec.Unit() * diPhoton_[vid].Vect().Unit();
		double val = tkPtVec.Mod();
		if ( cosTk < -0.5 )	{
			sumawy_[vid] += val;
		} else if ( cosTk > 0.5 ){
			sumtwd_[vid] += val;
		} else {
			sumtrv_[vid] += val;
		}
		sumweight_[vid] += e.tkWeight(i);
		vtxP_[vid] += tkPVec;
		tksPt_[vid].push_back(tkPt);

		Double_t p[3] = {0.,0.,0.};
		tkPVec.GetXYZ(p);
		for(int j=3; j--;){
			for(int k=j+1; k--;){
				(sphers_[vid])[j][k] += pow(tkPVec.Mag(),spherPwr_-2.) * p[j]*p[k];
			}
		}
		sumpr_[vid] += pow(tkPVec.Mag(),spherPwr_);
	}
	
	preselection_.clear();
	// filling loop over vertexes
	for(int i=0; i<e.nvtx(); ++i) {
		
		sphers_[i] *= 1./sumpr_[i];
		
		TVectorD eigVals(3);
		eigVals = TMatrixDSymEigen(sphers_[i]).GetEigenValues();

		spher_[i] = 1.5 * (eigVals[1]+eigVals[2]);
		tspher_[i] = 2. * eigVals[1] / (eigVals[0]+eigVals[1]);
		aplan_[i] = 1.5 * eigVals[2];

		threejetC_[i] = 3. * (eigVals[0]*eigVals[1] + eigVals[0]*eigVals[2] + eigVals[1]*eigVals[2]);
		fourjetD_[i] = 27. * eigVals[0]*eigVals[1]*eigVals[2];


		diPhotonPt_[i] = diPhoton_[i].Vect().XYvector();
		vtxPt_[i]      = vtxP_[i].XYvector();
		diPhotonPz_[i]   = diPhoton_[i].Vect().Pz();
		
 		sort(tksPt_[i].begin(), tksPt_[i].end(), greater<double>());
		
		acosA_[i] =  	acos(vtxPt_[i].Unit() * diPhotonPt_[i].Unit());
		ptasym_[i] = 	(vtxPt_[i].Mod() - diPhotonPt_[i].Mod())/(vtxPt_[i].Mod() + diPhotonPt_[i].Mod());

		ptmax_ [i] = 	tksPt_[i][0];
		ptmax3_[i] = 	accumulate(tksPt_[i].begin(),tksPt_[i].begin() + min(tksPt_[i].size(),(size_t)3), 0.0) ;
		thrust_[i] = 	ptbal_[i]/sumpt_[i];
		
		ptratio_[i] = 	vtxPt_[i].Mod()/diPhotonPt_[i].Mod();
		pzasym_[i] = 	fabs( (vtxP_[i].Pz() - diPhotonPz_[i])/(vtxP_[i].Pz() + diPhotonPz_[i]) );
		
		awytwdasym_[i] = (sumawy_[i]-sumtwd_[i])/(sumawy_[i]+sumtwd_[i]);
		
		preselection_.push_back(i);
	}
}


// -------------------------------------------------------------------------------------------------------------------------------------------------------------
VertexInfoAdapter::~VertexInfoAdapter()
{}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
TupleVertexInfo::~TupleVertexInfo()
{}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
TupleVertexInfo::TupleVertexInfo(int nvtx, float * vtxx, float * vtxy, float * vtxz, 
				 int ntracks, float * tkpx, float * tkpy, float * tkpz,
				 float * tkPtErr, int * tkVtxId, float * tkWeight,
				 float * tkd0, float * tkd0Err, float * tkdz, float * tkdzErr,
				 bool * tkIsHighPurity
) :
	nvtx_(nvtx),
	vtxx_(vtxx),
	vtxy_(vtxy),
	vtxz_(vtxz),
	
	ntracks_(ntracks),
	tkpx_(tkpx),
	tkpy_(tkpy),
	tkpz_(tkpz),
	tkPtErr_(tkPtErr),
	tkVtxId_(tkVtxId),
	tkWeight_(tkWeight),
	tkd0_(tkd0),
	tkd0Err_(tkd0Err),
	tkdz_(tkdz),
	tkdzErr_(tkdzErr),

	tkIsHighPurity_(tkIsHighPurity)
{
}

// Local Variables:
// mode: c++       
// mode: sensitive 
// c-basic-offset: 8
// End:             
