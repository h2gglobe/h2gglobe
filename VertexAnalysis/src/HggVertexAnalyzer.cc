#include "../interface/HggVertexAnalyzer.h"
#include "../interface/PhotonInfo.h"

#include "stdio.h"
#include "math.h"

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
vector<float> HggVertexAnalyzer::vars_;
vector<HggVertexAnalyzer::getter_t> HggVertexAnalyzer::varmeths_;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
HggVertexAnalyzer::dict_t & HggVertexAnalyzer::dictionary() 
{
	static HggVertexAnalyzer::dict_t dictionary;
	if( dictionary.empty() ) { fillDictionary(dictionary); }
	return dictionary;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void HggVertexAnalyzer::fillDictionary(HggVertexAnalyzer::dict_t& dictionary)
{
	//// dictionary_["mva"] =   make_pair(&HggVertexAnalyzer::mva,false);
	
	dictionary["diphopt"]   = make_pair(&HggVertexAnalyzer::diphopt,false);
	dictionary["nch"]   = make_pair(&HggVertexAnalyzer::nch,false);
	dictionary["ptmax"] = make_pair(&HggVertexAnalyzer::ptmax,false);
	dictionary["sumpt"] = make_pair(&HggVertexAnalyzer::sumpt,false);
	dictionary["ptvtx"] = make_pair(&HggVertexAnalyzer::ptvtx,false);
	dictionary["acosA"] = make_pair(&HggVertexAnalyzer::acosA,false);
	dictionary["ptasym"] = make_pair(&HggVertexAnalyzer::ptasym,false);
	dictionary["ptbal"] = make_pair(&HggVertexAnalyzer::ptbal,false);
	
	dictionary["nchthr"] = make_pair(&HggVertexAnalyzer::nchthr,false);
	dictionary["ptmax3"] = make_pair(&HggVertexAnalyzer::ptmax3,false);
	dictionary["thrust"] = make_pair(&HggVertexAnalyzer::thrust,false);
	
	dictionary["sumweight"] = make_pair(&HggVertexAnalyzer::sumweight,false);
	dictionary["logsumpt2"] = make_pair(&HggVertexAnalyzer::logsumpt2,false);
	dictionary["log(sumPt2)"] = make_pair(&HggVertexAnalyzer::logsumpt2,false);
	dictionary["ptratio"] =   make_pair(&HggVertexAnalyzer::ptratio,false);
	dictionary["pzasym"] =    make_pair(&HggVertexAnalyzer::pzasym,false);
	
	dictionary["spher"] = make_pair(&HggVertexAnalyzer::spher,false);
	dictionary["aplan"] = make_pair(&HggVertexAnalyzer::aplan,false);
	dictionary["sumpr"] = make_pair(&HggVertexAnalyzer::sumpr,false);
	
	dictionary["sumawy"] = make_pair(&HggVertexAnalyzer::sumawy,false);
	dictionary["sumtrv"] = make_pair(&HggVertexAnalyzer::sumtrv,false);
	dictionary["sumtwd"] = make_pair(&HggVertexAnalyzer::sumtwd,false);
	dictionary["awytwdasym"] = make_pair( &HggVertexAnalyzer::awytwdasym,false);
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
	vars_.resize(order.size(),0.);
	varmeths_.resize(order.size(),0);
	for(size_t ivar=0; ivar<order.size(); ++ivar) {
		reader.AddVariable( order[ivar], &vars_[ivar] );
		varmeths_[ivar] = dictionary()[order[ivar]].first;
	}
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void HggVertexAnalyzer::bookSpectators(TMVA::Reader & reader,const std::vector<std::string> & order)
{
	for(size_t ivar=0; ivar<order.size(); ++ivar) {
		vars_.push_back(0.); 
		reader.AddSpectator( order[ivar], &*vars_.rbegin() );
		varmeths_.push_back( dictionary()[order[ivar]].first );
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
	dict_t::value_type::second_type m = dictionary()[method];
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
		meths.push_back(dictionary()[*ivar]);
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
		meths.push_back(dictionary()[*ivar]);
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
		meths.push_back(dictionary()[*ivar]);
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
		meths.push_back(dictionary()[*ivar]);
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
		meths.push_back(dictionary()[*ivar]);
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
	mva_.clear(); mva_.resize(nvtx,0.);
	ptbal_.clear(); ptbal_.resize(nvtx,0.);
	thrust_.clear(); thrust_.resize(nvtx,0.);
	sumpt_.clear(); sumpt_.resize(nvtx,0.);
	sumpt2_.clear(); sumpt2_.resize(nvtx,0.);
	sumawy_.clear(); sumawy_.resize(nvtx,0.);
	sumtwd_.clear(); sumtwd_.resize(nvtx,0.);
	sumtrv_.clear(); sumtrv_.resize(nvtx,0.);
	sumweight_.clear(); sumweight_.resize(nvtx,0.);
	ptmax_.clear(); ptmax_.resize(nvtx,0.);
	nchthr_.clear(); nchthr_.resize(nvtx,0.);
	nch_.clear(); nch_.resize(nvtx,0.);
	vtxP_.clear(); vtxP_.resize(nvtx,0.);
	tksPt_.clear(); tksPt_.resize(nvtx, vector<double>(1));
	sphers_.clear(); sphers_.resize(nvtx,TMatrixDSym(3));
	sumpr_.clear(); sumpr_.resize(nvtx,0.);
	spher_.clear(); spher_.resize(nvtx,0.);
	tspher_.clear(); tspher_.resize(nvtx,0.);
	aplan_.clear(); aplan_.resize(nvtx,0.);
	threejetC_.clear(); threejetC_.resize(nvtx,0.);
	fourjetD_.clear(); fourjetD_.resize(nvtx,0.);
	
	diPhotonPt_.clear(); diPhotonPt_.resize(nvtx);
	vtxPt_.clear(); vtxPt_.resize(nvtx);
	diPhotonPz_.clear(); diPhotonPz_.resize(nvtx);
	
	acosA_.clear(); acosA_.resize(nvtx);
	ptasym_.clear(); ptasym_.resize(nvtx);
	
	ptmax3_.clear(); ptmax3_.resize(nvtx);
	thrust_.clear(); thrust_.resize(nvtx);
	
	ptratio_.clear(); ptratio_.resize(nvtx);
	pzasym_.clear(); pzasym_.resize(nvtx);
	
	awytwdasym_.clear(); awytwdasym_.resize(nvtx);
	
	diPhoton_.clear(); diPhoton_.resize(nvtx);
	for(int i=0; i<nvtx; ++i) {
		diPhoton_[i] = 
			p1.p4(e.vtxx(i),e.vtxy(i),e.vtxz(i)) +
			p2.p4(e.vtxx(i),e.vtxy(i),e.vtxz(i));
	}
	
	std::vector<unsigned short> vtxTracksBuf;
	std::vector<int> vtxTracksSizeBuf;
	if( ! e.hasVtxTracks() ) {
		vtxTracksBuf.resize(nvtx*e.ntracks());
		vtxTracksSizeBuf.resize(nvtx,0);
		for(int it=0; it<e.ntracks(); ++it) {
			int vid = e.tkVtxId(it);
			int & ntks = vtxTracksSizeBuf[vid];
			vtxTracksBuf[ vid*e.ntracks() + ntks ] = it;
			++ntks;
		}
	}

	///////// //calculating loop over tracks
	///////// for(int i=0; i<e.ntracks(); ++i) {
	///////// 	
	///////// 	int vid = e.tkVtxId(i);
	///////// 	
	///////// 	if( ( params_.highPurityOnly && !e.tkIsHighPurity(i)  )
	///////// 	    || fabs(e.tkd0(i,vid)/e.tkd0Err(i,vid)) > params_.maxD0Signif 
	///////// 	    || fabs(e.tkd0(i,vid)/e.tkdzErr(i,vid)) > params_.maxDzSignif ) { 
	///////// 		continue; 
	///////// 	}
	///////// 	
	///////// 	const TVector3 tkPVec(e.tkpx(i),e.tkpy(i),e.tkpz(i));
	///////// 	assert(vid >= 0 && vid < nvtx);
	///////// 
	///////// 	TVector2 tkPtVec = tkPVec.XYvector();
	///////// 	double tkPt = tkPtVec.Mod();
	///////// 	/// cout << "tkPt "<< tkPt << " tkPtErr " << i << " " << e.tkPtErr(i) << endl;
	///////// 	const double modpt = tkPt > e.tkPtErr(i) ? tkPt - e.tkPtErr(i)  : 0.;
	///////// 	if( modpt == 0. ) { continue; }
	///////// 	
	///////// 	// correct track pt a la POG
	///////// 	if( params_.rescaleTkPtByError ) {
	///////// 		const double ptcorr = modpt/tkPt;
	///////// 		tkPtVec *= ptcorr;
	///////// 		tkPt = modpt;
	///////// 	}
	///////// 	
	///////// 	if(tkPt > params_.trackCountThr) nchthr_[vid] += 1;
	///////// 	nch_[vid] += 1;
	///////// 	ptbal_[vid] -= tkPtVec * diPhoton_[vid].Vect().XYvector().Unit();
	///////// 	sumpt_[vid] += tkPt;
	///////// 	sumpt2_[vid] += tkPtVec.Mod2();
	///////// 	double cosTk = tkPVec.Unit() * diPhoton_[vid].Vect().Unit();
	///////// 	double val = tkPtVec.Mod();
	///////// 	if ( cosTk < -0.5 )	{
	///////// 		sumawy_[vid] += val;
	///////// 	} else if ( cosTk > 0.5 ){
	///////// 		sumtwd_[vid] += val;
	///////// 	} else {
	///////// 		sumtrv_[vid] += val;
	///////// 	}
	///////// 	sumweight_[vid] += e.tkWeight(i,vid);
	///////// 	vtxP_[vid] += tkPVec;
	///////// 	tksPt_[vid].push_back(tkPt);
	///////// 
	///////// 	Double_t p[3] = {0.,0.,0.};
	///////// 	tkPVec.GetXYZ(p);
	///////// 	for(int j=3; j--;){
	///////// 		for(int k=j+1; k--;){
	///////// 			(sphers_[vid])[j][k] += pow(tkPVec.Mag(),spherPwr_-2.) * p[j]*p[k];
	///////// 		}
	///////// 	}
	///////// 	sumpr_[vid] += pow(tkPVec.Mag(),spherPwr_);
	///////// }
	
	preselection_.clear();
	// filling loop over vertexes
	for(int vid=0; vid<e.nvtx(); ++vid) {
		
		const unsigned short * vtxTracks = e.hasVtxTracks() ? e.vtxTracks(vid) : &vtxTracksBuf[ vid*e.ntracks() ];
		int ntracks = e.hasVtxTracks() ? e.vtxNTracks(vid) : vtxTracksSizeBuf[ vid ];

		//calculating loop over tracks
		for(int it=0; it<ntracks; ++it) {
			
			int tid = vtxTracks[it];
			float tkWeight = e.tkWeight(tid,vid);			
			
			if( ( params_.highPurityOnly && !e.tkIsHighPurity(tid)  )
			    || fabs(e.tkd0(tid,vid)/e.tkd0Err(tid,vid)) > params_.maxD0Signif 
			    || fabs(e.tkd0(tid,vid)/e.tkdzErr(tid,vid)) > params_.maxDzSignif ) { 
				continue; 
			}
			
			const TVector3 tkPVec(e.tkpx(tid),e.tkpy(tid),e.tkpz(tid));
			assert(vid >= 0 && vid < nvtx);
		
			
			// remove tracks in a cone around the photon direction
			if ( params_.removeTracksInCone ) {
			  float dr1 = tkPVec.DeltaR(p1.p4(e.vtxx(vid),e.vtxy(vid),e.vtxz(vid)).Vect()); 
			  float dr2 = tkPVec.DeltaR(p2.p4(e.vtxx(vid),e.vtxy(vid),e.vtxz(vid)).Vect()); 
			  if ( dr1 < params_.coneSize  || dr2 < params_.coneSize)
			    continue;
			}
		
	
			TVector2 tkPtVec = tkPVec.XYvector();
			double tkPt = tkPtVec.Mod();
			const double modpt = tkPt > e.tkPtErr(tid) ? tkPt - e.tkPtErr(tid)  : 0.;
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
			sumweight_[vid] += tkWeight;
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
		
		sphers_[vid] *= 1./sumpr_[vid];
		
		TVectorD eigVals(3);
		eigVals = TMatrixDSymEigen(sphers_[vid]).GetEigenValues();

		spher_[vid] = 1.5 * (eigVals[1]+eigVals[2]);
		tspher_[vid] = 2. * eigVals[1] / (eigVals[0]+eigVals[1]);
		aplan_[vid] = 1.5 * eigVals[2];

		threejetC_[vid] = 3. * (eigVals[0]*eigVals[1] + eigVals[0]*eigVals[2] + eigVals[1]*eigVals[2]);
		fourjetD_[vid] = 27. * eigVals[0]*eigVals[1]*eigVals[2];


		diPhotonPt_[vid] = diPhoton_[vid].Vect().XYvector();
		vtxPt_[vid]      = vtxP_[vid].XYvector();
		diPhotonPz_[vid]   = diPhoton_[vid].Vect().Pz();
		
 		sort(tksPt_[vid].begin(), tksPt_[vid].end(), greater<double>());
		
		acosA_[vid] =  	acos(vtxPt_[vid].Unit() * diPhotonPt_[vid].Unit());
		ptasym_[vid] = 	(vtxPt_[vid].Mod() - diPhotonPt_[vid].Mod())/(vtxPt_[vid].Mod() + diPhotonPt_[vid].Mod());

		ptmax_ [vid] = 	tksPt_[vid][0];
		ptmax3_[vid] = 	accumulate(tksPt_[vid].begin(),tksPt_[vid].begin() + min(tksPt_[vid].size(),(size_t)3), 0.0) ;
		thrust_[vid] = 	ptbal_[vid]/sumpt_[vid];
		
		ptratio_[vid] = 	vtxPt_[vid].Mod()/diPhotonPt_[vid].Mod();
		pzasym_[vid] = 	fabs( (vtxP_[vid].Pz() - diPhotonPz_[vid])/(vtxP_[vid].Pz() + diPhotonPz_[vid]) );
		
		awytwdasym_[vid] = (sumawy_[vid]-sumtwd_[vid])/(sumawy_[vid]+sumtwd_[vid]);
		
		preselection_.push_back(vid);
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
