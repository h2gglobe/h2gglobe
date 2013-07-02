#ifndef _FunctionHelpers_h_
#define _FunctionHelpers_h_

#include "Math/IFunction.h"
#include "TF1.h"
#include "TSpline.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

#include <list>
#include <set>
#include <algorithm>


// -----------------------------------------------------------------------------------------------
class HistoConverter : public ROOT::Math::IBaseFunctionMultiDim {
public:
	HistoConverter() : sp_(0), hist_(0), g_(0) {};
	
	~HistoConverter() {
		if( sp_ )   { delete sp_; }
		if( hist_ ) { delete hist_; }
		if( g_ )    { delete g_; }
	};

	virtual double operator() (double *x, double *p) = 0;

	double DoEval(const double * x) const { 
		std::vector<double> xv(3,0.);
	}
		
	ROOT::Math::IBaseFunctionMultiDim * Clone() const { return clone(); }
	
	double eval(double x) { return (*this)(&x,0); };
	virtual HistoConverter * clone() const = 0;

protected:
	TSpline * sp_;
	TH1  * hist_;
	TGraph * g_;
	
};

// -----------------------------------------------------------------------------------------------
class TF1ToFunctor : public ROOT::Math::IBaseFunctionMultiDim
{
public:
	TF1ToFunctor( TF1 * tf1 ) : tf1_(tf1) {};
	
	double DoEval(const double * x) const { 
		std::vector<double> xv(3,0.);
		for(int ii=0; ii<NDim(); ++ii) { xv[ii] = x[ii]; }
		return  tf1_->Eval(x[0],x[1],x[3]);
	}; 
	
	ROOT::Math::IBaseFunctionMultiDim * Clone() const { return new TF1ToFunctor(*this); }

	unsigned int NDim() const { return tf1_->GetNdim(); };
	
private:
	TF1 * tf1_;
};



// -----------------------------------------------------------------------------------------------
template<class IntegratorT> class HistoIntegrator
{
public:
	HistoIntegrator(HistoConverter * integrand, double tolerance=0.01) : 
		integrand_(integrand)
	{
		integrator_.SetFunction(*integrand_);
		integrator_.SetRelTolerance(tolerance);
	};
	
	double operator() (double *x, double *y) {
		return integrator_.Integral(x,y);
	};

	const HistoConverter * integrand() { return integrand_; };
	IntegratorT integrator() { return integrator_; };
	
private:
	IntegratorT integrator_;
	HistoConverter * integrand_;
};

// -----------------------------------------------------------------------------------------------
class SparseToTF1 : public HistoConverter 
{
public:
	SparseToTF1(THnSparse * sparse) : hsparse_(sparse)
	{};
		
	double operator() (double *x, double *p) {
		return hsparse_->GetBinContent( hsparse_->GetBin(x) );
	};
	
	unsigned int NDim() const { return hsparse_->GetNdimensions(); };
	
	~SparseToTF1() {
		delete hsparse_;
	};
	
	const THnSparse * getSparseHisto() { return hsparse_; };
	
	HistoConverter * clone() const { new SparseToTF1(*this); };
	
protected:
	THnSparse * hsparse_;

};

// -----------------------------------------------------------------------------------------------
template<class IntegratorT> class HistoMultiDimCdf : public SparseToTF1
{
public:
	typedef HistoIntegrator<IntegratorT> histo_integrator_t;
	
	HistoMultiDimCdf(THnSparse * integrand, const double * endpoint, 
			 double tolerance=0.01) : 
		SparseToTF1(0),
		integrand_(integrand),
		endpoint_(endpoint,endpoint+integrand->GetNdimensions()),
		histo_integrator_(&integrand_,tolerance) {
		
		hsparse_ = (THnSparse*)integrand->Clone(Form("%s_cdf", integrand->GetName()));
		hsparse_->Reset();
		std::vector<int> idx(integrand->GetNdimensions());
		std::vector<double> coord(integrand->GetNdimensions());
		double norm = integrand->GetWeightSum();
		for(int ii=0; ii<integrand->GetNbins(); ++ii) {
			for(int idim=0; idim<idx.size(); ++idim) {
				TAxis * iaxis = integrand->GetAxis(idim);
				coord[idim] = iaxis->GetBinLowEdge(idx[idim]);
			}
			hsparse_->Fill( &coord[0], histo_integrator_(&coord[0],&endpoint_[0])/norm );
		}
	};

private:
	SparseToTF1 integrand_;
	std::vector<double> endpoint_;
	std::vector<double> cache_;
	double cacheval_;
	histo_integrator_t histo_integrator_;
};


// -----------------------------------------------------------------------------------------------
class  HistoToTF1 : public HistoConverter {
public:
	HistoToTF1( TString name, TH1 * g ) { sp_ = new TSpline3(g,name); hist_ = g; };
	double operator() (double *x, double *p) {
		return sp_->Eval( x[0] );
	};
	
	unsigned int NDim() const { return 1; }
	
	HistoConverter * clone() const { new HistoToTF1(*this); };
};

// -----------------------------------------------------------------------------------------------
class  GraphToTF1 : public HistoConverter {
public:
	GraphToTF1( TString name, TGraph * g) : capped_(false)
		{ sp_ = new TSpline3(name,g); };
	GraphToTF1( TString name, TGraph * g, double xmin, double valmin, double xmax, double valmax ) :
		capped_(true), xmin_(xmin), valmin_(valmin), xmax_(xmax), valmax_(valmax)
		{ sp_ = new TSpline3(name,g); };
	double operator() (double *x, double *p) {
		if( capped_ ) {
			if( x[0] <= xmin_ ) { return valmin_; }
			if( x[0] >= xmax_ ) { return valmax_; }
		}
		return sp_->Eval( x[0] );
	};

	unsigned int NDim() const { return 1; }
	HistoConverter * clone() const { new GraphToTF1(*this); };

private:
	bool capped_;
	double xmin_, valmin_, xmax_, valmax_;
	
};


// -----------------------------------------------------------------------------------------------
class  LinGraphToTF1 : public HistoConverter {
public:
	LinGraphToTF1( TString name, TGraph * g) : capped_(false)
		{ g_ = (TGraph*)g->Clone(name); };
	LinGraphToTF1( TString name, TGraph * g, double xmin, double valmin, double xmax, double valmax ) :
		capped_(true), xmin_(xmin), valmin_(valmin), xmax_(xmax), valmax_(valmax)
		{ g_ = (TGraph*)g->Clone(name); };
	double operator() (double *x, double *p) {
		double val = g_->Eval( x[0] );
		if( capped_ ) {
			if( x[0] <= xmin_ || val < valmin_ ) { return valmin_; }
			if( x[0] >= xmax_ || val > valmax_ ) { return valmax_; }
		}
		return val;
	};

	unsigned int NDim() const { return 1; }
	HistoConverter * clone() const { new LinGraphToTF1(*this); };
	
private:
	bool capped_;
	double xmin_, valmin_, xmax_, valmax_;
	
};

// ------------------------------------------------------------------------------------------------
class  SimpleHistoToTF2 : public HistoConverter {
public:
	SimpleHistoToTF2( TString name, TH2 * g ) { hist_ = g; };
	double operator() (double *x, double *p) {
		return ((TH2*)hist_)->GetBinContent( ((TH2*)hist_)->FindBin(x[0],x[1]) );
	};
	
	unsigned int NDim() const { return 2; }
	HistoConverter * clone() const { new SimpleHistoToTF2(*this); };
};


// ------------------------------------------------------------------------------------------------
TH1 * integrate1D(TH1 * h, bool normalize=true);
TH2 * integrate2D(TH2 * h, bool normalize=true);
HistoConverter * cdfInv(TH1 * h, double min, double max);
HistoConverter * cdf(TH1 * h, double min, double max);


// -----------------------------------------------------------------------------------------------
class IntegrationNode 
{
public:
	IntegrationNode(int id, std::vector<double> coord, double w) : 
		id_(id), coord_(coord), weight_(w), sum_(weight_), hasSum_(false)
		{};
	
	class weakLess {
	public:
		bool operator() (const IntegrationNode * a, const IntegrationNode *b) {
			return  this->operator()(*a,*b);
		};

		bool operator() (const IntegrationNode & a, const IntegrationNode &b) {
			for(size_t icoord=0; icoord<a.coord_.size(); ++icoord) {
				if ( a.coord_[icoord] < b.coord_[icoord] ) {
					return true;
				} else if ( a.coord_[icoord] > b.coord_[icoord] ) {
					return false;
				} 
			}
			return false;
		};
		
	};
	
	class strictLess {
	public:
		bool operator() (const IntegrationNode * a, const IntegrationNode *b) {
			return  this->operator()(*a,*b);
		};

		bool operator() (const IntegrationNode & a, const IntegrationNode &b) {
			for(size_t icoord=0; icoord<a.coord_.size(); ++icoord) {
				if ( a.coord_[icoord] > b.coord_[icoord] ) {
					return false;
				}
			}
			return true;
		};
	};
	
	void fill(double w) { weight_+=w; sum_+=w; };
	
	void print(std::ostream & out) const {
		out << "IntegrationNode " << id_ << " " << weight_;
		for(int ic=0; ic<coord_.size(); ++ic) {
			out << " " << coord_[ic];
		}
		out << "\n " << sum_ ;
		for( std::list<IntegrationNode *>::const_iterator ichild=children_.begin(); ichild!=children_.end(); ++ichild) {
			out << " " << (*ichild)->id_ ;
			/// std::vector<double> & ccoord = (*ichild)->coord_;
			/// out << "(";
			/// for(int jc=0; jc<coord_.size(); ++jc) {
			/// 	out << " " << ccoord[jc];
			/// }
			/// out << ") ";
		}
		out << "\n";
	};

	void addChild(IntegrationNode * node) { 
		//// children_.push_back(node);
		sum_ += node->weight_;
	};
	
	double sumEntries() { 
		//// if( hasSum_ ) { return sum_; }
		//// for( std::list<IntegrationNode *>::iterator ichild=children_.begin(); ichild!=children_.end(); ++ichild) {
		//// 	sum_ += (*ichild)->sumEntries();
		//// }
		//// hasSum_ = true;
		return sum_;
	};
	
	int id() { return id_; };
	const std::vector<double> coord() { return coord_; };
private:
	int id_;
	std::vector<double> coord_;
	double weight_;
	double sum_;
	bool hasSum_;
	
	std::list<IntegrationNode *> children_;
	
};

// -----------------------------------------------------------------------------------------------
class IntegrationWeb : public std::set<IntegrationNode*,IntegrationNode::weakLess>
{
public:
	IntegrationWeb() : scale_(1.) {};

	~IntegrationWeb() {
		for(iterator it=begin(); it!=end(); ++it) {
			delete *it;
		}
	};
	
	void link() {
		for(reverse_iterator it=rbegin(); it!=rend(); ++it) {
			reverse_iterator jt=it; ++jt;
			while(jt!=rend()) {
				if( grid_.size() == 1 || IntegrationNode::strictLess()(**jt,**it) ) {
					(*jt)->addChild(*it);
				}
				++jt;
			}
		}
	};
	
	int fill(const double* coord, double w) {
		IntegrationNode * node = get(coord,false);
		node->fill(w);
		return node->id();
	};

	void integrate() {
		for(reverse_iterator it=rbegin(); it!=rend(); ++it) {
			(*it)->sumEntries();
		}
	};
	
	double getIntegral(const double* coord) {
		return get(coord)->sumEntries()*scale_;
	};
	
	double getIntegral(const double* a, const double* b) {
		return getIntegral(b) - getIntegral(a);
	};
	
	void print( std::ostream & out ) const {
		out << "Integration grid:\n";
		for (size_t d=0; d <grid_.size(); ++d) {
			out << "[" << d << "]: ";
			for(std::set<double>::iterator it=grid_[d].begin(); it!=grid_[d].end(); ++it) {
				out << *it << " ";
			}
			out << "\n";
		}
		

		for(iterator it=begin(); it!=end(); ++it) {
			(*it)->print(out);
			out << std::endl;
		}
	};
	
	void scale(double x) { scale_ = x; };

protected:
	double scale_;
	std::vector<std::set<double> > grid_;
	
	void volume_coordinates(std::vector<double> & point) {
		for(size_t idim=0; idim<grid_.size(); ++idim) {
			/// std::cout << idim << ":" << point[idim];
			std::set<double>::iterator ibound = grid_[idim].lower_bound(point[idim]);
			if( ibound == grid_[idim].end() ) {
				point[idim] = ( point[idim] <= *(grid_[idim].begin()) ? *(grid_[idim].begin()) : *(grid_[idim].rbegin()) );
			} else if( point[idim] != *ibound ) {
				point[idim] = *(--ibound);
			}
			/// std::cout << "->" << point[idim] << " ";
		}
		/// std::cout << std::endl;
	};
	
	IntegrationNode * get(const double* coord, bool link=true) {
		std::vector<double> volume(coord,coord+grid_.size());
		volume_coordinates(volume);
		IntegrationNode tmp(-1,volume,0.);
		iterator inode = lower_bound(&tmp);
		if( inode == end() || key_comp()(**inode,tmp) || key_comp()(tmp,**inode) ) {
			inode = insert( inode, new IntegrationNode(size(),volume,0.) );
			iterator jnode = inode;
			++jnode;
			if( link ) { 
				while( jnode != end() ) { 
					if( IntegrationNode::strictLess()(**inode,**jnode) ) {
						(*inode)->addChild(*jnode);
					}
					++jnode;
				}
			}
		}
		//// (*inode)->print(std::cout);
		return *inode;
	};
	
	
};

// -----------------------------------------------------------------------------------------------
class SparseIntegrator : public IntegrationWeb, public HistoConverter
{
public:
	SparseIntegrator(THnSparse * integrand,double scale=1.) : hsp_(integrand) {
		/// reserve(integrand->GetNbins());
		scale_ = scale;
		Int_t dim = integrand->GetNdimensions();
		std::vector<int> bins(dim);
		std::vector<double> coord(dim);
		grid_.resize(dim);
		for (Int_t d = 0; d < dim; ++d) {
			const TAxis * axis = integrand->GetAxis(d);
			for(int ibin=0; ibin<=axis->GetNbins(); ++ibin) {
				grid_[d].insert( axis->GetBinLowEdge(ibin) );
			}
			//// for(std::set<double>::iterator it=grid_[d].begin(); it!=grid_[d].end(); ++it) {
			//// 	std::cout << *it << " ";
			//// }
			//// std::cout << std::endl;
		}
		for (Long64_t i = 0; i < integrand->GetNbins(); ++i) {
			double weight = integrand->GetBinContent(i, &bins[0]);
			for (Int_t d = 0; d < dim; ++d) {
				coord[d] = integrand->GetAxis(d)->GetBinLowEdge(bins[d]);
			}
			insert( new IntegrationNode(i,coord,weight) );
		}
		link();
	};

	~SparseIntegrator() {
		/// delete hsp_;
	};

	THnSparse * getIntegrand() { return hsp_; };
	
	double operator() (double *x, double *p) {
		return getIntegral(x);
	};
	
	unsigned int NDim() const { return grid_.size(); };

	HistoConverter * clone() const { return new SparseIntegrator(*this); };

private:
	THnSparse * hsp_;

};


#endif 

