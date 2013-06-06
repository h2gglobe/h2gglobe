#ifndef _FunctionHelpers_h_
#define _FunctionHelpers_h_

#include "Math/IFunction.h"
#include "TF1.h"
#include "TSpline.h"
#include "TH1.h"
#include "TH2.h"

// -----------------------------------------------------------------------------------------------
class HistoConverter {
public:
	HistoConverter() : sp_(0), hist_(0), g_(0) {};
	
	~HistoConverter() {
		if( sp_ )   { delete sp_; }
		if( hist_ ) { delete hist_; }
		if( g_ )    { delete g_; }
	};

	virtual double operator() (double *x, double *p) = 0;
	
	double eval(double x) { return (*this)(&x,0); };
	
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
class  HistoToTF1 : public HistoConverter {
public:
	HistoToTF1( TString name, TH1 * g ) { sp_ = new TSpline3(g,name); hist_ = g; };
	double operator() (double *x, double *p) {
		return sp_->Eval( x[0] );
	};
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
		if( capped_ ) {
			if( x[0] <= xmin_ ) { return valmin_; }
			if( x[0] >= xmax_ ) { return valmax_; }
		}
		return g_->Eval( x[0] );
	};

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
};


// ------------------------------------------------------------------------------------------------
TH1 * integrate1D(TH1 * h, bool normalize=true);
TH2 * integrate2D(TH2 * h, bool normalize=true);
HistoConverter * cdfInv(TH1 * h, double min, double max);
HistoConverter * cdf(TH1 * h, double min, double max);

#endif 

