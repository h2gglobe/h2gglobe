#include "TSpline.h"
#include "TH1.h"
#include "TF1.h"
#include "RooRealVar.h"
#include "Math/IFunction.h"

#include <algorithm>

class  HistoToTF1 {
public:
	HistoToTF1( TString name, TH1 * g ) { sp_ = new TSpline3(g,name); };
	double operator() (double *x, double *p) {
		return sp_->Eval( x[0] );
	};
	TSpline * sp_;
};


class  HistoIntegral {
public:
	HistoIntegral( TString name, TH1 * g ) { 
		func_ = new HistoToTF1(name, g);
		tf1_  = new TF1(name,func_,g->GetXaxis()->GetXmin(),g->GetXaxis()->GetXmax(),0,"HistoToTF1");
	};
	double operator() (double *x, double *p) {
		double min = x[0];
		double max = x[1];
		double & cutoff = p[0];
		// if( min > max ) { std::swap(min,max); }
		if( max < min+cutoff ) { max = min+cutoff; }
		return tf1_->Integral(min,max);
	};
	HistoToTF1 * func_;
	TF1 * tf1_;
};

class PieceWiseSignif : public ROOT::Math::IBaseFunctionMultiDim
{
public:
	PieceWiseSignif(int nbound, TF1 * sigInt, TF1 * bkgInt, double cutoff=2.5e-2) : 
		nbound_(nbound), sigInt_(sigInt), bkgInt_(bkgInt), cutoff_(cutoff) {};
	
	double operator() (double *x, double *p) const {
		/// std::sort(x,x+nbound_);
		
		double nsig2 = 0.;
		double lastb= x[0];		
		/// std::cout << "x[0] " << x[0] << std::endl;
		for(int ii=1; ii<nbound_; ++ii) {
			int jj=ii;
			double newb = lastb;
			
			for(; jj<nbound_; ++jj) {
				/// std::cout << "x[" << jj << "]" << x[jj] << std::endl;
				newb -= x[jj];
				if( lastb > newb + p[0] ) {
					break;
				}
			}
			if( newb < sigInt_->GetXmin() ) { newb = sigInt_->GetXmin(); }
			
			if( lastb > newb + p[0] ) {
				/// std::cout << "lastb " << lastb << " newb " << newb << " " << ii << " " << jj << std::endl;
				float csig = sigInt_->Integral(newb,lastb);
				csig *= csig;
				csig /= bkgInt_->Integral(newb,lastb);
				nsig2 -= csig;
			}
			ii = jj;
			lastb = newb;
		}
		
		//// std::sort(x,x+nbound_);
		//// for(int ii=nbound_-1; ii>=0; --ii) {
		//// 	int jj=ii-1;
		//// 	for(; jj>=0; --jj) {
		//// 		if( x[ii] >= x[jj] + p[0] ) {
		//// 			break;
		//// 		}
		//// 	}
		//// 	
		//// 	float csig = sigInt_->Integral(x[jj],x[ii]);
		//// 	csig *= csig;
		//// 	csig /= bkgInt_->Integral(x[jj],x[ii]);
		//// 	nsig2 -= csig;
		//// 	
		//// 	ii = jj;
		//// }
                return nsig2;
        };
	
	virtual double DoEval(const double * x) const { 
		std::vector<double> xv(x,x+nbound_);
		std::vector<double> pv(1,cutoff_);
		return this->operator()(&xv[0],&pv[0]); 
	}; 
	
	virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const { return new PieceWiseSignif(nbound_,sigInt_,bkgInt_,cutoff_); }

	virtual unsigned int NDim() const { return nbound_; }; 

	int nbound_;
	TF1 * sigInt_; 
	TF1 * bkgInt_;
	double cutoff_;
			      
};

///// class RooPieceWiseSignif : public RooAbsReal
///// {
///// public:
///// 	RooPieceWiseSignif(PieceWiseSignif * pwsig, RooArgList & vars, double cutoff=2.5e-2) 
///// 		: pwsig_(pwsig), vars_(vars), cutoff_(cutoff)
///// 	{
///// 		
///// 	}
///// 	
///// 	double pippo() { 
///// 		std::cout << "HERE" << std::endl;
///// 		return evaluate();
///// 	}
///// 	
///// 	virtual double evaluate() const {
///// 		std::vector<double> x(vars_.getSize()); 
///// 		double cutoff=cutoff_;
///// 		std::cout << "HERE" << std::endl;
///// 		for(int ii=0; ii<vars_.getSize(); ++ii) {
///// 			std::cout << dynamic_cast<RooAbsReal*>(vars_.at(ii))->getVal() << std::endl;
///// 			x[ii] = dynamic_cast<RooAbsReal*>(vars_.at(ii))->getVal();
///// 			std::cout << x[ii] << std::endl;
///// 		}
///// 		return (*pwsig_)(&x[0],&cutoff);
///// 	};
///// 
///// 	virtual TObject* clone(const char* newname) const { return new RooPieceWiseSignif(pwsig_,vars_,cutoff_); };
///// 	virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) { return false; }
///// 	virtual void writeToStream(std::ostream& os, Bool_t compact) const { return; };
///// 
///// private:
///// 	PieceWiseSignif * pwsig_;
///// 	RooArgList & vars_;
///// 	double cutoff_;
///// };

