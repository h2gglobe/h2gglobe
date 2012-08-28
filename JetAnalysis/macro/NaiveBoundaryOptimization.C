#include "TSpline.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "RooRealVar.h"
#include "Math/IFunction.h"

#include <algorithm>

// -----------------------------------------------------------------------------------------------------------------------------------------
class  HistoToTF1 {
public:
	HistoToTF1( TString name, TH1 * g ) { sp_ = new TSpline3(g,name); };
	double operator() (double *x, double *p) {
		return sp_->Eval( x[0] );
	};
	TSpline * sp_;
};

// -----------------------------------------------------------------------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------------------------------------------------------------------
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


// -----------------------------------------------------------------------------------------------------------------------------------------
class PieceWise2DSignif : public ROOT::Math::IBaseFunctionMultiDim
{
public:
	PieceWise2DSignif(int nbound, TF2 * sigInt, TF2 * bkg1Int, TF2 * bkg2Int, TF3 * fomden, double xcutoff=2.5e-2, 
			  double ycutoff=2.5e-2, double sigNorm=0., double bkg1Norm=0., double bkg2Norm=0.) : 
		nbound_(nbound), sigInt_(sigInt), bkg1Int_(bkg1Int), bkg2Int_(bkg2Int), fomden_(fomden) {
			cutoffs_[0] = xcutoff;
			cutoffs_[1] = ycutoff;
			double xmin =  sigInt_->GetXmin(), xmax =  sigInt_->GetXmax(), ymin =  sigInt_->GetYmin(), ymax =  sigInt_->GetYmax();
			std::cout << "PieceWise2DSignif "
				  << "xmin " << xmin << "\n"
				  << "xmax " << xmax << "\n"
				  << "ymin " << ymin << "\n" 
				  << "ymax " << ymax << "\n" 
				  << "xcutoff " << xcutoff << "\n" 
				  << "ycutoff " << ycutoff
				  << std::endl;
			sigNorm_  =( sigNorm == 0. ? sigInt_->Integral(xmin,xmax,ymin,ymax) : sigNorm );
			bkg1Norm_ =( bkg1Norm == 0. ? bkg1Int_->Integral(xmin,xmax,ymin,ymax) : bkg1Norm ); 
			bkg2Norm_ =( bkg2Norm == 0. ? bkg2Int_->Integral(xmin,xmax,ymin,ymax) : bkg2Norm );
			std::cout << "sigNorm " << sigNorm_  << "\n"
				  << "bkg1Norm " << bkg1Norm_ << "\n"
				  << "bkg2Norm " << bkg2Norm_
				  << std::endl;

		};
	
	double operator() (double *x, double *p) const {
		
		double nsig2 = 0.;
		double firstxb= x[0];		
		double firstyb= x[nbound_];
		double lastxb= firstxb;
		double lastyb= firstyb;
		double lastSigInt=0., lastBkg1Int=0., lastBkg2Int=0.;
		
		std::ostream_iterator< double > output( cout, "," );
		std::copy(&x[0],&x[2*nbound_-1],output); 
		std::cout << std::endl;
		
		double xmin =  sigInt_->GetXmin(), xmax =  sigInt_->GetXmax(), ymin =  sigInt_->GetYmin(), ymax =  sigInt_->GetYmax();
			
		for(int ii=1; ii<nbound_; ++ii) {
			int jj=ii;
			double newxb = lastxb;
			double newyb = lastyb;
			
			for(; jj<nbound_; ++jj) {
				newxb -= x[jj];
				newyb -= x[nbound_+jj];
				if( lastxb > newxb + p[0] && lastyb > newyb + p[1] ) {
					break;
				}
			}
			if( newxb < xmin ) { newxb = xmin; }
			if( newyb < ymin ) { newyb = ymin; }
			
			std::cout << "ibound " << jj << " (" << firstxb << "," << firstyb << ") "<< " (" << newxb << "," << newyb << ") " << nsig2;

			if( lastxb > newxb + p[0] && lastyb > newyb + p[1] ) {
				float sigInt  = sigInt_->Integral(lastxb,firstxb,lastyb,firstyb)  / sigNorm_;
				float bkg1Int = bkg1Int_->Integral(lastxb,firstxb,lastyb,firstyb) / bkg1Norm_;
				float bkg2Int = bkg2Int_->Integral(lastxb,firstxb,lastyb,firstyb) / bkg2Norm_;
				
				float sig  = sigInt  - lastSigInt;
				float bkg1 = bkg2Int - lastBkg1Int;
				float bkg2 = bkg2Int - lastBkg2Int;

				float num = sig*sig;
				float den = fomden_->Eval(sig,bkg1,bkg2);
				
				lastSigInt = sigInt;
				lastBkg1Int = bkg1Int;
				lastBkg2Int = bkg2Int;
				
				nsig2 -= num/den;
			}
			ii = jj;
			std::cout << " " << nsig2 << std::endl;
			lastxb = newxb;
			lastyb = newyb;
		}
		return nsig2;
        };

	virtual double DoEval(const double * x) const { 
		std::vector<double> xv(x,x+2*nbound_);
		std::vector<double> pv(cutoffs_,cutoffs_+2);
		return this->operator()(&xv[0],&pv[0]); 
	}; 
	
	virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const { return new PieceWise2DSignif(nbound_,sigInt_,bkg1Int_,bkg2Int_,fomden_,cutoffs_[0],cutoffs_[1]); }

	virtual unsigned int NDim() const { return 2*nbound_; }; 

	int nbound_;
	TF2 * sigInt_; 
	TF2 * bkg1Int_, * bkg2Int_;
	TF3 * fomden_;
	double cutoffs_[2];
	double sigNorm_, bkg1Norm_, bkg2Norm_;

};
