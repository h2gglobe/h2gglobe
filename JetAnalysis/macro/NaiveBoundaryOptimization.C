#include "TSpline.h"
#include "TH1.h"
#include "TH2.h"
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
class  SimpleHistoToTF2 {
public:
	SimpleHistoToTF2( TString name, TH2 * g ) { sp_ = g; };
	double operator() (double *x, double *p) {
		return sp_->GetBinContent( sp_->FindBin(x[0],x[1]) );
	};
	TH2 * sp_;
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
TH2 * integrate2D(TH2 * h) {
	TH2 * ret= (TH2*)h->Clone( Form("%s_cdf", h->GetName() ) );
	for(int xx=0; xx<ret->GetNbinsX()+1; ++xx) {
		for(int yy=1; yy<ret->GetNbinsY()+1; ++yy) {
			ret->SetBinContent( xx, yy, ret->GetBinContent(xx,yy) + ret->GetBinContent(xx,yy-1) );
		}
	}
	for(int yy=0; yy<ret->GetNbinsY()+1; ++yy) {
		for(int xx=1; xx<ret->GetNbinsX()+1; ++xx) {
			ret->SetBinContent( xx, yy, ret->GetBinContent(xx,yy) + ret->GetBinContent(xx-1,yy) );
		}
	}
	ret->Scale( 1./h->Integral() );
	return ret;
}

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
	PieceWise2DSignif(int nbound, bool isCdf, TF2 * sigInt, TF2 * bkg1Int, TF2 * bkg2Int, TF3 * fomden, double xcutoff=2.5e-2, 
			  double ycutoff=2.5e-2, double sigNorm=0., double bkg1Norm=0., double bkg2Norm=0.) : 
		nbound_(nbound), isCdf_(isCdf), sigInt_(sigInt), bkg1Int_(bkg1Int), bkg2Int_(bkg2Int), fomden_(fomden) {
			cutoffs_[0] = xcutoff;
			cutoffs_[1] = ycutoff;
			double xmin =  sigInt_->GetXmin(), xmax =  sigInt_->GetXmax(), ymin =  sigInt_->GetYmin(), ymax =  sigInt_->GetYmax();
			///// std::cout << "PieceWise2DSignif "
			///// 	  << "xmin " << xmin << "\n"
			///// 	  << "xmax " << xmax << "\n"
			///// 	  << "ymin " << ymin << "\n" 
			///// 	  << "ymax " << ymax << "\n" 
			///// 	  << "xcutoff " << xcutoff << "\n" 
			///// 	  << "ycutoff " << ycutoff
			///// 	  << std::endl;
			if( ! isCdf ) {
				sigNorm_  =( sigNorm == 0. ? sigInt_->Integral(xmin,xmax,ymin,ymax) : sigNorm );
				bkg1Norm_ =( bkg1Norm == 0. ? bkg1Int_->Integral(xmin,xmax,ymin,ymax) : bkg1Norm ); 
				bkg2Norm_ =( bkg2Norm == 0. ? bkg2Int_->Integral(xmin,xmax,ymin,ymax) : bkg2Norm );
			}
			//// std::cout << "sigNorm " << sigNorm_  << "\n"
			//// 	  << "bkg1Norm " << bkg1Norm_ << "\n"
			//// 	  << "bkg2Norm " << bkg2Norm_
			//// 	  << std::endl;

		};
	
	double operator() (double *x, double *p) const {
		
		double nsig2 = 0.;
		double firstxb= x[0];		
		double firstyb= x[nbound_];
		double lastxb= firstxb;
		double lastyb= firstyb;
		double xmin =  sigInt_->GetXmin(), ymin =  sigInt_->GetYmin();
		
		double lastSigInt,lastBkg1Int,lastBkg2Int;

		if( isCdf_ ) {
			lastSigInt=sigInt_->Eval(firstxb,firstyb) ;
			lastBkg1Int=bkg1Int_->Eval(firstxb,firstyb);
			lastBkg2Int=bkg2Int_->Eval(firstxb,firstyb);
		} else {
			lastSigInt=sigInt_->Integral(xmin,firstxb,ymin,firstyb) / sigNorm_; 
			lastBkg1Int=bkg1Int_->Integral(xmin,firstxb,ymin,firstyb) / bkg1Norm_; 
			lastBkg2Int=bkg2Int_->Integral(xmin,firstxb,ymin,firstyb) / bkg2Norm_;
		}
		/// std::ostream_iterator< double > output( cout, "," );
		/// std::copy(&x[0],&x[2*nbound_],output); 
		/// std::cout << std::endl;
		
		/// std::cout << xmin<< " " << firstxb<< " " << ymin<< " " << firstyb << " " << lastSigInt << " " << lastBkg1Int << " " <<  lastBkg2Int << std::endl;
		/// nsig2 = -lastSigInt*lastSigInt / fomden_->Eval(lastSigInt,lastBkg1Int,lastBkg2Int);
		
		for(int ii=1; ii<nbound_-1; ++ii) { // ignore the last boundary since Minuit does not float it
			int jj=ii;
			double newxb = lastxb;
			double newyb = lastyb;
			
			for(; jj<nbound_-1; ++jj) { // ignore the last boundary
				newxb -= x[jj];
				newyb -= x[nbound_+jj];
				if( lastxb > newxb + p[0] && lastyb > newyb + p[1] ) {
					break;
				}
			}
			if( newxb < xmin ) { newxb = xmin; }
			if( newyb < ymin ) { newyb = ymin; }
			
			/// std::cout << "ibound " << jj << " (" << firstxb << "," << firstyb << ") "<< " (" << newxb << "," << newyb << ") " << nsig2;

			if( lastxb > newxb + p[0] && lastyb > newyb + p[1] ) {
				float sigInt ;
				float bkg1Int;
				float bkg2Int;

				if( isCdf_ ) { 
					sigInt  = sigInt_->Eval(newxb,newyb) ;
					bkg1Int = bkg1Int_->Eval(newxb,newyb);
					bkg2Int = bkg2Int_->Eval(newxb,newyb);
				} else {
					sigInt  = sigInt_->Integral(xmin,newxb, ymin,newyb) / sigNorm_;
					bkg1Int = bkg1Int_->Integral(xmin,newxb,ymin,newyb) / bkg1Norm_;
					bkg2Int = bkg2Int_->Integral(xmin,newxb,ymin,newyb) / bkg2Norm_;
				}
				float sig  = lastSigInt  - sigInt  ;
				float bkg1 = lastBkg1Int - bkg1Int ;
				float bkg2 = lastBkg2Int - bkg2Int ;

				/// std::cout << " " << lastSigInt << " " << sigInt << sig 
				/// 	  << " " << lastBkg1Int << " "<< bkg1Int << bkg1 
				/// 	  << " " << lastBkg2Int << " "<< bkg2Int << bkg2
				/// 	;
				
				float num = sig*sig;
				float den = fomden_->Eval(sig,bkg1,bkg2);
				
				/// std::cout << " " << num << " " << den;

				lastSigInt = sigInt;
				lastBkg1Int = bkg1Int;
				lastBkg2Int = bkg2Int;
				
				nsig2 -= num/den;
			}
			ii = jj;
			//// std::cout << " " << nsig2 << std::endl;
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
	
	virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const { return new PieceWise2DSignif(nbound_,isCdf_,sigInt_,bkg1Int_,bkg2Int_,fomden_,cutoffs_[0],cutoffs_[1]); }

	virtual unsigned int NDim() const { return 2*nbound_; }; 

	int nbound_;
	bool isCdf_;
	TF2 * sigInt_; 
	TF2 * bkg1Int_, * bkg2Int_;
	TF3 * fomden_;
	double cutoffs_[2];
	double sigNorm_, bkg1Norm_, bkg2Norm_;

};
