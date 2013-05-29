#include "TMath.h"
#include "TSpline.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "RooRealVar.h"
#include "Math/IFunction.h"

#include <algorithm>

#include "../interface/CategoryOptimizer.h"
#include "../interface/NaiveCategoryOptimization.h"

class CutAndCountModelBuiler;


// ------------------------------------------------------------------------------------------------
double NaiveCutAndCountFomProvider::operator() ( std::vector<AbsModel *> sig, std::vector<AbsModel *> bkg) const
{
	double fom = 0.;
	int ncat = sig[0]->getNcat();
	for(int icat=0.; icat<ncat; ++icat) {
		/// std::cout << "NaiveCutAndCountFomProvider " << icat;
		double sig2 = 0.;
		double nbkg = 0.;
		std::vector<double> yields;
		for(std::vector<AbsModel *>::const_iterator isig = sig.begin(); isig!=sig.end(); ++isig) {
			double nisig = (*isig)->getCategoryYield(icat);
			//// std::cout << " " << nisig;
			yields.push_back(nisig);
			sig2 += nisig*nisig;
		}
		for(std::vector<AbsModel *>::const_iterator ibkg = bkg.begin(); ibkg!=bkg.end(); ++ibkg) {
			double nibkg = (*ibkg)->getCategoryYield(icat);  // FIXME width factor
			yields.push_back(nibkg);
			nbkg += nibkg;
			//// std::cout << " " << nibkg;
		}
		double den = ( denom_ ? (*denom_)(&yields[0]) : nbkg );
		fom -= sig2 / den;
		//// std::cout << " " << fom << std::endl;
	}
	return -sqrt(-fom);
};

// ------------------------------------------------------------------------------------------------
double PoissonCutAndCountFomProvider::operator() ( std::vector<AbsModel *> sig, std::vector<AbsModel *> bkg) const
{
	double qA = 0.;
	int ncat = sig[0]->getNcat();
	for(int icat=0.; icat<ncat; ++icat) {
		double nsig = 0.;
		double nbkg = 0.;
		for(std::vector<AbsModel *>::const_iterator isig = sig.begin(); isig!=sig.end(); ++isig) {
			nsig += (*isig)->getCategoryYield(icat);
		}
		for(std::vector<AbsModel *>::const_iterator ibkg = bkg.begin(); ibkg!=bkg.end(); ++ibkg) {
			nbkg += (*ibkg)->getCategoryYield(icat);  // FIXME width factor
		}
		qA -= 2.*(log( TMath::Poisson(nsig+nbkg, nsig+nbkg) ) - log( TMath::Poisson(nsig+nbkg, nbkg) ));
	}
	return -sqrt(-qA);
}

//////////// // -----------------------------------------------------------------------------------------------------------------------------------------
//////////// class  SignalWidthFactor {
//////////// public:
//////////// 	SignalWidthFactor( TString name, TGraph * sigW, TH1 * sigCdf) : 
//////////// 		sigWFunc_(new GraphToTF1(name,sigW)), sigCdf_(new HistoToTF1(name,sigCdf) ) {
//////////// 		sigWTF1_ = new TF1("sigTF1",sigWFunc_,-1.,1,0,"HistoToTF1");
//////////// 	};
//////////// 	double operator() (double a, double b) {
//////////// 		double x[] = {a,b};
//////////// 		return sigWTF1_->Integral(a,b) / ((*sigCdf_)(&a,0) - (*sigCdf_)(&b,0));
//////////// 	};
//////////// 	
//////////// 	GraphToTF1 * sigWFunc_;
//////////// 	TF1 * sigWTF1_;
//////////// 	HistoToTF1 * sigCdf_;
//////////// };
//////////// 
