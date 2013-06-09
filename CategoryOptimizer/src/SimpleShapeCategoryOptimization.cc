#include "TMath.h"
#include "TSpline.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooArgSet.h"
#include "RooProduct.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooMinimizer.h"
#include "RooAddition.h"
#include "Math/IFunction.h"
#include "TCanvas.h"

#include <algorithm>

#include "../interface/CategoryOptimizer.h"
#include "../interface/SimpleShapeCategoryOptimization.h"


// ------------------------------------------------------------------------------------------------
SecondOrderModel::SecondOrderModel(std::string name,
				   RooRealVar * x, AbsModel::type_t type, RooRealVar * mu, shape_t shape) : 
	name_(name),
	x_(x), mu_(mu),
	shape_(shape),
	likeg_(0) {
	type_ = type; 
	if( shape_ == automatic ) {
		shape_ = ( type_ == AbsModel::sig ? gaus : expo );
	}
	if( shape_ == expo ) {
		likeg_ = new TF1(Form("likeg_%s",name.c_str()),"[0]-1./x+[1]*exp(-[1]*x)/(1.-exp(-[1]*x))",0.,100.);
	}
};

// ------------------------------------------------------------------------------------------------
SecondOrderModel::~SecondOrderModel()
{
	if( likeg_ ) { delete likeg_; }
}

// ------------------------------------------------------------------------------------------------
TH1 * SecondOrderModelBuilder::getPdf(int idim)
{
	if( ranges_.size() == 1 ) {
		return (TH1*)pdf_->Clone();
	} else if( ranges_.size() == 2 ) {
		return ( idim == 0 ? ((TH2*)pdf_)->ProjectionX() : ((TH2*)pdf_)->ProjectionY() );
	} else if( ranges_.size() == 3 ) {
		return ( idim == 0 ? ((TH3*)pdf_)->ProjectionX() : 
			 ( idim == 1 ? ((TH3*)pdf_)->ProjectionY() : ((TH3*)pdf_)->ProjectionZ() )
			);
	}
	assert(0);
}

// ------------------------------------------------------------------------------------------------
RooAbsPdf * SecondOrderModel::getCategoryPdf(int icat)
{
	return categoryPdfs_[icat];
}

// ------------------------------------------------------------------------------------------------
void SecondOrderModel::buildPdfs()
{
	std::vector<int> catToFix;
	double smallestYield=1e+30;
	double largestSigma=0.;
	
	for(size_t icat=0; icat<categoryYields_.size(); ++icat) {
		if( icat >= categoryPdfs_.size() ) {
			bookShape(icat);
		} else {
			RooAbsPdf * pdf = categoryPdfs_[icat];
			setShapeParams( icat );
		}
		if( categoryYields_[icat] == 0 || categoryRMSs_[icat] == 0 ) {
			catToFix.push_back(icat);
		} else {
			smallestYield = std::min(smallestYield,categoryYields_[icat] );
			largestSigma = std::max(largestSigma,categoryRMSs_[icat]);
		}
	}
	
	for(size_t ifix=0; ifix<catToFix.size(); ++ifix) {
		int icat = catToFix[ifix];
		categoryYields_[icat] = smallestYield/10.;
		categoryRMSs_[icat] = largestSigma*2.;
		setShapeParams( icat );
	}
	/// dump();
}

// ------------------------------------------------------------------------------------------------
void SecondOrderModel::bookShape(int icat) 
{	
	assert( icat == categoryPdfs_.size() );
	RooRealVar * norm = new RooRealVar(Form("%s_cat_model_norm0_%d",name_.c_str(),icat),
					   Form("%s_cat_model_norm0_%d",name_.c_str(),icat),categoryYields_[icat],
					   0.,1e+6);
	owned_.addOwned(*norm);
	
	RooAbsPdf * pdf;
	if( shape_ == gaus ) {
		RooRealVar * mean = new RooRealVar(Form("%s_cat_model_gaus_mean_%d",name_.c_str(),icat),
						   Form("%s_cat_model_gaus_mean_%d",name_.c_str(),icat),
						   categoryMeans_[icat]);
		mean->setConstant(true);
		RooRealVar * sigma = new RooRealVar(Form("%s_cat_model_gaus_sigma_%d",name_.c_str(),icat),
						    Form("%s_cat_model_gaus_sigma_%d",name_.c_str(),icat),
						    categoryRMSs_[icat]);
		sigma->setConstant(true);
		pdf = new RooGaussian( Form("%s_cat_model_gaus_%d",name_.c_str(),icat), Form("%s_cat_model_gaus_%d",name_.c_str(),icat), 
				       *x_, *mean, *sigma );
		
		owned_.addOwned(* mean ), owned_.addOwned(* sigma ); 		
	} else if( shape_ == expo ) {
		likeg_->SetParameters(categoryMeans_[icat]-x_->getMin(),x_->getMax()-x_->getMin());
		double lambdaval=-likeg_->GetX(0.,1./likeg_->GetParameter(0),0.1/likeg_->GetParameter(0));
		RooRealVar * lambda = new RooRealVar(Form("%s_cat_model_expo_lambda_%d",name_.c_str(),icat),
						     Form("%s_cat_model_expo_lambda_%d",name_.c_str(),icat),
						     lambdaval,-1e+30,0.); /// -1./categoryMeans_[icat]
		if( type_ == AbsModel::bkg ){
			lambda->setConstant(false);
			/// lambda->setConstant(true);
		} else {
			lambda->setConstant(true);
		}
		pdf = new RooExponential( Form("%s_cat_model_expo_%d",name_.c_str(),icat), Form("%s_cat_model_expo_%d",name_.c_str(),icat), 
					  *x_, *lambda );
		
		owned_.addOwned(* lambda );
	}
	
	RooExtendPdf * expdf = 0;
	if( mu_ != 0 ) {
		RooProduct * renorm = new RooProduct(Form("%s_cat_model_yield_%d",name_.c_str(),icat),
						     Form("%s_cat_model_yield_%d",name_.c_str(),icat),
						     RooArgList(*norm,*mu_));
		expdf = new RooExtendPdf(Form("%s_cat_model_ext_pdf_%d",name_.c_str(),icat),
					 Form("%s_cat_model_ext_pdf_%d",name_.c_str(),icat),
					 *pdf, *renorm);
		norm->setConstant(true);
		owned_.addOwned(*renorm);
	} else {
		norm->setConstant(false);
		expdf = new RooExtendPdf(Form("%s_cat_model_ext_pdf_%d",name_.c_str(),icat),
					 Form("%s_cat_model_ext_pdf_%d",name_.c_str(),icat),
					 *pdf, *norm);
	}
	categoryNorms_.push_back(norm); 
	categoryPdfs_.push_back(expdf); owned_.addOwned(*pdf); owned_.addOwned(*expdf);
}

// ------------------------------------------------------------------------------------------------
void SecondOrderModel::setShapeParams(int icat)
{
	categoryNorms_[icat]->setVal(categoryYields_[icat]);
	RooArgSet * params = categoryPdfs_[icat]->getParameters((RooArgSet*)0);
	if( shape_ == gaus ) {
		params->setRealValue(Form("%s_cat_model_gaus_mean_%d",name_.c_str(),icat),categoryMeans_[icat]);
		params->setRealValue(Form("%s_cat_model_gaus_sigma_%d",name_.c_str(),icat),categoryRMSs_[icat]);
	} else if( shape_ == expo ) {
		likeg_->SetParameters(categoryMeans_[icat]-x_->getMin(),x_->getMax()-x_->getMin());
		double lambdaval=-likeg_->GetX(0.,1./likeg_->GetParameter(0),0.1/likeg_->GetParameter(0));
		params->setRealValue(Form("%s_cat_model_expo_lambda_%d",name_.c_str(),icat),lambdaval);
	}
}


// ------------------------------------------------------------------------------------------------
void SecondOrderModel::dump() {
	for(int ii=0; ii<categoryYields_.size(); ++ii) {
		std::cout << ii << " " << categoryYields_[ii] << " " << categoryMeans_[ii] 
			  << " "<< categoryRMSs_[ii] << std::endl;
		categoryPdfs_[ii]->Print();
		categoryPdfs_[ii]->getParameters((RooArgSet*)0)->Print("v");
	}
}

// ------------------------------------------------------------------------------------------------
void throwAsimov( double nexp, RooDataHist *asimov, RooAbsPdf *pdf, RooRealVar *x)
{
	asimov->reset();       
	RooArgSet mset(*x);
	pdf->fillDataHist(asimov,&mset,1,false);
   
	for (int i=0 ; i<asimov->numEntries() ; i++) {
		asimov->get(i) ;
		
		// Expected data, multiply p.d.f by nEvents
		Double_t w=asimov->weight()*nexp;
		asimov->set(w,sqrt(w));
	}
   
	Double_t corr = nexp/asimov->sumEntries() ;
	for (int i=0 ; i<asimov->numEntries() ; i++) {
		RooArgSet theSet = *(asimov->get(i));
		asimov->set(asimov->weight()*corr,sqrt(asimov->weight()*corr));
	}

}

// ------------------------------------------------------------------------------------------------
void appendData(RooDataHist & dest, RooDataHist &src)
{
	for (int i=0 ; i<src.numEntries() ; i++) {
		RooArgSet theSet = *(src.get(i));
		dest.add( theSet, src.weight(), src.weightError() );
	}
}

// ------------------------------------------------------------------------------------------------
double SimpleShapeFomProvider::operator() ( std::vector<AbsModel *> sig, std::vector<AbsModel *> bkg) const
{
	static std::vector<RooDataHist> asimovs;
	static std::vector<RooAddPdf> pdfs;

	for(int ipoi=0; ipoi<pois_.size(); ++ipoi ) {
		pois_[ipoi]->setVal(1.);
	}
	for(int ireset=0; ireset<resets_.size(); ++ireset ) {
		resets_[ireset]->setVal(1.);
	}
	
	double fom = 0.;
	int ncat = sig[0]->getNcat();
	
	RooCategory roocat("SimpleShapeFomCat","");
	std::map<std::string,RooAbsData*> catData;
	std::vector<RooRealVar *> normsToFix;
	
	/// WARNING: assuming that signal and background pdfs objects don't change for a given category
	for(int icat=0; icat<ncat; ++icat) {
		if( icat >= asimovs.size() ) {
			asimovs.push_back(RooDataHist(Form("asimovhist_%d",icat),"",*(sig[0]->getX())) );
		}
		
		bool buildPdf = ( icat >= pdfs.size() );
		RooArgList lpdfs, bpdfs;
		
		double ntot = 0.;
		std::vector<double> yields;
		for(std::vector<AbsModel *>::const_iterator isig = sig.begin(); isig!=sig.end(); ++isig) {
			ntot += (*isig)->getCategoryYield(icat);
			//// std::cout << (*isig)->getCategoryPdf(icat)->expectedEvents(0) 
			//// 	  << " " << (*isig)->getCategoryYield(icat) << std::endl;
			if( buildPdf ) {
				lpdfs.add( *((*isig)->getCategoryPdf(icat)) ); 
			}
		}
		for(std::vector<AbsModel *>::const_iterator ibkg = bkg.begin(); ibkg!=bkg.end(); ++ibkg) {
			ntot += (*ibkg)->getCategoryYield(icat);
			//// std::cout << (*ibkg)->getCategoryPdf(icat)->expectedEvents(0) 
			//// 	  << " " << (*ibkg)->getCategoryYield(icat) << std::endl;
			if( buildPdf ) {
				lpdfs.add( *((*ibkg)->getCategoryPdf(icat)) ); 
			}
		}
		
		if( buildPdf ) {
			pdfs.push_back(RooAddPdf(Form("sbpdf_%d",icat),"",lpdfs));
		}
		
		//// std::cout << Form("cat_%d",icat) << std::endl;
		assert( ! roocat.defineType(Form("cat_%d",icat)) );
		throwAsimov( ntot, &asimovs[icat], &pdfs[icat], sig[0]->getX());
		
		/// catData[Form("cat_%d",icat)] = &asimovs[icat];
	}
	
	RooAbsReal * nll;
	std::vector<RooAbsReal*> garbageColl;
	std::vector<RooAbsData*> garbageData;
	if( useRooSimultaneous_ ) {
		RooSimultaneous * roosim = new RooSimultaneous("SimpleShapeFomFit","",roocat);
		garbageColl.push_back(roosim);
		//// RooRealVar * weight = new RooRealVar("weight","weight",1.); 
		//// garbageColl.push_back(weight);
		/// RooDataSet * combData = new RooDataSet("combData","combData",RooArgList(*(sig[0]->getX()),roocat,*weight),"weight");
		RooDataHist * combData = new RooDataHist("combData","combData",RooArgList(*(sig[0]->getX()),roocat));
		garbageData.push_back(combData);
		for(int icat=0; icat<ncat; ++icat) {
			roocat.setLabel(Form("cat_%d",icat));
			//// roocat.Print();
			appendData( *combData, asimovs[icat] );
			//// std::cout << pdfs[icat].expectedEvents(0) << std::endl; 
			roosim->addPdf( pdfs[icat], Form("cat_%d",icat) );
		}
		
	
		nll = roosim->createNLL( *combData, RooFit::Extended(), RooFit::NumCPU(ncpu_) );
		garbageColl.push_back(nll);
	} else { 
		RooArgSet nlls;
		for(int icat=0; icat<ncat; ++icat) {
			RooAbsReal *inll = pdfs[icat].createNLL( asimovs[icat], RooFit::Extended() );
			nlls.add(*inll);
			garbageColl.push_back(inll);
		}
		nll = new RooAddition("nll","nll",nlls);
		garbageColl.push_back(nll);
	}
	
	// S+B fit
	for(int ipoi=0; ipoi<pois_.size(); ++ipoi ) {
		pois_[ipoi]->setVal(1.);
		pois_[ipoi]->setConstant(false);
	}
	RooMinimizer minimsb(*nll);
	minimsb.setMinimizerType(minimizer_.c_str());
	minimsb.setPrintLevel(-1);
	for(int ii=minStrategy_; ii<3; ++ii) {
		minimsb.setStrategy(ii);
		if( ! minimsb.migrad() ) { break; }
	}
	double minNllsb = nll->getVal();

	std::vector<RooPlot *> frames;
	if( debug_ ) {
		for(int icat=0; icat<ncat; ++icat) {
			RooPlot* frame = sig[0]->getX()->frame(RooFit::Title(Form("Category %d/%d",icat,ncat)));
			asimovs[icat].plotOn(frame);
			pdfs[icat].plotOn(frame);
			frames.push_back(frame);
		}
	}
	
	// B-only fit
	for(int ipoi=0; ipoi<pois_.size(); ++ipoi ) {
		pois_[ipoi]->setVal(0.);
		pois_[ipoi]->setConstant(true);
	}
	RooMinimizer minimb(*nll);
	minimb.setMinimizerType(minimizer_.c_str());
	minimb.setPrintLevel(-1);
	for(int ii=minStrategy_; ii<3; ++ii) {
		minimb.setStrategy(ii);
		if( ! minimb.migrad() ) { break; }
	}
	double minNllb = nll->getVal();     
	
	double qA = -2.*(minNllb - minNllsb);
	
	for(size_t ii=0; ii<garbageColl.size(); ++ii) {
		delete garbageColl[ii];
	}
	for(size_t ii=0; ii<garbageData.size(); ++ii) {
		delete garbageData[ii];
	}
	
	if( debug_ ) {
		for(int icat=0; icat<ncat; ++icat) {
			RooPlot* frame = frames[icat];
			pdfs[icat].plotOn(frame,RooFit::LineColor(kRed));
			TCanvas * canvas = new TCanvas(Form("cat_%d_%d",ncat,icat),Form("cat_%d_%d",ncat,icat));
			canvas->cd();
			frame->Draw();
			canvas->SaveAs(Form("cat_%d_%d.png",ncat,icat));
		}
	}

	return -sqrt(-qA);
}


