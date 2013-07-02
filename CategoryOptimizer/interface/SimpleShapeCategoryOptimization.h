#ifndef _SimpleShapeCategoryOptimization_h_
#define _SimpleShapeCategoryOptimization_h_

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Math/IFunction.h"
#include "RooRealVar.h"

#include <algorithm>

#include "CategoryOptimizer/interface/CategoryOptimizer.h"
#include "CategoryOptimizer/interface/FunctionHelpers.h"

// ------------------------------------------------------------------------------------------------
class SecondOrderModel : public AbsModel 
{
public:
	enum shape_t { automatic=0, gaus, expo };
	
	SecondOrderModel(std::string name, RooRealVar * x, AbsModel::type_t type=AbsModel::sig,
			 RooRealVar * mu=0, shape_t shape=automatic);
	~SecondOrderModel();

	RooRealVar * getX() { return x_; };
	void setMu(RooRealVar *mu) { mu_ = mu; };

	virtual void clear() { 
		categoryYields_.clear(); 
		categoryMeans_.clear();
		categoryRMSs_.clear();
	};
	
	void addCategory(double counts, double mean, double rms) { 
		categoryYields_.push_back(counts); 
		categoryMeans_.push_back(mean);
		categoryRMSs_.push_back(rms);
	};
	
	RooAbsPdf * getCategoryPdf(int icat);
	
	void dump();
	
	void buildPdfs();
	
	std::string name() { return name_; };

	void setShape(shape_t x);
	shape_t getShape() { return shape_; };

	void minEvents(double x) { minEvents_ = x; };
	double minEvents() { return minEvents_; };

private:
	void bookShape(int icat);
	void setShapeParams(int icat);
	
	std::string name_;
	RooRealVar *x_, *mu_;
	shape_t shape_;
	std::vector<double> categoryMeans_, categoryRMSs_;
	std::vector<RooAbsPdf *> categoryPdfs_;
	std::vector<RooRealVar *> categoryNorms_;
	RooArgSet owned_;
	TF1 * likeg_;
	double minEvents_;
};

// ------------------------------------------------------------------------------------------------
class SecondOrderModelBuilder : public AbsModelBuilder
{
public:
	SecondOrderModelBuilder(AbsModel::type_t type, std::string name,
				RooRealVar * x, 
				TH1 * pdf, TH1 * var, TH1 * var2, 
				double norm, double min,  double max) :
		ndim_(1),norm_(norm),
		pdf_(pdf), hsparse_(0),
		converterN_(new HistoToTF1(Form("%s_integrator",pdf->GetName()),integrate1D(pdf))),
		converterX_(new HistoToTF1(Form("%s_integrator",var->GetName()),integrate1D(var,false))),
		converterX2_(new HistoToTF1(Form("%s_integrator",var2->GetName()),integrate1D(var2,false))),
		model_(name,x,type) {
		ranges_.push_back(std::make_pair(min,max));
	};
	
	SecondOrderModelBuilder(AbsModel::type_t type, RooRealVar * x,
				std::string name,
				TH2 * pdf, TH2 * var, TH2 *var2, 
				double norm, double xmin, double xmax, 
				double ymin, double ymax)  : 
		ndim_(2),norm_(norm),
		pdf_(pdf), hsparse_(0),
		converterN_(new SimpleHistoToTF2(Form("%s_integrator",pdf->GetName()),integrate2D(pdf))),
		converterX_(new SimpleHistoToTF2(Form("%s_integrator",var->GetName()),integrate2D(var,false))),
		converterX2_(new SimpleHistoToTF2(Form("%s_integrator",var2->GetName()),integrate2D(var2,false))),
		model_(name,x,type) {
		ranges_.push_back(std::make_pair(xmin,xmax));
		ranges_.push_back(std::make_pair(ymin,ymax));
	};

	SecondOrderModelBuilder(AbsModel::type_t type, std::string name, RooRealVar * x, 
				TTree * tree, const RooArgList * varlist, const RooArgList * sellist,
				const char * weightBr);

	~SecondOrderModelBuilder() {
		if( converterN_ ) { delete converterN_; }
		if( converterX_ ) { delete converterX_; }
		if( converterX2_ ) { delete converterX2_; }
		if( hsparse_ ) { delete hsparse_; }
	};
	
	TTree * getTree();

	AbsModel * getModel() { return &model_; };
	void beginIntegration(double * boundaries) { 
		std::vector<double> extboundaries(ndim_+selectionCuts_.size());
		std::copy(boundaries,boundaries+ndim_,extboundaries.begin());
		std::copy(selectionCutsBegin_.begin(),selectionCutsBegin_.end(),extboundaries.begin()+ndim_);
		lastIntegral_ = (*converterN_)(&extboundaries[0],0);
		lastSumX_     = (*converterX_)(&extboundaries[0],0);
		lastSumX2_    = (*converterX2_)(&extboundaries[0],0);
		
		model_.clear(); 
	};
	
	void endIntegration() {
		model_.buildPdfs();
	};

	void setOrthoCuts(double * cuts) { std::copy(cuts,cuts+selectionCuts_.size(),selectionCuts_.begin()); };
	
	bool addBoundary(double * boundaries) {
		bool ret = true;
		std::vector<double> extboundaries(ndim_+selectionCuts_.size());
		std::copy(boundaries,boundaries+ndim_,extboundaries.begin());
		std::copy(selectionCuts_.begin(),selectionCuts_.end(),extboundaries.begin()+ndim_);
		double integral = (*converterN_)(&extboundaries[0],0);
		double sumX     = (*converterX_)(&extboundaries[0],0);
		double sumX2    = (*converterX2_)(&extboundaries[0],0);
		double norm     = norm_*(integral-lastIntegral_); 
		double mean     = (sumX - lastSumX_)/norm;
		double rms      = (sumX2 - lastSumX2_)/norm;
		///// std::copy( extboundaries.begin(), extboundaries.end(), std::ostream_iterator<double>(std::cout, ",") );
		///// std::cout << " " << integral << std::endl;
		if( mean*mean - rms > 0. ) {
			rms = sqrt( mean*mean - rms );
		} else {
			rms = 1.e-2*mean;
			/// norm = 0.;
			///// if( model_.getShape() == SecondOrderModel::gaus ) { 
			///// 	ret = false;
			///// 	penalty_ = 0.;
			///// }
		}
		if( norm <= model_.minEvents()*1.02 ) { 
			std::cout << " too few events  " << norm << " " << model_.minEvents() <<std::endl;
			penalty_ = norm/model_.minEvents(); 
			ret = false; 
		}
		model_.addCategory(norm, mean, rms);
		lastIntegral_ = integral;
		lastSumX_  = sumX;
		lastSumX2_ = sumX2;
		return ret;
	};
	
	double getMin(int idim) { return ranges_[idim].first;  };
	double getMax(int idim) { return ranges_[idim].second; };
	
	HistoConverter * getInputModelN() { return converterN_; };
	HistoConverter * getInputModelX() { return converterX_; };
	HistoConverter * getInputModelX2() { return converterX2_; };

	TF1 * getTF1N()  { return new TF1(Form("tf1N%s",model_.name().c_str()),converterN_ ,getMin(0),getMax(0),0); };
	TF1 * getTF1X()  { return new TF1(Form("tf1X%s",model_.name().c_str()),converterX_ ,getMin(0),getMax(0),0); };
	TF1 * getTF1X2() { return new TF1(Form("tf1X2%s",model_.name().c_str()),converterX2_,getMin(0),getMax(0),0); };
	
	TH1 * getPdf(int idim);
	
	double getPenalty() { return penalty_; };
	
private:
	SecondOrderModel model_;
	int ndim_;
	double norm_, lastIntegral_, lastSumX_, lastSumX2_, penalty_;
	std::vector<std::pair<double,double> > ranges_;
	std::vector<double> selectionCuts_, selectionCutsBegin_;
	
	TH1 * pdf_;
	THnSparse * hsparse_;
	HistoConverter *converterN_, *converterX_, *converterX2_;
};


// ------------------------------------------------------------------------------------------------
class SimpleShapeFomProvider : public AbsFomProvider 
{
public:
	SimpleShapeFomProvider(int nSubcats=1,RooRealVar *poi=0, int ncpu=4, const char * minimizer="Minuit2", 
			       int minStrategy=2) : 
		nSubcats_(nSubcats), ncpu_(ncpu), minimizer_(minimizer),minStrategy_(minStrategy), useRooSimultaneous_(false)
		{ 
			if( poi ) { addPOI(poi); }
			assert(minStrategy_<3); 
		};
		
	double operator()  ( std::vector<AbsModel *> sig, std::vector<AbsModel *> bkg) const;
	void addPOI(RooRealVar * poi) { pois_.push_back(poi); };
	void addResetP(RooRealVar * p) { resets_.push_back(p); };
	
	void minStrategy(int x) { minStrategy_=x; };
	void minimizer(const char * x) { minimizer_=x; };
	void useRooSimultaneous(bool x=true) { useRooSimultaneous_=x; };
	void nSubcats(int x) { nSubcats_ = x; };
	
	void addNuisance(RooRealVar * x, RooAbsReal *p=0) { 
		resets_.push_back(x); 
		if(p!=0) { 
			constrained_.add(*x);
			constraints_.add(*p,false); 
		} 
	};

private:
	int nSubcats_;
	int ncpu_;
	std::vector<RooRealVar *> pois_;
	std::vector<RooRealVar *> resets_;
	RooArgSet constraints_, constrained_;
	std::string minimizer_;
	int minStrategy_;
	bool useRooSimultaneous_;

};


#endif
