#ifndef _SimpleShapeCategoryOptimization_h_
#define _SimpleShapeCategoryOptimization_h_

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Math/IFunction.h"

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
};

// ------------------------------------------------------------------------------------------------
class SecondOrderModelBuilder : public AbsModelBuilder
{
public:
	SecondOrderModelBuilder(AbsModel::type_t type, std::string name,
				RooRealVar * x, 
				TH1 * pdf, TH1 * var, TH1 * var2, 
				double norm, double min,  double max) :
		norm_(norm),
		pdf_(pdf),
		converterN_(new HistoToTF1(Form("%s_integrator",pdf->GetName()),integrate1D(pdf))),
		converterX_(new HistoToTF1(Form("%s_integrator",var->GetName()),integrate1D(var,false))),
		converterX2_(new HistoToTF1(Form("%s_integrator",var2->GetName()),integrate1D(var2,false))),
		model_(name,x,type) {
		std::cout << pdf_ << " "  << converterN_ << " " << converterX_ << " " << converterX2_ << std::endl;
		ranges_.push_back(std::make_pair(min,max));
	};
	
	SecondOrderModelBuilder(AbsModel::type_t type, RooRealVar * x,
				std::string name,
				TH2 * pdf, TH2 * var, TH2 *var2, 
				double norm, double xmin, double xmax, 
				double ymin, double ymax)  : 
		norm_(norm),
		pdf_(pdf),
		converterN_(new SimpleHistoToTF2(Form("%s_integrator",pdf->GetName()),integrate2D(pdf))),
		converterX_(new SimpleHistoToTF2(Form("%s_integrator",var->GetName()),integrate2D(var,false))),
		converterX2_(new SimpleHistoToTF2(Form("%s_integrator",var2->GetName()),integrate2D(var2,false))),
		model_(name,x,type) {
		ranges_.push_back(std::make_pair(xmin,xmax));
		ranges_.push_back(std::make_pair(ymin,ymax));
	};

	~SecondOrderModelBuilder() {
		delete converterN_;
		delete converterX_;
		delete converterX2_;
	};
	
	AbsModel * getModel() { return &model_; };
	void beginIntegration(double * boundaries) { 
		lastIntegral_ = (*converterN_)(boundaries,0);
		lastSumX_     = (*converterX_)(boundaries,0);
		lastSumX2_    = (*converterX2_)(boundaries,0);
		
		model_.clear(); 
	};
	
	void endIntegration() {
		model_.buildPdfs();
	};

	void addBoundary(double * boundaries) {
		double integral = (*converterN_)(boundaries,0);
		double sumX     = (*converterX_)(boundaries,0);
		double sumX2    = (*converterX2_)(boundaries,0);
		double norm     = norm_*(integral-lastIntegral_); 
		double mean     = (sumX - lastSumX_)/norm;
		double rms      = (sumX2 - lastSumX2_)/norm;
		if( mean*mean - rms > 0. ) {
			rms = sqrt( mean*mean - rms );
		} else {
			rms = 0.;
			norm = 0.;
		}
		model_.addCategory(norm, mean, rms);
		lastIntegral_ = integral;
		lastSumX_  = sumX;
		lastSumX2_ = sumX2;
	};
	
	double getMin(int idim) { return ranges_[idim].first;  };
	double getMax(int idim) { return ranges_[idim].second; };
	
	HistoConverter * getInputModelN() { return converterN_; };
	HistoConverter * getInputModelX() { return converterX_; };
	HistoConverter * getInputModelX2() { return converterX2_; };
	
	TH1 * getPdf(int idim);
	
private:
	SecondOrderModel model_;
	double norm_, lastIntegral_, lastSumX_, lastSumX2_;
	std::vector<std::pair<double,double> > ranges_;
	
	TH1 * pdf_;
	HistoConverter *converterN_, *converterX_, *converterX2_;
};


// ------------------------------------------------------------------------------------------------
class SimpleShapeFomProvider : public AbsFomProvider 
{
public:
	SimpleShapeFomProvider(RooRealVar *poi=0, int ncpu=4, const char * minimizer="Minuit2", 
			       int minStrategy=2) : 
		ncpu_(ncpu), minimizer_(minimizer),minStrategy_(minStrategy), useRooSimultaneous_(false)
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
	
	
private:
	int ncpu_;
	std::vector<RooRealVar *> pois_;
	std::vector<RooRealVar *> resets_;
	std::string minimizer_;
	int minStrategy_;
	bool useRooSimultaneous_;

};


#endif
